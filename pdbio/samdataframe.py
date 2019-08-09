#!/usr/bin/env python
"""
Pandas-based SAM Data Handler.
https://github.com/dceoy/pdbio
"""

import logging
import os
import re
import subprocess
import sys
from collections import OrderedDict
from multiprocessing import cpu_count

import numpy as np
import pandas as pd

from .biodataframe import BaseBioDataFrame


class SamDataFrame(BaseBioDataFrame):
    """SAM DataFrame handler."""

    def __init__(self, path=None, samtools=None, n_thread=None, rname=None,
                 startpos=None, endpos=None, load=True):
        self.__logger = logging.getLogger(__name__)
        self.__samtools = samtools
        self.__n_thread = n_thread
        if rname and startpos and endpos:
            self.__region = '{0}:{1:d}-{2:d}'.format(rname, startpos, endpos)
        elif rname and startpos:
            self.__region = '{0}:{1:d}'.format(rname, startpos)
        elif rname:
            self.__region = rname
        else:
            self.__region = None
        self.regions = [self.__region] if self.__region else None
        self.__col_dtypes = OrderedDict([
            ('QNAME', str), ('FLAG', int), ('RNAME', str), ('POS', int),
            ('MAPQ', int), ('CIGAR', str), ('RNEXT', str), ('PNEXT', int),
            ('TLEN', int), ('SEQ', str), ('QUAL', str), ('OPT', str)
        ])
        self.__cols = list(self.__col_dtypes.keys())
        self.__n_cols = len(self.__cols)
        super().__init__(
            path=path, format_name='SAM', delimiter='\t', column_header=False,
            chrom_column='RNAME', pos_columns=['POS'],
            txt_file_exts=['.sam', '.txt', '.tsv'],
            bin_file_exts=['.bam', '.cram'], load=load
        )

    def load_table(self):
        if self.path.endswith(('.bam', '.cram')) or self.regions:
            self.df = self.convert_lines_to_df(
                lines=list(self.view(options=['-h'], regions=self.regions))
            )
        else:
            with self.open_readable_file(path=self.path) as f:
                self.df = self.convert_lines_to_df(lines=list(f))
        return self

    def parse_line(self, string, into_ordereddict=False):
        if re.search(r'^@[A-Z]{2}', string):
            self.header.append(string.strip())
        elif string.strip():
            items = string.strip().split('\t', maxsplit=(self.__n_cols - 1))
            df = pd.DataFrame(
                [[*items, ''] if len(items) < self.__n_cols else items],
                columns=self.__cols
            ).astype(dtype=self.__col_dtypes)
            if into_ordereddict:
                return df.iloc[0].to_dict(into=OrderedDict)
            else:
                return df

    def view(self, options=None, regions=None):
        args = [
            (self.__samtools or self.fetch_executable('samtools')), 'view',
            '-@', str(self.__n_thread or cpu_count()),
            *(options if options else list()), self.path,
            *(regions if regions else list())
        ]
        for s in self.run_and_parse_subprocess(args=args):
            yield s

    def write_body(self, path=None, mode='a'):
        if path:
            with open(path, mode=mode) as f:
                for s in self._generate_samline():
                    f.write(s)
        else:
            for s in self._generate_samline():
                sys.stdout.write(s)
                sys.stdout.flush()

    def _generate_samline(self):
        for s in self.df.apply(lambda r:
                               r.dropna().astype(str).str.cat(sep='\t'),
                               axis=1):
            yield (s.strip() + os.linesep)

    def tag_expanded_df(self, df=None, tag=None, cast_numeric_types=True,
                        drop=True):
        df_s = self.df if df is None else df
        if tag:
            tag_strs = df_s['OPT'].str.split('\t').apply(
                lambda t: ([s for s in t if s.startswith(tag)] or [np.nan])[0]
            )
            if tag_strs.isnull().all():
                raise RuntimeError('tag not found: {}'.format(tag))
            else:
                tag_col = tag_strs.dropna().iloc[0][:4]
                tag_vals = tag_strs.str.split(':', n=2).apply(
                    lambda t: (t[2] if t else t)
                )
                return df_s.assign(
                    **{
                        tag_col: (
                            tag_vals.astype(float) if tag_col[3] in 'if'
                            else tag_vals
                        )
                    }
                )
        else:
            return self._all_tagged_df(
                df=df_s, cast_numeric_types=cast_numeric_types, drop=drop
            )

    @staticmethod
    def _all_tagged_df(df, cast_numeric_types=True, drop=True):
        return pd.concat(
            [
                pd.DataFrame([{
                    'id': i, **({s[:4]: s[5:] for s in v} if v[0] else dict())
                }]) for i, v in df['OPT'].str.split('\t').items()
            ],
            ignore_index=True, sort=False
        ).set_index('id').pipe(
            lambda d: (
                d.astype(
                    dtype={c[:4]: float for c in d.columns if c[3] in 'if'}
                ) if cast_numeric_types else d
            )
        ).pipe(lambda d: df.join(d, how='left')).pipe(
            lambda d: (d.drop(columns='OPT') if drop else d)
        )

    def load_sam_by_region(self, rname, startpos, endpos,
                           completely_inclusion=True, options=None,
                           into_ordereddict=False):
        self.__logger.info('Load SAM file by region: {}'.format(self.path))
        if startpos > endpos:
            raise ValueError(
                'invalid range (startpos > endpos): {0}:{1}-{2}'.format(
                    rname, startpos, endpos
                )
            )
        elif startpos == endpos or not completely_inclusion:
            lines = list(
                self.view(
                    options=options,
                    regions=['{0}:{1:d}-{2:d}'.format(rname, startpos, endpos)]
                )
            )
        else:
            edges = [
                ['{0}:{1:d}-{1:d}'.format(rname, p)]
                for p in [startpos, endpos]
            ]
            lines = [
                s for s in self.view(options=options, regions=edges[0])
                if s in set(self.view(options=options, regions=edges[1]))
            ]
        self.header = list(self.view(options=['-H']))
        self.df = self.convert_lines_to_df(lines=lines, update_header=False)
        self.__logger.debug('self.df shape: {}'.format(self.df.shape))
        return self

    def median_depth(self, region=None, min_baseq=None, min_mapq=None):
        depths = self._depth_array(
            region=region, min_baseq=min_baseq, min_mapq=min_mapq
        )
        return (np.median(depths) if depths.size else 0)

    def describe_depth(self, region=None, min_baseq=None, min_mapq=None):
        return self._depth_array(
            region=region, min_baseq=min_baseq, min_mapq=min_mapq
        ).describe()

    def _depth_array(self, region=None, min_baseq=None, min_mapq=None):
        args = [
            (self.__samtools or self.fetch_executable('samtools')), 'depth',
            *(['-q', str(min_baseq)] if min_baseq is not None else list()),
            *(['-Q', str(min_mapq)] if min_mapq is not None else list()),
            *(['-r', region] if region else list()), self.path
        ]
        return pd.Series([
            s.strip().split('\t', maxsplit=2)[2]
            for s in self.run_and_parse_subprocess(args=args)
        ]).astype(np.int32)

    @staticmethod
    def cigar2reflen(cigar):
        """
        Args:
            cigar (str): CIGAR string of SAM

        Returns:
            int: length of a consumed reference bases

        Examples:
            >>> cigar2reflen(cigar='26S15M4D32M3I75M')
            126
        """
        return np.array([
            (s or '0')
            for s in re.split('M|D|N|=|X', re.sub(r'[0-9]+[ISHP]', '', cigar))
        ]).astype(np.int16).sum()

    @staticmethod
    def cigar2qlen(cigar):
        """
        Args:
            cigar (str): CIGAR string of SAM

        Returns:
            int: length of a consumed query bases

        Examples:
            >>> cigar2reflen(cigar='26S15M4D32M3I75M')
            151
        """
        return np.array([
            (s or '0')
            for s in re.split('M|I|S|=|X', re.sub(r'[0-9]+[DNHP]', '', cigar))
        ]).astype(np.int16).sum()

    @staticmethod
    def cigar2oplen(cigar, sum=True):
        """
        Args:
            cigar (str): CIGAR string of SAM
            sum (bool): sum lengths by character

        Returns:
            dict: total lengths of CIGAR operations

        Examples:
            >>> cigar2oplen(cigar='26S15M4D32M3I75M')
            {'M': 122, 'I': 3, 'D': 4, 'N': 0, 'S': 26, 'H': 0, 'P': 0, '=': 0,
             'X': 0}
        """
        ops = [
            (k, np.int16(v)) for k, v in zip(
                re.sub(r'[0-9]', '', cigar),
                re.split('M|I|D|N|S|H|P|=|X', cigar)[:-1]
            )
        ]
        if sum:
            return (lambda o: {k: (o.get(k) or 0) for k in 'MIDNSHP=X'})(
                o=pd.DataFrame(
                    ops, columns=['op', 'len']
                ).groupby('op')['len'].sum().to_dict()
            )
        else:
            return ops

    @staticmethod
    def cigar2alnrange(cigar):
        """
        Args:
            cigar (str): CIGAR string of SAM

        Returns:
            tuple: aligned sequence range of read

        Examples:
            >>> cigar2alnrange(cigar='26S15M4D32M3I75M')
            (26, 150)
        """
        ops = [
            (k, np.int16(v)) for k, v in zip(
                re.sub(r'[0-9]', '', cigar),
                re.split('M|I|D|N|S|H|P|=|X', cigar)[:-1]
            )
        ]
        seq_range = [0, sum([v for k, v in ops if k in 'MIS=X']) - 1]
        for k, v in ops:
            if k in 'ISHP':
                seq_range[0] += v
            else:
                break
        for k, v in reversed(ops):
            if k in 'ISHP':
                seq_range[1] -= v
            else:
                break
        return tuple(seq_range)

    @staticmethod
    def cigar2chars(cigar, only_aligned=False):
        """
        Args:
            cigar (str): CIGAR string of SAM
            only_aligned (bool): cut unmapped characters on edges

        Returns:
            str: CIGAR character sequence

        Examples:
            >>> cigar2reflen(cigar='6S5M4D12M3I5M')
            'SSSSSSMMMMMDDDDMMMMMMMMMMMMIIIMMMMM'
        """
        chars = ''.join([
            (k * int(v)) for k, v in zip(
                re.sub(r'[0-9]', '', cigar),
                re.split('M|I|D|N|S|H|P|=|X', cigar)[:-1]
            )
        ])
        if only_aligned:
            return re.sub(r'^[ISHP]+', '', re.sub(r'[ISHP]+$', '', chars))
        else:
            return chars

    @staticmethod
    def cigar2matchchars(cigar, md):
        """
        Args:
            cigar (str): CIGAR string of SAM
            md (str): MD tag value of SAM

        Returns:
            str: CIGAR character sequence
                (match: '=', unmatch: 'X', deletion: 'D')

        Examples:
            >>> cigar2matchchars(cigar='11M1D16M', md='5C5^T16')
            '=====X=====D================'
        """
        aln_cigar_list = list(
            re.sub(
                r'^[ISHP]+', '',
                re.sub(
                    r'[ISHP]+$', '',
                    ''.join([
                        (k * int(v)) for k, v in zip(
                            re.sub(r'[0-9]', '', cigar),
                            re.split('M|I|D|N|S|H|P|=|X', cigar)[:-1]
                        )
                    ])
                )
            )
        )
        not_match_md = [
            (('D' * len(s[1:])) if s.startswith('^') else ('X' if s else s))
            for s in re.split(r'[0-9]+', md)
        ]
        n_not_match_md = len(not_match_md)
        match_md = [
            ('=' * int(s)) for s in re.split(r'[\^A-Z]', md) if s
        ]
        if n_not_match_md > len(match_md):
            match_md.append('')
        assert n_not_match_md == len(match_md)
        md_cigars = ''.join([(u + m) for m, u in zip(match_md, not_match_md)])
        ri = 0
        for i, c in enumerate(aln_cigar_list):
            if c == 'M':
                aln_cigar_list[i] = md_cigars[ri]
                ri += 1
            elif c in 'DN=X':
                ri += 1
        return ''.join(aln_cigar_list)


class BamFileWriter(object):
    def __init__(self, out_bam_path, samtools, n_thread=1):
        self.__logger = logging.getLogger(__name__)
        self.__bam_path = out_bam_path
        self.__logger.debug('Write STDIN into BAM: {}'.format(self.__bam_path))
        self.__samtools = samtools
        args = [
            self.__samtools, 'view', '-@', str(n_thread), '-bS', '-', '-o',
            out_bam_path
        ]
        self.__logger.debug('STDIN => `{}`'.format(' '.join(args)))
        self.proc = subprocess.Popen(
            args=args, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )

    def write(self, *args, **kwargs):
        self.proc.stdin.write(*args, **kwargs)
        self.proc.stdin.flush()

    def close(self):
        self.__logger.debug('Finish writing BAM: {}'.format(self.__bam_path))
        self.proc.stdin.close()
        self.proc.wait()
        if self.proc.returncode == 0:
            subprocess.run([self.__samtools, 'quickcheck', self.__bam_path])
        else:
            self.__logger.error(
                'STDERR from subprocess `{0}`:{1}{2}'.format(
                    self.proc.args, os.linesep,
                    self.proc.stderr.decode('utf-8')
                )
            )
            raise subprocess.CalledProcessError(
                returncode=self.proc.returncode, cmd=self.proc.args,
                output=self.proc.stdout, stderr=self.proc.stderr
            )

    def __enter__(self):
        return self

    def __exit__(self, ex_type, ex_value, trace):
        args = {'ex_type': ex_type, 'ex_value': ex_value, 'trace': trace}
        self.__logger.debug('with-statement exit: {}'.format(args))
        self.close()

    def __del__(self):
        self.close()
