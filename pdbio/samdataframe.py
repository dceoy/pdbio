#!/usr/bin/env python
"""
Pandas-based SAM Data Handler.
https://github.com/dceoy/pdbio
"""

import logging
import os
import re
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
