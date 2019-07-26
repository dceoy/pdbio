#!/usr/bin/env python
"""
Pandas-based SAM Data Handler.
https://github.com/dceoy/pdbio
"""

import io
import logging
import re
from collections import OrderedDict
from multiprocessing import cpu_count

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
        self.__fixed_cols = [
            'QNAME', 'FLAG', 'RNAME', 'POS', 'MAPQ', 'CIGAR', 'RNEXT', 'PNEXT',
            'TLEN', 'SEQ', 'QUAL'
        ]
        self.__fixed_col_dtypes = {
            'QNAME': str, 'FLAG': int, 'RNAME': str, 'POS': int, 'MAPQ': int,
            'CIGAR': str, 'RNEXT': str, 'PNEXT': int, 'TLEN': int, 'SEQ': str,
            'QUAL': str
        }
        self.__detected_cols = list()
        self.__detected_col_dtypes = dict()
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
        if re.match(r'@[A-Z]{1}', string):
            self.header.append(string.strip())
        else:
            if not self.__detected_cols:
                n_fixed_cols = len(self.__fixed_cols)
                n_detected_cols = string.count('\t') + 1
                self.__detected_cols = [
                    *self.__fixed_cols,
                    *[
                        'OPT{}'.format(i)
                        for i in range(max(n_detected_cols - n_fixed_cols, 0))
                    ]
                ]
                self.__detected_col_dtypes = {
                    k: (self.__fixed_col_dtypes.get(k) or str)
                    for k in self.__detected_cols
                }
            df = pd.read_csv(
                io.StringIO(string), sep='\t', header=None,
                names=self.__detected_cols, dtype=self.__detected_col_dtypes
            )
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

    def load_sam_by_region(self, rname, startpos, endpos,
                           completely_inclusion=True, options=None,
                           into_ordereddict=False):
        self.__logger.info('Load SAM file by region: {}'.format(self.path))
        if completely_inclusion:
            edges = [
                ['{0}:{1:d}-{1:d}'.format(rname, p)]
                for p in [startpos, endpos]
            ]
            lines = [
                s for s in self.view(options=options, regions=edges[0])
                if s in set(self.view(options=options, regions=edges[1]))
            ]
        else:
            lines = list(
                self.view(
                    options=options,
                    regions=['{0}:{1:d}-{2:d}'.format(rname, startpos, endpos)]
                )
            )
        self.df = self.convert_lines_to_df(lines=lines)
        self.__logger.debug('self.df shape: {}'.format(self.df.shape))
        return self
