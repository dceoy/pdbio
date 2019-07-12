#!/usr/bin/env python
"""
Pandas-based SAM Data Handler.
https://github.com/dceoy/pdbio
"""

import io
import logging
import re
from multiprocessing import cpu_count

import pandas as pd

from .biodataframe import BaseBioDataFrame


class SamDataFrame(BaseBioDataFrame):
    """SAM DataFrame handler."""

    def __init__(self, path, samtools=None, n_thread=None):
        self.__logger = logging.getLogger(__name__)
        self.__samtools = samtools
        self.__n_thread = n_thread
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
            bin_file_exts=['.bam', '.cram']
        )

    def load(self):
        if self.path.endswith(('.bam', '.cram')):
            args = [
                (self.__samtools or self.fetch_executable('samtools')), 'view',
                '-@', str(self.__n_thread or cpu_count()), '-h', self.path
            ]
            for s in self.run_and_parse_subprocess(args=args):
                self._load_sam_line(string=s)
        else:
            with self.open_readable_file(path=self.path) as f:
                for s in f:
                    self._load_sam_line(string=s)

    def _load_sam_line(self, string):
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
            self.df = self.df.append(
                pd.read_csv(
                    io.StringIO(string), sep='\t', header=None,
                    names=self.__detected_cols,
                    dtype=self.__detected_col_dtypes
                ),
                ignore_index=True
            )
