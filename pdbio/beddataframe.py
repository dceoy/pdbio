#!/usr/bin/env python
"""
Pandas-based BED Data Handler.
https://github.com/dceoy/pdbio
"""

import io
import logging

import pandas as pd

from .biodataframe import BaseBioDataFrame


class BedDataFrame(BaseBioDataFrame):
    """BED DataFrame handler."""

    def __init__(self, path, opt_cols=None, load=True):
        self.__logger = logging.getLogger(__name__)
        self.__fixed_cols = ['chrom', 'chromStart', 'chromEnd']
        self.__opt_cols = opt_cols or [
            'name', 'score', 'strand', 'thickStart', 'thickEnd', 'itemRgb',
            'blockCount', 'blockSizes', 'blockStarts'
        ]
        self.__fixed_col_dtypes = {
            'chrom': str, 'chromStart': int, 'chromEnd': int, 'name': str,
            'score': int, 'strand': str, 'thickStart': int, 'thickEnd': int,
            'itemRgb': str, 'blockCount': int, 'blockSizes': int,
            'blockStarts': int
        }
        self.__detected_cols = list()
        self.__detected_col_dtypes = dict()
        super().__init__(
            path=path, format_name='BED', delimiter='\t', column_header=False,
            chrom_column='chrom', pos_columns=['chromStart', 'chromEnd'],
            txt_file_exts=['.bed', '.txt', '.tsv'], load=load
        )

    def load(self):
        with self.open_readable_file(path=self.path) as f:
            for s in f:
                self._load_bed_line(string=s)

    def _load_bed_line(self, string):
        if string.startswith(('browser', 'track')):
            self.header.append(string.strip())
        else:
            if not self.__detected_cols:
                self.__detected_cols = [
                    *self.__fixed_cols, *self.__opt_cols
                ][:(string.count('\t') + 1)]
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
