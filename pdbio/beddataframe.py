#!/usr/bin/env python
"""
Pandas-based BED Data Handler.
https://github.com/dceoy/pdbio
"""

import logging
import sys
from collections import OrderedDict

import pandas as pd

from .biodataframe import BaseBioDataFrame


class BedDataFrame(BaseBioDataFrame):
    """BED DataFrame handler."""

    def __init__(self, path=None, opt_cols=None, load=True):
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

    def load_table(self):
        with self.open_readable_file(path=self.path) as f:
            self.df = self.convert_lines_to_df(lines=list(f))
        return self

    def parse_line(self, string, into_ordereddict=False):
        if string.startswith(('browser', 'track')):
            self.header.append(string.strip())
        elif string:
            if not self.__detected_cols:
                self.__detected_cols = [
                    *self.__fixed_cols, *self.__opt_cols
                ][:(string.count('\t') + 1)]
                self.__logger.debug(
                    'self.__detected_cols: {}'.format(self.__detected_cols)
                )
                self.__detected_col_dtypes = {
                    k: (self.__fixed_col_dtypes.get(k) or str)
                    for k in self.__detected_cols
                }
                self.__logger.debug(
                    'self.__detected_col_dtypes: {}'.format(
                        self.__detected_col_dtypes
                    )
                )
            df = pd.DataFrame(
                [string.strip().split('\t')], columns=self.__detected_cols
            ).astype(dtype=self.__detected_col_dtypes)
            if into_ordereddict:
                return df.iloc[0].to_dict(into=OrderedDict)
            else:
                return df

    def write_body(self, path=None, mode='a', **kwargs):
        self.df.to_csv(
            (path or sys.stdout), mode=mode, index=False, sep='\t', **kwargs
        )
