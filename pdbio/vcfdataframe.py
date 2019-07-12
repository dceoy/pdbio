#!/usr/bin/env python
"""
Pandas-based VCF Data Handler.
https://github.com/dceoy/pdbio
"""

import io
import logging
from itertools import chain
from multiprocessing import cpu_count

import pandas as pd

from .biodataframe import BaseBioDataFrame


class VcfDataFrame(BaseBioDataFrame):
    """VCF DataFrame handler."""

    def __init__(self, path=None, bcftools=None, n_thread=1):
        self.__logger = logging.getLogger(__name__)
        self.__bcftools = bcftools
        self.__n_thread = n_thread
        self.__fixed_cols = [
            '#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO'
        ]
        self.__opt_cols = ['FORMAT']
        self.__fixed_col_dtypes = {
            '#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
            'QUAL': str, 'FILTER': str, 'INFO': str
        }
        self.__detected_cols = list()
        self.__detected_col_dtypes = dict()
        self.samples = list()
        super().__init__(
            path=path, format_name='VCF', delimiter='\t', column_header=True,
            chrom_column='#CHROM', txt_file_exts=['.vcf', '.txt', '.tsv'],
            bin_file_exts=['.bcf']
        )

    def load(self):
        if self.path.endswith('.bcf'):
            args = [
                (self.__bcftools or self.fetch_executable('bcftools')), 'view',
                '--threads', str(self.__n_thread or cpu_count()), self.path
            ]
            for s in self.run_and_parse_subprocess(args=args):
                self._load_vcf_line(string=s)
        else:
            with self.open_readable_file(path=self.path) as f:
                for s in f:
                    self._load_vcf_line(string=s)

    def _load_vcf_line(self, string):
        if string.startswith('##'):
            self.header.append(string.strip())
        elif string.startswith('#CHROM'):
            items = string.strip().split('\t')
            self.__detected_cols = items
            self.__detected_col_dtypes = {
                k: (self.__fixed_col_dtypes.get(k) or str) for k in items
            }
            for i, c in enumerate(items[:len(self.__fixed_cols)]):
                assert c == self.__fixed_cols[i], (
                    'invalid VCF column: {}'.format(c)
                )
            self.samples = [
                s for s in items
                if s not in [*self.__fixed_cols, *self.__opt_cols]
            ]
        else:
            self.df = self.df.append(
                pd.read_csv(
                    io.StringIO(string), sep='\t', header=None,
                    names=self.__detected_cols,
                    dtype=self.__detected_col_dtypes
                ),
                ignore_index=True
            )

    def expand_info_col(self, df=None):
        self.__logger.info('Expand the INFO column.')
        return (df or self.df).pipe(
            lambda d: d.drop(columns='INFO').join(
                self._parse_info(df=d), how='left'
            )
        )

    @staticmethod
    def _parse_info(df):
        return pd.DataFrame([
            dict([
                ('index', id),
                *[tuple(s.split('=', maxsplit=1)) for s in string.split(':')]
            ]) for id, string in df['INFO'].iteritems()
        ]).set_index('index')

    def expand_samples_cols(self, df=None):
        self.__logger.info('Expand the columns of samples.')
        return (df or self.df).pipe(
            lambda d: d[self.__fixed_cols[:8]].join(
                self._parse_format_and_samples(df=d), how='left'
            )
        )

    @staticmethod
    def _parse_format_and_samples(df):
        format_key_list = df['FORMAT'].str.split(':').tolist()
        return pd.DataFrame([
            dict([
                ('index', id),
                *chain.from_iterable([
                    [
                        ('{0}_{1}'.format(n, k), v)
                        for k, v in zip(format_key_list[id], s.split(':'))
                    ] for n, s in row.items()
                ])
            ]) for id, row in df[df.columns[9:]].iterrows()
        ]).set_index('index')
