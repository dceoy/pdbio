#!/usr/bin/env python
"""
Pandas-based VCF Data Handler.
https://github.com/dceoy/pdbio
"""

import io
import logging
from collections import OrderedDict
from itertools import chain, product
from multiprocessing import cpu_count

import pandas as pd

from .biodataframe import BaseBioDataFrame


class VcfDataFrame(BaseBioDataFrame):
    """VCF DataFrame handler."""

    def __init__(self, path, return_df=False, bcftools=None, n_thread=1):
        super().__init__(
            path=path,
            supported_exts=[
                *[
                    (e + c) for e, c in
                    product(['.vcf', '.txt', '.tsv'], ['', '.gz', '.bz2'])
                ], '.bcf'
            ]
        )
        self.__logger = logging.getLogger(__name__)
        self.__bcftools = bcftools
        self.__n_thread = n_thread
        self.__fixed_cols = [
            '#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO',
            'FORMAT'
        ]
        self.__fixed_col_dtypes = {
            '#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
            'QUAL': str, 'FILTER': str, 'INFO': str
        }
        self.__detected_cols = list()
        self.__detected_col_dtypes = dict()
        self.sample_dict = OrderedDict()
        if return_df:
            self.load_and_output_df()
        else:
            self.load()

    def load(self):
        self.__logger.info('Load a VCF file: {}'.format(self.path))
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
        self.df = self.df.reset_index(drop=True)

    def _load_vcf_line(self, string):
        if string.startswith('##'):
            self.header.append(string.strip())
        elif string.startswith('#CHROM'):
            items = string.strip().split('\t')
            if items[:len(self.__fixed_cols)] == self.__fixed_cols:
                samples = [s for s in items if s not in self.__fixed_cols]
                self.sample_dict = OrderedDict([
                    ('SAMPLE_{}'.format(i), n) for i, n in enumerate(samples)
                ])
                n_fixed_cols = len(self.__fixed_cols)
                n_detected_cols = len(items)
                self.__detected_cols = [
                    *self.__fixed_cols,
                    *[
                        'SAMPLE_{}'.format(i)
                        for i in range(max(n_detected_cols - n_fixed_cols, 0))
                    ]
                ]
                self.__detected_col_dtypes = {
                    k: (self.__fixed_col_dtypes.get(k) or str)
                    for k in self.__detected_cols
                }
            else:
                raise ValueError('invalid VCF columns')
        else:
            self.df = self.df.append(
                pd.read_csv(
                    io.StringIO(string), sep='\t', header=None,
                    names=self.__detected_cols,
                    dtype=self.__detected_col_dtypes
                )
            )

    def sort_df(self):
        self.__logger.info('Sort the VCF dataframe.')
        self.sort_df_by_chrom(chrom_col='#CHROM', add_cols=self.__fixed_cols)

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

    def write_vcf(self, path):
        self.__logger.info('Write a VCF file: {}'.format(path))
        self.write_header(path=path)
        self.df.rename(self.sample_dict).to_csv(
            path, mode=('a' if self.header else 'w'), sep='\t', index=False
        )
