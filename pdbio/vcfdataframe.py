#!/usr/bin/env python
"""
Pandas-based VCF Data Handler.
https://github.com/dceoy/pdbio
"""

import io
import logging
import re
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
        self.sample_dict = dict()
        super().__init__(
            path=path, format_name='VCF', delimiter='\t', column_header=True,
            chrom_column='#CHROM', pos_columns=['POS'],
            txt_file_exts=['.vcf', '.txt', '.tsv'], bin_file_exts=['.bcf']
        )

    def rename_samples_cols(self, prefix='SAMPLE_', sample_dict=None):
        self.__logger.info('Rename columns of samples')
        new_sample_dict = (
            sample_dict
            or {s: (prefix + str(i)) for i, s in enumerate(self.samples)}
        )
        self.__logger.debug('new_sample_dict: {}'.format(new_sample_dict))
        renaming_dict = (
            {v: new_sample_dict[k] for k, v in self.sample_dict.items()}
            if self.sample_dict else new_sample_dict
        )
        self.sample_dict = new_sample_dict
        self.df = self.df.rename(columns=renaming_dict)

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

    def expanded_df(self, df=None, by_info=True, by_samples=True, drop=True):
        df_x = (self.df if df is None else df)
        if by_info:
            df_x = self._expand_info_col(df=df_x, drop=drop)
        if by_samples:
            df_x = self._expand_samples_cols(df=df_x, drop=drop)
        return df_x

    def _expand_info_col(self, df=None, drop=True):
        self.__logger.info('Expand the INFO column.')
        return (self.df if df is None else df).pipe(
            lambda d: (
                d.drop(columns='INFO') if drop else d
            ).join(self._parse_info(df=d), how='left')
        )

    @staticmethod
    def _parse_info(df):
        return pd.DataFrame([
            dict([
                ('index', id),
                *[
                    (s.split('=', maxsplit=1)[0], re.sub(r'^[^=]+=?', '', s))
                    for s in val.split(';')
                ]
            ]) for id, val in df['INFO'].iteritems()
        ]).set_index('index').pipe(
            lambda d: d.rename(columns={k: 'INFO_' + k for k in d.columns})
        )

    def _expand_samples_cols(self, df=None, drop=True):
        self.__logger.info('Expand the columns of samples.')
        return (self.df if df is None else df).pipe(
            lambda d: (
                d.drop(columns=['FORMAT', *self.samples]) if drop else d
            ).join(self._parse_samples(df=d, samples=self.samples), how='left')
        )

    @staticmethod
    def _parse_samples(df, samples):
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
            ]) for id, row in df[samples].iterrows()
        ]).set_index('index')
