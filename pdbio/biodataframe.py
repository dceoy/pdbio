#!/usr/bin/env python
"""
Pandas-based Data Handler for VCF and BED Files.
https://github.com/dceoy/pdbio
"""

import bz2
import gzip
import io
import logging
import os
from abc import ABCMeta, abstractmethod
from collections import OrderedDict
from itertools import product

import pandas as pd


class BaseBioDataFrame(object, metaclass=ABCMeta):
    """Base DataFrame handler."""

    def __init__(self, path, supported_exts=None):
        if os.path.isfile(path):
            self.path = os.path.abspath(
                os.path.expanduser(os.path.expandvars(path))
            )
        else:
            raise FileNotFoundError('file not found: {}'.format(path))
        if (isinstance(supported_exts, (list, tuple))
                and not [x for x in supported_exts if path.endswith(x)]):
            raise ValueError('invalid file extension: {}'.format(path))
        else:
            pass
        self.header = list()
        self.df = pd.DataFrame()

    @abstractmethod
    def load(self):
        pass

    def load_and_output_df(self):
        self.load()
        return self.df

    @staticmethod
    def open_readable_file(path):
        if path.endswith('.gz'):
            return gzip.open(path, 'rt')
        elif path.endswith('.bz2'):
            return bz2.open(path, 'rt')
        else:
            return open(path, 'r')

    def write_header(self, path):
        if self.header:
            with open(path, mode='w') as f:
                for h in self.header:
                    f.write(h + os.linesep)


class VcfDataFrame(BaseBioDataFrame):
    """VCF DataFrame handler."""

    def __init__(self, path, return_df=False):
        super().__init__(
            path=path,
            supported_exts=[
                (e + c) for e, c in
                product(['.vcf', '.txt', '.tsv'], ['', '.gz', '.bz2'])
            ]
        )
        self.__logger = logging.getLogger(__name__)
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
                    ('SAMPLE{}'.format(i), n) for i, n in enumerate(samples)
                ])
                n_fixed_cols = len(self.__fixed_cols)
                n_detected_cols = len(items)
                self.__detected_cols = self.__fixed_cols + (
                    [
                        'SAMPLE{}'.format(i)
                        for i in range(n_detected_cols - n_fixed_cols)
                    ] if n_detected_cols > n_fixed_cols else list()
                )
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

    def write_vcf(self, path):
        self.__logger.info('Write a VCF file: {}'.format(path))
        self.write_header(path=path)
        self.df.pipe(lambda d: d.rename(self.sample_dict)).to_csv(
            path, mode=('a' if self.header else 'w'), sep='\t', index=False
        )


class BedDataFrame(BaseBioDataFrame):
    """BED DataFrame handler."""

    def __init__(self, path, opt_cols=None, return_df=False):
        super().__init__(
            path=path,
            supported_exts=[
                (e + c) for e, c in
                product(['.bed', '.txt', '.tsv'], ['', '.gz', '.bz2'])
            ]
        )
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
        if return_df:
            self.load_and_output_df()
        else:
            self.load()

    def load(self):
        self.__logger.info('Load a BED file: {}'.format(self.path))
        with self.open_readable_file(path=self.path) as f:
            for s in f:
                self._load_bed_line(string=s)
        self.df = self.df.reset_index(drop=True)

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
                )
            )

    def write_bed(self, path):
        self.__logger.info('Write a BED file: {}'.format(path))
        self.write_header(path=path)
        self.df.to_csv(
            path, header=False, mode=('a' if self.header else 'w'), sep='\t',
            index=False
        )
