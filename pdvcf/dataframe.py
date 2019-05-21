#!/usr/bin/env python
#
# Pandas-based Data Handler for VCF Files
# https://github.com/dceoy/pdvcf

from abc import ABCMeta, abstractmethod
import bz2
import gzip
import io
import os
import pandas as pd


def fetch_abspath(path):
    return os.path.abspath(os.path.expanduser(os.path.expandvars(path)))


class BaseBioDataFrame(object, metaclass=ABCMeta):
    def __init__(self, path, supported_exts=[]):
        abspath = fetch_abspath(path=path)
        if os.path.isfile(abspath):
            self.path = abspath
        else:
            raise FileNotFoundError('file not found: {}'.format(abspath))
        hit_exts = [x for x in supported_exts if abspath.endswith(x)]
        if supported_exts and not hit_exts:
            raise ValueError('invalid file extension: {}'.format(abspath))
        self.df = pd.DataFrame()

    @abstractmethod
    def load(self):
        pass

    def load_and_output_df(self):
        self.load()
        return self.df

    def write_df(self, path, mode='w', **kwargs):
        if self.header:
            with open(path, mode=mode) as f:
                for h in self.header:
                    f.write(h + os.linesep)
        self.df.to_csv(path, mode=('a' if self.header else 'w'), **kwargs)


class VcfDataFrame(BaseBioDataFrame):
    def __init__(self, path):
        super().__init__(
            path=path, supported_exts=['.vcf', '.vcf.gz', '.vcf.bz2']
        )
        self.__fixed_cols = [
            '#CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO',
            'FORMAT'
        ]
        self.__fixed_col_dtypes = {
            '#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
            'QUAL': str, 'FILTER': str, 'INFO': str
        }
        self.__detected_cols = []
        self.__detected_col_dtypes = {}
        self.header = []
        self.samples = []

    def load(self):
        if self.path.endswith('.gz'):
            f = gzip.open(self.path, 'rt')
        elif self.path.endswith('.bz2'):
            f = bz2.open(self.path, 'rt')
        else:
            f = open(self.path, 'r')
        for s in f:
            self._load_vcf_line(string=s)
        f.close()
        self.df = self.df.reset_index(drop=True)

    def _load_vcf_line(self, string):
        if string.startswith('##'):
            self.header.append(string.strip())
        elif string.startswith('#CHROM'):
            items = string.strip().split('\t')
            if items[:len(self.__fixed_cols)] == self.__fixed_cols:
                self.samples = [s for s in items if s not in self.__fixed_cols]
                n_fixed_cols = len(self.__fixed_cols)
                n_detected_cols = len(items)
                self.__detected_cols = self.__fixed_cols + (
                    [
                        'SAMPLE{}'.format(i)
                        for i in range(n_detected_cols - n_fixed_cols)
                    ] if n_detected_cols > n_fixed_cols else []
                )
                self.__detected_col_dtypes = {
                    k: (self.__fixed_col_dtypes.get(k) or str)
                    for k in self.__detected_cols
                }
            else:
                raise ValueError('invalid VCF columns')
        else:
            self.df = self.df.append(
                pd.read_table(
                    io.StringIO(string), header=None,
                    names=self.__detected_cols,
                    dtype=self.__detected_col_dtypes
                )
            )
