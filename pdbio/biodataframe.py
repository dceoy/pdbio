#!/usr/bin/env python
"""
Pandas-based Table Data Handler
https://github.com/dceoy/pdbio
"""

import bz2
import gzip
import logging
import os
import subprocess
import sys
from abc import ABCMeta, abstractmethod
from collections import OrderedDict
from itertools import product

import pandas as pd


class BaseBioDataFrame(object, metaclass=ABCMeta):
    """Base DataFrame handler for Table Files."""

    def __init__(self, path=None, format_name='TSV', delimiter='\t',
                 column_header=True, chrom_column=None, pos_columns=None,
                 txt_file_exts=None, bin_file_exts=None, load=True):
        for a in [pos_columns, txt_file_exts, bin_file_exts]:
            assert type(a) is not str
        self.__logger = logging.getLogger(__name__)
        self.__format_name = format_name
        self.__column_header = column_header
        self.__delimiter = delimiter
        self.__chrom_column = chrom_column
        self.__pos_columns = pos_columns
        self.__file_exts = [
            *(
                [e + c for e, c in product(txt_file_exts, ['', '.gz', '.bz2'])]
                if txt_file_exts else list()
            ),
            *(bin_file_exts or list())
        ]
        self.path = path
        self.header = list()
        self.df = pd.DataFrame()
        if path and load:
            self.load_table(path=path)

    def load_table(self, path):
        self._update_path(path=path)
        self.__logger.info(
            'Load {0} file: {1}'.format(self.__format_name, self.path)
        )
        self.header = list()
        self.df = pd.DataFrame()
        self.load()
        self.__logger.debug('self.df shape: {}'.format(self.df.shape))
        self.__logger.debug('self.df:{0}{1}'.format(os.linesep, self.df))
        return self

    def _update_path(self, path):
        abspath = self.normalize_path(path=path)
        if not os.path.isfile(abspath):
            raise FileNotFoundError('file not found: {}'.format(abspath))
        elif (self.__file_exts and
              not [x for x in self.__file_exts if abspath.endswith(x)]):
            raise ValueError('invalid file extension: {}'.format(abspath))
        elif self.path != abspath:
            self.__logger.debug('abspath: {}'.format(abspath))
            self.path = abspath

    @staticmethod
    def normalize_path(path):
        return os.path.abspath(os.path.expanduser(os.path.expandvars(path)))

    @abstractmethod
    def load(self):
        self.df = pd.read_csv(
            self.path, header=self.__column_header, sep=self.__delimiter
        )

    @staticmethod
    def open_readable_file(path):
        if path.endswith('.gz'):
            return gzip.open(path, mode='rt')
        elif path.endswith('.bz2'):
            return bz2.open(path, mode='rt')
        else:
            return open(path, mode='r')

    def output_table(self, path=None, **kwargs):
        if path:
            abspath = self.normalize_path(path=path)
            self.__logger.info(
                'Write {0} file: {1}'.format(self.__format_name, abspath)
            )
        else:
            abspath = None
            self.__logger.info('Print {}'.format(self.__format_name))
        if self.header:
            self.output_header(path=abspath)
        self.df.to_csv(
            (abspath or sys.stdout), mode=('a' if self.header else 'w'),
            index=False, header=self.__column_header, sep=self.__delimiter,
            **kwargs
        )

    def output_header(self, path=None):
        if path:
            with open(path, mode='w') as f:
                for h in self.header:
                    f.write(h + os.linesep)
        else:
            for h in self.header:
                print(h, flush=True)

    def sort(self, **kwargs):
        self.df = self._sort_df(**kwargs)
        return self

    def sorted_df(self, **kwargs):
        return self._sort_df(**kwargs)

    def _sort_df(self, df=None, by_chrom=True, **kwargs):
        self.__logger.info('Sort the dataframe.')
        return (self.df if df is None else df).pipe(
            lambda d: (
                self._sort_by_chrom_and_pos(df=d, **kwargs)
                if self.__chrom_column and by_chrom else
                d.sort_values(by=d.columns, **kwargs)
            )
        )

    def _sort_by_chrom_and_pos(self, df, **kwargs):
        ci = self.__chrom_column + '_sort_index'
        pis = self.__pos_columns or list()
        return df.assign(**{
            ci: (lambda d: d[self.__chrom_column].apply(self._chrom2int))
        }).sort_values(
            by=list(
                OrderedDict.fromkeys([
                    ci, self.__chrom_column, *pis, *df.columns
                ])
            ), **kwargs
        ).drop(columns=ci)

    @staticmethod
    def _chrom2int(chrom):
        id = chrom[3:] if chrom.lower().startswith('chr') else chrom
        return (
            int(id) if id.isdigit() else
            (({'X': 0, 'Y': 1, 'M': 2, 'MT': 2}.get(id.upper()) or 3) + 1000)
        )

    def fetch_executable(self, cmd):
        executables = [
            x for x in
            [os.path.join(p, cmd) for p in str.split(os.environ['PATH'], ':')]
            if os.access(x, os.X_OK)
        ]
        if executables:
            self.__logger.debug('path to {0}: {1}'.format(cmd, executables[0]))
            return executables[0]
        else:
            raise RuntimeError('command not found:: {}'.format(cmd))

    def run_and_parse_subprocess(self, args, stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE, **kwargs):
        with subprocess.Popen(args=args, stdout=stdout, stderr=stderr,
                              **kwargs) as p:
            for line in p.stdout:
                yield line.decode('utf-8')
            outs, errs = p.communicate()
            if p.returncode == 0:
                pass
            else:
                self.__logger.error(
                    'STDERR from subprocess `{0}`:{1}{2}'.format(
                        p.args, os.linesep, errs.decode('utf-8')
                    )
                )
                raise subprocess.CalledProcessError(
                    returncode=p.returncode, cmd=p.args, output=outs,
                    stderr=errs
                )
