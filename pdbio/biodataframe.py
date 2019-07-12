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
from abc import ABCMeta, abstractmethod

import pandas as pd


class BaseBioDataFrame(object, metaclass=ABCMeta):
    """Base DataFrame handler for Table Files."""

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
        self.__logger = logging.getLogger(__name__)
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

    def sort_df_by_chrom(self, chrom_col, add_cols=None):
        self.df = self.df.assign(
            chrom_int=lambda d: d[chrom_col].apply(self._chrom2int)
        ).sort_values(
            ['chrom_int', *(add_cols if add_cols else list())]
        ).drop(columns='chrom_int')

    @staticmethod
    def _chrom2int(chrom):
        id = chrom[3:] if chrom.lower().startswith('chr') else chrom
        return (
            int(id) if id.isdigit() else
            (({'X': 0, 'Y': 1, 'M': 2, 'MT': 2}.get(id.upper()) or 3) + 1000)
        )
