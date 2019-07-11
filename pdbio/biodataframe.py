#!/usr/bin/env python
"""
Pandas-based Data Handler for VCF, BED, and SAM Files.
https://github.com/dceoy/pdbio
"""

import bz2
import gzip
import io
import logging
import os
import re
import subprocess
from abc import ABCMeta, abstractmethod
from collections import OrderedDict
from itertools import product
from multiprocessing import cpu_count

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
                    ('SAMPLE{}'.format(i), n) for i, n in enumerate(samples)
                ])
                n_fixed_cols = len(self.__fixed_cols)
                n_detected_cols = len(items)
                self.__detected_cols = [
                    *self.__fixed_cols,
                    *[
                        'SAMPLE{}'.format(i)
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
        self.sort_df_by_chrom(chrom_col='#CHROM', add_cols=self.__fixed_cols)

    def write_vcf(self, path):
        self.__logger.info('Write a VCF file: {}'.format(path))
        self.write_header(path=path)
        self.df.rename(self.sample_dict).to_csv(
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

    def sort_df(self):
        self.sort_df_by_chrom(chrom_col='chrom', add_cols=self.__fixed_cols)

    def write_bed(self, path):
        self.__logger.info('Write a BED file: {}'.format(path))
        self.write_header(path=path)
        self.df.to_csv(
            path, header=False, mode=('a' if self.header else 'w'), sep='\t',
            index=False
        )


class SamDataFrame(BaseBioDataFrame):
    """SAM DataFrame handler."""

    def __init__(self, path, return_df=False, samtools=None, n_thread=None):
        super().__init__(
            path=path,
            supported_exts=[
                *[
                    (e + c) for e, c in
                    product(['.sam', '.txt', '.tsv'], ['', '.gz', '.bz2'])
                ], '.bam', '.cram'
            ]
        )
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
        if return_df:
            self.load_and_output_df()
        else:
            self.load()

    def load(self):
        self.__logger.info('Load a SAM file: {}'.format(self.path))
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
        self.df = self.df.reset_index(drop=True)

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
                )
            )

    def sort_df(self):
        self.sort_df_by_chrom(chrom_col='RNAME', add_cols=self.__fixed_cols)

    def write_sam(self, path):
        self.__logger.info('Write a SAM file: {}'.format(path))
        self.write_header(path=path)
        self.df.to_csv(
            path, header=False, mode=('a' if self.header else 'w'), sep='\t',
            index=False
        )
