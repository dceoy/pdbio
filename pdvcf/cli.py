#!/usr/bin/env python
"""
Pandas-based Data Handler for VCF Files

Usage:
    pdvcf csv [--debug] [--header=<txt>] <vcf> <csv>
    pdvcf --version
    pdvcf -h|--help

Options:
    --debug, --info     Execute a command with debug|info messages
    --header=<txt>      Write a VCF header into a text file
    --tsv               Use tab instead of comma for a field delimiter
    --version           Print version and exit
    -h, --help          Print help and exit

Commands:
    csv                 Convert a VCF file to a CSV file

Arguments:
    <vcf>               Path to a VCF file
    <csv>               Path to a CSV/TSV file
"""

import logging
import os

from docopt import docopt

from . import __version__
from .dataframe import VcfDataFrame


def main():
    args = docopt(__doc__, version='fract {}'.format(__version__))
    _set_log_config(debug=args['--debug'], info=args['--info'])
    logger = logging.getLogger(__name__)
    logger.debug('args:{0}{1}'.format(os.linesep, args))
    if args['tsv']:
        _convert_vcf_to_csv(
            vcf_path=args['<vcf>'], csv_path=args['<csv>'],
            header_txt_path=args['--header'],
            sep=('\t' if args['--tsv'] else ',')
        )


def _set_log_config(debug=None, info=None):
    if debug:
        lv = logging.DEBUG
    elif info:
        lv = logging.INFO
    else:
        lv = logging.WARNING
    logging.basicConfig(
        format='%(asctime)s %(levelname)-8s %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S', level=lv
    )


def _convert_vcf_to_csv(vcf_path, csv_path, header_txt_path=None, sep=','):
    vcfdf = VcfDataFrame(path=_fetch_abspath(vcf_path))
    vcfdf.load()
    vcfdf.df.rename(
        {'SAMPLE{}'.format(i): n for i, n in enumerate(vcfdf.df.samples)}
    ).to_csv(
        _fetch_abspath(csv_path), sep=sep, index=False
    )
    if header_txt_path:
        with open(_fetch_abspath(header_txt_path), 'r') as f:
            for h in vcfdf.header:
                f.write(h + os.linesep)


def _fetch_abspath(path):
    return os.path.abspath(os.path.expanduser(os.path.expandvars(path)))
