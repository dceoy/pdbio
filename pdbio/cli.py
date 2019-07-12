#!/usr/bin/env python
"""
Pandas-based Data Handler for VCF, BED, and SAM Files.

Usage:
    pdbio vcf2csv [--debug|--info] [--sort] [--tsv] [--header=<txt>] <src>
                  <dst>
    pdbio bed2csv [--debug|--info] [--sort] [--tsv] [--header=<txt>] <src>
                  <dst>
    pdbio sam2csv [--debug|--info] [--sort] [--tsv] [--header=<txt>] <src>
                  <dst>
    pdbio --version
    pdbio -h|--help

Options:
    --debug, --info     Execute a command with debug|info messages
    --sort              Sort a dataframe
    --tsv               Use tab instead of comma for a field delimiter
    --header=<txt>      Write a VCF header into a text file
    --version           Print version and exit
    -h, --help          Print help and exit

Commands:
    vcf2csv             Convert a VCF/BCF file to a CSV file
                        (BCF files require `bcftools` command)
    bed2csv             Convert a BED file to a CSV file
    sam2csv             Convert a SAM/BAM/CRAM file to a CSV file
                        (BAM/CRAM files require `samtools` command)

Arguments:
    <src>               Path to an input file
    <dst>               Path to an output file
"""

import logging
import os

from docopt import docopt

from . import __version__
from .beddataframe import BedDataFrame
from .samdataframe import SamDataFrame
from .vcfdataframe import VcfDataFrame


def main():
    args = docopt(__doc__, version='fract {}'.format(__version__))
    _set_log_config(debug=args['--debug'], info=args['--info'])
    logger = logging.getLogger(__name__)
    logger.debug('args:{0}{1}'.format(os.linesep, args))
    csv_convert = [k for k in ['vcf2csv', 'bed2csv', 'sam2csv'] if args[k]]
    if csv_convert:
        _convert_file_to_csv(
            src_path=args['<src>'], csv_dst_path=args['<dst>'],
            sort=args['--sort'], sep=('\t' if args['--tsv'] else ','),
            header_txt_dst_path=args['--header'],
            file_format=csv_convert[0].split('2')[0]
        )


def _convert_file_to_csv(src_path, csv_dst_path, sort=False, sep=',',
                         header_txt_dst_path=None, file_format='vcf'):
    if file_format == 'vcf':
        biodf = VcfDataFrame(path=src_path)
    elif file_format == 'bed':
        biodf = BedDataFrame(path=src_path)
    elif file_format == 'sam':
        biodf = SamDataFrame(path=src_path)
    else:
        raise ValueError('invalid file format: {}'.format(file_format))
    if sort:
        biodf.sort()
    biodf.df.to_csv(biodf.normalize_path(csv_dst_path), sep=sep, index=False)
    if header_txt_dst_path and biodf.header:
        biodf.write_header(path=biodf.normalize_path(header_txt_dst_path))


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
