#!/usr/bin/env python
"""
Pandas-based Data Handler for VCF, BED, and SAM Files.

Usage:
    pdbio vcf2csv [--debug] [--tsv] [--header=<txt>] <vcf> <csv>
    pdbio bed2csv [--debug] [--tsv] [--header=<txt>] <bed> <csv>
    pdbio sam2csv [--debug] [--tsv] [--header=<txt>] <sam> <csv>
    pdbio --version
    pdbio -h|--help

Options:
    --debug, --info     Execute a command with debug|info messages
    --tsv               Use tab instead of comma for a field delimiter
    --header=<txt>      Write a VCF header into a text file
    --version           Print version and exit
    -h, --help          Print help and exit

Commands:
    vcf2csv             Convert a VCF/BCF file to a CSV file
    bed2csv             Convert a BED file to a CSV file
    sam2csv             Convert a SAM/BAM/CRAM file to a CSV file

Arguments:
    <vcf>               Path to a VCF/BCF file
                        (BCF files require `bcftools` command)
    <bed>               Path to a BED file
    <sam>               Path to a SAM/BAM/CRAM file
                        (BAM/CRAM files require `samtools` command)
    <csv>               Path to a CSV/TSV file
"""

import logging
import os

from docopt import docopt

from . import __version__
from .biodataframe import BedDataFrame, SamDataFrame, VcfDataFrame


def main():
    args = docopt(__doc__, version='fract {}'.format(__version__))
    _set_log_config(debug=args['--debug'], info=args['--info'])
    logger = logging.getLogger(__name__)
    logger.debug('args:{0}{1}'.format(os.linesep, args))
    if args['vcf2csv']:
        _convert_file_to_csv(
            src_path=args['<vcf>'], csv_dst_path=args['<csv>'],
            header_txt_dst_path=args['--header'],
            sep=('\t' if args['--tsv'] else ',')
        )
    elif args['bed2csv']:
        _convert_file_to_csv(
            src_path=args['<bed>'], csv_dst_path=args['<csv>'],
            header_txt_dst_path=args['--header'],
            sep=('\t' if args['--tsv'] else ',')
        )
    elif args['sam2csv']:
        _convert_file_to_csv(
            src_path=args['<sam>'], csv_dst_path=args['<csv>'],
            header_txt_dst_path=args['--header'],
            sep=('\t' if args['--tsv'] else ',')
        )


def _convert_file_to_csv(src_path, csv_dst_path, header_txt_dst_path=None,
                         sep=',', file_format='vcf'):
    if file_format == 'vcf':
        biodf = VcfDataFrame(path=src_path)
        df = biodf.df.pipe(lambda d: d.rename(biodf.sample_dict))
    elif file_format == 'bed':
        biodf = BedDataFrame(path=src_path)
        df = biodf.df
    elif file_format == 'sam':
        biodf = SamDataFrame(path=src_path)
        df = biodf.df
    else:
        raise ValueError('invalid file format: {}'.format(file_format))
    df.to_csv(_fetch_abspath(csv_dst_path), sep=sep, index=False)
    if header_txt_dst_path and biodf.header:
        with open(_fetch_abspath(header_txt_dst_path), 'w') as f:
            for h in biodf.header:
                f.write(h + os.linesep)


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


def _fetch_abspath(path):
    return os.path.abspath(os.path.expanduser(os.path.expandvars(path)))
