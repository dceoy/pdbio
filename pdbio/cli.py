#!/usr/bin/env python
"""
Pandas-based Data Handler for VCF, BED, and SAM Files.

Usage:
    pdbio vcf2csv [--debug|--info] [--sort] [--tsv] [--expand-info]
                  [--expand-samples] [--header=<path>] <src> [<dst>]
    pdbio bed2csv [--debug|--info] [--sort] [--tsv] [--header=<path>]
                  <src> [<dst>]
    pdbio sam2csv [--debug|--info] [--sort] [--tsv] [--header=<path>]
                  <src> [<dst>]
    pdbio vcfsort [--debug|--info] <src> [<dst>]
    pdbio bedsort [--debug|--info] <src> [<dst>]
    pdbio samsort [--debug|--info] <src> [<dst>]
    pdbio idvars [--debug|--info] <fa> <chrom> <variant>...
    pdbio --version
    pdbio -h|--help

Options:
    --debug, --info     Execute a command with debug|info messages
    --sort              Sort a dataframe
    --tsv               Use tab instead of comma for a field delimiter
    --expand-info       Expand the INFO column in a VCF file
    --expand-samples    Expand columns of samples in a VCF file
    --header=<path>     Write a header into a text file
    --version           Print version and exit
    -h, --help          Print help and exit

Commands:
    vcf2csv             Convert a VCF/BCF file to a CSV file
    bed2csv             Convert a BED file to a CSV file
    sam2csv             Convert a SAM/BAM/CRAM file to a CSV file
    vcfsort             Sort a VCF/BCF file
    bedsort             Sort a BED file
    samsort             Sort a SAM/BAM/CRAM file
    idvars              Identify if variants are identical

Arguments:
    <src>               Path to an input file
    <dst>               Path to an output file
    <fa>                Path to a reference genome FASTA file
    <chrom>             Chromosome name (e.g., chr8)
    <variant>           Variant expression (e.g., 32529588:C_CGG)
"""

import logging
import os
import signal
import sys

from docopt import docopt

from . import __version__
from .beddataframe import BedDataFrame
from .identifier import identify_variants
from .samdataframe import SamDataFrame
from .vcfdataframe import VcfDataFrame


def main():
    args = docopt(__doc__, version='pdbio {}'.format(__version__))
    _set_log_config(debug=args['--debug'], info=args['--info'])
    logger = logging.getLogger(__name__)
    logger.debug('args:{0}{1}'.format(os.linesep, args))
    csv_convert = [k for k in ['vcf', 'bed', 'sam'] if args[k + '2csv']]
    chrom_sort = [k for k in ['vcf', 'bed', 'sam'] if args[k + 'sort']]
    if args['idvars']:
        identify_variants(
            fa_path=args['<fa>'], chrom=args['<chrom>'],
            variants=args['<variant>']
        )
    elif csv_convert:
        _convert_file_to_csv(
            src_path=args['<src>'], dst_path=args['<dst>'],
            sort=args['--sort'], sep=('\t' if args['--tsv'] else ','),
            header_dst_path=args['--header'], file_format=csv_convert[0],
            expand_info=args['--expand-info'],
            expand_samples=args['--expand-samples']
        )
    elif chrom_sort:
        _sort_by_chrom(
            src_path=args['<src>'], dst_path=args['<dst>'],
            file_format=chrom_sort[0]
        )


def _sort_by_chrom(src_path, dst_path=None, file_format='vcf'):
    signal.signal(signal.SIGINT, signal.SIG_DFL)
    if file_format == 'vcf':
        biodf = VcfDataFrame(path=src_path)
    elif file_format == 'bed':
        biodf = BedDataFrame(path=src_path)
    elif file_format == 'sam':
        biodf = SamDataFrame(path=src_path)
    else:
        raise ValueError('invalid file format: {}'.format(file_format))
    biodf.sort().write_table(path=dst_path)


def _convert_file_to_csv(src_path, dst_path=None, sort=False, sep=',',
                         header_dst_path=None, file_format='vcf',
                         expand_info=False, expand_samples=False):
    signal.signal(signal.SIGINT, signal.SIG_DFL)
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
    df = biodf.df
    if file_format == 'vcf':
        df = biodf.expanded_df(
            df=df, by_info=expand_info, by_samples=expand_samples
        )
    df.to_csv(
        (biodf.normalize_path(dst_path) if dst_path else sys.stdout),
        sep=sep, index=False
    )
    if header_dst_path and biodf.header:
        biodf.write_header(path=header_dst_path)


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
