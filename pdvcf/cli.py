#!/usr/bin/env python
"""
Pandas-based Data Handler for VCF Files

Usage:
    pdvcf tsv [--debug] [--header=<tsv>] [--variant=<tsv>] [--expand] [--quiet]
              <vcf>...
    pdvcf --version
    pdvcf -h|--help

Options:
    --debug, --info     Execute a command with debug|info messages
    --header=<tsv>      Set an output path for headers
    --variant=<tsv>     Set an output path for variants
    --expand            Parse and expand INFO and SAMPLE columns
    --quiet             Suppress record printing
    --version           Print version and exit
    -h, --help          Print help and exit

Commands:
    tsv                 Write records in VCF files into a single TSV file

Arguments:
    <vcf>               Path to a VCF file for 1 sample
"""

import logging
import os
from docopt import docopt
from .util import set_log_config
from . import __version__


def main():
    args = docopt(__doc__, version='fract {}'.format(__version__))
    set_log_config(debug=args['--debug'], info=args['--info'])
    logger = logging.getLogger(__name__)
    logger.debug('args:{0}{1}'.format(os.linesep, args))
    if args['tsv']:
        _write_variant_tsv()


def _write_variant_tsv(vcf_paths, variant_tsv):
    pass
