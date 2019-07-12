pdbio
=====

Pandas-based Data Handler for VCF, BED, and SAM Files

[![wercker status](https://app.wercker.com/status/fe5472a8f890edd3169e7bae5b648bac/s/master "wercker status")](https://app.wercker.com/project/byKey/fe5472a8f890edd3169e7bae5b648bac)

Installation
------------

```sh
$ pip install -U https://github.com/dceoy/pdbio/archive/master.tar.gz
```

Python API
----------

Example of API call

```py
from pprint import pprint
from pdbio.biodataframe import VcfDataFrame

vcf_path = '/path/to/input/vcf'
vcfdf = VcfDataFrame(path=vcf_path)

pprint(vcfdf.header)
pprint(vcfdf.samples)
print(vcfdf.df)
```

Command-line interface
----------------------

```sh
$ pdbio --help
Pandas-based Data Handler for VCF, BED, and SAM Files.

Usage:
    pdbio vcf2csv [--debug|--info] [--sort] [--tsv] [--header=<path>]
                  [--dst=<path>] <src>
    pdbio bed2csv [--debug|--info] [--sort] [--tsv] [--header=<path>]
                  [--dst=<path>] <src>
    pdbio sam2csv [--debug|--info] [--sort] [--tsv] [--header=<path>]
                  [--dst=<path>] <src>
    pdbio vcfsort [--debug|--info] [--dst=<path>] <src>
    pdbio bedsort [--debug|--info] [--dst=<path>] <src>
    pdbio samsort [--debug|--info] [--dst=<path>] <src>
    pdbio --version
    pdbio -h|--help

Options:
    --debug, --info     Execute a command with debug|info messages
    --sort              Sort a dataframe
    --tsv               Use tab instead of comma for a field delimiter
    --header=<path>     Write a header into a text file
    --dst=<path>        Write results into a text file
    --version           Print version and exit
    -h, --help          Print help and exit

Commands:
    vcf2csv             Convert a VCF/BCF file to a CSV file
    bed2csv             Convert a BED file to a CSV file
    sam2csv             Convert a SAM/BAM/CRAM file to a CSV file
    vcfsort             Sort a VCF/BCF file
    bedsort             Sort a BED file
    samsort             Sort a SAM/BAM/CRAM file

Arguments:
    <src>               Path to an input file
```
