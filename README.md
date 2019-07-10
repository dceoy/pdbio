pdbio
=====

Pandas-based Data Handler for VCF and BED Files

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
pprint(vcfdf.sample_dict)
print(vcfdf.df)
```

Command-line interface
----------------------

```sh
$ pdbio --help
Pandas-based Data Handler for VCF and BED Files

Usage:
    pdbio vcf2csv [--debug] [--tsv] [--header=<txt>] <vcf> <csv>
    pdbio bed2csv [--debug] [--tsv] [--header=<txt>] <bed> <csv>
    pdbio --version
    pdbio -h|--help

Options:
    --debug, --info     Execute a command with debug|info messages
    --tsv               Use tab instead of comma for a field delimiter
    --header=<txt>      Write a VCF header into a text file
    --version           Print version and exit
    -h, --help          Print help and exit

Commands:
    vcf2csv             Convert a VCF file to a CSV file
    bed2csv             Convert a BED file to a CSV file

Arguments:
    <vcf>               Path to a VCF file
    <bed>               Path to a BED file
    <csv>               Path to a CSV/TSV file
```
