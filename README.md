pdvcf
=====

Pandas-based Data Handler for VCF Files

[![wercker status](https://app.wercker.com/status/fe5472a8f890edd3169e7bae5b648bac/s/master "wercker status")](https://app.wercker.com/project/byKey/fe5472a8f890edd3169e7bae5b648bac)

Installation
------------

```sh
$ pip install -U https://github.com/dceoy/pdvcf/archive/master.tar.gz
```

Python API
----------

Example of API call

```py
from pprint import pprint
from pdvcf.dataframe import VcfDataFrame

vcf_path = '/path/to/input/vcf'
vcfdf = VcfDataFrame(path=vcf_path)
vcfdf.load()

pprint(vcfdf.samples)
pprint(vcfdf.header)
print(vcfdf.df)
```

Command-line interface
----------------------

```sh
$ pdvcf --help
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
```
