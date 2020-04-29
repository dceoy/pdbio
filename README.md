pdbio
=====

Pandas-based Data Handler for VCF, BED, and SAM Files

[![wercker status](https://app.wercker.com/status/fe5472a8f890edd3169e7bae5b648bac/s/master "wercker status")](https://app.wercker.com/project/byKey/fe5472a8f890edd3169e7bae5b648bac)
![Upload Python Package](https://github.com/dceoy/pdbio/workflows/Upload%20Python%20Package/badge.svg)

Installation
------------

```sh
$ pip install -U pdbio
```

Python API
----------

Example of API call

```py
from pprint import pprint
from pdbio.vcfdataframe import VcfDataFrame

vcf_path = 'test/example.vcf'
vcfdf = VcfDataFrame(path=vcf_path)

pprint(vcfdf.header)      # list of header
pprint(vcfdf.samples)     # list of samples
print(vcfdf.df)           # VCF dataframe

vcfdf.sort()              # sort by CHROM, POS, and the other
print(vcfdf.df)           # sorted dataframe
```

Command-line interface
----------------------

Example of commands

```sh
# Convert VCF data into sorted TSV data
$ pdbio vcf2csv --sort --tsv test/example.vcf

# Convert VCF data into expanded CSV data
$ pdbio vcf2csv --expand-info --expand-samples test/example.vcf

# Sort VCF data by CHROM, POS, and the other
$ pdbio vcfsort test/example.vcf
```

Run `pdbio --help` for more information.
