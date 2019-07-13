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

vcf_path = 'test/example.vcf'
vcfdf = VcfDataFrame(path=vcf_path)

pprint(vcfdf.header)      # list of header
pprint(vcfdf.samples)     # list of samples
print(vcfdf.df)           # VCF dataframe

vcffdf.sort()             # sort by CHROM, POS, and the other
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
