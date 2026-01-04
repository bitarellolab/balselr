# Input data to run ncd1 example

An R object in data.table format

## Usage

``` r
ncd1_input
```

## Format

A data.table with 838 rows and 6 columns

- CHR:

  Chromsome, single value - 1

- POS:

  Position along the chromosome - first is 92, last is 29992

- REF:

  Reference allele - A, C, G, T

- ALT:

  Alternate allele - A, C, G, T

- tx_1:

  Number of alternate allele copies - min=0, max=216

- tn_1:

  Total number of alleles, single value - 216

## Source

created from a VCF (Variant Call Format) input file generated with the
SLiM simulator by reading it in with parse_vcf(type="ncd1")
