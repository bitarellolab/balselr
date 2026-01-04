# Input data to run ncd2 example

An R object in data.table format

## Usage

``` r
ncd2_input
```

## Format

A data.table with 838 rows and 8 columns

- CHR:

  Chromsome, single value - 1

- POS:

  Position along the chromosome - first is 92, last is 29992

- REF:

  Reference allele - A, C, G, T

- ALT:

  Alternate allele - A, C, G, T

- tx_1:

  Number of alternate allele copies in the ingroup - min=0, max=216

- tn_1:

  Total number of alleles in the ingroup, single value - 216

- tx_2:

  Number of alternate allele copies in the outgroup - min=0, max=4

- tn_2:

  Total number of alleles in the outgroup, single value - 216

## Source

created from a VCF (Variant Call Format) input file generated with the
SLiM simulator by reading it in with parse_vcf(type="ncd2")
