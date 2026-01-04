# Make input file for NCD2 from VCF

Make input file for NCD2 from VCF

## Usage

``` r
.vcf_ncd2(x, nind = nind, index.col = index.col, verbose = T)
```

## Arguments

- x:

  A data.table object generated from reading a VCF file with read_vcf.

- nind:

  A vector containing the number of diploid individuals from each
  population to be used. The length of the vector should be equal to the
  number of populations being sampled.

- index.col:

  First genotype column in VCF file.

- verbose:

  Logical. Printed updates. Default is T.
