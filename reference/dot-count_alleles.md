# Count derived or ancestral alleles for a given vcf position

Count derived or ancestral alleles for a given vcf position

## Usage

``` r
.count_alleles(x, split = "|")
```

## Arguments

- x:

  A vector containing 0s and 1s corresponding to allele states. 0/0:
  homozygous reference; 0/1: heterozygous; 1/1: homozygous alternate

- split:

  Pattern to look for. Inherited from stringr::str_split
