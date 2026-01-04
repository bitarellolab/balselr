# Parse vcf

Parse vcf

## Usage

``` r
parse_vcf(
  infile = NULL,
  vcf_data = NULL,
  n0 = NULL,
  n1 = NULL,
  type = NULL,
  verbose = F
)
```

## Arguments

- infile:

  The path and name of a vcf file - single chromosome. Can be NULL if
  vcf_data is provided.

- vcf_data:

  A data.table object containing VCF data. Can be NULL if infile is
  provided.

- n0:

  Number of individuals in pop0. Must be \>=2 for both ncd1 and ncd2.

- n1:

  Number of individuals in pop1. Must be \>=1 for ncd2. Must be NULL for
  ncd1. NULL by default.

- type:

  Which input format will be generated. Options: ncd1, ncd2.

- verbose:

  Logical. Printed updates. Default is F.

## Value

Returns a data table object with relevant columns based on type

## Examples

``` r
parse_vcf(vcf_data=vcf_obj, n0=108, n1=2,type="ncd2")
#>        CHR   POS    REF    ALT  tx_1  tn_1  tx_2  tn_2
#>      <int> <int> <char> <char> <int> <int> <int> <int>
#>   1:     1    92      A      C   216   216     0     4
#>   2:     1   177      C      A   216   216     0     4
#>   3:     1   283      T      C     0   216     0     4
#>   4:     1   289      C      G     0   216     4     4
#>   5:     1   327      G      T     0   216     2     4
#>  ---                                                  
#> 834:     1 29734      A      T    53   216     0     4
#> 835:     1 29761      T      C     0   216     0     4
#> 836:     1 29795      A      T     1   216     0     4
#> 837:     1 29915      C      A     0   216     0     4
#> 838:     1 29992      G      T    11   216     0     4
parse_vcf(infile=system.file(package = "balselr", "example.vcf"), n0=108, type="ncd1")
#>        CHR   POS    REF    ALT  tx_1  tn_1
#>      <int> <int> <char> <char> <int> <int>
#>   1:     1    92      A      C   216   216
#>   2:     1   177      C      A   216   216
#>   3:     1   283      T      C     0   216
#>   4:     1   289      C      G     0   216
#>   5:     1   327      G      T     0   216
#>  ---                                      
#> 834:     1 29734      A      T    53   216
#> 835:     1 29761      T      C     0   216
#> 836:     1 29795      A      T     1   216
#> 837:     1 29915      C      A     0   216
#> 838:     1 29992      G      T    11   216
```
