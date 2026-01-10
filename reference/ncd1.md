# Calculate non-central statistic (NCD1)

Calculate non-central statistic (NCD1)

## Usage

``` r
ncd1(x = x, tf = 0.5, w = 3000, ncores = 2, minIS = 2, by = "POS")
```

## Arguments

- x:

  A data.table object with columns: CHR, POS, REF, ALT, tx_1 (number of
  alternate allele copies), tn_1 (total number of alleles).

- tf:

  Target frequency.

- w:

  Window size in bp. Default is 3000.

- ncores:

  Number of cores. Increasing this can speed things up for you. Default
  is 2.

- minIS:

  Minimum number of informative sites. Default is 2. Windows with less
  informative sites than this threshold are discarded.

- by:

  Define how to scan the genome. "POS" (default) defined sliding windows
  based on w. Future implementation: "IS" defined windows around each
  informative site.

## Value

A data.table object with columns: Win.ID, S (sites), IS (informative
sites), tf (target frequency), ncd1

## References

Bitarello et al. Genome Biology and Evolution, Volume 10, Issue 3, March
2018, Pages 939–955, https://doi.org/10.1093/gbe/evy054

## Examples

``` r
ncd1(x=ncd1_input, tf=0.5, w=3000, ncores=2, minIS=2)
#>            Win.ID     S    IS    tf      ncd1
#>            <char> <int> <int> <num>     <num>
#>  1:     1_92_3092    20    20   0.5 0.4286047
#>  2:   1_1592_4592    21    21   0.5 0.4392331
#>  3:   1_3092_6092    15    15   0.5 0.4046166
#>  4:   1_4592_7592    12    12   0.5 0.3897881
#>  5:   1_6092_9092    19    19   0.5 0.3793086
#>  6:  1_7592_10592    28    28   0.5 0.3641235
#>  7:  1_9092_12092    25    25   0.5 0.3797493
#>  8: 1_10592_13592    15    15   0.5 0.4063328
#>  9: 1_12092_15092    25    25   0.5 0.3935370
#> 10: 1_13592_16592    30    30   0.5 0.3914497
#> 11: 1_15092_18092    22    22   0.5 0.4010361
#> 12: 1_16592_19592    18    18   0.5 0.4038215
#> 13: 1_18092_21092    11    11   0.5 0.4385437
#> 14: 1_19592_22592    15    15   0.5 0.4623715
#> 15: 1_21092_24092    20    20   0.5 0.4395955
#> 16: 1_22592_25592    16    16   0.5 0.4394827
#> 17: 1_24092_27092    17    17   0.5 0.4659706
#> 18: 1_25592_28592    20    20   0.5 0.4454260
#> 19: 1_27092_30092    22    22   0.5 0.4174539
#> 20: 1_28592_31592    12    12   0.5 0.4136724
#>            Win.ID     S    IS    tf      ncd1
#>            <char> <int> <int> <num>     <num>
```
