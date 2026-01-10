# Calculate non-central statistic (NCD2)

Calculate non-central statistic (NCD2)

## Usage

``` r
ncd2(x = x, tf = 0.5, w = 3000, ncores = 2, minIS = 2, by = "POS")
```

## Arguments

- x:

  A data.table object with columns: CHR, POS, REF, ALT, tx_1 (number of
  alternate allele copies), tn_1 (total number of alleles), tx_2, and
  tn_2.

- tf:

  Target frequency. Any value between 0 and 0.5. Default is 0.5.

- w:

  Window size in bp. Default is 1000.

- ncores:

  Number of cores. Increasing this can speed things up for you. Default
  is 2.

- minIS:

  Minimum number of informative sites. Default is 2. Windows with less
  informative sites than this threshold are discarded.

- by:

  Define how to scan the genome. "POS" (default) defined sliding windows
  based on w. "IS" defined windows around each informative site.

## Value

A data.table object with columns: Win.ID, S (sites), FD (fixed
differences),IS (informative sites), tf (target frequency), ncd2

## References

Bitarello et al. Genome Biology and Evolution, Volume 10, Issue 3, March
2018, Pages 939–955, https://doi.org/10.1093/gbe/evy054

## Examples

``` r
ncd2(x=ncd2_input, tf=0.5, w=3000, ncores=2, minIS=2)
#>            Win.ID      ncd2     S    FD    IS    tf
#>            <char>     <num> <int> <int> <int> <num>
#>  1:     1_92_3092 0.4721648    20    29    49   0.5
#>  2:   1_1592_4592 0.4749101    21    28    49   0.5
#>  3:   1_3092_6092 0.4637019    15    22    37   0.5
#>  4:   1_4592_7592 0.4671135    12    25    37   0.5
#>  5:   1_6092_9092 0.4595906    19    33    52   0.5
#>  6:  1_7592_10592 0.4361486    28    27    55   0.5
#>  7:  1_9092_12092 0.4462502    25    27    52   0.5
#>  8: 1_10592_13592 0.4743746    15    36    51   0.5
#>  9: 1_12092_15092 0.4631460    25    42    67   0.5
#> 10: 1_13592_16592 0.4507518    30    32    62   0.5
#> 11: 1_15092_18092 0.4553608    22    24    46   0.5
#> 12: 1_16592_19592 0.4612429    18    24    42   0.5
#> 13: 1_18092_21092 0.4797632    11    21    32   0.5
#> 14: 1_19592_22592 0.4832595    15    18    33   0.5
#> 15: 1_21092_24092 0.4746826    20    26    46   0.5
#> 16: 1_22592_25592 0.4825311    16    37    53   0.5
#> 17: 1_24092_27092 0.4889201    17    34    51   0.5
#> 18: 1_25592_28592 0.4812203    20    36    56   0.5
#> 19: 1_27092_30092 0.4709149    22    37    59   0.5
#> 20: 1_28592_31592 0.4662210    12    17    29   0.5
#>            Win.ID      ncd2     S    FD    IS    tf
#>            <char>     <num> <int> <int> <int> <num>
```
