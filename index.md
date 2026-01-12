# balselr

## What is `balselr`?

Balancing selection with R allows you to run **NCD statistics** to
detect **long-term balancing selection** in genomic datasets as
described in the paper: [Bitarello et
al. (2018)](https://academic.oup.com/gbe/article/10/3/939/4938688)

## How to cite us

For now, please cite us as follows:

> Hansell, D and Bitarello, B (2026). balselr: Balancing Selection in R
> Zenodo. <https://doi.org/10.5281/zenodo.18148850>

Thank you!

## Installation

You can install the development version of balselr from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
library(devtools)
devtools::install_github("bitarellolab/balselr")
library(balselr)
```

## Example

This is a basic example which shows you how to read in a vcf file:

``` r
read_vcf(x = system.file(package = "balselr", "example.vcf"))
```

This is an example which shows how to parse a vcf file and output an
input file for `ncd1`:

``` r
parse_vcf(
        infile = system.file(package = "balselr", "example.vcf"),
        n0 = 108,
        type = "ncd1"
)
```

This is an example which shows how to parse a vcf file and output an
input file for `ncd2`:

``` r
parse_vcf(
        infile = system.file(package = "balselr", "example.vcf"),
        n0 = 108,
        n1 = 2,
        type = "ncd2"
)
```

Run `ncd1 (tf=0.5)` with a 3000 basepair window and a minimum of 8
informative sites per window using 2 cores:

``` r
data(ncd1_input)
ncd1(x = ncd1_input, tf = 0.5, w = 3000, ncores = 2, minIS=8)
```

Run `ncd2 (tf=0.5)` with a 3000 basepair window and a minimum of 2
informative sites per window using 2 cores:

``` r
data(ncd2_input)
ncd2(x = ncd2_input, tf = 0.5, w = 3000, ncores = 2, minIS=2)
```
