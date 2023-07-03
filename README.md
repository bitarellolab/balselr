
<!-- README.md is generated from README.Rmd. Please edit that file -->

# balselr

This package is under construction. If you want to run NCD1 or NCD2,
please use this [other
repo](https://github.com/bbitarello/NCD-Statistics) for now.

<!-- badges: start -->

<img src="images/balselr2.png" style="width:40.0%;height:40.0%" />
<!-- badges: end -->

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

This is a basic example which shows you how to read in a vcf file and
keep all samples/individuals

``` r
read_vcf(x="inst/example.vcf", only.bi=T)
```

This is an example which shows how to parse a vcf file and output an
input file for ncd1 or ncd2:

``` r
parse_vcf(infile="inst/example.vcf",n0=108, n1=1, type="ncd2", outfile="example_parse_ncd2.out")
```

This is an example which shows how to parse a vcf file and output an
input file for ncd1 or ncd1:

``` r
#upcoming
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this. You could also
use GitHub Actions to re-render `README.Rmd` every time you push. An
example workflow can be found here:
<https://github.com/r-lib/actions/tree/v1/examples>.
