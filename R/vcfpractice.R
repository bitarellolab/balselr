library("data.table")


read_vcf("example_vcf")


read_vcf("/Users/dhansell/Documents/GitHub//chrom22.vcf")
parsed_vcf22<-parse_vcf("/Users/dhansell/Documents/GitHub//chrom22.vcf", type = "ncd1")
stat22<-ncd1(parsed_vcf22)



read_vcf("/Users/dhansell/Documents/GitHub//chrom13.vcf")
parsed_vcf13<-parse_vcf("/Users/dhansell/Documents/GitHub//chrom13.vcf", type = "ncd1")
stat13<-ncd1(parsed_vcf13)

devtools::install_github("bitarellolab/balselr")
library(balselr)
