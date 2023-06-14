library("data.table")


read_vcf("example_vcf")


read_vcf("/Users/dhansell/Documents/GitHub//chrom22.vcf")
parsed_vcf22<-parse_vcf("/Users/dhansell/Documents/GitHub//chrom22.vcf", type = "ncd1")

stat22<-ncd1(parsed_vcf22)
