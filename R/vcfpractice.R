library("data.table")


read_vcf("example_vcf")

read_vcf("chrom22.vcf")
parse_vcf("chrom22.vcf", type = "ncd1")
