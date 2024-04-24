library(devtools)
devtools::install_github("bitarellolab/balselr")
library(balselr)


test_read<-read_vcf("test7.vcf")
head(test_read)

test_parse1<-parse_vcf("test7.vcf", type="ncd1", n0=27)
head(test_parse1)
test_parse2<-parse_vcf("test7.vcf", type="ncd2", n0=27, n1=4)
head(test_parse2)

test_ncd1<-ncd1(x=test_parse1, tf=0.5, w=3000, ncores=2, minIS=2)
head(test_ncd1)
test_ncd2<-ncd2(x=test_parse2, tf=0.5, w=3000, ncores=2, minIS=2)
head(test_ncd2)?ncd2
