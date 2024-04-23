library(data.table)
vcf_obj<-data.table::fread("data-raw/example.vcf", skip = "##", header = TRUE)
devtools::use_data(vcf_obj, overwrite = TRUE)
