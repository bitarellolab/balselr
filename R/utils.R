now <- Sys.time()
timestamp <- function(time) format(time, "%Y-%B-%d_%H-%M-%S")

#' Create output file with timestamp
#' @param infile Name of input file to modify
#' @examples outfile_path("inst/example.vcf")
#' @export
outfile_path <- function(infile = "example.vcf") {

        sub("\\..*", paste0("_",timestamp(now),".out"), infile)
}

.get_fds<-function(x){ #x: number of FDS, y=tf
        rep(0, x)
}


utils::globalVariables(c(
        "all_of", "across", "Mid", "temp2", "CenMaf", "MidMaf", "Method",
        "temp4", "IS", ".", "FDs", "AF", "AF2", "ALT", "CHR", "FD", "ID",
        "MAF", "NCD1", "NCD2", "POS", "REF", "SNP", "Win.ID", "tn_1", "tn_2",
        "tx_1", "tx_2", "is","end","start", "fd", "S","POS2"
))
