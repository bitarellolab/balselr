now <- Sys.time()
timestamp <- function(time) format(time, "%Y-%B-%d_%H-%M-%S")

#' Create output file with timestamp
#' @param infile Name of input file to modify
#' @examples outfile_path("inst/example.vcf")
#' @export
outfile_path <- function(infile = "example.vcf") {

        sub("\\..*", paste0("_",timestamp(now),".out"), infile)
}


.even_or_odd<-function(num){
        if((num %% 2) == 0) {
        return(0)
} else {
        return(1)
}
}

