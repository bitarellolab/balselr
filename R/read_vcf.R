#' Read vcf
#'
#' @param x
#'
#' @return
#' @export
#'
#' @examples
#' @import data.table
#' @importFrom data.table ":="
read_vcf <- function(x = "inst/example.vcf",
                     only.bi = T,
                     inds = "all") {
  inp <- data.table::fread(x, skip = "##", header = T)
  data.table::setnames(inp, "#CHROM", "CHR")
  if (only.bi == T) {
    inp <- inp[REF %in% c("A", "C", "T", "G")][ALT %in% c("A", "C", "T", "G")]
  } else {
    return(inp)
  }

  formatcol <- which(colnames(inp) == "FORMAT")

  if (inds != "all") {
    cbind(inp[, 1:formatcol], x2[, ...inds])
  } else {
    return(inp)
  }
}
