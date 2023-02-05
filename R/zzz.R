.onAttach <- function(libname, pkgname) {
  # output the specific message when package is loaded!
  RFver <- read.dcf(file=system.file("DESCRIPTION", package=pkgname),
                    fields=c("Description"))
  packageStartupMessage(paste(pkgname,': ', RFver))
}
