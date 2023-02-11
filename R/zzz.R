.onAttach <- function(libname, pkgname) {
  # output the specific message when package is loaded!
  RFver <- read.dcf(file=system.file("DESCRIPTION", package=pkgname),
                    fields=c("Description"))
  website <-  "Check out our Package website (https://feiyoung.github.io/PRECAST/index.html) for a more complete description of the methods and analyses"
  packageStartupMessage(paste(pkgname,': ', RFver, " ", website))
}
