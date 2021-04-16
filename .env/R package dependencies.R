# System
options("repos" = c(CRAN = "https://cloud.r-project.org"))
library("devtools")
library("BiocManager")

install_github("jalvesaq/colorout", ref = "7ea9440", upgrade = "never")

# QTLseqr and its dependencies
install_version("modeest",    version = "2.4.0",   upgrade = "never")
install_version("ggplot2",    version = "3.3.3",   upgrade = "never")
install_version("readr",      version = "1.4.0",   upgrade = "never")
install_version("Rcpp",       version = "1.0.1",   upgrade = "never")
install_version("pillar",     version = "1.3.1",   upgrade = "never")
install_version("pkgconfig",  version = "2.0.2",   upgrade = "never")
install_version("tibble",     version = "2.1.1",   upgrade = "never")
install_version("purrr",      version = "0.3.2",   upgrade = "never")
install_version("tidyselect", version = "1.1.0",   upgrade = "never")
install_version("dplyr",      version = "0.8.2",   upgrade = "never")
install_version("tidyr",      version = "1.1.2",   upgrade = "never")
install_version("locfit",     version = "1.5-9.4", upgrade = "never")
install_github("bmansfeld/QTLseqr", ref = "5e76137", upgrade = "never")

#install_bioc("3.8/rtracklayer", upgrade="never")