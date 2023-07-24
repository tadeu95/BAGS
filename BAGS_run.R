check.packages <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
}
packages<-c("bold","data.table","worms","stringr","readr","fingerprint","dplyr","ggplot2","shiny","shinyWidgets","snakecase","httr","plyr")
check.packages(packages)
if (!"seqRFLP" %in% installed.packages()[, "Package"]) {
  url <- "http://cran.r-project.org/src/contrib/Archive/seqRFLP/seqRFLP_1.0.1.tar.gz"
  pkgFile <- "seqRFLP_1.0.1.tar.gz"
  download.file(url = url, destfile = pkgFile)
  install.packages(pkgs=pkgFile, type="source", repos=NULL)
  }
if (!"worms" %in% installed.packages()[, "Package"]) {
  url <- "http://cran.r-project.org/src/contrib/Archive/worms/worms_0.2.2.tar.gz"
  pkgFile <- "worms_0.2.2.tar.gz"
  download.file(url = url, destfile = pkgFile)
  install.packages(pkgs=pkgFile, type="source", repos=NULL)
  }
library(shiny)
runGitHub("BAGS", "tadeu95")
