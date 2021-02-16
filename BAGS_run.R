check.packages <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
}
packages<-c("seqRFLP","bold","data.table","worms","stringr","readr","fingerprint","dplyr","ggplot2","shiny","shinyWidgets","snakecase")
check.packages(packages)
library(shiny)
runGitHub("BAGS", "tadeu95")