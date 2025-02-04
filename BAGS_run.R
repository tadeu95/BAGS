# Function to check and install CRAN packages
check.packages <- function(pkg) {
  new.pkg <- pkg[!(pkg %in% installed.packages()[, "Package"])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
}

# List of required CRAN packages (excluding archived ones)
packages <- c("data.table", "stringr", "readr", "fingerprint", "dplyr", 
              "ggplot2", "shiny", "shinyWidgets", "snakecase", "httr", "plyr")

# Install missing CRAN packages
check.packages(packages)

# Install archived seqRFLP package if missing
if (!"seqRFLP" %in% installed.packages()[, "Package"]) {
  url <- "http://cran.r-project.org/src/contrib/Archive/seqRFLP/seqRFLP_1.0.1.tar.gz"
  pkgFile <- "seqRFLP_1.0.1.tar.gz"
  download.file(url = url, destfile = pkgFile)
  install.packages(pkgs = pkgFile, type = "source", repos = NULL)
}

# Install archived worms package if missing
if (!"worms" %in% installed.packages()[, "Package"]) {
  url <- "http://cran.r-project.org/src/contrib/Archive/worms/worms_0.2.2.tar.gz"
  pkgFile <- "worms_0.2.2.tar.gz"
  download.file(url = url, destfile = pkgFile)
  install.packages(pkgs = pkgFile, type = "source", repos = NULL)
}

# Install archived bold package if missing
if (!"bold" %in% installed.packages()[, "Package"]) {
  url <- "http://cran.r-project.org/src/contrib/Archive/bold/bold_1.3.0.tar.gz"
  pkgFile <- "bold_1.3.0.tar.gz"
  download.file(url = url, destfile = pkgFile)
  install.packages(pkgs = pkgFile, type = "source", repos = NULL)
}

# Load shiny and run GitHub app
library(shiny)
runGitHub("BAGS", "tadeu95")

