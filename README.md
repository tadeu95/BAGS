# *BAGS: Barcode, Audit & Grade System*

### WARNING: the information needed for grade assignment has not been updated in the weblinks since July 2023 due to a recent version of "rsconnect". If you run BAGS from your R environment, then the information is up to date.
### Available links for direct access to BAGS:
1. [Link1](https://tadeu-apps.shinyapps.io/bags)
2. [Link2](https://tadeu-apps.shinyapps.io/bags2)
3. [Link3](https://tadeu-apps.shinyapps.io/bags3)
4. [Link4](https://tadeu-apps.shinyapps.io/bags4)
5. [Link5](https://tadeu-apps.shinyapps.io/bags5)

### NOTE: the web links have limited capacity for the retrieval or large taxonomic groups, and they are limited to one user at a time. To run BAGS without limitations, host the application locally in your R environment

## **How to host BAGS locally in your R environment**

### Running BAGS from GitHub:

1. Download and install the most suited version of [R for your operating system.](https://www.r-project.org/)

2. Optionally you can download and install [RStudio.](https://rstudio.com/products/rstudio/download/)

3. Open RGui or Rstudio and run this line of code to install the necessary packages:
```
install.packages(c("seqRFLP","bold","data.table","plyr","httr","stringr","readr","fingerprint","dplyr","ggplot2","shiny","shinyWidgets","snakecase"))
```
4. The packages "seqRFLP" and "worms" have been removed from the CRAN repository. To install the package "seqRFLP" from the CRAN archive, run these lines of code:
```
url <- "http://cran.r-project.org/src/contrib/Archive/seqRFLP/seqRFLP_1.0.1.tar.gz"
pkgFile <- "seqRFLP_1.0.1.tar.gz"
download.file(url = url, destfile = pkgFile)
install.packages(pkgs=pkgFile, type="source", repos=NULL)
```
5. To install the package "worms" from the CRAN archive, run these  lines of code:
```
url <- "http://cran.r-project.org/src/contrib/Archive/worms/worms_0.2.2.tar.gz"
pkgFile <- "worms_0.2.2.tar.gz"
download.file(url = url, destfile = pkgFile)
install.packages(pkgs=pkgFile, type="source", repos=NULL)
```
6. Load the "shiny" package by running:
```
library(shiny)
```
7. Run the app from GitHub:
```
runGitHub("BAGS", "tadeu95")
```

### Running BAGS by manually downloading the app file to your computer:

1. Download and install the most suited version of [R for your operating system.](https://www.r-project.org/)

2. Optionally you can download and install [RStudio.](https://rstudio.com/products/rstudio/download/)

3. Open RGui or Rstudio and install the necessary packages:
```
install.packages(c("bold","data.table","stringr","readr","fingerprint","dplyr","ggplot2","shiny","shinyWidgets","snakecase"))
url <- "http://cran.r-project.org/src/contrib/Archive/seqRFLP/seqRFLP_1.0.1.tar.gz"
pkgFile <- "seqRFLP_1.0.1.tar.gz"
download.file(url = url, destfile = pkgFile)
install.packages(pkgs=pkgFile, type="source", repos=NULL)
url <- "http://cran.r-project.org/src/contrib/Archive/worms/worms_0.2.2.tar.gz"
pkgFile <- "worms_0.2.2.tar.gz"
download.file(url = url, destfile = pkgFile)
install.packages(pkgs=pkgFile, type="source", repos=NULL)
```
4. Go to "File" on the top left corner, click "Open file" / "Open script" and choose the "app" file.

5. In the top right corner of the opened app script file, click "run app" if you're on RStudio. Another alternative is to run the script line by line using ctrl+ENTER for each chunk of code. If you're on RGui, select all the text in the script file (Ctrl+A) then click the right button of the mouse and choose "Run line or selection" (Ctrl+R).

6. After running the script, a window with the app will pop up, and you can click "open in browser" to run the app in your default browser if you're on RStudio. 
If you're on RGui, the app will automatically open in your default browser.



