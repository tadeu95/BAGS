# *BAGS: Barcode, Audit & Grade System*

<div style="border: 2px solid red; padding: 15px; background-color: #ffe6e6; font-size: 18px;">
  <h1 style="color:red; text-align:center;">⚠️ Important Notice ⚠️</h1>
  <p>
The BAGS output may be affected due to the removal of the <code>bold</code> R package from the CRAN repository. While earlier versions of the <code>bold</code> package remain available in the archive, the removal from CRAN may impact its accessibility and usability. BAGS was designed to use the BOLD version 4 API, and with the ongoing transition to BOLD version 5, its current functionality and compatibility remain uncertain. Additionally, the reference files used for grade assignments are currently static and cannot be updated. As a result, some grade assignments may include minor inaccuracies or outdated information. We are actively working on a solution to restore full functionality.
  </p>
  <p>
    <strong>Please be cautious when interpreting results and aware of potential limitations or inaccuracies.</strong>
  </p>
</div>




### Available links for direct access to BAGS:
1. [Link1](https://tadeu-apps.shinyapps.io/bags)
2. [Link2](https://tadeu-apps.shinyapps.io/bags2)
3. [Link3](https://tadeu-apps.shinyapps.io/bags3)
4. [Link4](https://tadeu-apps.shinyapps.io/bags4)
5. [Link5 (new column with problematic BINs for grade E species in testing)](https://tadeu-apps.shinyapps.io/bags5)

### NOTE: the web links have limited capacity for the retrieval or large taxonomic groups, and they are limited to one user at a time. To run BAGS without limitations, host the application locally in your R environment

## **How to host BAGS locally in your R environment**

### Running BAGS from GitHub:

1. Download and install the most suited version of [R for your operating system.](https://www.r-project.org/)

2. Optionally you can download and install [RStudio.](https://rstudio.com/products/rstudio/download/)

3. Open RGui or Rstudio and run this line of code to install the necessary packages:
```
install.packages(c("data.table","plyr","httr","stringr","readr","fingerprint","dplyr","ggplot2","shiny","shinyWidgets","snakecase"))
```
4. The packages "bold", "seqRFLP", and "worms" have been removed from the CRAN repository. To install the package "bold" from the CRAN archive, run these lines of code:

```
install.packages("crul")
url <- "http://cran.r-project.org/src/contrib/Archive/bold/bold_1.3.0.tar.gz"
pkgFile <- "bold_1.3.0.tar.gz"
download.file(url = url, destfile = pkgFile)
install.packages(pkgs=pkgFile, type="source", repos=NULL)
```
5. To install the package "seqRFLP" from the CRAN archive, run these lines of code:

```
url <- "http://cran.r-project.org/src/contrib/Archive/seqRFLP/seqRFLP_1.0.1.tar.gz"
pkgFile <- "seqRFLP_1.0.1.tar.gz"
download.file(url = url, destfile = pkgFile)
install.packages(pkgs=pkgFile, type="source", repos=NULL)
```
6. To install the package "worms" from the CRAN archive, run these  lines of code:
```
url <- "http://cran.r-project.org/src/contrib/Archive/worms/worms_0.2.2.tar.gz"
pkgFile <- "worms_0.2.2.tar.gz"
download.file(url = url, destfile = pkgFile)
install.packages(pkgs=pkgFile, type="source", repos=NULL)
```
7. Load the "shiny" package by running:
```
library(shiny)
```
8. Run the app from GitHub:
```
runGitHub("BAGS", "tadeu95")
```

### Running BAGS by manually downloading the app file to your computer:

1. Download and install the most suited version of [R for your operating system.](https://www.r-project.org/)

2. Optionally you can download and install [RStudio.](https://rstudio.com/products/rstudio/download/)

3. Open RGui or Rstudio and install the necessary packages:
```
install.packages(c("data.table","stringr","readr","fingerprint","dplyr","ggplot2","shiny","shinyWidgets","snakecase"))
url <- "http://cran.r-project.org/src/contrib/Archive/bold/bold_1.3.0.tar.gz"
pkgFile <- "bold_1.3.0.tar.gz"
download.file(url = url, destfile = pkgFile)
install.packages(pkgs=pkgFile, type="source", repos=NULL)
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



