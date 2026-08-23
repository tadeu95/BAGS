# *BAGS: Barcode, Audit & Grade System*

<div style="border: 2px solid red; padding: 15px; background-color: #ffe6e6; font-size: 18px;">
  <h1 style="color:red; text-align:center;">⚠️ Important Notice ⚠️</h1>

  <p>
    <strong>BAGS is currently not operational.</strong>
  </p>

  <p>
    BAGS was originally developed to retrieve barcode records from BOLD Systems using the BOLD v4 API and the R package <code>bold</code>. Following the transition of BOLD Systems to version 5 and its new API infrastructure, the v4 services on which BAGS relies no longer provide the data access required by the application. Consequently, BAGS is currently unable to retrieve the BOLD records needed to perform its analyses and grade assignments.
  </p>

  <p>
    In addition, the R package <code>bold</code>, which BAGS uses for data retrieval, was removed from the CRAN repository in August 2024. Although archived versions of the package remain available, they were developed for the BOLD v4 API and therefore do not restore the current functionality of BAGS.
  </p>

  <p>
    <strong>We are currently assessing the possibility of updating BAGS to work with the new BOLD v5 API.</strong> This would require adapting the data-retrieval workflow and parts of the application to the structure and endpoints of the new API.
  </p>

  <p>
    In the meantime, users interested in DNA barcode reference library assessment and curation may wish to explore the <strong>Biodiversity Genomics Europe (BGE) Reference Library Pipeline and Curation Tool</strong>. This tool was developed by partners of the Biodiversity Genomics Europe (BGE) project and incorporates a BAGS species-level assessment alongside additional quality-control criteria, phylogenetic analyses, and manual curation functionalities.
  </p>

  <p>
    <strong>More information about the BGE Reference Library Pipeline and Curation Tool is available here:</strong><br>
    <strong><a href="https://iboleurope.org/bge-reference-library-pipeline-and-curation-tool/">https://iboleurope.org/bge-reference-library-pipeline-and-curation-tool/</a></strong>
  </p>

  <p>
    We hope to provide further information regarding the future development of BAGS as the transition to the BOLD v5 infrastructure progresses.
  </p>
</div>




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



