This is a small tutorial on how to host and use the app in your own R environment

How to run the app from github:

1-Download and install the most suited version of R for your computer at https://www.r-project.org/

2-(Optional) Download and install Rstudio at https://rstudio.com/products/rstudio/download/

3- Open RGui or RStudio

4-Install the packages by running the following line of code:

install.packages(c("seqRFLP","bold","data.table","worms","stringr","readr","fingerprint","dplyr","ggplot2","shiny","shinyWidgets"))

5- Run the following line of code to load the package "shiny":

library(shiny)

6- Launch the app by running:

runGitHub("BAGS", "tadeu95")

7-You can now start using BAGS.


How to launch the app by manually downloading the app file to your computer:

1-Download and install the most suited version of R for your computer at https://www.r-project.org/

2-(Optional) Download and install Rstudio at https://rstudio.com/products/rstudio/download/

3-Open RGui or Rstudio and run the following line of code to install the required packages to launch the app:

install.packages(c("seqRFLP","bold","data.table","worms","stringr","readr","fingerprint","dplyr","ggplot2","shiny","shinyWidgets"))

4-If you're on RStudio, go to "File" on the top left corner, click "Open file" and choose the "app" file.

5-In the top right corner of the opened script, click "run app" if you're on RStudio. Another alternative is to run the script line by line using ctrl+ENTER for each chunk of code. If you're on RGui, just copy and paste the entire script to the command line and click the "enter" key.

6-After running the script, a window with the app will pop up, and you can click "open in browser" to run the app in your default browser if you're on RStudio. 
If you're on RGui, the app will automatically open in your default browser.

7-You can now start using BAGS.

