This is a small tutorial on how to host and use the app in your own R environment

How to run from GitHub:

1-Download and install the most suited version of R for your computer at https://www.r-project.org/

2-(Optional) Download and install Rstudio at https://rstudio.com/products/rstudio/download/

3- Open RGui or Rstudio and run the following line of code to install the required packages to launch the app:

install.packages(c("seqRFLP","bold","data.table","worms","stringr","readr","fingerprint","dplyr","ggplot2","shiny","shinyWidgets","snakecase"))

4- Run the following line of code to load the "shiny" package:

library(shiny)

5- Launch the app by running:

runGitHub("BAGS", "tadeu95")


How to launch the app by manually downloading the app file to your computer:

1-Download and install the most suited version of R for your computer at https://www.r-project.org/

2-(Optional) Download and install Rstudio at https://rstudio.com/products/rstudio/download/

3-Open RGui or Rstudio and run the following line of code to install the required packages to launch the app:

install.packages(c("seqRFLP","bold","data.table","worms","stringr","readr","fingerprint","dplyr","ggplot2","shiny","shinyWidgets"))

4- Go to "File" on the top left corner, click "Open file" / "Open script" and choose the "app" file.

5-In the top right corner of the opened app script file, click "run app" if you're on RStudio. Another alternative is to run the script line by line using ctrl+ENTER for each chunk of code. If you're on RGui, select all the text in the script file (Ctrl+A) then click the right button of the mouse and choose "Run line or selection" (Ctrl+R).

6-After running the script, a window with the app will pop up, and you can click "open in browser" to run the app in your default browser if you're on RStudio. 
If you're on RGui, the app will automatically open in your default browser.
