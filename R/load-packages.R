# Here are the names of all the packages needed
needed_packages <- c(
  "plyr",
  "dplyr",
  "stringr",
  "reshape2",
  "ggplot2",
  "RColorBrewer",
  "Rrunstruct",
  "lazyeval",
  "maps",
  "mapdata",
  "maptools",
  "grid"
)


### Here is a snippet of code that will install all the necessary packages.
## Run the code in the if(FALSE) block below to install all the needed packages.
if(FALSE) {
  ## first install all the regular packages you can get from CRAN
  # be warned that this will update any existing packages on your system
  reg_packages <- needed_packages[!needed_packages == "Rrunstruct"]
  install.packages(needed_packages, dependencies = TRUE)
  
  ## then get Rrunstruct from github
  
  # get devtools if needed
  if(any(rownames(installed.packages()) == "devtools") == FALSE) {
    install.packages("devtools", dependencies = TRUE)
  }
  # then get Rrunstruct
  devtools::install_github("eriqande/Rrunstruct", ref = "01a18b1df")
  
}  # end if(FALSE) block







### this loads all the packages, and throws an error if any of them fail
dump <- lapply(needed_packages, function(x) {
  success <- require(x, character.only = TRUE)
  if(!success) {
    stop("FAILURE!  Unable to load package ", x, " you might need to install it.  See the source code of ./R/load-packages for a code snippet that will download an install all needed packages.")
  }
})






