# sturgeon_fluidigm


## Steps To Reproduce the Analysis

### Preliminaries and Data Sources
All the data and scripts to reproduce our analyses are included herein, but
before you begin, you should study this carefully.

1. __Preliminaries__: First, you have to get the repository.  Do it like this:

    ```sh
    git clone https://github.com/ngthomas/sturgeon_fluidigm
    ```
  
  Then you should ensure that you have R (version > 3.2) and all the R libraries that we use:
  
    ```r
    install.packages(
    c("MASS", 
      "dplyr", 
      "ggplot2", 
      "grid", 
      "gridExtra", 
      "plyr", 
      "reshape2", 
      "scales", 
      "stringr"
    ))
    ```
    
2. __Metadata__: Each of the genetic samples was received with a variety of information with them.
  * `data/meta/acipenser.mer` holds the standard meta data attached to the samples in the SWFSC repository.
  * `data/meta/private.csv` does not live with the repository.  These are data such as geographic locations of bycatch that are considered confidential.  So, if you want to reproduce _all_ of the work here, you need to get that file from Carlos and put it in the correct location in the repository (and be sure not to commit  it --- note that it is gitignored.)
  * `data/meta/summarized.csv` holds the confidential data summarized sufficiently that it is no longer considered confidential, and is used in analyses here.
3. __Genetic Data__: The genetic data that appear in this repo are of two different varieties.  
  * The first are the totally "raw" data off the Fluidigm machine.  There are four chips of these sorts of data that were used for training our scoring procedure. These are "Detailed Table" results from Fluidigm and they include the automatic genotype calls (which are intended for diploids, and hence not very reliable) and the intensities of the different dye-sets. They are located in `data/four_training_chips`.  Specifically they are:
    - `data/four_training_chips/1381905043.csv`
    - `data/four_training_chips/1381935339.csv`
    - `data/four_training_chips/1381992037.csv`
    - `data/four_training_chips/1381992302.csv`
  * After those four training chips were used to develop the tetraploid scoring methodology, we scored those four plates using that methodology.  The results of these are in the directory `data/four_scored_chips`.  They are in a similar format to the data files in `four_training_chips` but they __NEED TO KNOW FROM THOMAS HOW THESE ARE DIFFERENT__.  Specifically these files are:
    - `data/four_scored_chips/recall1_1381905043.csv`
    - `data/four_scored_chips/recall1_1381935339.csv`
    - `data/four_scored_chips/recall1_1381992037.csv`
    - `data/four_scored_chips/recall1_1381992302.csv`


### Scripts to Run the analyses.

The analyses are done in a series of R scripts.  We have all the necessary scripts set up to run in 
order in the directory `./R-main/`.  Each of them should be run with the current
working directory in R set to be the top level of the repository. (i.e. the directory
that includes the subdirectory `R-main`). Here we briefly describe what each script does:

1. _01-develop-chip-scoring.R_  describe, describe
2. _02-process-scored-training-data-chips.R_ describe describe
