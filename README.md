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
  * The first are the relatively "raw" data off the Fluidigm machine. These data are intensity data that have not been normalized using the negative control values.  In order to access these raw values in a typical analysis, you have to contact Fluidigm and ask for an access code that lets you access the raw data. There are four chips of these sorts of data that were used for training our scoring procedure. These are "Detailed Table" results from Fluidigm and they include the automatic genotype calls (which are intended for diploids, and hence not very reliable) and the intensities of the different dye-sets. They are located in `data/training_chips`.  Specifically they are:
    - `data/training_chips/1381905043_raw.csv`
    - `data/training_chips/1381935339_raw.csv`
    - `data/training_chips/1381992037_raw.csv`
    - `data/training_chips/1381992302_raw.csv`
  * After those four training chips were used to develop the tetraploid scoring methodology, we scored those four plates using that methodology.  The resulting scored genotypes for the different training individuals are in `data/training_chips/training_chips_scored.csv`. This file contains the genotypes that can be assigned from the Fluidigm genotyping software.  However, the Fluidigm software is designed for only the 3 genotypes of a diploid.  Hence, at loci at which we can reliably score more than three genotypes, we overload some of the genotype clusters in that we might call a cluster at the lower right of the image as an "XX" and then we also will score a cluster at the upper left as an "XX" and then, during post-processing, we correct the calls in the upper left to refer to either the 4th or 5th cluster.  This is taken care of in `R-main/01-develop-chip-scoring.R`.


### Scripts to Run the analyses.

The analyses are done in a series of R scripts.  We have all the necessary scripts set up to run in 
order in the directory `./R-main/`.  Each of them should be run with the current
working directory in R set to be the top level of the repository. (i.e. the directory
that includes the subdirectory `R-main`). Here we briefly describe what each script does:

1. `R-main/01-develop-chip-scoring.R`:  This script operates on the data that are found in the `*_raw.csv` files in `data/training_chips`. The script reads the data in and then creates 4 pages of plots.  Each page has 24 panels and each panel has the data from one locus across the four chips. Line segments connect the same individual typed on different plates.  This helps us identify reliable clusters. The outputs from this script are the PDF files of the figures:
  - `outputs/plate_x_y1.pdf`
  - `outputs/plate_x_y2.pdf`
  - `outputs/plate_x_y3.pdf`
  - `outputs/plate_x_y4.pdf`

 After those outputs came out, they were used to develop the scoring methodology and we scored those training chips accordingly. We have provided images of all the intensities of the points as in the plate_x_y outputs above, but coloring things by genotype.  This was done both for all plates combined and for each plate separately in files with names of the forms:
  - `outputs/plate_x_y_by_final_genotype_plate_1381905043_only1.pdf`
  - `outputs/plate_x_y_by_final_genotype_plate_1381905043_only2.pdf`
  -  . . . 
  - `outputs/late_x_y_by_final_genotype_plate_1381935339_only1.pdf`
  - `outputs/plate_x_y_by_final_genotype_plate_1381935339_only2.pdf`
  -  . . .
  -  `outputs/plate_x_y_by_final_genotype_plates_combined1.pdf`
  -  `outputs/plate_x_y_by_final_genotype_plates_combined2.pdf`
  -  . . .



2. `R-main/02-process-scored-training-data-chips.R` describe describe
