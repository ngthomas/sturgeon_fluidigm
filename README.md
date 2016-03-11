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
  
  Get all the R packages that you need.  Last I checked these were:
    ```r
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
    ```
    But, you should check the file `./R/load-packages.R` for the latest.  That file also has a code snippet in it that
    will let you download and install all the needed packages in a few lines of R code.
    
2. __Metadata__: Each of the genetic samples was received with a variety of information with them.
  * `data/meta/AM001_AM006.tab` is a tab-delimited dump of all the SWFSC repository information for all the samples that appear in the analysis. The data were exported with both marine and freshwater sample information into this file by Cassie Columbus.
  * `data/meta/sample_sheet.csv` is a CSV file that just shows which category each fish is in, keyed by
  its NMFS_DNA_ID.
  * Files starting with `data/meta/private*` do not live with the repository.  These are data such as geographic locations of bycatch that are considered confidential.  So, if you want to reproduce _all_ of the work here (specifically the plots of sampling locations), you won't be able to do it.  Those analyses have to be done by Carlos who has access to those data.  Note that `private-green-sturgeon-reconciled.xlsx` is the file that Carlos uses to make the map with the lat-long locations of bycaught fish.  You will see that filename in the code.  And now you will understand why that part of the code does not run for you.
  * `data/meta/bycatch_IDS.rds` holds the IDs of fish in the bycatch data. It is no longer considered confidential in this form because it ain't nothing but the IDs.  It is used in analyses and summaries here.  It contains the `NMFS_DNA_ID`s of the fish that were from the bycatch. It also includes information about which `NMFS_DNA_ID`s contain tissue that was sampled twice from the same fish.  The column of `NMFS_DNA_ID` is what should be used for most analyses, and the duplicates will be used for matching genotypes verifications.
3. __Genetic Data__: The genetic data that appear in this repo are of two different varieties.  
  * The first are the relatively "raw" data off the Fluidigm machine. These data are intensity data that have not been normalized using the negative control values.  In order to access these raw values in a typical analysis, you have to contact Fluidigm and ask for an access code that lets you access the raw data. There are four chips of these sorts of data that were used for training our scoring procedure. These are "Detailed Table" results from Fluidigm and they include the automatic genotype calls (which are intended for diploids, and hence not very reliable) and the intensities of the different dye-sets. They are located in `data/training_chips`.  Specifically they are:
    - `data/training_chips/1381905043_raw.csv`
    - `data/training_chips/1381935339_raw.csv`
    - `data/training_chips/1381992037_raw.csv`
    - `data/training_chips/1381992302_raw.csv`
  * After those four training chips were used to develop the tetraploid scoring methodology, we scored those four plates using that methodology.  The resulting scored genotypes for the different training individuals are in `data/training_chips/training_chips_scored.csv`. This file contains the genotypes that can be assigned from the Fluidigm genotyping software.  However, the Fluidigm software is designed for only the 3 genotypes of a diploid.  Hence, at loci at which we can reliably score more than three genotypes, we overload some of the genotype clusters in that we might call a cluster at the lower right of the image as an "XX" and then we also will score a cluster at the upper left as an "XX" and then, during post-processing, we correct the calls in the upper left to refer to either the 4th or 5th cluster.  This is taken care of in `R-main/02-develop-chip-scoring.R`.


### Scripts to Run the analyses.

The analyses are done in a series of R scripts.  We have all the necessary scripts set up to run in 
order in the directory `./R-main/`.  Each of them should be run with the current
working directory in R set to be the top level of the repository. (i.e. the directory
that includes the subdirectory `R-main`).  Many of these scripts write out plots or data
files that hold summarized analyses.  These all get written to a directory `./outputs` that
gets created if it does not exist.  All the contents of `./outputs` are .gititnored.

Here we briefly describe what each script does:

1. `R-main/01-meta-data-summaries.R`: Basically just count up how many individuals are in each
category of sample and write out a table.  This also loads the "sample sheet" and the bycatch 
IDs which are used later on, too.  There is some code inside this script that doesn't get run,
but it remains there to shows how the sample sheet was created.

1. `R-main/02-develop-chip-scoring.R`:  This script operates on the data that are found in the `*_raw.csv` files in `data/training_chips`. The script reads the data in and then creates 4 pages of plots.  Each page has 24 panels and each panel has the data from one locus across the four chips. Line segments connect the same individual typed on different plates.  This helps us identify reliable clusters. The outputs from this script include the PDF files of the figures:
  - `outputs/plate_x_y1.pdf`
  - `outputs/plate_x_y2.pdf`
  - `outputs/plate_x_y3.pdf`
  - `outputs/plate_x_y4.pdf`

After those outputs came out, they were used to develop the scoring methodology and we scored those training chips accordingly. We have provided images of all the intensities of the points as in the `plate_x_y*` outputs above, but coloring things by genotype.  This was done both for all plates combined and for each plate separately in files with names of the forms:
  - `outputs/plate_x_y_by_final_genotype_plate_1381905043_only1.pdf`
  - `outputs/plate_x_y_by_final_genotype_plate_1381905043_only2.pdf`
  -  . . . 
  - `outputs/late_x_y_by_final_genotype_plate_1381935339_only1.pdf`
  - `outputs/plate_x_y_by_final_genotype_plate_1381935339_only2.pdf`
  -  . . .
  -  `outputs/plate_x_y_by_final_genotype_plates_combined1.pdf`
  -  `outputs/plate_x_y_by_final_genotype_plates_combined2.pdf`
  -  . . .



2. `R-main/03-score-plate-five.R`: This file shows the sort of workflow that is adopted to score chips beyond the training chips (specifically looking at how we score chip #5).  It shows the workflow for plotting figures of the raw intensity data for the training chips and for chip 5 together on the same plots, then it reads in the scored version of chip 5.  This scored version of the plate is the output from the Fluidigm Software after you use the software to manually score the individuals at each locus.  The scored file is saved at `./data/more_chips/1382136064_scored.csv`.  The header on it looks like the following, in case that helps people to understand what format/type of output from the Fluidigm software it is:

  ```
  Chip Run Info,C:\Users\biopipe\Desktop\Sturgeon5rr_SNPtype_1382136064\ChipRun.bml,1382136064,96.96 (138x),,ROX,FAM-MGB : VIC-MGB,11/5/2015 3:10:49     PM,00:01:22,EP1-50017
  Application Version,3.1.3
  Application Build,20120816.1511
  Export Type,Detailed Table Results,Standard
  Number of Combined Chip Runs,1,9216
  Confidence Threshold,65.00
  Normalization Method,NTC Normalization
  Allele Probe Type Mapping,Allele X,FAM-MGB
  Allele Probe Type Mapping,Allele Y,VIC-MGB
  Allele Axis Mapping,Allele X,X
  Allele Axis Mapping,Allele Y,Y
  ```
  
  The rest of the script just converts those scores from the Fluidigm software into call of the genotype categories we have defined at these loci. Then it combines chip 5 with the training chips and saves the resulting data frame into `./outputs/genotype_from_five_chips.rds`.  Subsequently the script analyzes the patterns of missing data and creates the plot `outputs/successful-assay-histogram-crop.pdf`

1. `R-main/04-do-structure-and-gsi.R`:  This script creates a STRUCTURE input file from the data and then run STRUCTURE on it using the slg_pipe commands and programs in the directory StructureArea.  Then it runs CLUMP and DISTRUCT on the output and creates the distruct plots for the paper and writes out the output.  Then it runs `gsi-sim` on the bycatch data too.  

1. `R-main/05-indiv-id-sims.R`: Does the simulations of the log-likelihood ratio for individual identification, and then it computes those log-likelihood ratios for all pairs of fish genotypes.

1. `R-main/06-bycatch-map.R`: Creates the maps of the bycatch fish.  