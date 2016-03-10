source("R/load-packages.R")


if(!file.exists("outputs")) dir.create("outputs")

# read in the full repository information:
repo <- tbl_df(read.table("data/meta/AM001_AM006.tab", sep = "\t", header = T, stringsAsFactors = FALSE))






# read in the bycatch IDs info. is just a data frame of the NMFS_DNA_IDs
# that correspond to Bycatch fish from the region and also gives information
# about the duplicate tissues they sent us.
bycid <- readRDS("data/meta/bycatch_IDS.rds")

# read in the sheet of data about the samples we are using
samsheet <- read.csv("data/meta/sample_sheet.csv", stringsAsFactors = FALSE)

# so, in the future, we can left_join on samsheet to get a consistent set of 
# individuals to use.

#### Here is how I created sample categories and made a sample list.  Doesn't get re-run... ####
if(FALSE) {  # this doesn't typically get re-run because I made the file to use in the future
  sams <- repo %>% 
    filter(!(NMFS_DNA_ID %in% bycid$Duplicate_Tissue)) %>%  # toss out one of the duplicates of each duplicated one
    filter(!(NMFS_DNA_ID %in% c("AM000193"))) %>%   # bycatch individuals without any good meta-data in the Region's data base
    filter(!(NMFS_DNA_ID %in% c("AM000196"))) %>%  # third tissue sample from the same bycatch indiv.  Toss it.
    filter(!(NMFS_DNA_ID %in% c("AM000408"))) %>%  # duplicate tissue from a non-bycatch individual
    filter(NMFS_DNA_ID != "AM000474") %>%    # individual for which proper meta data are not available
    mutate(
      category = ifelse(WATERSHED %in% c("Sacramento River", "Klamath River"), 
                        "reference", "non-reference"),
      life_stage = ifelse(toupper(REPORTED_LIFE_STAGE) == "JUVENILE",
                          "juvenile", ifelse(str_detect(SAMPLE_COMMENTS, "^GST"), "egg", "non-juvenile")),
      collection_location = ifelse(NMFS_DNA_ID %in% bycid$NMFS_DNA_ID, "Bycatch", WATERSHED)
      ) %>% 
    select(NMFS_DNA_ID, category, life_stage, collection_location) %>%
    mutate(life_stage = ifelse(collection_location == "Klamath River", "juvenile", life_stage)) %>% # note that the klamath are all juvenile, but not recorded in meta data
    mutate(group_name = paste(
      str_sub(category, 1, 1), 
      str_sub(life_stage, 1, 1), 
      str_sub(collection_location, 1, 3),
      sep = ""
      )) %>%
    mutate(pipe_name = paste(group_name, str_sub(NMFS_DNA_ID, 6, 8), sep = ""))
  

  
  # write that out 
  write.csv(sams, file = "data/meta/sample_sheet.csv", row.names = FALSE)
}

if(FALSE) {
#### Summarize the Bycatch Data ####
# filter out things that aren't in the repo (there aren't any)
byc <- bycid %>%
  filter(NMFS_DNA_ID %in% repo$NMFS_DNA_ID)



# we want to make a LaTeX table with column and row totals
btab <- as.matrix(table(byc$Location, byc$Year)) %>%
  cbind(NA, Totals = rowSums(.)) %>%
  rbind(NA, Totals = colSums(.))

# now write that as a latex table
tmp <- cbind(Location = rownames(btab), btab) 
tmp[is.na(tmp)] <- " "
outf <- "outputs/bycatch_table.tex"
cat(paste(" &  \\multicolumn{", ncol(btab) - 2, "}{c}{\\underline{Year of Capture}} \\\\\n", sep = ""), file = outf)
cat(colnames(tmp), sep = " & ", file = outf, append = TRUE)
cat(" \\\\ \\hline \n", file = outf, append = TRUE)
write.table(tmp, sep = " & ", eol = "\\\\\n", 
            row.names = F, col.names = F, file = outf, 
            quote = FALSE, append = TRUE)

}  # took this out for data privacy reasons.  Code here if we want to go
                 # back to make such a table

#### Make a table of sample counts ####

stab <- samsheet %>%
  mutate(category = factor(category, levels = rev(unique(category)))) %>%
  mutate(collection_location = factor(collection_location, levels = rev(sort(unique(collection_location))))) %>%
  group_by(collection_location, category, life_stage, group_name) %>%
  tally() %>%
  setNames(c("Location", "Category", "Life Stage", "Group Short Name", "n"))


outf <- "outputs/sample_summary.tex"
cat(names(stab), sep = " & ", file = outf)
cat("\\\\ \\hline\n", file = outf, append = TRUE)
write.table(stab, sep = " & ", eol = "\\\\\n", col.names = F, 
            row.names = F, quote = F, file = outf, append = T)
         
         
         
         
         
         
         
         