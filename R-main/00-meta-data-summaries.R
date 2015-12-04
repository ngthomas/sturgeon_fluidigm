

# read in the full repository information:
repo <- tbl_df(read.table("data/meta/AM001_AM006.tab", sep = "\t", header = T, stringsAsFactors = FALSE))

# read in the bycatch IDs info:
bycid <- readRDS("data/meta/bycatch_IDS.rds")


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



#### Now summarize the freshwater samples ####

# toss out the ones we know are bycatch
nonbyc <- repo %>%
  filter(!(NMFS_DNA_ID %in% byc$NMFS_DNA_ID), !(NMFS_DNA_ID %in% byc$Duplicate_Tissue_Of))
         
         
         
         
         
         
         
         
         
         
         
         