library(tidyverse)
library(DESeq2)
library(future)
library(furrr)

dds <- readRDS("/data/rudensky/EYW/git_projects/SIG06/analysis_outs/dds_object_SIG07.rds")

# import feature names (used in function below)
featureNames <- read_csv("/data/rudensky/EYW/git_projects/SIG06/processing_outs/featureNames_SIG06.csv")
featureNames <- select(featureNames, -category)

# sample input can be single character or vector of characters
# control input must be single character, defaults to linker
iterative_DEG <- function(samples, control = "linker"){
  # empty tibble for full DEG list
  resAll <- tibble()
  
  for(ligand in samples){
    # print progress
    print(paste("processing",ligand))
  
    # perform IHW DEG analysis
    # don't perform independent filtering because using IHW
    res <- results(dds, contrast = c("treatment",ligand,control),
                 independentFiltering=F,
                 filterFun = ihw)
    
    # perform LFC shrinkage and make tidy
    # pass res to lfcShrink
    resTidy <- lfcShrink(dds, type = "ashr", res=res) %>%
      as_tibble(rownames = "ensembl_ID") %>%
      left_join(featureNames) %>%
      mutate(treatment = ligand) %>%
      arrange(padj)
    
    # add to full DEG list
    resAll <- bind_rows(resAll, resTidy) 
  }
  return(resAll)
}

# set options for parallel analysis
plan(multisession, workers=parallelly::availableCores())

# select viral treatments groups
viralGroups <- dds$treatment %>% unique()
viralGroups <- viralGroups[grep("recomb|none|linker", viralGroups, invert = T)]

resFull <- future_map_dfr(viralGroups, ~ iterative_DEG(.x))
write_csv(resFull, "/data/rudensky/EYW/git_projects/SIG06/analysis_outs/deg_full_SIG06.csv")

resSig <- resFull %>% filter(padj < 0.1)
write_csv(resSig, "/data/rudensky/EYW/git_projects/SIG06/analysis_outs/deg_sig_SIG06.csv")
