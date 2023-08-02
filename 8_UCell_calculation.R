# load libraries
library(escape)
library(Seurat)
library(rhdf5)
library(SeuratData)
library(SeuratDisk)
library(UCell)

# set working directory
setwd("/Volumes/Passport/CSDB/brain_scRNAseq/")

# convert the raw .h5ad file to seurat object
Convert('all_internship/data/processed/raw_merged_barcode_filtered.h5ad', dest = 'h5seurat', overwrite = TRUE)
raw_data <- LoadH5Seurat('all_internship/data/processed/raw_merged_barcode_filtered.h5seurat', meta.data = FALSE, misc = FALSE, assays = 'RNA')

# load in all expressed LINEs and LTRs subfamily names
LINEs <- readLines("all_internship/data/processed/LINEs.txt")

LTRs <- readLines("all_internship/data/processed/LTRs.txt")

# load in SenMayo gene set
SenMayo <- readRDS("Final/Utilities/out/gene_sets.rds")
SenMayo <- SenMayo$SenMayo

# specify gene sets
gene.sets <- list(LINEs = LINEs,
                  LTRs = LTRs,
                  SenMayo = SenMayo)

# add UCell scores
raw_data <- AddModuleScore_UCell(raw_data, features = gene.sets, name = NULL, ncores = 4)

# calculate top 10%
SenMayo_cutoff_top <- raw_data@meta.data$SenMayo
SenMayo_cutoff_top <- quantile(SenMayo_cutoff_top, probs = c(0.9))
raw_data@meta.data <- transform(raw_data@meta.data, senescent_top10_perc = ifelse(raw_data@meta.data$SenMayo > 0.09400847, "senescent", "other"))

# calculate top 10% in old condition
SenMayo_cutoff_top_old <- raw_data@meta.data[raw_data@meta.data$condition == 'old',]
SenMayo_cutoff_top_old <- quantile(SenMayo_cutoff_top_old$SenMayo, probs = c(0.9))
raw_data@meta.data <- transform(raw_data@meta.data, old_senescent_top10_perc = ifelse(raw_data@meta.data$SenMayo > 0.09482429 & raw_data@meta.data$condition == 'old', "senescent", "other"))

# write meta data to file
meta <- raw_data@meta.data
write.csv(meta, file = "all_internship/data/processed/metadata_ES.csv", row.names = TRUE)
