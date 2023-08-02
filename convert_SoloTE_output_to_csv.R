library(Seurat)
library(tidyverse)

options(timeout=100000000000)

# load the files from the sample folders in and assign to Seurat object named accordingly
input_dir = "/Volumes/Passport/CSDB/brain_scRNAseq/SoloTE_processedBam_output/"
sample_names = list()
for (sample_folder in list.files(input_dir)) {
  folder_dir = paste0(input_dir, sample_folder, '/out_SoloTE_output')
  file_names = list.files(folder_dir)
  sample_name <- substr(sample_folder, 1, 4)
  sample_names <- append(sample_names, sample_name)
  seurat_obj <- Read10X(data.dir = folder_dir)
  seurat_obj <- CreateSeuratObject(counts = seurat_obj,
                                   min.cells = 3,
                                   min.features = 200,
                                   project = sample_name)
  assign(sample_name, seurat_obj)
}

# subset the Seurat object into a Seurat object for the genes and one for the TEs
# create list with the names of the Seurat objects for TEs and genes
sample_genes_names = list()
sample_TEs_names = list()
for (sample in sample_names) {
  seurat_obj <- get(sample)
  genes <- rownames(seurat_obj)[!grepl('SoloTE', rownames(seurat_obj))]
  TEs <- rownames(seurat_obj)[grepl('SoloTE', rownames(seurat_obj))]
  seurat_genes <- subset(seurat_obj,
                         features = genes)
  seurat_TEs <- subset(seurat_obj,
                       features = TEs)
  name_genes <- paste0(sample, '_genes')
  sample_genes_names <- append(sample_genes_names, name_genes)
  name_TEs <- paste0(sample, '_TEs')
  sample_TEs_names <- append(sample_TEs_names, name_TEs)
  assign(name_genes, seurat_genes)
  assign(name_TEs, seurat_TEs)
}

# write the object for the genes to csv files for further processing in Python with ScanPy
output_dir = "/Volumes/Passport/CSDB/brain_scRNAseq/Final/Utilities/out/"
for (sample in sample_genes_names) {
  seurat_genes <- get(sample)
  folder_dir = paste0(output_dir, sample, '_raw_counts.csv')
  write.csv(seurat_genes@assays$RNA@counts, file = folder_dir)
}

rm(OX1X, OX1X_genes, OX2X, OX2X_genes, OX3X, OX3X_genes, OX4X, OX4X_genes, OX5X, OX5X_genes, OX6X, OX6X_genes, OX7X, OX7X_genes, OX8X, OX8X_genes,
   YX1L, YX1L_genes, YX2L, YX2L_genes, YX3R, YX3R_genes, YX4R, YX4R_genes, YX5R, YX5R_genes, YX6L, YX6L_genes, YX7R, YX7R_genes, YX8L, YX8L_genes)

rm(sample_genes_names, sample_names)
rm(seurat_genes, seurat_obj, seurat_TEs)

# write the object for the TEs to csv files
output_dir = "/Volumes/Passport/CSDB/brain_scRNAseq/Final/Utilities/out/"
for (sample in sample_TEs_names) {
  seurat_TEs <- get(sample)
  folder_dir = paste0(output_dir, sample, '_raw_counts.csv')
  write.csv(seurat_TEs@assays$RNA@counts, file = folder_dir)
}

