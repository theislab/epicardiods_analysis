.libPaths("/mnt/home/icb/laura.martens/miniconda3/envs/archr/lib/R/library")

library(Signac)
library(data.table)
library(stringr)
library(Seurat)

sampleFile = fread('/storage/groups/ml01/workspace/laura.martens/moretti_colab/signac/sample_file.csv')
inputFiles <- sampleFile$fragments


peaks <- CallPeaks(inputFiles)

setwd("/storage/groups/ml01/workspace/laura.martens/moretti_colab/transfer_data/")
frag_obj <- list()
for(frag_file in inputFiles){
    frag_obj[frag_file] <- CreateFragmentObject(frag_file, validate.fragments = FALSE )
}


counts <- FeatureMatrix(
  fragments = frag_obj,
  features = peaks
)


chrom_assay <- CreateChromatinAssay(
  counts,
  fragments = frag_obj,
  ranges=peaks,
genome = 'hg19'
)

adata <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
)
saveRDS(adata, file = "/storage/groups/ml01/workspace/laura.martens/moretti_colab/signac/atac.rds")