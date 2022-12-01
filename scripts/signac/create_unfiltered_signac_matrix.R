.libPaths("~/miniconda3/envs/archr/lib/R/library")

library(data.table)
library(stringr)
library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(EnsDb.Hsapiens.v75)
library(ggplot2)
library(patchwork)
set.seed(1234)

library(stringr)
library(dplyr)

# Load signac object

data_path <- '/lustre/groups/ml01/workspace/laura.martens/moretti_colab/cellranger/merged/outs/'

archr <- fread('/lustre/groups/ml01/workspace/laura.martens/moretti_colab/archr_final/Epicardiods/obs.csv')
archr$barcode = str_split(archr$V1, '#|-', simplify=TRUE)[,2]

atac <- readRDS(file ='/lustre/groups/ml01/workspace/laura.martens/moretti_colab/signac/atac_fltr.rds')

fragments <- CreateFragmentObject(paste0(data_path, "fragments.tsv.gz"))
counts <- FeatureMatrix(
  fragments = fragments,
  features = granges(atac)
)
print(dim(counts))
cellranger <- as.data.frame(str_split(colnames(counts), '-', simplify=TRUE))
colnames(cellranger) <- c('barcode', 'sample_id')
cellranger$colnames <- colnames(counts)

cellranger <- cellranger %>% mutate(Sample=recode(sample_id, "1" = "MUC26649_1234", "2"="MUC26650_1234", "3"="MUC26651_1234", "4"="MUC26652_1234", "5"="MUC26653_1234", "6"="MUC26654_1234", "7"="MUC26655_1234"))

metadata <- merge(cellranger, archr, by = c('barcode', 'Sample'), all.x = FALSE)
rownames(metadata) <- metadata$colnames

counts <- counts[, metadata$colnames]
print(dim(counts))

chrom_assay <- CreateChromatinAssay(
  counts = counts,
  sep = c(":", "-"),
  genome = 'hg19',
  fragments = paste0(data_path, "fragments.tsv.gz"),
  min.cells = 0,
  min.features = 0
)

atac <- CreateSeuratObject(
  counts = chrom_assay,
  assay = "peaks",
  meta.data = metadata
)


# extract gene annotations from EnsDb
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)

# change to UCSC style since the data was mapped to hg19
seqlevelsStyle(annotations) <- 'UCSC'

# add the gene information to the object
Annotation(atac) <- annotations

saveRDS(atac, file = "/lustre/groups/ml01/workspace/laura.martens/moretti_colab/signac/atac_archr.rds")