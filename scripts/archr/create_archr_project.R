library(ArchR)
library(data.table)
library(parallel)
sampleFile = fread('/storage/groups/ml01/workspace/laura.martens/moretti_colab/archr_new/sample_file.csv')
inputFiles <- sampleFile$fragments
names(inputFiles) <- sampleFile$library_id
addArchRGenome("hg19")

ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  filterTSS = 2, #Dont set this too high because you can always increase later
  filterFrags = 300, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE,
)

doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search.
  LSIMethod = 1
)

proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "Epicardiods",
  copyArrows = TRUE #This is recommened so that you maintain an unaltered copy for later usage.
)
proj <- filterDoublets(ArchRProj = proj)

proj <- addIterativeLSI(ArchRProj = proj, useMatrix = "TileMatrix", name = "IterativeLSI")
proj <- addClusters(input = proj, reducedDims = "IterativeLSI")
proj <- saveArchRProject(ArchRProj = proj)