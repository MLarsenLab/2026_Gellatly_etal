if (!require("pacman")) install.packages("pacman")
pacman::p_load(Matrix, Seurat, SeuratWrappers, scDblFinder, hdf5r)
install_github("chris-mcginnis-ucsf/DoubletFinder")
library(DoubletFinder)

#More space for globals 
options(future.globals.maxSize = 50 * 1024^3)

Read_CellBender_h5_Mat <- function(file_name, use.names = TRUE, unique.features = TRUE) {
  # https://samuel-marsh.github.io/scCustomize/reference/Read_CellBender_h5_Mat.html
  # Check hdf5r installed
  if (!requireNamespace('hdf5r', quietly = TRUE)) {
    cli_abort(message = c("Please install hdf5r to read HDF5 files",
                          "i" = "`install.packages('hdf5r')`")
    )
  }
  # Check file
  if (!file.exists(file_name)) {
    stop("File not found")
  }
  
  if (use.names) {
    feature_slot <- 'features/name'
  } else {
    feature_slot <- 'features/id'
  }
  
  # Read file
  infile <- hdf5r::H5File$new(filename = file_name, mode = "r")
  
  counts <- infile[["matrix/data"]]
  indices <- infile[["matrix/indices"]]
  indptr <- infile[["matrix/indptr"]]
  shp <- infile[["matrix/shape"]]
  features <- infile[[paste0("matrix/", feature_slot)]][]
  barcodes <- infile[["matrix/barcodes"]]
  
  
  sparse.mat <- sparseMatrix(
    i = indices[] + 1,
    p = indptr[],
    x = as.numeric(x = counts[]),
    dims = shp[],
    repr = "T"
  )
  
  if (unique.features) {
    features <- make.unique(names = features)
  }
  
  rownames(x = sparse.mat) <- features
  colnames(x = sparse.mat) <- barcodes[]
  sparse.mat <- as(object = sparse.mat, Class = "dgCMatrix")
  
  infile$close_all()
  
  return(sparse.mat)
}
removeDoublets <- function(seuratObject, pctDbl = 20){
  #uses scDblFinder and DoubletFinder to find droplets that may contain multiple cells
  #Returns a seruat object where cells that both scDblFinder and DoubletFinder agree on as doublets have been removed

  #scDblFinder needs a single cell experiment as input and outputs a score and class (doublet or singlet)
  #Assume 20% doublets, this is an over estimation of number of doublets as we are looking for agreement 
  aSCE <- as.SingleCellExperiment(seuratObject)
  aSCE <- scDblFinder(aSCE, dbr = (pctDbl/100))
  
  #copy the cell assignments to new metadata in the seurat object
  seuratObject$scDblFinder.score <- aSCE$scDblFinder.score
  seuratObject$scDblFinder.class <- aSCE$scDblFinder.class
  
  #DoubletFinder requires data to be normalized before it can identify doublets
  #The object passed in may not be ready to be normalized or may need to be normalized in a different way so
  #a copy is made and the copy is normalized for DoubletFinder and the cell assignments copied to the original
  dblSeurat <- seuratObject
  dblSeurat <- NormalizeData(dblSeurat)
  dblSeurat <- FindVariableFeatures(dblSeurat, selection.method = "vst", nfeatures = 2000)
  geneList <- rownames(dblSeurat)
  dblSeurat <- ScaleData(dblSeurat, feature = geneList)
  dblSeurat <- RunPCA(dblSeurat, features = VariableFeatures(object = dblSeurat))
  dblSeurat <- RunUMAP(dblSeurat, dims = 1:40)
  
  #DoubletFinder needs to know how many doublets to look for, default of using 20% of total number of cells
  dblCount <- as.integer(length(rownames(dblSeurat@meta.data))*(pctDbl/100))
  dblSeurat <- doubletFinder(dblSeurat, PCs = 1:40, pN = 0.25, pK = 0.1, nExp = dblCount, reuse.pANN = FALSE, sct = FALSE)
  #DoubletFinder returns the classification in a new metadata with the name based on the parameters used
  DF <- sprintf("DF.classifications_0.25_0.1_%i", dblCount)
  DFIndex <- which(colnames(dblSeurat@meta.data) == DF)
  
  #Copy cell assignment to the seurat object metadata
  seuratObject$DF.classifications <- dblSeurat@meta.data$DF
  
  #Find cells that are classified as a doublet by both scDblFinder and doubletFinder
  doublets1 <- subset(seuratObject, scDblFinder.class == "doublet")
  trueDoublets <- subset(doublets1, DF.classifications == "Doublet")
  seuratObject$True_Doublet <- ifelse(colnames(seuratObject)%in% colnames(trueDoublets), "Doublet", "Singlet")
  
  #Print the percentage of doublets found
  print(sprintf("%i doublets found in %i total cells (%.2f%%)", ncol(trueDoublets), ncol(seuratObject), 
                (ncol(trueDoublets)/ncol(seuratObject)*100)))
  
  #Make a subset of only singlets and return the cleaned object
  seuratObject <- subset(seuratObject, True_Doublet == "Singlet")
  return(seuratObject)
}
loadForProcess <- function(fileNamePath, cellBenderP = TRUE, projectName = "Project", pctDouble = 10){
  #Create a seurat object from a data file and do basic processing of the seurat object before returning the object
  #Accepts either cellranger output file or cellbender output file, defaults to cellbender output
  #debug assignments
  if(cellBenderP){
    aMatrix <- Read_CellBender_h5_Mat(fileNamePath)
    aName <- str_sub(basename(fileNamePath), 1, (str_length(basename(fileNamePath))-12))
    aSeurat <- CreateSeuratObject(aMatrix, project = aName, min.cells = 3, min.features = 200)
    aSeurat@misc["raw_cell_count"] <- ncol(aSeurat)
  }else{
    a10X <- Read10X(data.dir = fileNamePath)
    aName <- projectName
    aSeurat <- CreateSeuratObject(a10X, project = aName, min.cells = 3, min.features = 200)
    aSeurat@misc["raw_cell_count"] <- ncol(aSeurat)
  }
  #find mitochondria gene percentage and filter on to few (emtpy) to many (doubletts) features and over 25% mitochondria genes
  aSeurat[["percent.mt"]] <- PercentageFeatureSet(aSeurat, pattern = "^mt-")
  aSeurat <- subset(aSeurat, subset = nFeature_RNA > 200 & nFeature_RNA < 9000 & percent.mt < 25)
  aSeurat@misc["live_cell_count"] <- ncol(aSeurat)
  
  #remove doublets and return processed seurat object
  aSeurat <- removeDoublets(aSeurat, pctDbl = pctDouble)
  aSeurat@misc["cell_count"] <- ncol(aSeurat)
  return(aSeurat)
}
process10x <- function(seuratObject, projName = "project", percentMT = 25, clusterResolution = 1){
  #Input a seurat object
  #output normalized and clustered seurat object, also creates files for basic QC of data in projName folder in CWD
  dir.create(projName)
  
  seuratObject[["percent.mt"]] <- PercentageFeatureSet(seuratObject, pattern = "^mt-")
  
  pdf(paste(projName, "/Feature_Count_MT_Pre.pdf", sep = ""), width = 40, height = 30)
  print(VlnPlot(seuratObject, features = c("nFeature_RNA", "nCount_RNA", "percent.mt")) +
          theme(axis.title = element_text(size = 60), axis.text = element_text(size = 60), title = element_text(size = 60))
  )
  dev.off()
  
  seuratObject <- subset(seuratObject, subset = nFeature_RNA > 200 & nFeature_RNA < 9000 & percent.mt < percentMT)
  
  pdf(paste(projName, "/Feature_Count_MT_Post.pdf", sep = ""), width = 40, height = 30)
  print(VlnPlot(seuratObject, features = c("nFeature_RNA", "nCount_RNA", "percent.mt")) +
          theme(axis.title = element_text(size = 40), axis.text = element_text(size = 40)) 
  ) 
  dev.off()
  
  seuratObject <- SCTransform(seuratObject, vst.flavor = "v2", verbose = FALSE)
  
  seuratObject <- RunPCA(seuratObject, features = VariableFeatures(object = seuratObject))
  
  aplot <-ElbowPlot(seuratObject, ndims = 50)
  pdf(paste(projName, "/ElbowPlot.pdf", sep = ""), width = 40, height = 30)
  print(aplot)
  dev.off()
  
  seuratObject <- FindNeighbors(seuratObject, dims = 1:40)
  seuratTree <- FindClusters(seuratObject, resolution = (seq(0, (clusterResolution*2), by = (clusterResolution/5))))
  
  pdf(paste(projName, "/clusttree.pdf", sep = ""), width = 40, height = 30)
  print(clustree(seuratTree, prefix = "SCT_snn_res."))
  dev.off()
  
  seuratObject <- FindClusters(seuratObject, resolution = clusterResolution)
  
  seuratObject <- RunUMAP(seuratObject, dims = 1:40)
  
  pdf(paste(projName, "/umap.pdf", sep = ""), width = 40, height = 30)
  print(DimPlot(seuratObject, reduction = "umap", label = TRUE, label.size = 16, pt.size = 0.5, repel = TRUE) +
          NoLegend() +
          theme(axis.title = element_text(size = 40), axis.text = element_text(size = 40)))
  dev.off()
  
  #Acinar Marker genes and genes that tend to be contaminants in other cells
  featureList = c("Mucl2", "Muc19", "Aqp5", "Pecam1", "Cdh5", "Prol1", "Sox9", "Sox2", "Sox10", "Krt7", 
                    "Krt5", "Krt14", "Trp63", "Acta2", "Kcnip4", "Dcpp1", "Dcpp2", "Dcpp3")
  
  pdf(paste(projName, "/DotPlot.pdf", sep = ""), width = 40, height = 30)
  print(DotPlot(seuratObject, features = featureList, cols = c("lightblue", "red3"), 
                col.min = 0, dot.scale = 20) +
          theme(axis.title = element_text(size = 60), axis.text = element_text(size = 60), 
                legend.text = element_text(size = 60), legend.key.size = unit(3, 'cm'), 
                legend.title  = element_text(size = 60), legend.spacing = unit(3, 'cm')) +
          RotatedAxis())
  dev.off()
  
  seuratObject$dataset <- projName
  return(seuratObject)
}
mergeMeta <- function(x, y, ...){
  #input paramters to pass to seurats merge function
  #output returns merged seurat object with preserving metadata, metadata for each seurat object is saved under that objects 'project.name' in misc metadata of merged object
  retVal <- merge(x = x, y = y, ...)
  
  key <- x@project.name
  retVal@misc[[key]] <- x@misc
  
  if("list" == typeof(y)){
    for(i in 1:length(y)){
      key <- y[[i]]@project.name
      retVal@misc[[key]] <- y[[i]]@misc
    }
  } else {
    key <- y@project.name
    retVal@misc[[key]] <- y@misc
  }
  return(retVal)
}
formattedUMAP <- function(seuratObject){
  
  if(50000 > ncol(seuratObject)){
    aPtSize <- 0.2
  }
  else{
    aPtSize <- 0.05
  }
  
  aPlot <- DimPlot(seuratObject, reduction = "umap", label = FALSE, label.size = 7, pt.size = aPtSize, repel = TRUE, raster = FALSE) +
    NoLegend() +
    theme(axis.title = element_text(size = 15, face = "bold"), axis.text = element_text(size = 15, face = "bold"))
  aPlot <- LabelClusters(aPlot, id = "ident", fontface = "bold", color = "black", size = 7)
  return(aPlot)
}
highlightMeta <- function(seuratObject, path = ".", meta = "orig.ident", type = "PDF"){
  
  Idents(seuratObject) <- meta
  dataSetList <- mixedsort(unique(seuratObject@meta.data[[meta]]))
  colorList <- c(rep("grey", length(dataSetList)-1), "red")
  for(i in 1:length(unique(seuratObject@meta.data[[meta]]))){
    print(dataSetList)
    afileName <- sprintf("%s/UMAP_%s_highlight", path, dataSetList[1])
    
    if("PDF" == type)
    {
      pdf(sprintf("%s.pdf", afileName), width = 10, height = 7, bg = "white")
      print(DimPlot(seuratObject, reduction = "umap", label = FALSE, label.size = 16, pt.size = 0.5, repel = TRUE,
                    order = dataSetList, cols = colorList, raster = FALSE) +
              theme(axis.title = element_text(size = 15, face = "bold"), axis.text = element_text(size = 15, face = "bold")))
      dev.off()
    }else if ("TIFF" == type){
      tiff(file = sprintf("%s.tiff", afileName), units = "in", 
           res = 400, height = 7, width = 10, compression = "none")  
      print(DimPlot(seuratObject, reduction = "umap", label = FALSE, label.size = 16, pt.size = 0.5, repel = TRUE,
                    order = dataSetList, cols = colorList, raster = FALSE) +
              theme(axis.title = element_text(size = 15, face = "bold"), axis.text = element_text(size = 15, face = "bold")))
      dev.off()     
    }
    else{
      message(sprintf("Type %s not supported", type))
      return()
    }
    
    dataSetList <- c(dataSetList[2:length(dataSetList)], dataSetList[1])
  }
  
}
convert_human_to_mouse <- function(gene_list) {
  #Adapted from https://www.biostars.org/p/9567892/
  #Converts a list of human gene names to a list of mouse gene equivalents
  output = c()
  mouse_human_genes = read.csv("https://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",sep="\t")
  
  for(gene in gene_list) {
    class_key = (mouse_human_genes %>% filter(Symbol == gene & Common.Organism.Name == "human"))[['DB.Class.Key']]
    if( !identical(class_key, integer(0)) ) {
      human_genes = (mouse_human_genes %>% filter(DB.Class.Key == class_key & Common.Organism.Name=="mouse, laboratory"))[,"Symbol"]
      for(human_gene in human_genes) {
        output = rbind(c(gene, human_gene), output)
      }
    }
  }
  return (rev(output[,2]))
}
getAllDEGGeneNames <- function(aSeurat, aIdent2 = "", aMetadata = "", skipNormalization = FALSE){
  #aSeurat, seurat object to run the DEG from
  #aIdent2, Ident to compair all other Idents to, ie all vs Homeostatic, if blank, Ident vs all
  #aMetadata if not blank, which metadata field aIdent belongs to
  #skipNormailzation, skip PrepSCTFindMarkers and SCTransform, speeds up if already run
  if(!skipNormalization){
    aSeurat <- PrepSCTFindMarkers(aSeurat)
    aSeurat <- SCTransform(aSeurat)
  }
  if("" != aMetadata){
    Idents(aSeurat) <- aMetadata
  }
  geneList <- data.frame(tmp = "")
  if("" != aIdent2){
    for(aIdent in levels(aSeurat)){
      if(aIdent == aIdent2){next}
      aDF <- FindMarkers(aSeurat, ident.1 = aIdent, ident.2 = aIdent2, min.pct = 0.25,
                         logfc.threshold = 0.25, only.pos = TRUE)
      aDF <- aDF[aDF$p_val_adj <= 0.05, ]
      for(aGene in rownames(aDF)){
        geneList[aGene, aIdent] <- aDF[aGene, "avg_log2FC"]
      }
    }
  }else{
    aDF <- FindAllMarkers(aSeurat, min.pct = 0.25, logfc.threshold = 0.25, only.pos = TRUE)
    aDF <- aDF[aDF$p_val_adj <= 0.05, ]
    for(aIdent in unique(aDF$cluster)){
      for(aGene in aDF[aDF$cluster == aIdent, "gene"]){
        geneList[aGene, aIdent] <- aDF[aDF$cluster == aIdent & aDF$gene == aGene, "avg_log2FC"]
      }
    }
  }
  
  geneList$tmp <- NULL
  geneList <- geneList[-1, , drop = FALSE]
  geneList[is.na(geneList)] <- 0
  geneList$max_value <- apply(geneList, 1, max, na.rm = TRUE)
  geneList <- geneList[order(geneList$max_value, decreasing = TRUE), ]
  
  return(rownames(geneList))
}
DotPlotFromEnrichment <- function(aSeurat, FileName = "", outputPath = ".", outputPrefix = "Dotplot", title = "Title", topCount = -1, topIdent1 = "", topIdent2 = "", filterGenes = c()){
  if(!file.exists(FileName)){
    warning(sprintf("File %s not found.", FileName))
    return()
  }
  
  if(-1 != topCount & topIdent1 != "" & topIdent2 != ""){
    aDFMarkers <- FindMarkers(aSeurat, ident.1 = topIdent1, ident.2 = topIdent2, 
                              min.pct = 0.25, logfc.threshold = 0.25, only.pos = TRUE)
    aDFMarkers <- aDFMarkers[aDFMarkers$p_val_adj <= 0.05, ]
    aDFMarkers <- aDFSubset[order(aDFSubset$avg_log2FC, decreasing = TRUE), ]
    if(0 == length(filterGenes)){
      #No filter, order by Markers
      filterGenes <- rownames(aDFMarkers)
    }else{
      #FilterGenes is NOT sorted (topeIdent1/2 set)
      filterGenes <- filterGenes[order(match(filterGenes, rownames(aDFMarkers)))]
    }
  }
  #filterGenes should have ordered list of genes, either ordered when passed in (topIdent1/2 = ""), 
  #Passed in filter genes Ordered by avg_log2FC from FindMarkers, 
  #or set to the ordered list of FindMarkers
  
  aDF <- read.csv(file = FileName)
  if(colnames(aDF)[[1]] == "Enrichment.FDR"){
    #This is a ShinyGo file
    for(index in 1:nrow(aDF)){
      symbolList <- strsplit(aDF[index, "Genes"], split = "  ")[[1]]
      symbolList <- gsub(" ", "", symbolList)
      aDesc <- aDF[index, "Pathway"]
      
      DotPlotFromList(aSeurat, 
                      outputPath = outputPath, 
                      outputPrefix = sprintf("%s-%s", outputPrefix, index), 
                      title = title,
                      description = aDesc,
                      topCount = topCount,
                      features = symbolList, 
                      geneOrderFilter = filterGenes)
      
    }
  }else if(colnames(aDF)[[1]] == "GroupID"){
    #This is a Metascape enrichment tab saved as a .csv
    
    aDF_Filtered <- aDF[grepl("Summary", aDF$GroupID),]
    aDF_Filtered$GroupID <- sub("_.*", "", aDF_Filtered$GroupID)
    for(aDesc in aDF_Filtered$Description){
      symbolList <- strsplit(aDF_Filtered[aDF_Filtered$Description == aDesc, "Symbols"], split = ",")[[1]]
      index <- aDF_Filtered[aDF_Filtered$Description == aDesc, "GroupID"]
      
      DotPlotFromList(aSeurat, 
                      outputPath = outputPath, 
                      outputPrefix = sprintf("%s-%s", outputPrefix, index),
                      title = title,
                      description = aDesc,
                      topCount = topCount,
                      features = symbolList, 
                      geneOrderFilter = filterGenes)
    }
  }else{
    warning(sprintf("Format not supported for file %s", aFileName))
  }
}
DotPlotFromList <- function(aSeurat, outputPath = ".", outputPrefix = "Dotplot", title = "Title", description = "", topCount = -1, features = c(), geneOrderFilter = c()){
  
  #geneOrderFilter is an ordered list of gene names to filter by
  if(0 != length(geneOrderFilter)){
    features <- features[features %in% geneOrderFilter]
    features <- features[order(match(features, geneOrderFilter))]
  }
  if(-1 != topCount & topCount <= length(features)){
    features <- features[1:topCount]
  }
  if(0 == length(features)){
    warning(sprintf("No features found for %s", outputPrefix))
    return()
  }
  
  #Image size settings, if large number of terms, make the .pdf bigger
  if((length(features) > 30)){
    textSize <- 8 
    aWidth <- 20
    aHeight <- 7
  }else{
    textSize <- 12 
    aWidth <- 10
    aHeight <- 7
  }
  
  aFileName <- sprintf("%s/%s %s.pdf", outputPath, outputPrefix, gsub("[/:]", "-", description))
  pdf(file = aFileName, height = aHeight, width = aWidth)
  print(DotPlot(aSeurat, features = features, scale.min = 0, scale.max = 100) +
          theme(axis.text.x = element_text(size = textSize)) +
          #NoLegend() +
          theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
          ggtitle(sprintf("%s\n%s", title, description)))
  dev.off()
}
seuratIntegration <- function(seuratObject, referenceOI = c()){
  
  if("list" != typeof(seuratObject)){
    seuratList <- c(seuratObject)
  }else{
    seuratList <- seuratObject
  }
  seuratOIList <- c()
  print("Breaking seurat object into orig.ident")
  for(i in 1:length(seuratList)){
    Idents(seuratList[[i]]) <- "orig.ident"
    for(origIdent in levels(seuratList[[i]])){
      seuratOIList <- c(seuratOIList, subset(seuratList[[i]], idents = origIdent))
    }
  }
  
  print("SCTransform start")
  seuratOIList <- lapply(X = seuratOIList, FUN = function(x) {
    x <- SCTransform(x, vst.flavor = "v2", verbose = FALSE)
  })
  
  options(future.globals.maxSize = 50 * 1024^3)
  print("SelectIntegrationFeatures start")
  features <- SelectIntegrationFeatures(object.list = seuratOIList, nfeatures = 3000)
  seuratOIList <- PrepSCTIntegration(object.list = seuratOIList, anchor.features = features)
  print("FindIntegrationAnchors Start")
  if(length(referenceOI) >= 1){
    anchors <- FindIntegrationAnchors(object.list = seuratOIList, anchor.features = features, 
                                      normalization.method = "SCT", reference = referenceOI)
  }else{
    anchors <- FindIntegrationAnchors(object.list = seuratOIList, anchor.features = features, 
                                      normalization.method = "SCT")
  }
  rm(seuratOIList)
  rm(features)
  gc()
  
  print("IntegrateData start")
  aSeuratIntegrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
  rm(anchors)
  gc()
  print("Integration finished")
  
  DefaultAssay(aSeuratIntegrated) <- "integrated"
  aSeuratIntegrated <- SCTransform(aSeuratIntegrated, vst.flavor = "v2", verbose = FALSE)
  aSeuratIntegrated <- RunPCA(aSeuratIntegrated, npcs = 40, verbose = FALSE)
  aSeuratIntegrated <- RunUMAP(aSeuratIntegrated, reduction = "pca", dims = 1:40)
  aSeuratIntegrated <- FindNeighbors(aSeuratIntegrated, reduction = "pca", dims = 1:40)
  aSeuratIntegrated <- FindClusters(aSeuratIntegrated, resolution = 1)
  DefaultAssay(aSeuratIntegrated) <- "SCT"
  aSeuratIntegrated <- PrepSCTFindMarkers(aSeuratIntegrated)
  
  return(aSeuratIntegrated)
}

#Data loading, processing and integration--------------------------------------------------
#Load and merge Homeostatic, 2Week, 4Week, 8Week and 12Week IR cellbender files
#Merged timepoint seurat files are saved 
#input Cellbender processed files for 3 replicates for all timepoints in projectPath/cellbender, defaults to c:/Gellatly_2026/cellbender
#output saved RDS file of integrated and cleaned dataset

projectPath <- "C:/Gellatly_2026"
CellbenderFileBase <- sprintf("%s/cellbender", projectPath)

vDatasets <- c("Homeostatic", "2WeekIR", "4WeekIR", "8WeekIR", "12WeekIR")
for(aDataset in vDataSets){
  seuratList <- c()
  for(i in 1:3){
    afilePath <- sprintf("%s/%s_%i_cellbender_FPR_0.01_filtered.h5", CellbenderFileBase, aDataset, i)
    aSeurat <- loadForProcess(afilePath, cellBenderP = TRUE)
    seuratList <- c(seuratList, aSeurat)
  }
  aSeuratMerge <- mergeMeta(seuratList[[1]], y = seuratList[2:3], add.cell.ids = c("N1", "N2", "N3"), project = sprintf("%sMerge", aDataset))
  aSeurat <- process10x(aSeuratMerge, projName = aDataset, percentMT = 25, clusterResolution = 1)
  saveRDS(aSeurat, file = "seurat_%s_merge.rds", aDataset)
}

#Start integration
#Uses Homeostatic replicates as reference for integration
IntegrationSeuratList <- c()
for(aDataset in vDatasets){
  IntegrationSeuratList <- c(IntegrationSeuratList, readRDS(sprintf("seurat_%s_merge.rds", aDataset)))
}

seuratOIList <- c()
for(i in 1:length(IntegrationSeuratList)){
  Idents(IntegrationSeuratList[[i]]) <- "orig.ident"
  for(origIdent in levels(IntegrationSeuratList[[i]])){
    seuratOIList <- c(seuratOIList, subset(IntegrationSeuratList[[i]], idents = origIdent))
  }
}

for(i in 1:length(seuratOIList)){
  DefaultAssay(seuratOIList[[i]]) <- "RNA"
}

seuratOIList <- lapply(X = seuratOIList, FUN = function(x) {
  x <- SCTransform(x, vst.flavor = "v2", verbose = FALSE)
})

features <- SelectIntegrationFeatures(object.list = seuratOIList, nfeatures = 3000)
seuratOIList <- PrepSCTIntegration(object.list = seuratOIList, anchor.features = features)

anchors <- FindIntegrationAnchors(object.list = seuratOIList, anchor.features = features, 
                                  normalization.method = "SCT", reference = c(1, 2, 3))

aSeuratIntegrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
rm(anchors)

DefaultAssay(aSeuratIntegrated) <- "integrated"
aSeuratIntegrated <- RunPCA(aSeuratIntegrated, npcs = 40, verbose = FALSE)
aSeuratIntegrated <- RunUMAP(aSeuratIntegrated, reduction = "pca", dims = 1:40)
aSeuratIntegrated <- FindNeighbors(aSeuratIntegrated, reduction = "pca", dims = 1:40)
aSeuratIntegrated <- FindClusters(aSeuratIntegrated, resolution = 1)

DefaultAssay(aSeuratIntegrated) <- "SCT"
aSeuratIntegrated <- SCTransform(aSeuratIntegrated)
aSeuratIntegrated <- PrepSCTFindMarkers(aSeuratIntegrated)
saveRDS(aSeuratIntegrated, file = "seurat_Integrated.RDS")
rm(aSeuratIntegrated)

#Celltype ID and subsetting-----------------------------------------------

aSeuratC57IR <- readRDS(file = "seurat_Integrated.RDS")
aSeuratC57IR <- FindClusters(aSeuratC57IR, graph.name = "integrated_snn", resolution = 1.5)

aSeuratC57IR$ID <- "Missing"

aSeuratC57IR$ID[aSeuratC57IR$seurat_clusters == '0'] <- "Endothelial"

aSeuratC57IR$ID[aSeuratC57IR$seurat_clusters == '1'] <- "Endothelial"
aSeuratC57IR$ID[aSeuratC57IR$seurat_clusters == '2'] <- "Immune"
aSeuratC57IR$ID[aSeuratC57IR$seurat_clusters == '3'] <- "Epithelial"
aSeuratC57IR$ID[aSeuratC57IR$seurat_clusters == '4'] <- "Immune"
aSeuratC57IR$ID[aSeuratC57IR$seurat_clusters == '5'] <- "Epithelial"
aSeuratC57IR$ID[aSeuratC57IR$seurat_clusters == '6'] <- "Immune"
aSeuratC57IR$ID[aSeuratC57IR$seurat_clusters == '7'] <- "Epithelial"
aSeuratC57IR$ID[aSeuratC57IR$seurat_clusters == '8'] <- "Other"
aSeuratC57IR$ID[aSeuratC57IR$seurat_clusters == '9'] <- "Epithelial"
aSeuratC57IR$ID[aSeuratC57IR$seurat_clusters == '10'] <- "Fibroblast"

aSeuratC57IR$ID[aSeuratC57IR$seurat_clusters == '11'] <- "MNP"
aSeuratC57IR$ID[aSeuratC57IR$seurat_clusters == '12'] <- "Epithelial"
aSeuratC57IR$ID[aSeuratC57IR$seurat_clusters == '13'] <- "Epithelial"
aSeuratC57IR$ID[aSeuratC57IR$seurat_clusters == '14'] <- "Epithelial"
aSeuratC57IR$ID[aSeuratC57IR$seurat_clusters == '15'] <- "Epithelial"
aSeuratC57IR$ID[aSeuratC57IR$seurat_clusters == '16'] <- "Epithelial"
aSeuratC57IR$ID[aSeuratC57IR$seurat_clusters == '17'] <- "Other"
aSeuratC57IR$ID[aSeuratC57IR$seurat_clusters == '18'] <- "Epithelial"
aSeuratC57IR$ID[aSeuratC57IR$seurat_clusters == '19'] <- "Epithelial"
aSeuratC57IR$ID[aSeuratC57IR$seurat_clusters == '20'] <- "Myoepithelial"

aSeuratC57IR$ID[aSeuratC57IR$seurat_clusters == '21'] <- "Endothelial"
aSeuratC57IR$ID[aSeuratC57IR$seurat_clusters == '22'] <- "Immune"
aSeuratC57IR$ID[aSeuratC57IR$seurat_clusters == '23'] <- "Epithelial"
aSeuratC57IR$ID[aSeuratC57IR$seurat_clusters == '24'] <- "Immune"
aSeuratC57IR$ID[aSeuratC57IR$seurat_clusters == '25'] <- "Epithelial"
aSeuratC57IR$ID[aSeuratC57IR$seurat_clusters == '26'] <- "Epithelial"
aSeuratC57IR$ID[aSeuratC57IR$seurat_clusters == '27'] <- "Endothelial"
aSeuratC57IR$ID[aSeuratC57IR$seurat_clusters == '28'] <- "Epithelial"
aSeuratC57IR$ID[aSeuratC57IR$seurat_clusters == '29'] <- "Epithelial"
aSeuratC57IR$ID[aSeuratC57IR$seurat_clusters == '30'] <- "MNP"

aSeuratC57IR$ID[aSeuratC57IR$seurat_clusters == '31'] <- "Epithelial"
aSeuratC57IR$ID[aSeuratC57IR$seurat_clusters == '32'] <- "Epithelial"
aSeuratC57IR$ID[aSeuratC57IR$seurat_clusters == '33'] <- "Epithelial"
aSeuratC57IR$ID[aSeuratC57IR$seurat_clusters == '34'] <- "Endothelial"
aSeuratC57IR$ID[aSeuratC57IR$seurat_clusters == '35'] <- "Epithelial"
aSeuratC57IR$ID[aSeuratC57IR$seurat_clusters == '36'] <- "Epithelial"
aSeuratC57IR$ID[aSeuratC57IR$seurat_clusters == '37'] <- "Other"
aSeuratC57IR$ID[aSeuratC57IR$seurat_clusters == '38'] <- "Epithelial"
aSeuratC57IR$ID[aSeuratC57IR$seurat_clusters == '39'] <- "Epithelial"
aSeuratC57IR$ID[aSeuratC57IR$seurat_clusters == '40'] <- "Other"

aSeuratC57IR$ID[aSeuratC57IR$seurat_clusters == '41'] <- "MNP"
aSeuratC57IR$ID[aSeuratC57IR$seurat_clusters == '42'] <- "MNP"
aSeuratC57IR$ID[aSeuratC57IR$seurat_clusters == '43'] <- "Epithelial"
aSeuratC57IR$ID[aSeuratC57IR$seurat_clusters == '44'] <- "Dividing"
aSeuratC57IR$ID[aSeuratC57IR$seurat_clusters == '45'] <- "Epithelial"
aSeuratC57IR$ID[aSeuratC57IR$seurat_clusters == '46'] <- "Epithelial"
aSeuratC57IR$ID[aSeuratC57IR$seurat_clusters == '47'] <- "Other"
aSeuratC57IR$ID[aSeuratC57IR$seurat_clusters == '48'] <- "Epithelial"
aSeuratC57IR$ID[aSeuratC57IR$seurat_clusters == '49'] <- "Other"
aSeuratC57IR$ID[aSeuratC57IR$seurat_clusters == '50'] <- "Epithelial"

aSeuratC57IR$ID[aSeuratC57IR$seurat_clusters == '51'] <- "Other"
aSeuratC57IR$ID[aSeuratC57IR$seurat_clusters == '52'] <- "Fibroblast"
aSeuratC57IR$ID[aSeuratC57IR$seurat_clusters == '53'] <- "Epithelial"
aSeuratC57IR$ID[aSeuratC57IR$seurat_clusters == '54'] <- "Fibroblast"
aSeuratC57IR$ID[aSeuratC57IR$seurat_clusters == '55'] <- "Epithelial"
aSeuratC57IR$ID[aSeuratC57IR$seurat_clusters == '56'] <- "Other"
aSeuratC57IR$ID[aSeuratC57IR$seurat_clusters == '57'] <- "MNP"

Idents(aSeuratC57IR) <- "ID"

saveRDS(aSeuratC57IR, file = "seurat_Integrated.RDS")

#Figure 1B-----------------------------------------
aSeuratC57IR <- readRDS(file = "seurat_Integrated.RDS")
Idents(aSeuratC57IR) <- "ID"
localPath <- sprintf("%s/UMAP", projectPath)
tiff(file = sprintf("%s/UMAP_C57IR_ID.tiff", localPath), units = "in", 
     res = 400, height = 7, width = 10, compression = "none")  
formattedUMAP(aSeuratC57IR)
dev.off()

#Figure 1C-1G----------------------------

aSeuratC57IR <- readRDS(file = "seurat_Integrated.RDS")
localPath <- sprintf("%s/Highlight", projectPath)
dir.create(localPath)
highlightMeta(aSeuratC57IR, path = localPath, meta = "dataset", type = "TIFF")

#Create Epi Subset------------------------------

aSeuratC57IR <- readRDS(file = "seurat_Integrated.RDS")
Idents(aSeuratIntegrated) <- "seurat_clusters"

EpiSubset <- subset(aSeuratIntegrated, idents = c("3","5","7","9","12","13","14",
                                                  "15","16","18","19","23","25",
                                                  "26","28","29","31","32","33",
                                                  "35","36","37","38","39","43",
                                                  "45","46","48","50","53","55"))

seuratList<-c()

Idents(EpiSubset) <- "orig.ident"
for(origIdent in levels(EpiSubset)){
  seuratList <- c(seuratList, subset(EpiSubset, idents = origIdent))
}

#SCTransform has the same purpose as NormalizeData(), FindVariableFeatures() and ScaleData() 
seuratList <- lapply(X = seuratList, FUN = function(x) {
  x <- SCTransform(x, vst.flavor = "v2", verbose = FALSE)
})

#PrepSCTIntegration
seuratList <- PrepSCTIntegration(object.list = seuratList, anchor.features = features)
#anchors <- FindIntegrationAnchors(object.list = seuratList, anchor.features = features, normalization.method = "SCT")

#Normal method didn't work, so using control replicates as reference for integration
anchors <- FindIntegrationAnchors(object.list = seuratList, anchor.features = features, 
                                  normalization.method = "SCT", reference = c(1, 2, 3))

#Integrate the EpiSUbset using anchors derived from homeostatic reference
aSeuratC57Epi <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

DefaultAssay(aSeuratC57Epi) <- "integrated"
aSeuratC57Epi <- RunPCA(aSeuratC57Epi, npcs = 40, verbose = FALSE)
aSeuratC57Epi <- RunUMAP(aSeuratC57Epi, reduction = "pca", dims = 1:40)
aSeuratC57Epi <- FindNeighbors(aSeuratC57Epi, reduction = "pca", dims = 1:40)
aSeuratC57Epi <- FindClusters(aSeuratC57Epi, resolution = 1)

#Note the default assay will be set to 'integrated', you want this for clustering, change 
#default assay to 'SCT' before doing further data analysis such as positive gene averages
DefaultAssay(aSeuratC57Epi) <- "SCT"

#Set project.name to something meaningful
aSeuratC57Epi@project.name <- "C57 Integration Epi Subset"

#Normalize for read depth across all samples
aSeuratC57Epi <- PrepSCTFindMarkers(aSeuratC57Epi)

#Save the integrated dataset
saveRDS(aSeuratC57Epi, file = "seurat_Episubset.RDS")




#Figure 1I-------------------------------

aSeuratC57Epi <- readRDS(file = "seurat_Episubset.RDS")
aSeuratC57Epi$ID_dataset <- sprintf("%s_%s", aSeuratC57Epi$ID, aSeuratC57Epi$dataset)
Idents(aSeuratC57Epi) <- "ID_dataset"
table(aSeuratC57Epi$ID_dataset)

#Figure 3A-3B------------------------
#The DEG lists will need to be run through Metascape
#Species M. musculus, Analysis as species M. musculus
#Custom Analysis, Only select GO Biological Processes and KEGG Pathway
#Enrichment Analysis -> Analysis Report Page -> Gene List Report Excel Sheets
#Open the downloaded file and save the file in ./SMG/Metascape
#save the 'Enrichment' tab as a new .csv with the name '<timepoint> vs Homeostatic Metascape_enrichment.csv'

aSeuratC57Epi <- readRDS(file = "seurat_Episubset.RDS")

#Generate DEG for Metascape
localPath <- sprintf("%s/SMG/DEG", projectPath)
dir.create(localPath)

aSeurat <- subset(aSeuratC57Epi, ID == "SMG Acinar")
aSeurat <- SCTransform(aSeurat, return.only.var.genes = FALSE)
aSeurat <- PrepSCTFindMarkers(aSeurat)
Idents(aSeurat) <- "dataset"

for(aDataset in unique(aSeurat$dataset)){
  if(aDataset == "Homeostatic"){next} #Skip over homeostatic, DEG is relative to homeostatic
  aDFMarkers <- FindMarkers(aSeurat, ident.1 = aDataset, ident.2 = "Homeostatic", min.pct = 0.25, 
                            logfc.threshold = 0.25, only.pos = TRUE)
  aDFMarkers <- aDFMarkers[aDFMarkers$p_val_adj <= 0.05, ] #Filter for less than p_val_adj <- 0.5
  write.csv(aDFMarkers, file = sprintf("%s/%s v Homeostatic.csv", localPath, aDataset))
}

#Figure 3C-----------------------------

#List of stress genes
#https://www.sciencedirect.com/science/article/pii/S3050620425000417?via%3Dihub
# PMID: 40766395
# Unfolded Protein Response (UPR)
# DNA Damage Response (DDR)
# Oxidative Stress Response  (OSR)
# Heat Shock Response (HSR)
# Peroxisome Proliferator Activated Receptor Alpha (PPA)
# Apoptosis (APO)
# Autophagy (AUT)
# Cell Cycle Arrest (CCA)

stressList <- list(x = "")
stressList[["UPR"]] <- convert_human_to_mouse(c("ATF6", "XBP1", "HSPA5", "DDIT3", "ERN1", "EIF2AK3", "ATF4"))
stressList[["DDR"]] <- convert_human_to_mouse(c("TP53", "ATM", "CHK2", "BRCA1", "RAD51", "BRCA2", "MDM2", "CHK1", "NBN"))
stressList[["OSR"]] <- convert_human_to_mouse(c("SOD1", "GPX1", "CAT", "NFE2L2", "HMOX1", "GSR", "NQO1", "PRDX1", "TXN", "GCLC"))
stressList[["HSR"]] <- convert_human_to_mouse(c("HSPA1A", "HSP90AA1", "DNAJB1", "HSPB1", "HSPH1", "HSPA8", "DNAJC3", "CRYAB"))
stressList[["PPA"]] <- convert_human_to_mouse(c("ACAA1", "ACADM", "ACADVL", "ACOX1", "ANGPTL4", "CPT1A", "CPT2", "FABP4", "LPL", "PDK4", "PLIN2", "PPARGC1A", "CD36"))
stressList[["APO"]] <- convert_human_to_mouse(c("CASP3", "BAX", "BCL2", "CASP8", "CASP9", "FAS", "FADD", "CYCS", "BID"))
stressList[["AUT"]] <- convert_human_to_mouse(c("ATG5", "ATG7", "ATG12", "SQSTM1", "BECN1", "ULK1", "AMBRA1", "LAMP2", "BCL2"))
stressList[["CCA"]] <- convert_human_to_mouse(c("CDKN1A", "TP53", "RB1", "GADD45A", "CDK2", "CCNE1", "CDK4", "CHEK1", "CHEK2"))
stressList[["x"]] <- NULL

aSeuratC57Epi <- readRDS(file = "seurat_Episubset.RDS")
aSeuratSMG <- subset(aSeuratC57Epi, ID == "SMG Acinar")
aSeuratSMG <- PrepSCTFindMarkers(aSeuratSMG)
aSeuratSMG <- SCTransform(aSeuratSMG)
Idents(aSeuratSMG) <- "dataset"
masterDEGList <- list(tmp = 1)
for(aID in c("SMG Acinar")){
  aSeurat <- subset(aSeuratSMG, ID == aID)
  for(aDataset in c("2WeekIR", "4WeekIR", "8WeekIR", "12WeekIR")){
    id1 <- sprintf("%s", aDataset)
    id2 <- sprintf("Homeostatic")
    dataName <- sprintf("%s_%s_vs_%s", aID, id1, id2)
    aDF <- FindMarkers(aSeurat, ident.1 = id1, ident.2 = id2, min.pct = 0.25, logfc.threshold = 0.25, only.pos = TRUE)
    aDF <- aDF[aDF$p_val_adj <= 0.05, ]
    masterDEGList[[dataName]] <-  rownames(aDF)
    
    dataName <- sprintf("%s_%s_vs_%s", aID, id2, id1)
    aDF <- FindMarkers(aSeurat, ident.1 = id2, ident.2 = id1, min.pct = 0.25, logfc.threshold = 0.25, only.pos = TRUE)
    aDF <- aDF[aDF$p_val_adj <= 0.05, ]
    masterDEGList[[dataName]] <-  rownames(aDF)
  }
}
masterDEGList$tmp <- NULL

fullList <- c()
for(aProcess in names(stressList)){
  fullList <- c(fullList, stressList[[aProcess]])
}
fullList <- unique(fullList)

foundList <- c()
for(aData in names(masterDEGList)){
  aList <- fullList[fullList %in% masterDEGList[[aData]]]
  foundList <- c(foundList, aList)
}
foundList <- unique(foundList)

#foundList was manualy sorted based on type
sortedList <- c("Hspa5", "Xbp1", "Txn1", "Prdx1", "Gpx1", "Dnajc3", "Hspa8", "Cycs")

tiff(file = sprintf("%s/SMG_Stress_Markers_Founds.tiff", localPath), units = "in", 
     res = 400, height = 7, width = 10, compression = "none")
print(formattedDotPlot(aSeuratSMG, features = sortedList, textScale.pct = 150, dot.scale = 14))
dev.off()

#Figure 3D, 3F------------------------------------
aSeuratC57Epi <- readRDS(file = "seurat_Episubset.RDS")
aSeuratSMG <- subset(aSeuratC57Epi, ID == "SMG Acinar")
aSeuratSMG <- PrepSCTFindMarkers(aSeuratSMG)
aSeuratSMG <- SCTransform(aSeuratSMG)

#Generate a list of genes are are differently expressed at any timepoint compared to homeostatic
#Order by the largest fold change at any timepoint

geneList <- getAllDEGGeneNames(aSeuratSMG, aIdent2 = "Homeostatic", aMetadata = "dataset")

localPath <- sprintf("%s/SMG", projectPath)
dir.create(sprintf("%s/Plots", localPath))
Idents(aSeuratSMG) <- "dataset"
for(aDataset in unique(aSeuratSMG$dataset)){
  if(aDataset == "Homeostatic") {next}
  
  #Read the metascape .csv, for all the genes in the 'Symbols' cell of the 'Summary' row, generate a dotplot ordered by highest fold change at any timepoint
  DotPlotFromEnrichment(aSeuratSMG, 
                        FileName = sprintf("%s/Metascape/%s vs Homeostatic metascape_enrichment.csv", localPath, aDataset),
                        outputPath = sprintf("%s/Plots", localPath), 
                        outputPrefix = sprintf("Dotplot %s vs Homeostatic Filtered", aDataset),
                        title = sprintf("%s vs Homeostatic", aDataset),
                        filterGenes = geneList)
  #Same as above, but only the top 10 genes for better readability
  DotPlotFromEnrichment(aSeuratSMG, 
                        FileName = sprintf("%s/Metascape/%s vs Homeostatic metascape_enrichment.csv", localPath, aDataset),
                        outputPath = sprintf("%s/Plots", localPath), 
                        outputPrefix = sprintf("Dotplot %s vs Homeostatic Filtered Top10", aDataset),
                        title = sprintf("%s vs Homeostatic", aDataset),
                        filterGenes = geneList,
                        topCount = 10)
}

#Gene lists copied from dotplots created above
localPath <- sprintf("%s/SMG/Plots", projectPath)

translation_2Week <- c("Lars2", "Eif4g3", "Cald1", "Eif3j1", "Rpl23a", "Eif5b", "Jun", "mt-Rnr1", "Rpl13a", "Rps29")
translation_12Week <- c("Lars2", "Cald1", "Eif5b", "mt-Rnr1", "Rpl13a", "Rps29", "Pik3r1", "Rpl38", "Sorbs2", "Rps28")
tiff(file = sprintf("%s/DotPlot_SMG_Translation.tiff", localPath), units = "in", 
     res = 400, height = 7, width = 10, compression = "none")  
print(formattedDotPlot(aSeuratSMG, features = unique(c(translation_2Week, translation_12Week)), textScale.pct = 150, dot.scale = 12))
dev.off()

oxPhos_2Week <- c("Gphn", "Ephx1", "Hexb", "Dnah11", "ND3", "Jun", "ND1", "Atp5me", "Pde4b", "Tpr")
oxPhos_12Week <- c("Gphn", "Hexb", "Dnah11", "Atp5me", "Pde4b", "Tubb4b", "Gne", "Polr2l", "Pik3r1", "Cox6a1")
tiff(file = sprintf("%s/DotPlot_SMG_OxPhos.tiff", localPath), units = "in", 
     res = 400, height = 7, width = 10, compression = "none")  
print(formattedDotPlot(aSeuratSMG, features = unique(c(oxPhos_2Week, oxPhos_12Week)), textScale.pct = 150, dot.scale = 12))
dev.off()


#Figure 3E--------------------------------------------
#Data used to generate graph

aSeuratC57Epi <- readRDS(file = "seurat_Episubset.RDS")
aSeuratSMG <- subset(aSeuratC57Epi, ID == "SMG Acinar")
aSeuratSMG <- PrepSCTFindMarkers(aSeuratSMG)
aSeuratSMG <- SCTransform(aSeuratSMG)

#Find All Ribosome genes
list40S <- c("Rpsa", "Rps2", "Rps3", "Rps3a", "Rps4x", "Rps4y1", "Rps4y2", "Rps5", "Rps6", "Rps7", "Rps8", 
             "Rps9", "Rps10", "Rps11", "Rps12", "Rps13", "Rps14", "Rps15", "Rps15a", "Rps16", "Rps17", "Rps18", 
             "Rps19", "Rps20", "Rps21", "Rps23", "Rps24", "Rps25", "Rps26", "Rps27", "Rps27a", "Rps28", "Rps29", 
             "Rps30", "Rack1")
list60S <- c("Rpl3", "Rpl3l", "Rpl4", "Rpl5", "Rpl6", "Rpl7", "Rpl7a", "Rpl8", "Rpl9", "Rpl10", "Rpl10a", "Rpl11", 
             "Rpl12", "Rpl13", "Rpl13a", "Rpl14", "Rpl15", "Rpl17", "Rpl18", "Rpl18a", "Rpl19", "Rpl21", "Rpl22", 
             "Rpl23", "Rpl23a", "Rpl24", "Rpl26", "Rpl27", "Rpl27a", "Rpl28", "Rpl29", "Rpl30", "Rpl31", "Rpl32", 
             "Rpl34", "Rpl35", "Rpl35a", "Rpl36", "Rpl36a", "Rpl36al", "Rpl37", "Rpl37a", "Rpl38", "Rpl39", "Rpl40", 
             "Rpl41", "Rplp0", "Rplp1", "Rplp2")
listAllS <- unique(c(list40S, list60S))

for(aTimepoint in c("2WeekIR", "4WeekIR", "8WeekIR", "12WeekIR")){
  id1 <- aTimepoint
  id2 <- "Homeostatic"
  dataName <- sprintf("%s_vs_%s", id1, id2)
  aDF <- FindMarkers(aSeuratSMG, ident.1 = id1, ident.2 = id2, min.pct = 0, logfc.threshold = 0, only.pos = FALSE)
  write.csv(aDF, file = sprintf("%s/%s_All_Genes.csv", localPath, dataName))
  aDFFiltered <- aDF[listAllS, ]
  write.csv(aDFFiltered, file = sprintf("%s/%s_AllRibosomeGenes.csv", localPath, dataName))
  
  dataName <- sprintf("%s_vs_%s", id2, id1)
  aDF <- FindMarkers(aSeuratSMG, ident.1 = id2, ident.2 = id1, min.pct = 0, logfc.threshold = 0, only.pos = FALSE)
  write.csv(aDF, file = sprintf("%s/%s_All_Genes.csv", localPath, dataName))
  aDFFiltered <- aDF[listAllS, ]
  write.csv(aDFFiltered, file = sprintf("%s/%s_AllRibosomeGenes.csv", localPath, dataName))
}
#Create SMG Acinar subset---------------------------
#The Acinar subset was created in steps, first step was acinar + dividing cells from epi subset

aSeuratC57Epi <- readRDS(file = "seurat_Episubset.RDS")
aSeuratEpi@misc[[1]] <- "EpiSubset"

aSeuratSMGAcinarDividing <- subset(aSeuratEpi, ID %in% c("SMG Acinar", "Dividing"))
aSeuratSMGAcinarDividing@misc[[1]] <- "SMGAcinarDividingSubset"
seuratList <- c(aSeuratSMGAcinarDividing)

seuratOIList <- c()
for(i in 1:length(seuratList)){
  Idents(seuratList[[i]]) <- "orig.ident"
  for(origIdent in levels(seuratList[[i]])){
    seuratOIList <- c(seuratOIList, subset(seuratList[[i]], idents = origIdent))
  }
}

seuratOIList <- lapply(X = seuratOIList, FUN = function(x) {
  x <- SCTransform(x, vst.flavor = "v2", verbose = FALSE)
})

features <- SelectIntegrationFeatures(object.list = seuratOIList, nfeatures = 3000)
seuratOIList <- PrepSCTIntegration(object.list = seuratOIList, anchor.features = features)
anchors <- FindIntegrationAnchors(object.list = seuratOIList, anchor.features = features, normalization.method = "SCT",
                                  reference = c(1, 2, 3))
aSeuratAcinar <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

DefaultAssay(aSeuratAcinar) <- "integrated"
aSeuratAcinar <- SCTransform(aSeuratAcinar, vst.flavor = "v2", verbose = FALSE)
aSeuratAcinar <- RunPCA(aSeuratAcinar, npcs = 40, verbose = FALSE)
aSeuratAcinar <- RunUMAP(aSeuratAcinar, reduction = "pca", dims = 1:40)
aSeuratAcinar <- FindNeighbors(aSeuratAcinar, reduction = "pca", dims = 1:40)
aSeuratAcinar <- FindClusters(aSeuratAcinar, resolution = 1)
DefaultAssay(aSeuratAcinar) <- "SCT"
aSeuratAcinar <- PrepSCTFindMarkers(aSeuratAcinar)

#Dcpp3 is a SL gene, cells high in Dcpp3 are likely SL cells that co-clustered with SMG Acinar cells
aGeneList <- c("Dcpp3")
aSeuratAcinar$Dcpp <- 0
for(aGene in aGeneList){
  print(aGene)
  
  aSeuratAcinar$Dcpp <- ifelse((aSeuratAcinar@assays$SCT@counts[aGene, ] >= 20),
                               1,
                               aSeuratAcinar$Dcpp)
}
aSeuratAcinar <- subset(aSeuratAcinar, Dcpp == "0")

aSeuratAcinar <- seuratIntegration(aSeuratAcinar, referenceOI = c(1, 2, 3))

aSeuratAcinar <-  FindClusters(aSeuratAcinar, resolution = 0.12)
aSeuratAcinar <- subset(aSeuratAcinar, seurat_clusters == "5", invert = TRUE)
aSeuratAcinar <- SCTransform(aSeuratAcinar)
aSeuratAcinar <- PrepSCTFindMarkers(aSeuratAcinar)

aSeuratAcinar$seurat_clusters2[aSeuratAcinar$seurat_clusters == '0'] <- "0"
aSeuratAcinar$seurat_clusters2[aSeuratAcinar$seurat_clusters == '1'] <- "1"
aSeuratAcinar$seurat_clusters2[aSeuratAcinar$seurat_clusters == '2'] <- "2"
aSeuratAcinar$seurat_clusters2[aSeuratAcinar$seurat_clusters == '3'] <- "3"
aSeuratAcinar$seurat_clusters2[aSeuratAcinar$seurat_clusters == '4'] <- "4"
aSeuratAcinar$seurat_clusters2[aSeuratAcinar$seurat_clusters == '6'] <- "5"

aSeuratAcinar$Subpop_Metascape <- "Missing"
aSeuratAcinar$Subpop_Metascape[aSeuratAcinar$seurat_clusters2 == '0'] <- "Active"
aSeuratAcinar$Subpop_Metascape[aSeuratAcinar$seurat_clusters2 == '1'] <- "Baseline"
aSeuratAcinar$Subpop_Metascape[aSeuratAcinar$seurat_clusters2 == '2'] <- "Stressed"
aSeuratAcinar$Subpop_Metascape[aSeuratAcinar$seurat_clusters2 == '3'] <- "High Metabolic_1"
aSeuratAcinar$Subpop_Metascape[aSeuratAcinar$seurat_clusters2 == '4'] <- "High Metabolic_2"
aSeuratAcinar$Subpop_Metascape[aSeuratAcinar$seurat_clusters2 == '5'] <- "Proliferating"

saveRDS(aSeuratAcinar, file = "seurat_Acinar.RDS")

#Figure 4A-----------------------------------
aSeuratAcinar <- readRDS(file = "seurat_Acinar.RDS")
Idents(aSeuratAcinar) <- "seurat_clusters2"
localPath <- sprintf("%s/SMG_Subcluster", projectPath)
dir.create(localPath)
tiff(file = sprintf("%s/UMAP_SMG_Subcluster_seurat_clusters.tiff", localPath), units = "in", 
     res = 400, height = 7, width = 10, compression = "none")  
print(formattedUMAP(aSeuratAcinar))
dev.off()

#Figure 4B-----------------------------------
#Create the data used to create the graph in figure 4B
aSeuratAcinar <- readRDS(file = "seurat_Acinar.RDS")
aDF <- data.frame(tmp = 0)
for(aDataset in unique(aSeuratAcinar$dataset)){
  nTotalCells <- length(WhichCells(aSeuratAcinar, expression = dataset == aDataset))
  for(aCluster in unique(aSeuratAcinar$seurat_clusters2)){
    nCells <- length(WhichCells(aSeuratAcinar, expression = dataset == aDataset & seurat_clusters2 == aCluster))
    aDF[aDataset, aCluster] <- nCells/nTotalCells
  }
}
aDF$tmp <- NULL
aDF <- aDF[-1,]
write.csv(aDF, file = sprintf("%s/SMG_Subcluster/clusterCellCounts.csv", projectPath))

#Figure 4C-------------------------------
aSeuratAcinar <- readRDS(file = "seurat_Acinar.RDS")
Idents(aSeuratAcinar) <- "Subpop_Metascape"
localPath <- sprintf("%s/SMG_Subcluster", projectPath)
dir.create(localPath)
tiff(file = sprintf("%s/UMAP_SMG_Subcluster_subpop.tiff", localPath), units = "in", 
     res = 400, height = 7, width = 10, compression = "none")  
print(formattedUMAP(aSeuratAcinar))
dev.off()

#Figure 4D-4I and Appendix------------------------------------------------
aSeuratAcinar <- readRDS(file = "seurat_Acinar.RDS")

geneList <- getAllDEGGeneNames(aSeuratAcinar, aIdent2 = "Homeostatic", aMetadata = "dataset")

localPath <- sprintf("%s/SMG_Subcluster/DEG/Metascape_Custom", projectPath)
Idents(aSeuratAcinar) <- "dataset"
for(aCluster in unique(aSeuratAcinar$seurat_clusters2)){
  DotPlotFromEnrichment(aSeuratAcinar, 
                        FileName = sprintf("%s/Cluster %s metascape_enrichment.csv", localPath, aCluster),
                        outputPath = sprintf("%s/Plots/Dataset", localPath), 
                        outputPrefix = sprintf("Dotplot Cluster %s by dataset Filtered", aCluster),
                        title = sprintf("Cluster %s", aCluster),
                        filterGenes = geneList)
  
  DotPlotFromEnrichment(aSeuratAcinar, 
                        FileName = sprintf("%s/Cluster %s metascape_enrichment.csv", localPath, aCluster),
                        outputPath = sprintf("%s/Plots/Dataset", localPath), 
                        outputPrefix = sprintf("Dotplot Cluster %s by dataset Filtered Top10", aCluster),
                        title = sprintf("Cluster %s", aCluster),
                        filterGenes = geneList,
                        topCount = 10)
}
#Gene lists copied from dotplots created above

Idents(aSeuratAcinar) <- "dataset"
cluster_0_processing <- c("Rrbp1", "Sec62", "Dnajc1", "Canx", "P4hb", "Nucb2", "Sec63", "Nktr", "Vcp", "Eif2ak3")
localPath <- sprintf("%s/SMG_Subcluster/DEG/Metascape_Custom/Plots/Dataset", projectPath)
tiff(file = sprintf("%s/Dotplot_SMG_Subcluster_Cluster0_processing.tiff", localPath), units = "in", 
     res = 400, height = 7, width = 10, compression = "none")  
formattedDotPlot(aSeuratAcinar, features = cluster_0_processing, textScale.pct = 200, dot.scale = 15) 
dev.off()

Idents(aSeuratAcinar) <- "dataset"
cluster_1_mitochondrial <- c("Rps26", "Mrps21")
localPath <- sprintf("%s/SMG_Subcluster/DEG/Metascape_Custom/Plots/Dataset", projectPath)
tiff(file = sprintf("%s/Dotplot_SMG_Subcluster_Cluster1_mitoochondrial.tiff", localPath), units = "in", 
     res = 400, height = 7, width = 10, compression = "none")  
formattedDotPlot(aSeuratAcinar, features = cluster_1_mitochondrial, textScale.pct = 200, dot.scale = 15) 
dev.off()

Idents(aSeuratAcinar) <- "dataset"
cluster_2_Stress <- c("Trp53inp1", "Ephx1", "Jun", "Txnip", "Ralbp1", "Fos", "Dusp1", "Mt1", "Errfi1", "Cd2ap")
localPath <- sprintf("%s/SMG_Subcluster/DEG/Metascape_Custom/Plots/Dataset", projectPath)
tiff(file = sprintf("%s/Dotplot_SMG_Subcluster_Cluster2_stress.tiff", localPath), units = "in", 
     res = 400, height = 7, width = 10, compression = "none")  
formattedDotPlot(aSeuratAcinar, features = cluster_2_Stress, textScale.pct = 200, dot.scale = 15) 
dev.off()

Idents(aSeuratAcinar) <- "dataset"
cluster_3_translation <- c("Rpl23a", "Rps27", "Rpl31", "Mrpl52", "Rps25", "Rpl12", "Myl6", "Mrps21", "Rps20", "Rps6")
localPath <- sprintf("%s/SMG_Subcluster/DEG/Metascape_Custom/Plots/Dataset", projectPath)
tiff(file = sprintf("%s/Dotplot_SMG_Subcluster_Cluster3_translation.tiff", localPath), units = "in", 
     res = 400, height = 7, width = 10, compression = "none")  
formattedDotPlot(aSeuratAcinar, features = cluster_3_translation, textScale.pct = 200, dot.scale = 15) 
dev.off()

Idents(aSeuratAcinar) <- "dataset"
cluster_4_ribosome <- c("Rpl23a", "Rps29", "Rps28", "Rps21", "Rps12", "Rps27", "Rpl31", "Rpl41", "Mrpl52", "Rps8")
localPath <- sprintf("%s/SMG_Subcluster/DEG/Metascape_Custom/Plots/Dataset", projectPath)
tiff(file = sprintf("%s/Dotplot_SMG_Subcluster_Cluster4_ribosome.tiff", localPath), units = "in", 
     res = 400, height = 7, width = 10, compression = "none")  
formattedDotPlot(aSeuratAcinar, features = cluster_4_ribosome, textScale.pct = 200, dot.scale = 14) 
dev.off()

Idents(aSeuratAcinar) <- "dataset"
cluster_5_division <- c("Ccng1", "Pafah1b1", "Tubb4b", "Hmgb1", "Lmna", "Hsp90ab1", "Ncor1", "Hnrnpa2b1", "Hnrnpu", "Anapc13")
localPath <- sprintf("%s/SMG_Subcluster/DEG/Metascape_Custom/Plots/Dataset", projectPath)
tiff(file = sprintf("%s/Dotplot_SMG_Subcluster_Cluster5_division.tiff", localPath), units = "in", 
     res = 400, height = 7, width = 10, compression = "none")  
formattedDotPlot(aSeuratAcinar, features = cluster_5_division, textScale.pct = 200, dot.scale = 14) 
dev.off()
