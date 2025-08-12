.libPaths(c("/xtdisk/jiangl_group/wangqifei/hu/envTool/Rlibs"))
####### Requirements #######
packages <- c("data.table", "tidyverse", "Seurat", "DoubletFinder")
# Import packages
invisible(lapply(packages, function(p){
suppressPackageStartupMessages(library(p, character.only=T, quietly=T))
}))

#### Arg Set ####
library(argparse)
parser <- ArgumentParser(description='Do cell filter befor scTE')
parser$add_argument('-i', '--input', type="character", help='input matrix path for seurat read')
parser$add_argument('-s', '--sample', type="character", help='sample name ident')
parser$add_argument('-b', '--bc', type="integer", default=1, help='what is bc')
# parser$add_argument('-p', '--platform', type="character", help='sample platform [10x-v1, 10x-v2, 10x-v3]')
parser$add_argument('-o', '--output', type="character", help='output dir path')
parser$add_argument('-t', '--type', type="character", default="STAR", help='input data construction type')
parser$add_argument('-c', '--core', type="integer", default=5, help='cpu core max counts')
parser$add_argument('-gs', '--gsize', type="integer", default=50, help='max future core global memory size(Gb)')
parser$add_argument('-e', '--env', type="character", default=NULL, help='use conda env to transfer Seurat > h5ad')
args <- parser$parse_args()   # 获取参数

#### Parallelization in Seurat ####
library(future)
plan(multicore, workers = args$core)
# Set total 50G RAM
options(future.globals.maxSize = args$gsize*1024^3)


####### Functions #######
# 1.mat2seuob: Create a seurat object from raw csv and Pre-QC
# csv: csv realpath
# bc: barcode list (optional)
csv2seuob <- function(csv, bc=NA) {
   ## Check requirements
   pkg <- c("data.table", "tidyverse", "Seurat")
   load_pkg <- sapply(pkg, require, character.only=T, quietly=T)
   if(!all(load_pkg)){
      return(paste("Error: No packages ", names(load_pkg[!load_pkg]), sep=","))
   }

   EC <- fread(csv, data.table = FALSE, header = TRUE) %>% 
      column_to_rownames("V1")
   EC <- t(EC)
   ## Pre-QC
   EC <- CreateSeuratObject(EC, min.cells = 3, min.features = 200)
   EC[["percent.mt"]] <- PercentageFeatureSet(EC, pattern = "^MT-")
   EC <- subset(EC, subset = nFeature_RNA > 500 & nCount_RNA > 1000 & percent.mt < 10)

   if(!is.na(bc)) {
      cell.call <- readLines(bc)
      EC <- subset(EC, cells = intersect(cell.call, colnames(EC)))
   }
   return(EC)
}


# 2.find_doub: Find doublets by DoubletFinder pipeline
# seuob: processed seurat object
Find_doublet <- function(seuob, dim.usage=30) {
   ## Check requirements
   pkg <- c("Seurat", "DoubletFinder")
   load_pkg <- sapply(pkg, require, character.only=T, quietly=T)
   if(!all(load_pkg)){
      return(paste("Error: No packages ", names(load_pkg[!load_pkg]), sep=","))
   }
   if (packageVersion("DoubletFinder")>2) {
      paramSweep_v3 = paramSweep
      doubletFinder_v3 = doubletFinder
   }
   sweep.res.list <- paramSweep_v3(seuob, 
      PCs=1:dim.usage, sct=F, num.cores=args$core)
   sweep.stats <- summarizeSweep(sweep.res.list, GT=F)
   bcmvn <- find.pK(sweep.stats)
   nExp_poi <- round(0.05 * ncol(seuob))
   p <- as.numeric(as.vector(bcmvn[bcmvn$MeanBC == max(bcmvn$MeanBC),]$pK))
   seuob <- doubletFinder_v3(seuob, PCs=1:dim.usage, pN=0.25, pK=p, nExp=nExp_poi, reuse.pANN=F, sct=F)
   colnames(seuob@meta.data)[ncol(seuob@meta.data)] <- "doublet.info"
   return(seuob)
}


# 3.rm_doub: Remove doublets, include csv2seuob(), find_doub()
# nsam: sample name
# inpath: csv realpath
rm_doub <- function(nsam, inpath, bc=NA, dim.usage=30, type="csv") {
   ## Check requirements
   pkg <- c("Seurat")
   load_pkg <- sapply(pkg, require, character.only=T, quietly=T)
   if(!all(load_pkg)){
      return(paste("Error: No packages ", names(load_pkg[!load_pkg]), sep=","))
   }

   ## Create a seurat object
   if (type == "csv") {
      EC <- csv2seuob(inpath, bc)
   } else if (type == "10x") {
      EC <- Read10X(inpath)
      EC <- CreateSeuratObject(EC, min.cells = 3, min.features = 200)
      EC[["percent.mt"]] <- PercentageFeatureSet(EC, pattern = "^MT-")
      EC <- subset(EC, subset = nFeature_RNA > 500 & nCount_RNA > 1000 & percent.mt < 10)
   } else if (type == "STAR") {
      EC <- ReadSTARsolo(inpath)
      print("Read:")
      print(dim(EC))
      EC <- CreateSeuratObject(EC, min.cells = 3, min.features = 200)
      print("Create:")
      print(dim(EC))
      print(str(EC@meta.data))
      EC[["percent.mt"]] <- PercentageFeatureSet(EC, pattern = "^MT-")
      EC <- subset(EC, subset = nFeature_RNA > 500 & nCount_RNA > 1000 & percent.mt < 10)
   }


   ## Process
   EC <- NormalizeData(EC)
   EC <- FindVariableFeatures(EC, selection.method="vst", nfeatures=3000)
   EC <- ScaleData(EC)
   EC <- RunPCA(EC)
   EC <- RunUMAP(EC, dims=1:dim.usage)

   ## Remove doublets
   EC <- Find_doublet(EC)
   EC <- subset(EC, doublet.info == "Singlet")

   ## Remove columns "^pANN_" in meta.data
   c <- grep("pANN_", colnames(EC@meta.data))
   EC@meta.data <- EC@meta.data[, -c]
   EC@meta.data$orig.ident <- nsam

   ## Output ncell of the sample
   print(paste0(nsam," cells: ", length(EC@meta.data$orig.ident), " is read!", 
         " Time:", format(Sys.time(), "%Y%m%d %X"), sep = " "))
   return(EC)
}


####### Main #######
# nsam <- commandArgs(trailingOnly = TRUE)[1]
# inpath <- commandArgs(trailingOnly = TRUE)[2]
# outpath <- commandArgs(trailingOnly = TRUE)[3]
# bc.path <- commandArgs(trailingOnly = TRUE)[4]  # optional

obj <- rm_doub(args$sample, args$input, bc=args$bc, type="STAR")
bc <- colnames(obj)
# Save
writeLines(bc, con=paste0(args$output, args$sample, ".cells.csv"))
saveRDS(obj, paste0(args$output, args$sample, ".rds"))
# rds -> h5ad
library(sceasy)
library(reticulate)
use_condaenv(args$env)

# V5 -> V4
if (class(obj$RNA)=='Assay5') {
   print("Transfer V5 to V4")
   obj[["RNA"]] <- as(obj[["RNA"]], "Assay")
}

sceasy::convertFormat(obj, from="seurat", to="anndata",
                       outFile=paste0(args$output, args$sample, ".h5ad"))

