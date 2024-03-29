library(ggpubr)
library(Seurat)
library(dplyr)
library(patchwork)
library(ggplot2)
library(cowplot)
library(stringr)
library(ensembldb)
library(AnnotationHub)
library(dplyr)

getDistribution = function( SO, genes, out_dir, COND ){
  
  ## get the raw read counts for genes in the list 'genes'
  # SO=hct_auxi_T120_filt; genes=all_useful_genes;i=1; COND='DMSO'
  ## out_dir=''
  ind = match( genes, rownames(SO) )
  ind = ind[!is.na(ind)]
  rrc = SO@assays$RNA@counts
  rrc = rrc[ind,]
 
  for( i in 1:nrow( rrc )) {
    x = rrc[i,]
    gene=rownames(rrc)[i]
    tp = data.frame( hist(x, seq( -1, max(x) ), plot = FALSE)$counts,
                     seq(0, max(x) ) )
    write.table( tp, file=paste0(out_dir,gene,'_',COND,'.txt' ), 
                 quote=FALSE, row.names=FALSE, col.names=FALSE,sep='\t' )
  }}

getDistribution2 = function( SO, genes, out_dir, COND ){
  
  ## get the raw read counts for genes in the list 'genes'
  # SO=hct_auxi_T120_filt; genes=all_useful_genes;i=1; COND='DMSO'
  ## out_dir=''
  ind = match( genes, rownames(SO) )
  ind = ind[!is.na(ind)]
  rrc = SO@assays$RNA@counts
  rrc = rrc[ind,]
  rrc = rrc + 5
 
  for( i in 1:nrow( rrc )) {
    x = rrc[i,]
    gene=rownames(rrc)[i]
    tp = data.frame( hist(x, seq( -1, max(x) ), plot = FALSE)$counts,
                     seq(0, max(x) ) )
    write.table( tp, file=paste0(out_dir,gene,'_',COND,'.txt' ), 
                 quote=FALSE, row.names=FALSE, col.names=FALSE,sep='\t' )
  }}

options(future.globals.maxSize = 4000 * 1024^2)
load("cell_cycle_genes.RData")

### CH12 WT samples process ###
T1.data <- Read10X(data.dir="WT_r1/filtered_feature_bc_matrix/")
data = CreateSeuratObject(counts = T1.data, min.cells = 3, min.features = 200, project = "WTr1")
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern="^mt-")
Wdata1<- data

T1.data <- Read10X(data.dir="WT_r2/filtered_feature_bc_matrix/")
data = CreateSeuratObject(counts = T1.data, min.cells = 3, min.features = 200, project = "WTr2")
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern="^mt-")
Wdata2 <- data

Wdata1_filt <- subset(Wdata1, subset = nFeature_RNA > 4000 & nFeature_RNA < 7500 & percent.mt < 5 & nCount_RNA > 2000 & nCount_RNA < 70000)
Wdata2_filt <- subset(Wdata2, subset = nFeature_RNA > 2000 & nFeature_RNA < 4500 & percent.mt < 5 & nCount_RNA > 2000 & nCount_RNA < 20000)

dataBM <- merge( x= Wdata1_filt, y=Wdata2_filt, add.cell.id = c("sWTr1","sWTr2"))

dataBM$log10GenesPerUMI <- log10(dataBM$nFeature_RNA) / log10(dataBM$nCount_RNA)
metadata <-dataBM@meta.data
metadata$cells <- rownames(metadata)
# Rename columns
metadata <- metadata %>%
        dplyr::rename(seq_folder = orig.ident,
                      nUMI = nCount_RNA,
                      nGene = nFeature_RNA)
metadata$sample <- NA
metadata$sample[which(str_detect(metadata$cells, "^sWTr1_"))] <- "sWTr1"
metadata$sample[which(str_detect(metadata$cells, "^sWTr2_"))] <- "sWTr2"

dataBM$sample <- metadata$sample

counts <- GetAssayData(object = dataBM, slot = "counts")
nonzero <- counts > 0
keep_genes <- Matrix::rowSums(nonzero) >= 10
filtered_counts <- counts[keep_genes, ]


dataBM <- CreateSeuratObject(filtered_counts, meta.data = dataBM@meta.data)


seurat_phase <- NormalizeData(dataBM)
seurat_phase <- CellCycleScoring(seurat_phase, 
                                 g2m.features = g2m_genes, 
                                 s.features = s_genes)
seurat_phase <- FindVariableFeatures(seurat_phase, 
                     selection.method = "vst",
                     nfeatures = 2000, 
                     verbose = FALSE)
		     
# Scale the counts
seurat_phase <- ScaleData(seurat_phase)
seurat_phase <- RunPCA(seurat_phase)

# make histogram  
data_filt <- subset(seurat_phase, sample=="sWTr1")
write.table(data_filt@assays$RNA@counts, file="CH12_sWT_r1.txt", quote=F, sep="\t")

genes = rownames(data_filt)
length(genes)
getDistribution(data_filt, genes=genes, out_dir="sWT_r1_Sfiltered/", COND='sWTr1')

data_filt <- subset(seurat_phase, sample=="sWTr2")
write.table(data_filt@assays$RNA@counts, file="CH12_sWT_r2.txt", quote=F, sep="\t")
genes = rownames(data_filt)
length(genes)
getDistribution(data_filt, genes=genes, out_dir="sWT_r2_Sfiltered/", COND='sWTr2')

### ESC samples process ###
## CTCF degron
T1.data <- Read10X(data.dir="ESC_CDMSO/filtered_feature_bc_matrix/")
data = CreateSeuratObject(counts = T1.data, min.cells = 3, min.features = 200, project = "CDMSO")
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern="^mt-")
Ddata1<- data

T1.data <- Read10X(data.dir="ESC_CAUXIN/filtered_feature_bc_matrix/")
data = CreateSeuratObject(counts = T1.data, min.cells = 3, min.features = 200, project = "CAUXIN")
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern="^mt-")
Adata1 <- data

Ddata1_filt <- subset(Ddata1, subset = nFeature_RNA > 3000 & nFeature_RNA < 6500 & percent.mt < 5 & nCount_RNA > 10000 & nCount_RNA < 35000)
Adata1_filt <- subset(Adata1, subset = nFeature_RNA > 3000 & nFeature_RNA < 6500 & percent.mt < 5 & nCount_RNA > 10000 & nCount_RNA < 35000)

## Rad21 degron
T1.data <- Read10X(data.dir="ESC_RDMSO/filtered_feature_bc_matrix/")
data = CreateSeuratObject(counts = T1.data, min.cells = 3, min.features = 200, project = "RDMSO")
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern="^mt-")
Ddata2<- data

T1.data <- Read10X(data.dir="ESC_RAUXIN/filtered_feature_bc_matrix/")
data = CreateSeuratObject(counts = T1.data, min.cells = 3, min.features = 200, project = "RAUXIN")
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern="^mt-")
Adata2 <- data

Ddata2_filt <- subset(Ddata2, subset = nFeature_RNA > 4000 & nFeature_RNA < 7500 & percent.mt < 5 & nCount_RNA > 10000 & nCount_RNA < 50000)
Adata2_filt <- subset(Adata2, subset = nFeature_RNA > 4000 & nFeature_RNA < 7500 & percent.mt < 5 & nCount_RNA > 10000 & nCount_RNA < 50000)

dataB <- merge( x= Ddata1_filt, y=c(Adata1_filt,Ddata2_filt,Adata2_filt), add.cell.id = c("CDMSO","CAUXIN","RDMSO","RAUXIN"))

dataB$log10GenesPerUMI <- log10(dataB$nFeature_RNA) / log10(dataB$nCount_RNA)
metadata <-dataB@meta.data
metadata$cells <- rownames(metadata)
# Rename columns
metadata <- metadata %>%
        dplyr::rename(seq_folder = orig.ident,
                      nUMI = nCount_RNA,
                      nGene = nFeature_RNA)
metadata$sample <- NA
metadata$sample[which(str_detect(metadata$cells, "^CDMSO_"))] <- "CDMSO"
metadata$sample[which(str_detect(metadata$cells, "^CAUXIN_"))] <- "CAUXIN"
metadata$sample[which(str_detect(metadata$cells, "^RDMSO_"))] <- "RDMSO"
metadata$sample[which(str_detect(metadata$cells, "^RAUXIN_"))] <- "RAUXIN"

dataB$sample <- metadata$sample

counts <- GetAssayData(object = dataB, slot = "counts")
nonzero <- counts > 0
keep_genes <- Matrix::rowSums(nonzero) >= 10
filtered_counts <- counts[keep_genes, ]


dataB <- CreateSeuratObject(filtered_counts, meta.data = dataB@meta.data)


seurat_phase <- NormalizeData(dataB)
seurat_phase <- CellCycleScoring(seurat_phase, 
                                 g2m.features = g2m_genes, 
                                 s.features = s_genes)
seurat_phase <- FindVariableFeatures(seurat_phase, 
                     selection.method = "vst",
                     nfeatures = 2000, 
                     verbose = FALSE)
		     
# Scale the counts
seurat_phase <- ScaleData(seurat_phase)
seurat_phase <- RunPCA(seurat_phase)

# make histogram for CDMSO samples 
data_filt <- subset(seurat_phase, sample=="CDMSO")
sum(data_filt[["RNA"]]@counts)

genes = rownames(data_filt)
length(genes)
getDistribution(data_filt, genes=genes, out_dir="ESC_CDMSO_Sfiltered/", COND='ESCCDMSO')
write.table(data_filt@assays$RNA@counts, file="ESC_CDMSO.txt", quote=F, sep="\t")

data_filt <- subset(seurat_phase, sample=="CAUXIN")
sum(data_filt[["RNA"]]@counts)

genes = rownames(data_filt)
length(genes)
getDistribution(data_filt, genes=genes, out_dir="ESC_CAUXIN_Sfiltered/", COND='ESCCAUXIN')
write.table(data_filt@assays$RNA@counts, file="ESC_CAUXIN.txt", quote=F, sep="\t")

# make histogram for RDMSO samples 
data_filt <- subset(seurat_phase, sample=="RDMSO")
sum(data_filt[["RNA"]]@counts)

genes = rownames(data_filt)
length(genes)
getDistribution(data_filt, genes=genes, out_dir="ESC_RDMSO_Sfiltered/", COND='ESCRDMSO')
write.table(data_filt@assays$RNA@counts, file="ESC_RDMSO.txt", quote=F, sep="\t")

data_filt <- subset(seurat_phase, sample=="RAUXIN")
sum(data_filt[["RNA"]]@counts)

genes = rownames(data_filt)
length(genes)
getDistribution(data_filt, genes=genes, out_dir="ESC_RAUXIN_Sfiltered/", COND='ESCRAUXIN')
write.table(data_filt@assays$RNA@counts, file="ESC_RAUXIN.txt", quote=F, sep="\t")

### HCT cohesion degron time course ###
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

T1.data <- Read10X(data.dir="HCT_DMSO0h/filtered_gene_bc_matrices/")
data = CreateSeuratObject(counts = T1.data, min.cells = 3, min.features = 200, project = "DMSO0h")
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern="^MT-")
Ddata0h<- data

T1.data <- Read10X(data.dir="HCT_AUXIN0h/filtered_gene_bc_matrices/")
data = CreateSeuratObject(counts = T1.data, min.cells = 3, min.features = 200, project = "AUXIN0h")
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern="^MT-")
Adata0h<- data

Adata0h_filt <- subset(Adata0h, subset = nFeature_RNA > 1000 & nFeature_RNA < 5000 & percent.mt < 7.5 & nCount_RNA > 2000 & nCount_RNA < 35000)
Ddata0h_filt <- subset(Ddata0h, subset = nFeature_RNA > 1000 & nFeature_RNA < 5000 & percent.mt < 7.5 & nCount_RNA > 2000 & nCount_RNA < 35000)

T1.data <- Read10X(data.dir="HCT_DMSO30m/filtered_gene_bc_matrices/")
data = CreateSeuratObject(counts = T1.data, min.cells = 3, min.features = 200, project = "DMSO30m")
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern="^MT-")
Ddata30m<- data

T1.data <- Read10X(data.dir="HCT_AUXIN30m/filtered_gene_bc_matrices/")
data = CreateSeuratObject(counts = T1.data, min.cells = 3, min.features = 200, project = "AUXIN30m")
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern="^MT-")
Adata30m<- data

Adata30m_filt <- subset(Adata30m, subset = nFeature_RNA > 2000 & nFeature_RNA < 6000 & percent.mt < 7.5 & nCount_RNA > 2000 & nCount_RNA < 45000)
Ddata30m_filt <- subset(Ddata30m, subset = nFeature_RNA > 2000 & nFeature_RNA < 6000 & percent.mt < 7.5 & nCount_RNA > 2000 & nCount_RNA < 45000)

T1.data <- Read10X(data.dir="HCT_DMSO2h/filtered_gene_bc_matrices/")
data = CreateSeuratObject(counts = T1.data, min.cells = 3, min.features = 200, project = "DMSO2h")
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern="^MT-")
Ddata2h<- data

T1.data <- Read10X(data.dir="HCT_AUXIN2h/filtered_gene_bc_matrices/")
data = CreateSeuratObject(counts = T1.data, min.cells = 3, min.features = 200, project = "AUXIN2h")
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern="^MT-")
Adata2h<- data

Adata2h_filt <- subset(Adata2h, subset = nFeature_RNA > 2500 & nFeature_RNA < 6000 & percent.mt < 7.5 & nCount_RNA > 2000 & nCount_RNA < 60000)
Ddata2h_filt <- subset(Ddata2h, subset = nFeature_RNA > 2500 & nFeature_RNA < 6000 & percent.mt < 7.5 & nCount_RNA > 2000 & nCount_RNA < 60000)

Adata0h_filt@meta.data$time <- "0h"
Adata30m_filt@meta.data$time <- "30m"
Adata2h_filt@meta.data$time <- "2h"
Ddata0h_filt@meta.data$time <- "0h"
Ddata30m_filt@meta.data$time <- "30m"
Ddata2h_filt@meta.data$time <- "2h"


dataB <- merge( x= Ddata0h_filt, y=c(Ddata30m_filt, Ddata2h_filt,  Adata0h_filt, Adata30m_filt, Adata2h_filt), add.cell.id = c("D0","D1","D2","A0","A1","A2"))

dataB$log10GenesPerUMI <- log10(dataB$nFeature_RNA) / log10(dataB$nCount_RNA)
metadata <-dataB@meta.data
metadata$cells <- rownames(metadata)
# Rename columns
metadata <- metadata %>%
        dplyr::rename(seq_folder = orig.ident,
                      nUMI = nCount_RNA,
                      nGene = nFeature_RNA)
metadata$sample <- NA
metadata$sample[which(str_detect(metadata$cells, "^D0_"))] <- "D0"
metadata$sample[which(str_detect(metadata$cells, "^D1_"))] <- "D1"
metadata$sample[which(str_detect(metadata$cells, "^D2_"))] <- "D2"
metadata$sample[which(str_detect(metadata$cells, "^A0_"))] <- "A0"
metadata$sample[which(str_detect(metadata$cells, "^A1_"))] <- "A1"
metadata$sample[which(str_detect(metadata$cells, "^A2_"))] <- "A2"

dataB$sample <- metadata$sample

counts <- GetAssayData(object = dataB, slot = "counts")
nonzero <- counts > 0
keep_genes <- Matrix::rowSums(nonzero) >= 10
filtered_counts <- counts[keep_genes, ]

dataB <- CreateSeuratObject(filtered_counts, meta.data = dataB@meta.data)

dataB$sample2 <- NA
dataB$sample2[which(str_detect(dataB$sample, "^D"))] <- "DMSO"
dataB$sample2[which(str_detect(dataB$sample, "^A"))] <- "AUXIN"

seurat_phase <- NormalizeData(dataB)
seurat_phase <- CellCycleScoring(seurat_phase, 
                                 g2m.features = g2m.genes, 
                                 s.features = s.genes)
seurat_phase <- FindVariableFeatures(seurat_phase, 
                     selection.method = "vst",
                     nfeatures = 2000, 
                     verbose = FALSE)
		   
# Scale the counts
seurat_phase <- ScaleData(seurat_phase)
seurat_phase <- RunPCA(seurat_phase)

# make histogram for G1 samples and all samples
data_filt <- subset(seurat_phase, sample=="D0")
sum(data_filt[["RNA"]]@counts)

genes = rownames(data_filt)
length(genes)
getDistribution(data_filt, genes=genes, out_dir="HCT_DMSO0h_Sfiltered/", COND='HCTDMSO0h')
write.table(data_filt@assays$RNA@counts, file="HCT_DMSO0h.txt", quote=F, sep="\t")

data_filt <- subset(seurat_phase, sample=="D1")
sum(data_filt[["RNA"]]@counts)

genes = rownames(data_filt)
length(genes)
getDistribution(data_filt, genes=genes, out_dir="HCT_DMSO30m_Sfiltered/", COND='HCTDMSO30m')
write.table(data_filt@assays$RNA@counts, file="HCT_DMSO30m.txt", quote=F, sep="\t")

data_filt <- subset(seurat_phase, sample=="D2")
sum(data_filt[["RNA"]]@counts)

genes = rownames(data_filt)
length(genes)
getDistribution(data_filt, genes=genes, out_dir="HCT_DMSO2h_Sfiltered/", COND='HCTDMSO2h')
write.table(data_filt@assays$RNA@counts, file="HCT_DMSO2h.txt", quote=F, sep="\t")

data_filt <- subset(seurat_phase, sample=="A0")
sum(data_filt[["RNA"]]@counts)

genes = rownames(data_filt)
length(genes)
getDistribution(data_filt, genes=genes, out_dir="HCT_AUXIN0h_Sfiltered/", COND='HCTAUXIN0h')
write.table(data_filt@assays$RNA@counts, file="HCT_AUXIN0h.txt", quote=F, sep="\t")

data_filt <- subset(seurat_phase, sample=="A1")
sum(data_filt[["RNA"]]@counts)

genes = rownames(data_filt)
length(genes)
getDistribution(data_filt, genes=genes, out_dir="HCT_AUXIN30m_Sfiltered/", COND='HCTAUXIN30m')
write.table(data_filt@assays$RNA@counts, file="HCT_AUXIN30m.txt", quote=F, sep="\t")

data_filt <- subset(seurat_phase, sample=="A2")
sum(data_filt[["RNA"]]@counts)

genes = rownames(data_filt)
length(genes)
getDistribution(data_filt, genes=genes, out_dir="HCT_AUXIN2h_Sfiltered/", COND='HCTAUXIN2h')
write.table(data_filt@assays$RNA@counts, file="HCT_AUXIN2h.txt", quote=F, sep="\t")

### mouse spleen resting and activation B cell data analysis ###
load("cell_cycle_genes.RData")

T1.data <- Read10X(data.dir="Ctrl_SprB/filtered_feature_bc_matrix/")
data = CreateSeuratObject(counts = T1.data, min.cells = 3, min.features = 200, project = "Ctrl_SprB")
mito.genes_T1 = grep(pattern = "^mt-", x = rownames(x = data@assays[["RNA"]]), value = TRUE)
percent.mito_T1 = Matrix::colSums(data@assays[["RNA"]][mito.genes_T1, ]) / Matrix::colSums(data@assays[["RNA"]])
data = AddMetaData(object = data, metadata = percent.mito_T1, col.name = "percent.mito")
rB_filt = subset(x = data, subset = nFeature_RNA > 1000 & nFeature_RNA < 3000 & percent.mito < 0.05 & nCount_RNA < 12000 & nCount_RNA > 2000 )

T1.data <- Read10X(data.dir="Ctrl_SpaB/filtered_feature_bc_matrix/")
data = CreateSeuratObject(counts = T1.data, min.cells = 3, min.features = 200, project = "Ctrl_SpaB")
mito.genes_T1 = grep(pattern = "^mt-", x = rownames(x = data@assays[["RNA"]]), value = TRUE)
percent.mito_T1 = Matrix::colSums(data@assays[["RNA"]][mito.genes_T1, ]) / Matrix::colSums(data@assays[["RNA"]])
data = AddMetaData(object = data, metadata = percent.mito_T1, col.name = "percent.mito")
aB_filt = subset(x = data, subset = nFeature_RNA > 1500 & nFeature_RNA < 4500 & percent.mito < 0.05 & nCount_RNA < 25000 & nCount_RNA > 2000 )


dataB <- merge( x= rB_filt, y=aB_filt, add.cell.id = c("SprB","SpaB"))

dataB$log10GenesPerUMI <- log10(dataB$nFeature_RNA) / log10(dataB$nCount_RNA)
metadata <-dataB@meta.data
metadata$cells <- rownames(metadata)
# Rename columns
metadata <- metadata %>%
        dplyr::rename(seq_folder = orig.ident,
                      nUMI = nCount_RNA,
                      nGene = nFeature_RNA)
metadata$sample <- NA
metadata$sample[which(str_detect(metadata$cells, "^SprB_"))] <- "rB"
metadata$sample[which(str_detect(metadata$cells, "^SpaB_"))] <- "aB"

dataB$sample <- metadata$sample

counts <- GetAssayData(object = dataB, slot = "counts")
nonzero <- counts > 0
keep_genes <- Matrix::rowSums(nonzero) >= 10
filtered_counts <- counts[keep_genes, ]

dataB <- CreateSeuratObject(filtered_counts, meta.data = dataB@meta.data)


seurat_phase <- NormalizeData(dataB)
seurat_phase <- CellCycleScoring(seurat_phase, 
                                 g2m.features = g2m_genes, 
                                 s.features = s_genes)
seurat_phase <- FindVariableFeatures(seurat_phase, 
                     selection.method = "vst",
                     nfeatures = 2000, 
                     verbose = FALSE)
		     
# Scale the counts
seurat_phase <- ScaleData(seurat_phase)
seurat_phase <- RunPCA(seurat_phase)

# make histogram 
data_filt <- subset(seurat_phase, sample=="rB")
sum(data_filt[["RNA"]]@counts)

genes = rownames(data_filt)
length(genes)
getDistribution(data_filt, genes=genes, out_dir="CtrlSprB_Sfiltered/", COND='CtrlSprB')
write.table(data_filt@assays$RNA@counts, file="CtrlSprB.txt", quote=F, sep="\t")

data_filt <- subset(seurat_phase, sample=="aB")
sum(data_filt[["RNA"]]@counts)

genes = rownames(data_filt)
length(genes)
getDistribution(data_filt, genes=genes, out_dir="CtrlSpaB_Sfiltered/", COND='CtrlSpaB')
write.table(data_filt@assays$RNA@counts, file="CtrlSpaB.txt", quote=F, sep="\t")

## Myc KO mouse data analysis
T1.data <- Read10X(data.dir="Mycm_r1/filtered_feature_bc_matrix/")
data = CreateSeuratObject(counts = T1.data, min.cells = 3, min.features = 200, project = "CH12_data")
mito.genes_T1 = grep(pattern = "^mt-", x = rownames(x = data@assays[["RNA"]]), value = TRUE)
percent.mito_T1 = Matrix::colSums(data@assays[["RNA"]][mito.genes_T1, ]) / Matrix::colSums(data@assays[["RNA"]])
data = AddMetaData(object = data, metadata = percent.mito_T1, col.name = "percent.mito")
M1data_filt = subset(x = data, subset = nFeature_RNA > 900 & nFeature_RNA < 2500 & percent.mito < 0.075 & nCount_RNA < 10000)

T1.data <- Read10X(data.dir="Mycm_r2/filtered_feature_bc_matrix/")
data = CreateSeuratObject(counts = T1.data, min.cells = 3, min.features = 200, project = "CH12_data")
mito.genes_T1 = grep(pattern = "^mt-", x = rownames(x = data@assays[["RNA"]]), value = TRUE)
percent.mito_T1 = Matrix::colSums(data@assays[["RNA"]][mito.genes_T1, ]) / Matrix::colSums(data@assays[["RNA"]])
data = AddMetaData(object = data, metadata = percent.mito_T1, col.name = "percent.mito")
M2data_filt = subset(x = data, subset = nFeature_RNA > 900 & nFeature_RNA < 2500 & percent.mito < 0.075 & nCount_RNA < 10000)

T1.data <- Read10X(data.dir="WT_r1/filtered_feature_bc_matrix/")
data = CreateSeuratObject(counts = T1.data, min.cells = 3, min.features = 200, project = "CH12_data")
mito.genes_T1 = grep(pattern = "^mt-", x = rownames(x = data@assays[["RNA"]]), value = TRUE)
percent.mito_T1 = Matrix::colSums(data@assays[["RNA"]][mito.genes_T1, ]) / Matrix::colSums(data@assays[["RNA"]])
data = AddMetaData(object = data, metadata = percent.mito_T1, col.name = "percent.mito")
W1data_filt = subset(x = data, subset = nFeature_RNA > 900 & nFeature_RNA < 2500 & percent.mito < 0.075 & nCount_RNA < 10000)

T1.data <- Read10X(data.dir="WT_r2/filtered_feature_bc_matrix/")
data = CreateSeuratObject(counts = T1.data, min.cells = 3, min.features = 200, project = "CH12_data")
mito.genes_T1 = grep(pattern = "^mt-", x = rownames(x = data@assays[["RNA"]]), value = TRUE)
percent.mito_T1 = Matrix::colSums(data@assays[["RNA"]][mito.genes_T1, ]) / Matrix::colSums(data@assays[["RNA"]])
data = AddMetaData(object = data, metadata = percent.mito_T1, col.name = "percent.mito")
W2data_filt = subset(x = data, subset = nFeature_RNA > 900 & nFeature_RNA < 2500 & percent.mito < 0.075 & nCount_RNA < 10000)

M1data_filt@meta.data$rep <- "rep1"
M2data_filt@meta.data$rep <- "rep2"
W1data_filt@meta.data$rep <- "rep1"
W2data_filt@meta.data$rep <- "rep2"

dataB <- merge( x= W1data_filt, y=c(W2data_filt, M1data_filt,M2data_filt), add.cell.id = c("WTr1","WTr2","Mycmmr1","Mycmmr2"))

dataB$log10GenesPerUMI <- log10(dataB$nFeature_RNA) / log10(dataB$nCount_RNA)
metadata <-dataB@meta.data
metadata$cells <- rownames(metadata)
# Rename columns
metadata <- metadata %>%
        dplyr::rename(seq_folder = orig.ident,
                      nUMI = nCount_RNA,
                      nGene = nFeature_RNA)
metadata$sample <- NA
metadata$sample[which(str_detect(metadata$cells, "^WTr1_"))] <- "WTr1"
metadata$sample[which(str_detect(metadata$cells, "^Mycmmr1_"))] <- "Mycmmr1"
metadata$sample[which(str_detect(metadata$cells, "^WTr2_"))] <- "WTr2"
metadata$sample[which(str_detect(metadata$cells, "^Mycmmr2_"))] <- "Mycmmr2"

dataB$sample <- metadata$sample

counts <- GetAssayData(object = dataB, slot = "counts")
nonzero <- counts > 0
keep_genes <- Matrix::rowSums(nonzero) >= 10
filtered_counts <- counts[keep_genes, ]


dataB <- CreateSeuratObject(filtered_counts, meta.data = dataB@meta.data)

dataB$sample2 <- NA
dataB$sample2[which(str_detect(dataB$sample, "^WT"))] <- "WT"
dataB$sample2[which(str_detect(dataB$sample, "^Mycmm"))] <- "Mycmm"

seurat_phase <- NormalizeData(dataB)
seurat_phase <- CellCycleScoring(seurat_phase, 
                                 g2m.features = g2m_genes, 
                                 s.features = s_genes)
seurat_phase <- FindVariableFeatures(seurat_phase, 
                     selection.method = "vst",
                     nfeatures = 2000, 
                     verbose = FALSE)
		     
# Scale the counts
seurat_phase <- ScaleData(seurat_phase)
seurat_phase <- RunPCA(seurat_phase)

data_filt <- subset(seurat_phase, sample=="Mycmmr1")
sum(data_filt[["RNA"]]@counts)

genes = rownames(data_filt)
length(genes)
getDistribution(data_filt, genes=genes, out_dir="Mycm_r1_Sfiltered/", COND='Mycmr1')
write.table(data_filt@assays$RNA@counts, file="Mycmr1.txt", quote=F, sep="\t")

data_filt <- subset(seurat_phase, sample=="Mycmmr2")
sum(data_filt[["RNA"]]@counts)

genes = rownames(data_filt)
length(genes)
getDistribution(data_filt, genes=genes, out_dir="Mycm_r2_Sfiltered/", COND='Mycmr2')
write.table(data_filt@assays$RNA@counts, file="Mycmr2.txt", quote=F, sep="\t")

data_filt <- subset(seurat_phase, sample=="WTr1")
sum(data_filt[["RNA"]]@counts)

genes = rownames(data_filt)
length(genes)
getDistribution(data_filt, genes=genes, out_dir="xWT_r1_Sfiltered/", COND='xWTr1')
write.table(data_filt@assays$RNA@counts, file="xWTr1.txt", quote=F, sep="\t")

data_filt <- subset(seurat_phase, sample=="WTr2")
sum(data_filt[["RNA"]]@counts)

genes = rownames(data_filt)
length(genes)
getDistribution(data_filt, genes=genes, out_dir="xWT_r2_Sfiltered/", COND='xWTr2')
write.table(data_filt@assays$RNA@counts, file="xWTr2.txt", quote=F, sep="\t")

### LNK, SKIN, mast, BoneMarrow cell process ###
T1.data <- Read10X(data.dir="LNK/filtered_feature_bc_matrix/")
data = CreateSeuratObject(counts = T1.data, min.cells = 3, min.features = 200, project = "NK")
mito.genes_T1 = grep(pattern = "^mt-", x = rownames(x = data@assays[["RNA"]]), value = TRUE)
percent.mito_T1 = Matrix::colSums(data@assays[["RNA"]][mito.genes_T1, ]) / Matrix::colSums(data@assays[["RNA"]])
data = AddMetaData(object = data, metadata = percent.mito_T1, col.name = "percent.mito")

data_filt = subset(x = data, subset = nFeature_RNA > 3000 & nFeature_RNA < 6000 & percent.mito < 0.05 & nCount_RNA < 50000 & nCount_RNA > 10000)

dataB <- data_filt

dataB$log10GenesPerUMI <- log10(dataB$nFeature_RNA) / log10(dataB$nCount_RNA)
metadata <-dataB@meta.data
metadata$cells <- rownames(metadata)
# Rename columns
metadata <- metadata %>%
        dplyr::rename(seq_folder = orig.ident,
                      nUMI = nCount_RNA,
                      nGene = nFeature_RNA)

counts <- GetAssayData(object = dataB, slot = "counts")
nonzero <- counts > 0
keep_genes <- Matrix::rowSums(nonzero) >= 10
filtered_counts <- counts[keep_genes, ]


dataB <- CreateSeuratObject(filtered_counts, meta.data = dataB@meta.data)


seurat_phase <- NormalizeData(dataB)
seurat_phase <- CellCycleScoring(seurat_phase, 
                                 g2m.features = g2m_genes, 
                                 s.features = s_genes)
seurat_phase <- FindVariableFeatures(seurat_phase, 
                     selection.method = "vst",
                     nfeatures = 2000, 
                     verbose = FALSE)
		     
# Scale the counts
seurat_phase <- ScaleData(seurat_phase)
seurat_phase <- RunPCA(seurat_phase)
# make histogram
data_filt <- seurat_phase
genes = rownames(data_filt)
length(genes)
getDistribution(data_filt, genes=genes, out_dir="NK_Sfiltered/", COND='NK')
write.table(data_filt@assays$RNA@counts, file="NK.txt", quote=F, sep="\t")

T1.data <- Read10X(data.dir="Mast/filtered_feature_bc_matrix/")
data = CreateSeuratObject(counts = T1.data, min.cells = 3, min.features = 200, project = "MAST")
mito.genes_T1 = grep(pattern = "^mt-", x = rownames(x = data@assays[["RNA"]]), value = TRUE)
percent.mito_T1 = Matrix::colSums(data@assays[["RNA"]][mito.genes_T1, ]) / Matrix::colSums(data@assays[["RNA"]])
data = AddMetaData(object = data, metadata = percent.mito_T1, col.name = "percent.mito")

data_filt = subset(x = data, subset = nFeature_RNA > 1500 & nFeature_RNA < 3000 & percent.mito < 0.05 & nCount_RNA < 10000 & nCount_RNA > 2000)

dataB <- data_filt

dataB$log10GenesPerUMI <- log10(dataB$nFeature_RNA) / log10(dataB$nCount_RNA)
metadata <-dataB@meta.data
metadata$cells <- rownames(metadata)
# Rename columns
metadata <- metadata %>%
        dplyr::rename(seq_folder = orig.ident,
                      nUMI = nCount_RNA,
                      nGene = nFeature_RNA)

counts <- GetAssayData(object = dataB, slot = "counts")
nonzero <- counts > 0
keep_genes <- Matrix::rowSums(nonzero) >= 10
filtered_counts <- counts[keep_genes, ]


dataB <- CreateSeuratObject(filtered_counts, meta.data = dataB@meta.data)


seurat_phase <- NormalizeData(dataB)
seurat_phase <- CellCycleScoring(seurat_phase, 
                                 g2m.features = g2m_genes, 
                                 s.features = s_genes)
seurat_phase <- FindVariableFeatures(seurat_phase, 
                     selection.method = "vst",
                     nfeatures = 2000, 
                     verbose = FALSE)
		     
# Scale the counts
seurat_phase <- ScaleData(seurat_phase)
seurat_phase <- RunPCA(seurat_phase)

# make histogram 
data_filt <- seurat_phase
sum(data_filt[["RNA"]]@counts)

genes = rownames(data_filt)
length(genes)
getDistribution(data_filt, genes=genes, out_dir="MAST_Sfiltered/", COND='MAST')
write.table(data_filt@assays$RNA@counts, file="MAST.txt", quote=F, sep="\t")

T1.data <- Read10X(data.dir="B220BM_rep1/filtered_feature_bc_matrix/")
data = CreateSeuratObject(counts = T1.data, min.cells = 3, min.features = 200, project = "BM_r1")
mito.genes_T1 = grep(pattern = "^mt-", x = rownames(x = data@assays[["RNA"]]), value = TRUE)
percent.mito_T1 = Matrix::colSums(data@assays[["RNA"]][mito.genes_T1, ]) / Matrix::colSums(data@assays[["RNA"]])
data = AddMetaData(object = data, metadata = percent.mito_T1, col.name = "percent.mito")

data1_filt = subset(x = data, subset = nFeature_RNA > 900 & nFeature_RNA < 2500 & percent.mito < 0.05 & nCount_RNA < 10000)

T1.data <- Read10X(data.dir="B220BM_rep2/filtered_feature_bc_matrix/")
data = CreateSeuratObject(counts = T1.data, min.cells = 3, min.features = 200, project = "BM_r2")
mito.genes_T1 = grep(pattern = "^mt-", x = rownames(x = data@assays[["RNA"]]), value = TRUE)
percent.mito_T1 = Matrix::colSums(data@assays[["RNA"]][mito.genes_T1, ]) / Matrix::colSums(data@assays[["RNA"]])
data = AddMetaData(object = data, metadata = percent.mito_T1, col.name = "percent.mito")

data2_filt = subset(x = data, subset = nFeature_RNA > 900 & nFeature_RNA < 2500 & percent.mito < 0.05 & nCount_RNA < 10000)

dataB <- merge( x= data1_filt, y=data2_filt, add.cell.id = c("BMr1","BMr2"))

dataB$log10GenesPerUMI <- log10(dataB$nFeature_RNA) / log10(dataB$nCount_RNA)
metadata <-dataB@meta.data
metadata$cells <- rownames(metadata)
# Rename columns
metadata <- metadata %>%
        dplyr::rename(seq_folder = orig.ident,
                      nUMI = nCount_RNA,
                      nGene = nFeature_RNA)
metadata$sample <- NA
metadata$sample[which(str_detect(metadata$cells, "^BMr1_"))] <- "r1"
metadata$sample[which(str_detect(metadata$cells, "^BMr2_"))] <- "r2"

dataB$sample <- metadata$sample

counts <- GetAssayData(object = dataB, slot = "counts")
nonzero <- counts > 0
keep_genes <- Matrix::rowSums(nonzero) >= 10
filtered_counts <- counts[keep_genes, ]

dataB <- CreateSeuratObject(filtered_counts, meta.data = dataB@meta.data)

seurat_phase <- NormalizeData(dataB)
seurat_phase <- CellCycleScoring(seurat_phase,
                                 g2m.features = g2m_genes,
                                 s.features = s_genes)
seurat_phase <- FindVariableFeatures(seurat_phase,
                     selection.method = "vst",
                     nfeatures = 2000,
                     verbose = FALSE)

# Scale the counts
seurat_phase <- ScaleData(seurat_phase)
seurat_phase <- RunPCA(seurat_phase)

#histogram
data_filt <- subset(seurat_phase, sample == "r1")

genes = rownames(data_filt)
length(genes)
getDistribution(data_filt, genes=genes, out_dir="WTB220BM_rep1_Sfiltered/", COND='WTB220BMr1')
write.table(data_filt@assays$RNA@counts, file="WTB220BM_rep1.txt", quote=F, sep="\t")

data_filt <- subset(seurat_phase, sample == "r2")

genes = rownames(data_filt)
length(genes)
getDistribution(data_filt, genes=genes, out_dir="WTB220BM_rep2_Sfiltered/", COND='WTB220BMr2')
write.table(data_filt@assays$RNA@counts, file="WTB220BM_rep2.txt", quote=F, sep="\t")

T1.data <- Read10X(data.dir="Skin_r1/filtered_feature_bc_matrix/")
data = CreateSeuratObject(counts = T1.data, min.cells = 3, min.features = 200, project = "Skin_r1")
mito.genes_T1 = grep(pattern = "^mt-", x = rownames(x = data@assays[["RNA"]]), value = TRUE)
percent.mito_T1 = Matrix::colSums(data@assays[["RNA"]][mito.genes_T1, ]) / Matrix::colSums(data@assays[["RNA"]])
data = AddMetaData(object = data, metadata = percent.mito_T1, col.name = "percent.mito")

data1_filt = subset(x = data, subset = nFeature_RNA > 1000 & nFeature_RNA < 8000 & percent.mito < 0.05 & nCount_RNA < 50000)

T1.data <- Read10X(data.dir="Skin_r2/filtered_feature_bc_matrix/")
data = CreateSeuratObject(counts = T1.data, min.cells = 3, min.features = 200, project = "Skin_r2")
mito.genes_T1 = grep(pattern = "^mt-", x = rownames(x = data@assays[["RNA"]]), value = TRUE)
percent.mito_T1 = Matrix::colSums(data@assays[["RNA"]][mito.genes_T1, ]) / Matrix::colSums(data@assays[["RNA"]])
data = AddMetaData(object = data, metadata = percent.mito_T1, col.name = "percent.mito")

data2_filt = subset(x = data, subset = nFeature_RNA > 1000 & nFeature_RNA < 8000 & percent.mito < 0.05 & nCount_RNA < 50000)

dataB <- merge( x= data1_filt, y=data2_filt, add.cell.id = c("Sr1","Sr2"))

dataB$log10GenesPerUMI <- log10(dataB$nFeature_RNA) / log10(dataB$nCount_RNA)
metadata <-dataB@meta.data
metadata$cells <- rownames(metadata)
# Rename columns
metadata <- metadata %>%
        dplyr::rename(seq_folder = orig.ident,
                      nUMI = nCount_RNA,
                      nGene = nFeature_RNA)
metadata$sample <- NA
metadata$sample[which(str_detect(metadata$cells, "^Sr1_"))] <- "r1"
metadata$sample[which(str_detect(metadata$cells, "^Sr2_"))] <- "r2"

dataB$sample <- metadata$sample

counts <- GetAssayData(object = dataB, slot = "counts")
nonzero <- counts > 0
keep_genes <- Matrix::rowSums(nonzero) >= 10
filtered_counts <- counts[keep_genes, ]

dataB <- CreateSeuratObject(filtered_counts, meta.data = dataB@meta.data)

seurat_phase <- NormalizeData(dataB)
seurat_phase <- CellCycleScoring(seurat_phase,
                                 g2m.features = g2m_genes,
                                 s.features = s_genes)
seurat_phase <- FindVariableFeatures(seurat_phase,
                     selection.method = "vst",
                     nfeatures = 2000,
                     verbose = FALSE)

# Scale the counts
seurat_phase <- ScaleData(seurat_phase)
seurat_phase <- RunPCA(seurat_phase)
#histogram
data_filt <- subset(seurat_phase, sample == "r1")

genes = rownames(data_filt)
length(genes)
getDistribution(data_filt, genes=genes, out_dir="SKIN_rep1_Sfiltered/", COND='SKINr1')
write.table(data_filt@assays$RNA@counts, file="SKIN_rep1.txt", quote=F, sep="\t")

data_filt <- subset(seurat_phase, sample == "r2")

genes = rownames(data_filt)
length(genes)
getDistribution(data_filt, genes=genes, out_dir="SKIN_rep2_Sfiltered/", COND='SKINr2')
write.table(data_filt@assays$RNA@counts, file="SKIN_rep2.txt", quote=F, sep="\t")

### HCT non-starved data process ###
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

T1.data <- Read10X(data.dir="HCT_DMSO2h_rep1m/filtered_feature_bc_matrix/")
data = CreateSeuratObject(counts = T1.data, min.cells = 3, min.features = 200, project = "Drep1")
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern="^MT-")
Ddata1 <- data

Ddata1_filt <- subset(Ddata1, subset = nFeature_RNA > 4000 & nFeature_RNA < 8000 & percent.mt < 5 & nCount_RNA > 20000 & nCount_RNA < 70000)


T1.data <- Read10X(data.dir="HCT_DMSO2h_rep2m/filtered_feature_bc_matrix/")
data = CreateSeuratObject(counts = T1.data, min.cells = 3, min.features = 200, project = "Drep2")
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern="^MT-")
Ddata2 <- data

Ddata2_filt <- subset(Ddata2, subset = nFeature_RNA > 4000 & nFeature_RNA < 8000 & percent.mt < 5 & nCount_RNA > 20000 & nCount_RNA < 70000)

Ddata1_filt@meta.data$rep <- "rep1"
Ddata2_filt@meta.data$rep <- "rep2"

dataB <- merge( x= Ddata1_filt, y=Ddata2_filt, add.cell.id = c("Dr1","Dr2"))

dataB$log10GenesPerUMI <- log10(dataB$nFeature_RNA) / log10(dataB$nCount_RNA)
metadata <-dataB@meta.data
metadata$cells <- rownames(metadata)
# Rename columns
metadata <- metadata %>%
        dplyr::rename(seq_folder = orig.ident,
                      nUMI = nCount_RNA,
                      nGene = nFeature_RNA)
metadata$sample <- NA
metadata$sample[which(str_detect(metadata$cells, "^Dr1_"))] <- "Dr1"
metadata$sample[which(str_detect(metadata$cells, "^Dr2_"))] <- "Dr2"

dataB$sample <- metadata$sample

counts <- GetAssayData(object = dataB, slot = "counts")
nonzero <- counts > 0
keep_genes <- Matrix::rowSums(nonzero) >= 10
filtered_counts <- counts[keep_genes, ]

dataB <- CreateSeuratObject(filtered_counts, meta.data = dataB@meta.data)

dataB$sample2 <- NA
dataB$sample2[which(str_detect(dataB$sample, "^D"))] <- "DMSO"
dataB$sample2[which(str_detect(dataB$sample, "^J"))] <- "JQ1"

seurat_phase <- NormalizeData(dataB)
seurat_phase <- CellCycleScoring(seurat_phase, 
                                 g2m.features = g2m.genes, 
                                 s.features = s.genes)
seurat_phase <- FindVariableFeatures(seurat_phase, 
                     selection.method = "vst",
                     nfeatures = 2000, 
                     verbose = FALSE)
		     
# Scale the counts
seurat_phase <- ScaleData(seurat_phase)
seurat_phase <- RunPCA(seurat_phase)

# make histogram for Dr1 samples 
data_filt <- subset(seurat_phase, sample=="Dr1")

genes = rownames(data_filt)
length(genes)
getDistribution(data_filt, genes=genes, out_dir="HCT_DMSOr1m_Sfiltered/", COND='HCTDMSOr1')
write.table(data_filt@assays$RNA@counts, file="HCT_DMSOr1mS.txt", quote=F, sep="\t")

# make histogram for Dr2 samples 
data_filt <- subset(seurat_phase, sample=="Dr2")

genes = rownames(data_filt)
length(genes)
getDistribution(data_filt, genes=genes, out_dir="HCT_DMSOr2m_Sfiltered/", COND='HCTDMSOr2')
write.table(data_filt@assays$RNA@counts, file="HCT_DMSOr2mS.txt", quote=F, sep="\t")


# select the same number of cells
data_filt <- subset(seurat_phase, sample=="Dr2" & Phase == "G1")
SampleN <- table(data_filt$Phase)

# make histogram for Dr1 samples 
data_filt <- subset(seurat_phase, sample=="Dr1" & Phase == "G1")
Nlist <- sample(table(data_filt$Phase), SampleN, replace = FALSE)
genes = rownames(data_filt)
getDistribution2(data_filt, genes=genes, out_dir="HCT_DMSOr1m_G1_Sfiltered2/", COND='HCTDMSOr1G1', Nlist)
write.table(data_filt@assays$RNA@counts[,Nlist], file="HCT_DMSOr1mS2_G1.txt", quote=F, sep="\t")

data_filt <- subset(seurat_phase, sample=="Dr1" & Phase == "G2M")
Nlist <- sample(table(data_filt$Phase), SampleN, replace = FALSE)
genes = rownames(data_filt)
getDistribution2(data_filt, genes=genes, out_dir="HCT_DMSOr1m_G2M_Sfiltered2/", COND='HCTDMSOr1G2M',Nlist)
write.table(data_filt@assays$RNA@counts[,Nlist], file="HCT_DMSOr1mS2_G2M.txt", quote=F, sep="\t")

data_filt <- subset(seurat_phase, sample=="Dr1" & Phase == "S")
Nlist <- sample(table(data_filt$Phase), SampleN, replace = FALSE)
genes = rownames(data_filt)
getDistribution2(data_filt, genes=genes, out_dir="HCT_DMSOr1m_S_Sfiltered2/", COND='HCTDMSOr1S',Nlist)
write.table(data_filt@assays$RNA@counts[,Nlist], file="HCT_DMSOr1mS2_S.txt", quote=F, sep="\t")

# make histogram for Dr2 samples 
data_filt <- subset(seurat_phase, sample=="Dr2" & Phase == "G2M")
Nlist <- sample(table(data_filt$Phase), SampleN, replace = FALSE)
genes = rownames(data_filt)
getDistribution2(data_filt, genes=genes, out_dir="HCT_DMSOr2m_G2M_Sfiltered2/", COND='HCTDMSOr2G2M',Nlist)
write.table(data_filt@assays$RNA@counts[,Nlist], file="HCT_DMSOr2mS2_G2M.txt", quote=F, sep="\t")

data_filt <- subset(seurat_phase, sample=="Dr2" & Phase == "S")
Nlist <- sample(table(data_filt$Phase), SampleN, replace = FALSE)
genes = rownames(data_filt)
getDistribution2(data_filt, genes=genes, out_dir="HCT_DMSOr2m_S_Sfiltered2/", COND='HCTDMSOr2S',Nlist)
write.table(data_filt@assays$RNA@counts[,Nlist], file="HCT_DMSOr2mS2_S.txt", quote=F, sep="\t")

### HCT JQ1 treatment data analysis ###
T1.data <- Read10X(data.dir="HCT_JQ1Sr1/filtered_feature_bc_matrix/")
data = CreateSeuratObject(counts = T1.data, min.cells = 3, min.features = 200, project = "Jrep1")
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern="^MT-")
Jdata1<- data

T1.data <- Read10X(data.dir="HCT_DMSOSr1/filtered_feature_bc_matrix/")
data = CreateSeuratObject(counts = T1.data, min.cells = 3, min.features = 200, project = "Drep1")
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern="^MT-")
Ddata1 <- data

Jdata1_filt <- subset(Jdata1, subset = nFeature_RNA > 4500 & nFeature_RNA < 10000 & percent.mt < 7.5 & nCount_RNA > 20000 & nCount_RNA < 70000)
Ddata1_filt <- subset(Ddata1, subset = nFeature_RNA > 4500 & nFeature_RNA < 10000 & percent.mt < 7.5 & nCount_RNA > 20000 & nCount_RNA < 70000)

T1.data <- Read10X(data.dir="HCT_JQ1Sr2/filtered_feature_bc_matrix/")
data = CreateSeuratObject(counts = T1.data, min.cells = 3, min.features = 200, project = "Jrep2")
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern="^MT-")
Jdata2<- data

T1.data <- Read10X(data.dir="HCT_DMSOSr2/filtered_feature_bc_matrix/")
data = CreateSeuratObject(counts = T1.data, min.cells = 3, min.features = 200, project = "Drep2")
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern="^MT-")
Ddata2 <- data

Jdata2_filt <- subset(Jdata2, subset = nFeature_RNA > 4500 & nFeature_RNA < 10000 & percent.mt < 7.5 & nCount_RNA > 20000 & nCount_RNA < 70000)
Ddata2_filt <- subset(Ddata2, subset = nFeature_RNA > 4500 & nFeature_RNA < 10000 & percent.mt < 7.5 & nCount_RNA > 20000 & nCount_RNA < 70000)

Jdata1_filt@meta.data$rep <- "rep1"
Jdata2_filt@meta.data$rep <- "rep2"
Ddata1_filt@meta.data$rep <- "rep1"
Ddata2_filt@meta.data$rep <- "rep2"

dataB <- merge( x= Ddata1_filt, y=c(Ddata2_filt, Jdata1_filt,Jdata2_filt), add.cell.id = c("Dr1","Dr2","Jr1","Jr2"))

dataB$log10GenesPerUMI <- log10(dataB$nFeature_RNA) / log10(dataB$nCount_RNA)
metadata <-dataB@meta.data
metadata$cells <- rownames(metadata)
# Rename columns
metadata <- metadata %>%
        dplyr::rename(seq_folder = orig.ident,
                      nUMI = nCount_RNA,
                      nGene = nFeature_RNA)
metadata$sample <- NA
metadata$sample[which(str_detect(metadata$cells, "^Dr1_"))] <- "Dr1"
metadata$sample[which(str_detect(metadata$cells, "^Jr1_"))] <- "Jr1"
metadata$sample[which(str_detect(metadata$cells, "^Dr2_"))] <- "Dr2"
metadata$sample[which(str_detect(metadata$cells, "^Jr2_"))] <- "Jr2"

dataB$sample <- metadata$sample

counts <- GetAssayData(object = dataB, slot = "counts")
nonzero <- counts > 0
keep_genes <- Matrix::rowSums(nonzero) >= 10
filtered_counts <- counts[keep_genes, ]

dataB <- CreateSeuratObject(filtered_counts, meta.data = dataB@meta.data)

dataB$sample2 <- NA
dataB$sample2[which(str_detect(dataB$sample, "^D"))] <- "DMSO"
dataB$sample2[which(str_detect(dataB$sample, "^J"))] <- "JQ1"

seurat_phase <- NormalizeData(dataB)
seurat_phase <- CellCycleScoring(seurat_phase, 
                                 g2m.features = g2m.genes, 
                                 s.features = s.genes)
seurat_phase <- FindVariableFeatures(seurat_phase, 
                     selection.method = "vst",
                     nfeatures = 2000, 
                     verbose = FALSE)
		     
# Scale the counts
seurat_phase <- ScaleData(seurat_phase)
seurat_phase <- RunPCA(seurat_phase)

# make histogram for Dr1 samples 
data_filt <- subset(seurat_phase, sample=="Dr1")
sum(data_filt[["RNA"]]@counts)

genes = rownames(data_filt)
length(genes)
getDistribution(data_filt, genes=genes, out_dir="HCT_DMSOSr1_Sfiltered/", COND='HCTDMSOr1')
write.table(data_filt@assays$RNA@counts, file="HCT_DMSOSr1.txt", quote=F, sep="\t")

# make histogram for Dr2 samples 
data_filt <- subset(seurat_phase, sample=="Dr2")
sum(data_filt[["RNA"]]@counts)

genes = rownames(data_filt)
length(genes)
getDistribution(data_filt, genes=genes, out_dir="HCT_DMSOSr2_Sfiltered/", COND='HCTDMSOr2')
write.table(data_filt@assays$RNA@counts, file="HCT_DMSOSr2.txt", quote=F, sep="\t")

# make histogram for Jr1 samples 
data_filt <- subset(seurat_phase, sample=="Jr1")
sum(data_filt[["RNA"]]@counts)

genes = rownames(data_filt)
length(genes)
getDistribution(data_filt, genes=genes, out_dir="HCT_JQ1Sr1_Sfiltered/", COND='HCTJQ1r1')
write.table(data_filt@assays$RNA@counts, file="HCT_JQ1Sr1.txt", quote=F, sep="\t")

# make histogram for Jr2 samples 
data_filt <- subset(seurat_phase, sample=="Jr2")
sum(data_filt[["RNA"]]@counts)

genes = rownames(data_filt)
length(genes)
getDistribution(data_filt, genes=genes, out_dir="HCT_JQ1Sr2_Sfiltered/", COND='HCTJQ1r2')
write.table(data_filt@assays$RNA@counts, file="HCT_JQ1Sr2.txt", quote=F, sep="\t")

### HCT MED26 degron data process ###
T1.data <- Read10X(data.dir="MED26D_IAAr1M/filtered_feature_bc_matrix/")
data = CreateSeuratObject(counts = T1.data, min.cells = 3, min.features = 200, project = "Irep1")
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern="^MT-")

Jdata1<- data

T1.data <- Read10X(data.dir="MED26D_DMSOr1M/filtered_feature_bc_matrix/")
data = CreateSeuratObject(counts = T1.data, min.cells = 3, min.features = 200, project = "Drep1")
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern="^MT-")

Ddata1 <- data

Jdata1_filt <- subset(Jdata1, subset = nFeature_RNA > 7500 & nFeature_RNA < 13000 & percent.mt < 5 & nCount_RNA > 20000 & nCount_RNA < 150000)
Ddata1_filt <- subset(Ddata1, subset = nFeature_RNA > 7000 & nFeature_RNA < 12000 & percent.mt < 5 & nCount_RNA > 20000 & nCount_RNA < 150000)

T1.data <- Read10X(data.dir="MED26D_IAAr2M/filtered_feature_bc_matrix/")
data = CreateSeuratObject(counts = T1.data, min.cells = 3, min.features = 200, project = "Irep2")
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern="^MT-")

Jdata2<- data

T1.data <- Read10X(data.dir="MED26D_DMSOr2M/filtered_feature_bc_matrix/")
data = CreateSeuratObject(counts = T1.data, min.cells = 3, min.features = 200, project = "Drep2")
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern="^MT-")

Ddata2 <- data

Jdata2_filt <- subset(Jdata2, subset = nFeature_RNA > 6000 & nFeature_RNA < 12000 & percent.mt < 5 & nCount_RNA > 20000 & nCount_RNA < 150000)
Ddata2_filt <- subset(Ddata2, subset = nFeature_RNA > 7000 & nFeature_RNA < 12000 & percent.mt < 5 & nCount_RNA > 20000 & nCount_RNA < 150000)

Jdata1_filt@meta.data$rep <- "rep1"
Jdata2_filt@meta.data$rep <- "rep2"
Ddata1_filt@meta.data$rep <- "rep1"
Ddata2_filt@meta.data$rep <- "rep2"

dataB <- merge( x= Ddata1_filt, y=c(Ddata2_filt, Jdata1_filt,Jdata2_filt), add.cell.id = c("Dr1","Dr2","Ir1","Ir2"))

dataB$log10GenesPerUMI <- log10(dataB$nFeature_RNA) / log10(dataB$nCount_RNA)
metadata <-dataB@meta.data
metadata$cells <- rownames(metadata)
# Rename columns
metadata <- metadata %>%
        dplyr::rename(seq_folder = orig.ident,
                      nUMI = nCount_RNA,
                      nGene = nFeature_RNA)
metadata$sample <- NA
metadata$sample[which(str_detect(metadata$cells, "^Dr1_"))] <- "Dr1"
metadata$sample[which(str_detect(metadata$cells, "^Ir1_"))] <- "Ir1"
metadata$sample[which(str_detect(metadata$cells, "^Dr2_"))] <- "Dr2"
metadata$sample[which(str_detect(metadata$cells, "^Ir2_"))] <- "Ir2"

dataB$sample <- metadata$sample

counts <- GetAssayData(object = dataB, slot = "counts")
nonzero <- counts > 0
keep_genes <- Matrix::rowSums(nonzero) >= 10
filtered_counts <- counts[keep_genes, ]

dataB <- CreateSeuratObject(filtered_counts, meta.data = dataB@meta.data)

dataB$sample2 <- NA
dataB$sample2[which(str_detect(dataB$sample, "^D"))] <- "DMSO"
dataB$sample2[which(str_detect(dataB$sample, "^I"))] <- "IAA"

seurat_phase <- NormalizeData(dataB)
seurat_phase <- CellCycleScoring(seurat_phase, 
                                 g2m.features = g2m.genes, 
                                 s.features = s.genes)
seurat_phase <- FindVariableFeatures(seurat_phase, 
                     selection.method = "vst",
                     nfeatures = 2000, 
                     verbose = FALSE)
		     
# Scale the counts
seurat_phase <- ScaleData(seurat_phase)
seurat_phase <- RunPCA(seurat_phase)

table(seurat_phase$sample)
table(seurat_phase$sample, seurat_phase$Phase)
# make histogram for Dr1 samples 
data_filt <- subset(seurat_phase, sample=="Dr1")
sum(data_filt[["RNA"]]@counts)

genes = rownames(data_filt)
length(genes)
getDistribution(data_filt, genes=genes, out_dir="HCTMED26D_DMSOr1M_Sfiltered/", COND='MED26DDMSOr1')
write.table(data_filt@assays$RNA@counts, file="HCTMED26D_DMSOr1M.txt", quote=F, sep="\t")

# make histogram for Dr2 samples 
data_filt <- subset(seurat_phase, sample=="Dr2")
sum(data_filt[["RNA"]]@counts)

genes = rownames(data_filt)
length(genes)
getDistribution(data_filt, genes=genes, out_dir="HCTMED26D_DMSOr2M_Sfiltered/", COND='MED26DDMSOr2')
write.table(data_filt@assays$RNA@counts, file="HCTMED26D_DMSOr2M.txt", quote=F, sep="\t")

# make histogram for Ir1 samples 
data_filt <- subset(seurat_phase, sample=="Ir1")
sum(data_filt[["RNA"]]@counts)

genes = rownames(data_filt)
length(genes)
getDistribution(data_filt, genes=genes, out_dir="HCTMED26D_IAAr1M_Sfiltered/", COND='MED26DIAAr1')
write.table(data_filt@assays$RNA@counts, file="HCTMED26D_IAAr1M.txt", quote=F, sep="\t")

# make histogram for Ir2 samples 
data_filt <- subset(seurat_phase, sample=="Ir2")
sum(data_filt[["RNA"]]@counts)

genes = rownames(data_filt)
length(genes)
getDistribution(data_filt, genes=genes, out_dir="HCTMED26D_IAAr2M_Sfiltered/", COND='MED26DIAAr2')
write.table(data_filt@assays$RNA@counts, file="HCTMED26D_IAAr2M.txt", quote=F, sep="\t")

