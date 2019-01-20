library(Matrix)
library(cellrangerRkit)
library(ggplot2)
library(scran)
library(scater)
library(edgeR)
library(limma)

## Read Data ##
temp.path <- '~/SingleCellData/single_cell_analysis/Marusyk_Joint3Aggr_ALK_redo'
genome <- 'GRCh38'
count <- load_cellranger_matrix(temp.path)
aggr.data <- read.csv('~/SingleCellData/single_cell_analysis/Marusyk_Joint3Aggr_ALK_redo/outs/aggregation_csv.csv')
targets <- as.character(aggr.data$library_id[as.numeric(gsub('^.+-', replacement = '', colnames(count)))])
count.mat <- as.matrix(count)

## Map Gene Symbols ##
ensembl <- useMart('ENSEMBL_MART_ENSEMBL', dataset = 'hsapiens_gene_ensembl')
mapping <- getBM(mart = ensembl,
                 attributes = c('ensembl_gene_id', 'hgnc_symbol'),
                 filters = 'ensembl_gene_id',
                 values = rownames(count.mat))
mapping <- na.omit(mapping[mapping$ensembl_gene_id != '' & mapping$hgnc_symbol != '',])
mapping <- mapping[match(rownames(count.mat), table = mapping$ensembl_gene_id),]
rownames(mapping) <- 1:nrow(mapping)
idx <- is.na(mapping$ensembl_gene_id) | is.na(mapping$hgnc_symbol)
count.mat <- count.mat[!idx,]
mapping <- mapping[!idx,]
rownames(mapping) <- 1:nrow(mapping)
idx <- order(Matrix::rowMeans(count.mat), decreasing = TRUE)
count.mat <- count.mat[idx,]
mapping <- mapping[idx,]
idx <- duplicated(mapping$hgnc_symbol)
count.mat <- count.mat[!idx,]
mapping <- mapping[!idx,]
rownames(count.mat) <- mapping$hgnc_symbol


## Filter ##
sc.data <- SingleCellExperiment(assays = list(counts = as(count.mat, 'dgCMatrix')))
sc.data <- calculateQCMetrics(sc.data, 
                              use_spikes = FALSE, 
                              detection_limit = 0)
plotExprsFreqVsMean(sc.data)
keep.total <- sc.data$total_counts > 6e3 ## Remove cells based on total counts
keep.n <- sc.data$total_features_by_counts > 1000 ## Remove cells based on number of genes expressed
sc.data <- sc.data[,keep.total & keep.n]
n.exprs <- nexprs(sc.data, 
                  byrow = TRUE, 
                  detection_limit = 0)
keep_feature <- n.exprs > 250
sc.data <- sc.data[keep_feature,]

## Normalize ##
targets <- as.character(aggr.data$library_id[as.numeric(gsub('^.+-', 
                                                             replacement = '', 
                                                             colnames(sc.data)))])
sc.data <- computeSumFactors(sc.data, 
                             cluster = targets)
saveRDS(sc.data, file = '~/sc_trajectory/data/processedData.rds')
