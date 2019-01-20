library(kohonen)
library(Matrix)
library(cellrangerRkit)
library(Rtsne)
library(irlba)
library(scran)
library(scater)
library(edgeR)
library(biomaRt)

## Read Data ##
temp.path <- '~/SingleCellData/single_cell_analysis/Marusyk_Joint3Aggr_ALK_redo'
genome <- 'GRCh38'
count <- load_cellranger_matrix(temp.path)
aggr.data <- read.csv('~/SingleCellData/single_cell_analysis/Marusyk_Joint3Aggr_ALK_redo/outs/aggregation_csv.csv')
targets <- as.character(aggr.data$library_id[as.numeric(gsub('^.+-', replacement = '', colnames(count)))])
count.mat <- as.matrix(count)

## Filter Drug ##
select <- c('H3122-1', 'H3122-2', 'H3122-Alec-4h', 
            'H3122-Alec-48h', 'H3122-Alec-3W', 'H3122-erAlec-1')
idx <- which(targets %in% select)
count.mat <- count.mat[,idx]

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
keep.total <- sc.data$total_counts > 6e3 ## Remove cells based on total counts
keep.n <- sc.data$total_features_by_counts > 1200 ## Remove cells based on number of genes expressed
sc.data <- sc.data[,keep.total & keep.n]
n.exprs <- nexprs(sc.data, 
                  byrow = TRUE, 
                  detection_limit = 0)
keep_feature <- n.exprs > 100
sc.data <- sc.data[keep_feature,]

## Normalize ##
targets <- as.character(aggr.data$library_id[as.numeric(gsub('^.+-', 
                                                             replacement = '', 
                                                             colnames(sc.data)))])
sc.data <- computeSumFactors(sc.data, cluster = targets)
dge.data <- convertTo(sc.data, type = 'edgeR')
cpm.data <- edgeR::cpm(dge.data, normalized.lib.sizes = TRUE, log = TRUE)
cpm.data <- scale(t(cpm.data))

## Run SOM ##
message('Running SOM')
som.res <- som(cpm.data,
               grid = somgrid(xdim = 18,
                              ydim = 18,
                              topo = 'hexagonal',
                              neighbourhood.fct = 'gaussian',
                              toroidal = TRUE), 
               rlen = 3600,
               alpha = c(1, 0.01))
saveRDS(som.res,
        file = '~/SingleCellSomResAlectinib.rds')