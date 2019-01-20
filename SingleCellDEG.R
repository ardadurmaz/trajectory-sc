library(Matrix)
library(cellrangerRkit)
library(ggplot2)
library(scran)
library(scater)
library(edgeR)
library(limma)
library(biomaRt)

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

## Subset Drug ##
select <- c('H3122-1', 'H3122-2', 'H3122-Alec-4h', 'H3122-Alec-48h', 'H3122-Alec-3W', 'H3122-erAlec-1', 'H3122-GFP-erAlec')
select <- c('H3122-1', 'H3122-2', 'H3122-erLor', 'H3122-erLor-1', 'H3122-Lor-3W', 'H3122-Lor-48h')
select <- c('H3122-1', 'H3122-2', 'H3122-Cz-48h', 'H3122-Cz-3W', 'H3122-erCz', 'H3122-erCz-3')

idx <- which(targets %in% select)
local.targets <- targets[idx]
local.count.mat <- count.mat[,idx]

## Filter Genes ##
sc.data <- SingleCellExperiment(assays = list(counts = as(local.count.mat, 'dgCMatrix')))
sc.data <- calculateQCMetrics(sc.data, 
                              use_spikes = FALSE, 
                              detection_limit = 0)
keep.total <- sc.data$total_counts > 6e3 ## Remove cells based on total counts
keep.n <- sc.data$total_features_by_counts > 1000 ## Remove cells based on number of genes expressed
sc.data <- sc.data[,keep.total & keep.n]
n.exprs <- nexprs(sc.data, byrow = TRUE, detection_limit = 0)
keep_feature <- n.exprs > 100
sc.data <- sc.data[keep_feature,]
local.targets <- as.character(aggr.data$library_id[as.numeric(gsub('^.+-', replacement = '', colnames(sc.data)))])

## Normalize Library Sizes ##
sc.data <- computeSumFactors(sc.data, cluster = local.targets, positive = TRUE)
dge.data <- convertTo(sc.data, type = 'edgeR')
local.targets <- gsub(pattern = '-', replacement = '', as.character(aggr.data$library_id[as.numeric(gsub('^.+-', replacement = '', colnames(dge.data)))]))
design.mat <- model.matrix(~ 0 + local.targets)

## Comparisons ##
all.combn <- combn(colnames(design.mat), 2)
cm.controls <- matrix(0, nrow = ncol(design.mat), ncol = ncol(all.combn), dimnames = list(colnames(design.mat), apply(all.combn, 2, function(x){paste(x, collapse = '_')})))

for(i in 1:ncol(cm.controls)){
  local.comp <- unlist(strsplit(colnames(cm.controls)[i], split = '_'))
  cm.controls[which(rownames(cm.controls) == local.comp[1]), i] <- 1
  cm.controls[which(rownames(cm.controls) == local.comp[2]), i] <- -1
}

## Fit ##
dge.data <- estimateDisp(dge.data, design = design.mat, robust = TRUE)
fit <- glmQLFit(dge.data, design.mat, robust = TRUE)

require(parallel)
cl <- makeCluster(12)
clusterExport(cl, varlist = c('fit', 'cm.controls'))
comp.res <- parSapply(cl, 1:ncol(cm.controls), simplify = FALSE, function(i){
  require(edgeR)
  tr <- glmTreat(fit, contrast = cm.controls[,i], lfc = log2(1.4))
  res.treat <- topTags(tr, adjust.method = 'bonferroni', n = Inf, sort.by = 'p.value')$table
  return(res.treat)
})
stopCluster(cl)
names(comp.res) <- colnames(cm.controls)

all.res <- list('Comparisons' = comp.res,
                'Data' = dge.data)
saveRDS(all.res, file = 'results/CrizatinibDEG.rds')
