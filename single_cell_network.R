library(Matrix)
library(cellrangerRkit)
library(ggplot2)
library(Rtsne)
library(irlba)
library(scran)
library(scater)
library(edgeR)
library(limma)
library(pheatmap)

sc.data <- readRDS('~/single_cell_deg/processedDataStrict.rds')
aggr.data <- read.csv('~/SingleCellData/single_cell_analysis/Marusyk_Joint3Aggr_ALK_redo/outs/aggregation_csv.csv')
targets <- as.character(aggr.data$library_id[as.numeric(gsub('^.+-', 
                                                             replacement = '', 
                                                             colnames(sc.data)))])
dge.data <- convertTo(sc.data, type = 'edgeR')

## Differential Expression ##
local.targets <- gsub(pattern = '-', 
                      replacement = '', 
                      as.character(aggr.data$library_id[as.numeric(gsub('^.+-', 
                                                                        replacement = '', 
                                                                        colnames(dge.data)))]))
design.mat <- model.matrix(~ 0 + local.targets)
cm.controls <- makeContrasts(AlectinibG1C1 = local.targetsH3122erAlec1 - local.targetsH31221,
                             AlectinibG1C2 = local.targetsH3122erAlec1 - local.targetsH31222,
                             LorlatinibG1C1 = local.targetsH3122erLor1 - local.targetsH31221,
                             LorlatinibG1C2 = local.targetsH3122erLor1 - local.targetsH31222,
                             LorlatinibG0C1 = local.targetsH3122erLor - local.targetsH31221,
                             LorlatinibG0C2 = local.targetsH3122erLor - local.targetsH31222,
                             CrizatinibG0C1 = local.targetsH3122erCz - local.targetsH31221,
                             CrizatinibG0C2 = local.targetsH3122erCz - local.targetsH31222,
                             CrizatinibG3C1 = local.targetsH3122erCz3 - local.targetsH31221,
                             CrizatinibG3C2 = local.targetsH3122erCz3 - local.targetsH31222,
                             levels = design.mat)
voom.data <- voom(dge.data, 
                  design = design.mat, 
                  normalize.method = 'none',
                  span = 0.2,
                  plot = TRUE)
fit <- lmFit(voom.data, 
             design = design.mat, 
             method = "robust")
cfit <- contrasts.fit(fit, contrasts = cm.controls)
efit <- eBayes(cfit, robust = TRUE)
all.deg.res <- sapply(1:ncol(cm.controls), simplify = FALSE, function(i){
  res <- topTable(efit, 
                  coef = i, 
                  adjust.method = 'bonferroni', 
                  number = Inf, 
                  lfc = 0, 
                  p.value = 1)
  return(res)
})
names(all.deg.res) <- colnames(cm.controls)



## Generate Mapping ##
require(biomaRt)
ensembl <- useMart('ENSEMBL_MART_ENSEMBL', 
                   dataset = 'hsapiens_gene_ensembl')
mapping <- getBM(mart = ensembl,
                 attributes = c('ensembl_gene_id', 'hgnc_symbol'),
                 filters = 'ensembl_gene_id',
                 values = rownames(all.deg.res[[1]]))

## Random Walk Restarts ##
ppi <- readRDS('~/resources/stringdb.medconf.rds')
ppi.n <- ppi %*% Matrix::Diagonal(x = Matrix::colSums(ppi)^-1)
all.seed.genes <- sapply(all.deg.res, simplify = FALSE, function(x){
  local.genes <- rownames(x)[abs(x$logFC) > log2(1.4) & x$adj.P.Val < 0.01]
  local.genes <- unique(na.omit(mapping$hgnc_symbol[match(local.genes, 
                                                          table = mapping$ensembl_gene_id)]))
  local.genes <- local.genes[local.genes != '']
  return(local.genes)
})

p0 <- Matrix(sapply(all.seed.genes, function(x){return(ifelse(colnames(ppi) %in% x, 1, 0))}))
p0 <- p0 %*% Matrix::Diagonal(x = Matrix::colSums(p0)^-1)
pt <- Matrix(1/ncol(ppi.n), ncol = ncol(p0), nrow = nrow(ppi.n))

r <- 0.75
delta <- 1
count <- 1
while(delta > 1e-16 && count < 100){
  px <- (1-r) * ppi.n %*% pt + r * p0
  delta <- sum(abs(px - pt))
  pt <- px
  count <- count + 1
}
colnames(pt) <- colnames(p0)
pt <- pt[!apply(pt, 1, function(x){all(x < 1e-4)}),]
pt <- as.matrix(pt)
pt.ft <- scale(t(log2(pt)))
rwr.pca <- prcomp(pt.ft, center = FALSE, scale. = FALSE)

plot.data <- data.frame('Coordinate.1' = rwr.pca$x[,1],
                        'Coordinate.2' = rwr.pca$x[,2],
                        'Coordinate.3' = rwr.pca$x[,3],
                        'Coordinate.4' = rwr.pca$x[,4],
                        'Label' = rownames(rwr.pca$x))
ggplot(plot.data, 
       aes(x = Coordinate.1, 
           y = Coordinate.2, 
           color = Label)) +
  geom_point(size = 6) +
  ggtitle('PPI Level DGE') +
  theme_minimal() +
  scale_color_brewer(palette = 'Paired')

temp.col.annot <- data.frame('Degree' = log2(Matrix::colSums(ppi)[match(colnames(pt.ft), 
                                                                         table = colnames(ppi))]),
                             row.names = colnames(pt.ft))
pheatmap(pt.ft, 
         scale = 'none', 
         fontsize_col = 0.5,
         annotation_col = temp.col.annot, 
         main = 'Biased Random Walk (DGE)')


## Plot LFC ##
all.genes <- Reduce(union, all.seed.genes)
lfc.ft <- sapply(all.deg.res, function(x){
  local.res <- setNames(x$logFC, 
                        nm = mapping$hgnc_symbol[match(rownames(x), table = mapping$ensembl_gene_id)])
  return(local.res[match(all.genes, table = names(local.res))])
})
lfc.ft <- scale(t(lfc.ft))
pca.res <- prcomp(lfc.ft, center = FALSE, scale. = FALSE)
plot.data <- data.frame('Coordinate.1' = pca.res$x[,1],
                        'Coordinate.2' = pca.res$x[,2],
                        'Coordinate.3' = pca.res$x[,3],
                        'Coordinate.4' = pca.res$x[,4],
                        'Label' = rownames(pca.res$x))
ggplot(plot.data, 
       aes(x = Coordinate.1, 
           y = Coordinate.2, 
           color = Label)) +
  geom_point(size = 6) +
  ggtitle('DGE') +
  theme_minimal() +
  scale_color_brewer(palette = 'Paired')



## Plot Network ##
require(network)
require(sna)
require(ggplot2)
require(ggnet)

## Control-1 ##
pt.local <- pt[,grep(pattern = 'C1$', names(all.seed.genes))]
pt.local <- pt.local[!apply(pt.local, 1, function(x){all(x < 1e-3)}),]
pt.local <- as.matrix(pt.local)
pt.local.ft <- scale(t(log10(pt.local)))
temp.idx <- colnames(ppi) %in% colnames(pt.local.ft)
ppi.local <- ppi[temp.idx,
                 temp.idx]

graph <- network(as.matrix(ppi.local), 
                 directed = FALSE, 
                 hyper = FALSE, 
                 loops = FALSE, 
                 bipartite = FALSE, 
                 matrix.type = 'adjacency')

## Get Coordinates ##
x <- gplot.layout.fruchtermanreingold(graph, NULL)
graph %v% 'x' <- x[,1]
graph %v% 'y' <- x[,2]

## Add logFC Colors ##
require(RColorBrewer)
breakList <- seq(-4, 4, 0.1)
man.colors <- colorRampPalette(c('navyblue', 'white', 'firebrick'))(length(breakList))

ret.colors <- function(x = NULL, col = man.colors, idx = breakList){
  mapped.idx <- sapply(x, function(t){
    if(is.na(t))
      return(NA)
    local.idx <- 1
    while(local.idx < length(breakList) && t > breakList[local.idx]){
      local.idx <- local.idx + 1
    }
    return(local.idx)
  })
  temp <- man.colors[mapped.idx]
  temp[is.na(temp)] <- 'gray65'
  return(temp)
}

## Add Colors ##
graph %v% 'Alectinib-G1' <- all.deg.res[['AlectinibG1C1']]$logFC[match(graph %v% 'vertex.names',
                                                                       table = mapping$hgnc_symbol[match(rownames(all.deg.res[[1]]),
                                                                                                         table = mapping$ensembl_gene_id)])]
graph %v% 'Lorlatinib-G1' <- all.deg.res[['LorlatinibG1C1']]$logFC[match(graph %v% 'vertex.names',
                                                                        table = mapping$hgnc_symbol[match(rownames(all.deg.res[[1]]),
                                                                                                          table = mapping$ensembl_gene_id)])]
graph %v% 'Lorlatinib-G0' <- all.deg.res[['LorlatinibG0C1']]$logFC[match(graph %v% 'vertex.names',
                                                                         table = mapping$hgnc_symbol[match(rownames(all.deg.res[[1]]),
                                                                                                           table = mapping$ensembl_gene_id)])]
graph %v% 'Crizatinib-G0' <- all.deg.res[['CrizatinibG0C1']]$logFC[match(graph %v% 'vertex.names',
                                                                         table = mapping$hgnc_symbol[match(rownames(all.deg.res[[1]]),
                                                                                                           table = mapping$ensembl_gene_id)])]
graph %v% 'Crizatinib-G3' <- all.deg.res[['CrizatinibG3C1']]$logFC[match(graph %v% 'vertex.names',
                                                                         table = mapping$hgnc_symbol[match(rownames(all.deg.res[[1]]),
                                                                                                           table = mapping$ensembl_gene_id)])]


## Add Filters ##
graph %v% 'Filter_Alectinib-G1' <- ifelse(pt.local[,'AlectinibG1C1'] >= 1e-3, 1, NA)[match(graph %v% 'vertex.names',
                                                                                           table = rownames(pt.local))]
graph %v% 'Filter_Lorlatinib-G1' <- ifelse(pt.local[,'LorlatinibG1C1'] >= 1e-3, 1, NA)[match(graph %v% 'vertex.names',
                                                                                           table = rownames(pt.local))]
graph %v% 'Filter_Lorlatinib-G0' <- ifelse(pt.local[,'LorlatinibG0C1'] >= 1e-3, 1, NA)[match(graph %v% 'vertex.names',
                                                                                             table = rownames(pt.local))]
graph %v% 'Filter_Crizatinib-G3' <- ifelse(pt.local[,'CrizatinibG3C1'] >= 1e-3, 1, NA)[match(graph %v% 'vertex.names',
                                                                                             table = rownames(pt.local))]
graph %v% 'Filter_Crizatinib-G0' <- ifelse(pt.local[,'CrizatinibG0C1'] >= 1e-3, 1, NA)[match(graph %v% 'vertex.names',
                                                                                             table = rownames(pt.local))]



## Plot ##
ggnet2(graph,
       label = TRUE,
       label.size = 4,
       size = 8,
       color = ret.colors(x = graph %v% 'Alectinib-G1', col = man.colors, idx = breakList),
       na.rm = 'Filter_Alectinib-G1',
       alpha = 0.9,
       legend.position = 'bottom',
       legend.size = 12,
       edge.size = 0.1)
