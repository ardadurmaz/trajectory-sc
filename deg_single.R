library(Matrix)
library(cellrangerRkit)
library(ggplot2)
library(Rtsne)
library(irlba)
library(scran)
library(scater)
library(edgeR)
library(limma)

sc.data <- readRDS('~/single_cell_deg/processedData.rds')
aggr.data <- read.csv('~/SingleCellData/single_cell_analysis/Marusyk_Joint3Aggr_ALK_redo/outs/aggregation_csv.csv')
targets <- as.character(aggr.data$library_id[as.numeric(gsub('^.+-', 
                                                             replacement = '', 
                                                             colnames(sc.data)))])


## Cluster ##
dge.data <- convertTo(sc.data, type = 'edgeR')
cpm.data <- edgeR::cpm(dge.data, normalized.lib.sizes = TRUE, log = TRUE)
cpm.data <- scale(t(cpm.data))

colors.manual <- c('antiquewhite3','antiquewhite4',
                   'aquamarine', 'aquamarine3', 'aquamarine4',
                   'azure4', 'beige', 
                   'blueviolet', 'burlywood4',
                   'cadetblue', 'cadetblue4', 'chartreuse', 
                   'chartreuse4', 'chocolate', 'cornflowerblue', 
                   'cyan', 'cyan4', 'darkgoldenrod1', 
                   'darkorchid4', 'darkred','hotpink1', 
                   'darksalmon', 'darkslategrey', 'greenyellow',
                   'lightpink4', 'seashell4', 'blue4',
                   'yellow', 'tan', 'slateblue4')


## Alectinib ##
cpm.data.local <- cpm.data[grep(pattern = 'alec|^h3122-1$|^h3122-2$',
                                ignore.case = TRUE,
                                targets),]

tsne.res <- Rtsne(cpm.data.local, 
                  perplexity = 30, 
                  theta = 0.2, 
                  partial_pca = TRUE, 
                  pca_center = FALSE, 
                  pca_scale = FALSE, 
                  normalize = FALSE,
                  num_threads = 16)

plot.data <- data.frame('Coordinate.1' = tsne.res$Y[,1],
                        'Coordinate.2' = tsne.res$Y[,2],
                        'Type' = as.character(aggr.data$library_id[as.numeric(gsub('^.+-', 
                                                                                   replacement = '', 
                                                                                   rownames(cpm.data.local)))]))
plot.colors <- setNames(colors.manual[1:length(unique(plot.data$Type))], nm = unique(plot.data$Type))

ggplot(plot.data, aes(x = Coordinate.1, 
                      y = Coordinate.2)) +
  geom_point(aes(color = Type), shape = 20) +
  theme_minimal() +
  theme(legend.key.size = unit(0.8, 'cm')) +
  stat_ellipse(aes(color = Type), alpha = 0.6) +
  guides(color = guide_legend(override.aes = list(size = 5)),
         fill = FALSE) +
  scale_color_manual(values = plot.colors)

## Create Dendrogram ##
require(fastcluster)
require(dendextend)


hc <- fastcluster::hclust(dist(tsne.res$Y, method = 'euclidean'), method = 'average')
dend <- as.dendrogram(hc)
local.order <- labels(dend)

dend %>% 
  set('branches_k_color', k = 15) %>%
  set('labels_cex', 0.01) %>%
  plot(axes = FALSE, main = 'Hierarchical Clustering')
colored_bars(colors = plot.colors[match(targets, table = names(plot.colors))], dend = dend, rowLabels = "Group")


## K-Means Cluster ##
k.cluster <- kmeans(tsne.res$Y, 
                    centers = 20,
                    iter.max = 1000)
plot.data$Cluster <- paste('Cluster-', k.cluster$cluster, sep = '')

ggplot(plot.data, aes(x = Coordinate.1, 
                      y = Coordinate.2)) +
  geom_point(aes(color = Type), shape = 20) +
  theme_minimal() +
  scale_color_manual(values = plot.colors) +
  theme(legend.key.size = unit(0.8, 'cm')) +
  guides(color = guide_legend(override.aes = list(size = 5))) +
  stat_ellipse(aes(fill = Cluster), alpha = 0.6)

## Density Clustering ##
require(densityClust)
dist.mat <- dist(tsne.res$Y, method = 'euclidean')
density.clust.res <- densityClust(distance = dist.mat, 
                                  gaussian = TRUE)
cluster.res <- findClusters(density.clust.res)





all.sig.genes <- Reduce(union, lapply(all.deg.res, function(x){rownames(x)[abs(x$logFC) > 1 & x$adj.P.Val < 0.01]}))
all.sig.mat <- sapply(all.deg.res, function(x){
  return(x$logFC[match(all.sig.genes, table = rownames(x))])
})
rownames(all.sig.mat) <- all.sig.genes
all.sig.genes.mapped <- mapping$hgnc_symbol[match(all.sig.genes,
                                                  table = mapping$ensembl_gene_id)]
idx <- duplicated(all.sig.genes.mapped)
all.sig.mat <- all.sig.mat[!idx,]
rownames(all.sig.mat) <- all.sig.genes.mapped[!idx]


## Plot Heatmap ##
require(pheatmap)

drug.annot <- data.frame('Drug.a' = c('erAlec1', 'GFPerAlec', 'erAlec1', 
                                      'GFPerAlec', 'erAlec1', 'GFPerAlec',
                                      'erAlec1', 'GFPerAlec', 'erLor',
                                      'erLor1', 'erLor', 'erLor1'),
                         'Drug.b' = c('erLor', 'erLor', 'erLor1', 
                                      'erLor1','erCz', 'erCz', 
                                      'erCz3', 'erCz3', 'erCz',
                                      'erCz', 'erCz3', 'erCz3'),
                         row.names = colnames(all.sig.mat))
column.annotation <- data.frame('Control' = c('1', '2', '1', '2', '1', '2','1', '2', '1', '2'),
                                'Resistance' = c('Alec-1', 'Alec-1', 
                                           'Lor-1', 'Lor-1', 'Lor-0', 'Lor-0',
                                           'Criz-0', 'Criz-0', 'Criz-3', 'Criz-3'),
                                row.names = colnames(all.sig.mat))
temp.row.annot <- data.frame('Degree' = log10(colSums(ppi)[match(rownames(all.sig.mat), table = colnames(ppi))]),
                             row.names = rownames(all.sig.mat))

pheatmap(t(all.sig.mat), 
         fontsize_row = 8, 
         fontsize_col = 2, 
         annotation_row = column.annotation, 
         annotation_col = temp.row.annot,
         show_colnames = FALSE, 
         show_rownames = FALSE,
         scale = 'none',
         main = 'Differential Gene Expression (Controls)')


## Plot Degree Distribution ##
require(ggplot2)
require(reshape2)

temp.res <- sapply(1:length(all.deg.res), simplify = FALSE, function(i){
  x <- all.deg.res[[i]]
  local.res <- subset(x, abs(x$logFC) > 0.48 & x$adj.P.Val < 0.05)
  local.res <- setNames(local.res$logFC, rownames(local.res))
  names(local.res) <- mapping$hgnc_symbol[match(names(local.res), mapping$ensembl_gene_id)]
  local.res <- sort(local.res, decreasing = TRUE)
  local.res <- na.omit(local.res[!duplicated(names(local.res))])
  local.res <- local.res[!is.na(names(local.res))]
  return(data.frame('Gene' = names(local.res),
                    'Type' = ifelse(local.res < 0, 'Down', 'Up'),
                    'Test' = rep(names(all.deg.res)[i], length(local.res)),
                    'Degree' = colSums(ppi)[match(names(local.res), table = colnames(ppi))],
                    stringsAsFactors = FALSE))
})
temp.res <- do.call('rbind', lapply(temp.res, function(x){return(x)}))
temp.res <- na.omit(temp.res)
rownames(temp.res) <- 1:nrow(temp.res)

ggplot(data = temp.res,
       aes(x = Test, y = log10(Degree))) +
  geom_boxplot(position = 'dodge', notch = TRUE, aes(fill = Type)) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))


pca.res <- prcomp(t(all.sig.mat), scale. = FALSE, center = FALSE)
plot.data <- data.frame('Coordinate.1' = pca.res$x[,1],
                        'Coordinate.2' = pca.res$x[,2],
                        'Coordinate.3' = pca.res$x[,3],
                        'Label' = rownames(pca.res$x))
ggplot(plot.data, 
       aes(x = Coordinate.1, 
           y = Coordinate.2, 
           color = Label)) +
  geom_point(size = 6) +
  ggtitle('Gene Level DGE') +
  theme_minimal() +
  scale_color_brewer(palette = 'Paired')

reactome <- read.table('~/Ensembl2Reactome_All_Levels.txt',
                       header = FALSE,
                       sep = '\t',
                       comment.char = '',
                       stringsAsFactors = FALSE,
                       quote = '')
reactome <- subset(reactome, reactome$V6 == 'Homo sapiens')
rownames(reactome) <- 1:nrow(reactome)
reactome.list <- split(reactome$V1, reactome$V4)
reactome.list <- lapply(reactome.list, function(x){
  return(unique(x[x %in% rownames(dge.data)]))
})
reactome.list <- reactome.list[sapply(reactome.list, length) >= 8]
reactome.list <- reactome.list[sapply(reactome.list, length) < 400]


## Pathway Hierarchy ##
hier <- read.table('~/resources/ReactomePathwaysRelation.txt',
                   header = FALSE,
                   sep = '\t',
                   stringsAsFactors = FALSE)
hier <- hier[intersect(grep(pattern = 'hsa', ignore.case = TRUE, hier$V1),
                       grep(pattern = 'hsa', ignore.case = TRUE, hier$V2)),]
root.nodes <- setdiff(unique(hier$V1),
                      unique(hier$V2))
non.root.nodes <- setdiff(union(hier$V1, hier$V2), root.nodes)

require(igraph)
local.graph <- simplify(graph_from_edgelist(as.matrix(hier), directed = TRUE))
temp.paths <- sapply(non.root.nodes,
                     simplify = FALSE,
                     function(v.x){
                       temp <- shortest_paths(local.graph, 
                                              v.x, 
                                              to = root.nodes, 
                                              mode = "in",
                                              weights = NULL, 
                                              output = "vpath",
                                              predecessors = FALSE, 
                                              inbound.edges = FALSE)$vpath
                       temp.size <- sapply(temp, function(i){V(local.graph)[i]$name})
                       names(temp.size) <- root.nodes
                       temp.size <- do.call('c', lapply(temp.size, function(x){return(length(x))}))
                       temp.size <- temp.size[temp.size != 0]
                       return(names(temp.size)[which.min(temp.size)])
                     })
root.mapping <- data.frame('root.path' = as.character(sapply(temp.paths, function(x){return(x)})),
                           'leaf.path' = names(temp.paths),
                           stringsAsFactors = FALSE)
root.mapping$root.path <- reactome$V4[match(root.mapping$root.path, table = reactome$V2)]
root.mapping$leaf.path <- reactome$V4[match(root.mapping$leaf.path, table = reactome$V2)]
root.mapping <- na.omit(root.mapping)
rownames(root.mapping) <- 1:nrow(root.mapping)

## Gene Set Enrichment Analysis ##
require(fgsea)
enr.res <- sapply(all.deg.res,
                  simplify = FALSE,
                  function(x){
                    local.res <- sort(setNames(x$logFC, nm = rownames(x)), decreasing = TRUE)
                    local.enr.res <- fgsea(pathways = reactome.list,
                                           stats = local.res,
                                           nperm = 100000,
                                           minSize = 8,
                                           nproc = 16)
                    return(setNames(local.enr.res$NES, local.enr.res$pathway))
                  })
enr.res.ft <- sapply(enr.res,
                     function(x){
                       return(x[match(names(reactome.list),
                                      table = names(x))])
                     })
reactome.pca.res <- prcomp(t(enr.res.ft), center = FALSE, scale. = FALSE)
plot.data <- data.frame('Coordinate.1' = reactome.pca.res$x[,1],
                        'Coordinate.2' = reactome.pca.res$x[,2],
                        'Coordinate.3' = reactome.pca.res$x[,3],
                        'Label' = rownames(reactome.pca.res$x))
ggplot(plot.data, 
       aes(x = Coordinate.2, 
           y = Coordinate.3, 
           color = Label)) +
  geom_point(size = 6) +
  ggtitle('Pathway Level DGE') +
  theme_minimal() +
  scale_color_brewer(palette = 'Paired')




enr.res <- sapply(all.deg.res,
                  simplify = FALSE,
                  function(local.genes){
                    local.res <- sapply(reactome.list, function(local.path){
                      test.mat <- matrix(c(sum(local.genes %in% local.path),
                                           sum(!(local.genes %in% local.path)),
                                           sum(rownames(dge.data) %in% local.path),
                                           sum(!(rownames(dge.data) %in% local.path))),
                                         byrow = FALSE,
                                         ncol = 2)
                      p.val <- fisher.test(test.mat, alternative = 'greater')$p.value
                      return(p.val)
                    })
                    local.res <- p.adjust(local.res, method = 'fdr')
                    return(local.res)
                  })
ft.data <- do.call('rbind', 
                   sapply(1:length(enr.res), 
                          simplify = FALSE, 
                          function(x){
                            return(data.frame('comparison' = rep(names(enr.res)[x], length(enr.res[[x]])),
                                              'pathway.id' = names(enr.res[[x]]),
                                              'root.pathway.id' = root.mapping$root.path[match(names(enr.res[[x]]), 
                                                                                               table = root.mapping$leaf.path)],
                                              'p.value' = enr.res[[x]],
                                              stringsAsFactors = FALSE))
                          }))

p.thr <- 0.01
ft.data <- subset(ft.data, ft.data$p.value < p.thr)
ft.data <- na.omit(ft.data)
rownames(ft.data) <- 1:nrow(ft.data)


require(ggplot2)
ggplot(data = ft.data,
       aes(x = comparison,
           y = pathway.id,
           fill = root.pathway.id)) +
  geom_point(aes(size = -log10(p.value),
                 color = root.pathway.id),
             shape = 20) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 4),
        axis.text.x = element_text(angle = 60, hjust = 1)) +
  guides(size = FALSE)







pheatmap(-log10(enr.res), fontsize_row = 4, fontsize_col = 8)



require(biomaRt)
ensembl <- useMart("ENSEMBL_MART_ENSEMBL", dataset = 'hsapiens_gene_ensembl')
mapping <- getBM(mart = ensembl,
                 attributes = c('ensembl_gene_id', 'hgnc_symbol'),
                 filters = 'ensembl_gene_id',
                 values = rownames(dge.data))

require(pheatmap)
temp <- edgeR::cpm(dge.data[all.res[,5] != 0, design.mat[,7] == 1 | design.mat[,1] == 1], noralized.lib.sizes = TRUE, log = TRUE)
pheatmap(temp, fontsize = 6)

