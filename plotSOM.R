library(ggplot2)
library(ggforce)
library(kohonen)
library(igraph)

## Read SOM Result ##
som.res <- readRDS('data/SingleCellSomRes.rds')

## Read Group Info ##
aggr.data <- read.csv('~/SingleCellData/single_cell_analysis/Marusyk_Joint3Aggr_ALK_redo/outs/aggregation_csv.csv')
targets <- as.character(aggr.data$library_id[as.numeric(gsub('^.+-', 
                                                             replacement = '', 
                                                             rownames(som.res$data[[1]])))])
my.data <- som.res$codes
clusters <- kmeans(my.data[[1]], centers = 15)

## Plot ##
grid.vals <- as.data.frame(som.res$grid$pts)
grid.vals$cluster.id <- factor(paste('Cluster-', clusters$cluster, sep = ''))
point.vals <- data.frame('id' = som.res$unit.classif,
                         'dist' = som.res$distances,
                         'Type' = targets)
point.vals$x <- grid.vals$x[point.vals$id]
point.vals$y <- grid.vals$y[point.vals$id]

ggplot(data = grid.vals, aes(x0 = x, y0 = y)) +
  geom_circle(aes(r = 0.5), color = 'grey80', n = 16) +
  theme_minimal() +
  ggtitle('SOM') +
  theme(panel.background = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        plot.title = element_text(hjust = 0.5, face = 'bold', size = 16)) +
  geom_jitter(data = point.vals, aes(x, y, col = Type)) +
  guides(color = guide_legend(override.aes = list(size = 4)))


## U-Matrix ##
all.dist <- as.matrix(object.distances(som.res, type = 'codes'))

## MST ##
require(igraph)
mst.res <- mst(graph_from_adjacency_matrix(all.dist, 
                                           mode = 'undirected', 
                                           weighted = TRUE, 
                                           diag = FALSE))
plot(mst.res,
     vertex.size = 6,
     vertex.label = NA,
     vertex.color = colors.man[match(grid.vals$cluster.id, table = names(colors.man))])

legend('topright',
       legend = names(colors.man), 
       fill = colors.man)
