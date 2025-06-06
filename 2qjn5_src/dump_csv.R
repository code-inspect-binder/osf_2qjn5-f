library(uwot)
library(dbscan)
library(corrplot)
library(xtable)
library(cluster)
library(phangorn)
library(Rtsne)

data.table <- read.csv('phoible_binarised.csv',
                       sep = '\t',
                       encoding = 'UTF-8',
                       h = F)
h <- scan('phoible_header.tsv', what = 'character', sep='\t', encoding = 'UTF-8')
colnames(data.table) <- h
row.names(data.table) <- data.table[,1]
data.table <- data.table[,-1]
data.table <- data.table[ , colSums(data.table) >= 30, ]
cor.mat <- cor(data.table)
cor.metric <- as.dist(acos(cor.mat))

## For separating Australian langs, which start at #
inv.ids = strtoi(unlist(lapply(as.data.frame(strsplit(rownames(data.table), '_'))[2,],
                               as.character)))

data.table.no.austr <- data.table[ inv.ids < 2629, ]
austr <- data.table[ inv.ids >= 2629, ]
cor.mat.no.austr <- cor(data.table.no.austr)
cor.metric.no.austr <- as.dist(acos(cor.mat.no.austr))

d.cor.umap <- umap(cor.metric,
                   min_dist = 1,
                   n_neighbors = 3,
                   )
## Use kNNdistplot(d.cor.umap, k = 3) to determine the value for eps
d.cor.dbscan <- dbscan(d.cor.umap, minPts = 3, eps = 1.8)


data.table.large.cluster <- data.table[,d.cor.dbscan$cluster == 1]
large.cluster.cor.metric <- as.dist(acos(cor(data.table.large.cluster)))

large.umap <- umap(
    large.cluster.cor.metric,
    min_dist = 1,
    n_neighbors = 3
)
large.dbscan <- dbscan(large.umap, minPts = 3, eps = 1)
large.cols <- rainbow(length(unique(large.dbscan$cluster)), start = 2/6)

write.csv(
    d.cor.umap,
    'd_acos_cor_fixed.csv',
    row.names = F
)

write.csv(
    d.cor.dbscan$cluster,
    'dbscan_clustering_fixed.csv'
)

write.csv(
    large.umap,
    'd_acos_cor_fixed_large_cluster.csv',
    row.names = F
)

write.csv(
    large.dbscan$cluster,
    'dbscan_clustering_of_large_cluster_fixed.csv'
)
