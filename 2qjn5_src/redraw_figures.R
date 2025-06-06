library(viridis)
library(extrafont)
library(ape)

font_import()
loadfonts()

get.chull.for.cluster <- function(df,
                                  clustering,
                                  cl.num) {
    df[clustering[,2]==cl.num,][chull(df[clustering[,2]==cl.num,]),]
}

data.table <- read.csv('phoible_binarised.csv',
                       sep = '\t',
                       encoding = 'UTF-8',
                       h = F)
h <- scan('phoible_header.tsv', what = 'character', sep='\t',
          encoding = 'UTF-8')

colnames(data.table)  <- h
row.names(data.table) <- data.table[,1]

data.table <- data.table[,-1]
data.table <- data.table[ , colSums(data.table) >= 30, ]

inv.ids = strtoi(unlist(lapply(as.data.frame(strsplit(rownames(data.table), '_'))[2,],
                               as.character)))
data.table.no.austr <- data.table[ inv.ids < 2629, ]
austr               <- data.table[ inv.ids >= 2629, ]
cor.mat.no.austr    <- cor(data.table.no.austr)
cor.metric.no.austr <- as.dist(acos(cor.mat.no.austr))

umap.large   <- read.csv('d_acos_cor_fixed_large_cluster.csv')
umap.all     <- read.csv('d_acos_cor_fixed.csv')
dbscan.all   <- read.csv('dbscan_clustering_fixed.csv')
dbscan.large <- read.csv('dbscan_clustering_of_large_cluster_fixed.csv')

cols.all   <- viridis(n = length(unique(dbscan.all[,2])), end = 0.8)
cols.large <- viridis(n = length(unique(dbscan.large[,2])), end = 0.8)

tmp_dir <- 'pdfs'
# UMAP + polygons for the whole dataset
cairo_pdf(
    sprintf('%s/umap_cor_metric_poly_new.pdf', tmp_dir),
    width = 10,
    height = 16,
    # units = 'in',
    # res = 300
)
par(family = 'Doulos SIL')
plot(umap.all[,1]~umap.all[,2],
     type = 'n',
     xlab = '',
     ylab = '',
     xaxt = 'n',
     yaxt = 'n')
text(
    umap.all[,2],
    umap.all[,1],
    labels = colnames(data.table),
    #col = cols.all[dbscan.all[,2]],
    col = cols.all[dbscan.all[,2]+1],
    cex = 1.5
)
for (i in unique(dbscan.all[,2])) {
    if (i == 0)
        next
    polygon(
        get.chull.for.cluster(umap.all[,c(2,1)],
                              dbscan.all,
                              i),
        border = NA,
        col = rgb(.8, .8, .8, .3)
    )
}
dev.off()

# UMAP + polygons for the large cluster

cairo_pdf(
    sprintf('%s/large-cluster-poly_new.pdf', tmp_dir),
    width = 10,
    height = 16,
    # units = 'in',
    # res = 300
)
par(family = 'Doulos SIL')
plot(umap.large[,1]~umap.large[,2],
     type = 'n',
     xlab = '',
     ylab = '',
     xaxt = 'n',
     yaxt = 'n')
text(umap.large[,2],
     umap.large[,1],
     labels = colnames(data.table)[dbscan.all[,2] == 1],
     col = cols.large[dbscan.large[,2]+1],
     cex = 1.1
     )
for (i in unique(dbscan.large[,2])) {
    if (i == 0)
        next
    polygon(
        get.chull.for.cluster(umap.large[,c(2,1)],
                              dbscan.large,
                              i),
        border = NA,
        col = rgb(.8, .8, .8, .4)
    )
}
dev.off()


# Corrplot

library(corrplot)

palatalised.segs <- colnames(data.table)[grepl('??', colnames(data.table))]
palatalised.cor <- cor(data.table[,palatalised.segs])

cairo_pdf('pdfs/corrplot.pdf',
    width = 8,
    height = 8
    # units = 'in',
    # res = 300
    )
par(family = 'Doulos SIL')                                                             
corrplot(palatalised.cor,
         order = 'hclust',
         diag = F,
         tl.srt = 0,
         tl.col = 'black',
         tl.cex = 1.4,
         tl.offset = 1,
         col = 'black',
         bg = 'gold2',
         addrect = 3,
         cl.pos = 'n')
dev.off()

# Dengrograms

bci = c('j', 'w', 'm', 'r', 'l', 'p', 'k', 't', 'n')
png('img_new/basic_inventory_single_linkage.png',
    width = 16,
    height = 10,
    units = 'in',
    res = 300)
par(mfrow = c(1,2), family = "Doulos SIL")
plot(as.dendrogram(
                   hclust(
                          as.dist(acos(cor(data.table[,bci]))),
                          method = 'single'), cex = 2),
     horiz = T,
     type = 'triangle',
     yaxt = 'n',
     main = "With Australia",
     nodePar = list(lab.cex = 2, pch = c(NA, NA)),
     cex.main = 2)
plot(as.dendrogram(
                   hclust(
                          as.dist(acos(cor(data.table.no.austr[,bci]))),
                          method = 'single'), cex = 2),
    horiz = T,
    type = 'triangle',
    yaxt = 'n',
    main = "Without Australia",
    nodePar = list(lab.cex = 2, pch = c(NA, NA)),
    cex.main = 2)
dev.off()

first.ext <- c('s', 'b', 'h', 'ɡ', 'd', 'f', 'ɲ',
               't̠ʃ', 'ʔ', 'ʃ', 'z', 'd̠ʒ', 'v', 'ʒ', 'ʎ')
png('img_new/first_extension_single_linkage.png',
    width = 16,
    height = 10,
    units = 'in',
    res = 300)
par(family = 'Doulos SIL')
plot(as.phylo(hclust(as.dist(acos(cor(data.table[,first.ext]))),
                     method = 'single')),
     yaxt = 'n',
     main = "",
     type = 'unrooted',
     cex = 2)
dev.off()

gutturals <- c('ħ', 'ʕ', 'ɢ', 'qʰ', 'q', 'ʁ', 'χ')
png('img_new/gutturals_single_linkage.png',
    width = 16,
    height = 10,
    units = 'in',
    res = 300)
par(family = 'Doulos SIL')
plot(as.phylo(hclust(as.dist(acos(cor(data.table[,gutturals]))),
                     method = 'single')),
     yaxt = 'n',
     main = "",
     type = 'cladogram',
     cex = 2)
dev.off()

# Frequency plot

top.freqs = sort(colSums(data.table), decreasing = T)[1:40]
cex.val = 1.6
cairo_pdf('pdfs/top_freq_segments.pdf',
    width = 16,
    height = 10)
par(family = "Doulos SIL")
plot(top.freqs ~ seq_along(top.freqs),
    xlab = 'Frequency rank',
    ylab = 'Frequency',
    xlim = c(1, 40),
    type = 'n',
    xaxt = 'n',
    cex.lab = cex.val,
    cex.axis = cex.val)
axis(1, at = seq_along(top.freqs),
    labels = seq_along(top.freqs),
    cex.axis = 1.1)
lines(seq_along(top.freqs),
    top.freqs,
    lty = 2, col = 'grey')
text(seq_along(top.freqs), top.freqs,
    labels = names(top.freqs),
    cex = 2)
dev.off()
