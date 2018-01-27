# |----------------------------------------------------------------------------------|
# | Project: Statistics for genomics data with no replicas                           |
# | Script: Analysis of db/db microarray data (public)                               |
# | Authors: Davit Sargsyan                                                          |
# | Created: 12/02/2017                                                              |
# |                                                                                  |
# | Dataset1: db/db model of obesity-induced type 2 diabetes: neuroretina            |
# | Data source: https://www.ncbi.nlm.nih.gov/sites/GDSbrowser#details               |
# | Dataset Number: GDS5041                                                          |
# | Analysis of neuroretina from 8-week C57BLKsJ-db/db diabetic mice.                |
# | These leptin-receptor-deficient db/db mice have been shown to exhibit            |
# | early features of diabetic retinopathy (DR). Results provide insight into        |
# | the molecular mechanisms underlying early diabetic retinal neurodegeneration.    |
# |----------------------------------------------------------------------------------|
# Header----
# Save consol output to a log file
# sink(file = "tmp/log_manuscript_dbdb_ncbi_analysis_v1.txt")

# source("https://bioconductor.org/biocLite.R")
# biocLite("pd.mogene.1.1.st.v1")
# biocLite("GenomicRanges")

require(data.table)
require(ggplot2)
require(oligo)
require(hugene20sttranscriptcluster.db)
require(affycoretools)
require(limma)
# require(edgeR)
# require(DNAMR) # Javier's package

# Part I: Data----
HOME <- paste(getwd(),
              "/data/dbdb_ncbi",
              sep = "")
lData <- dir(HOME)
lData

l1 <- oligo::read.celfiles(filenames = paste(HOME,
                                             lData,
                                             sep = "/"))
l1

# Annotation
l1@annotation
# ?pd.mogene.1.1.st.v1

# Sample names
trt.names <- c("DBDB 1",
               "DBDB 2",
               "DBDB 3",
               "DBDB 4",
               "DB+ 1",
               "DB+ 2",
               "DB+ 3",
               "DB+ 4")
sampleNames(l1) <- trt.names
l1

# Pseudo-image
for (i in 1:length(sampleNames(l1))) {
  tiff(filename = paste("tmp/",
                        sampleNames(l1)[i],
                        ".tiff",
                        sep = ""),
       height = 5,
       width = 5,
       units = 'in',
       res = 300,
       compression = "lzw+p")
  
  oligo::image(l1,
               which = i,
               transfo = rank)
  
  graphics.off()
}

# Raw data distribution----
dt0 <- data.table(l1@assayData$exprs)
dt0

set.seed(2017)
tmp <- melt.data.table(dt0[sample(x = 1:nrow(dt0),
                                  size = 10000)])

p1 <- ggplot(data = tmp) +
  geom_boxplot(aes(x = variable,
                   y = log2(value),
                   fill = variable)) +
  scale_x_discrete("Group") + 
  scale_y_continuous("Log2(Raw Gene Expressons)") + 
  ggtitle("Random Sample of 10,000 Raw Values") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none",
        axis.text.x = element_text(angle = 45,
                                   hjust = 1))
print(p1)

# Normalized data
genePS <- rma(l1,
              target = "core")
genePS
dt1 <- exprs(genePS)
dt1 <- data.table(PROBEID = rownames(dt1),
                  dt1)
dt1

# Annotate----
dt2 <- annotateEset(genePS,
                    pd.mogene.1.1.st.v1)
anno <- data.table(dt2@featureData@data)
anno$PROBEID <- as.character(anno$PROBEID)
anno
summary(anno$SYMBOL)
length(unique(anno$SYMBOL))
# 24,214 genes

# Remove NAs
anno <- droplevels(subset(anno,
                          !is.na(SYMBOL)))

# Merge expressions and annotation----
dt1 <- merge(anno, 
             dt1,
             by = "PROBEID")
dt1

dt1$SYMBOL <- as.character(dt1$SYMBOL)
length(unique(dt1$SYMBOL))

# Normilazed and annotated data distribution----
tmp <- melt.data.table(dt1,
                       measure.vars = 5:12)
tmp

p2 <- ggplot(data = tmp) +
  geom_boxplot(aes(x = variable,
                   y = value,
                   fill = variable)) +
  scale_x_discrete("Group") + 
  scale_y_continuous("Normalized Gene Expressions") + 
  ggtitle("All Normilazed Annotated Genes") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none",
        axis.text.x = element_text(angle = 45,
                                   hjust = 1))
print(p2)

gridExtra::grid.arrange(p1, p2, nrow = 1)

tiff(filename = "tmp/dbdb_raw_vs_rma.norm.tiff",
     height = 5,
     width = 10,
     units = 'in',
     res = 300,
     compression = "lzw+p")
gridExtra::grid.arrange(p1, p2, nrow = 1)
graphics.off()

# How do samples compare?
t1 <- cor(dt1[, 5:12])
t1
plot(dt1$`DBDB 1` ~ dt1$`DB+ 1`)

plot(log2(dt1$`DBDB 1`) ~ log2(dt1$`DB+ 1`))

# edgeR----
dt1[, 5:12]

dt2 <- DGEList(counts = as.matrix(dt1[, 5:12]),
               group = rep(c("DBDB",
                             "DB+"),
                             each = 4))
dt2
# Data normalization----
dt2 <- calcNormFactors(dt2)
dt2

# Explore the data----
# a. Multidimensional scaling plot of distances between 
#    gene expression profiles. 
#    This plots samples on a two-dimensional scatterplot so 
#    that distances on the plot approximate the typical 
#    log2 fold changes between the samples.
plotMDS(dt2)
#    Options: method = bcv: ??
plotMDS(dt2,
        method = "bcv")

# PCA----
m.pca <- prcomp(t(dt2$counts),
                center = TRUE,
                scale = TRUE)
summary(m.pca)
plot(m.pca)

# Keep only the most important cytokines (Javier)----
# Select PC-s to pliot (PC1 & PC2)
choices <- 1:2

nobs.factor <- sqrt(nrow(m.pca$x) - 1)
d <- m.pca$sdev
u <- m.pca$x
v <- m.pca$rotation

# Scores
df.u <- data.frame(u[, choices])
# Add grouping variables
df.u$smpl <- rownames(df.u)
df.u$grp <- rep(c("DBDB",
                  "DB+"),
                  each = 4)
df.u

# Directions
df.v <- as.data.frame(v[, choices])

# Annotate
df.v$Gene <- dt1$SYMBOL

# Separate top variables (largest arrows)
df.v$lgth <- sqrt(df.v$PC1^2 + df.v$PC2^2)
df.v <- df.v[order(df.v$lgth,
                   decreasing = TRUE), ]

df.v$Gene <- factor(df.v$Gene,
                    levels = unique(df.v$Gene))
df.v

p1 <- ggplot(df.v[1:200, ]) +
  geom_bar(aes(x = Gene,
               y = lgth),
           stat = "identity") +
  ggtitle("Genes") +
  scale_x_discrete("Gene") +
  scale_y_continuous("Axis Length") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 45,
                                   hjust = 1))
p1
tiff(filename = "tmp/mes13_methylseq_pca_axis.tiff",
     height = 5,
     width = 20,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p1)
graphics.off()

# Axis labels
u.axis.labs <- paste(colnames(df.v)[1:2], 
                     sprintf('(%0.1f%% explained var.)', 
                             100*m.pca$sdev[choices]^2/sum(m.pca$sdev^2)))

# Which cytokines to display
# var.keep.ndx <- 1:10
# # Genes with largest PC1 and PC1
# var.keep.ndx <- unique(c(which(order(abs(df.v$PC1)) %in% 1:10),
#                          which(order(abs(df.v$PC2)) %in% 1:10)))
# NEW(DS  01/24/2018): genes that are most influential in PC1 direction
dff <- df.v$PC1/df.v$PC2
var.keep.ndx <- which(order(dff) %in% c(1:10,
                                        (length(dff)-10):length(dff)))


xx <- 10000

p2 <- ggplot(data = df.v[var.keep.ndx,], 
             aes(x = PC1,
                 y = PC2)) +
  coord_equal() +
  geom_segment(aes(x = 0,
                   y = 0, 
                   xend = 0.95*xx*PC1,
                   yend = 0.95*xx*PC2),
               arrow = arrow(length = unit(1/2, 'picas')), 
               color = "black",
               size = 0.5) +
  geom_text(aes(x = xx*PC1,
                y = xx*PC2,
                label = df.v$Gene[var.keep.ndx]),
            size = 3,
            angle = 30,
            hjust = 0.5) +
  geom_point(data = df.u,
             aes(fill = smpl),
             shape = 21,
             size = 3,
             alpha = 0.8) +
  scale_x_continuous(u.axis.labs[1]) +
  scale_y_continuous(u.axis.labs[2]) +
  scale_fill_discrete(name = "Sample") +
  ggtitle("Biplot of genes with 20 largest and\n20 smallest PC1/PC2 ratios\nDBDB Microarray") +
  coord_equal(ratio = 1) +
  theme(plot.title = element_text(hjust = 0.5,
                                  size = 20))
p2

tiff(filename = "tmp/biplot_cytokines.tiff",
     height = 8,
     width = 8,
     units = 'in',
     res = 300,
     compression = "lzw+p")
print(p2)
graphics.off()

# Estimating dispersion----
dt2 <- estimateCommonDisp(dt2, verbose=T)
# dt2 <- estimateTagwiseDisp(dt2) ???
dt2
plotBCV(dt2)
# # Not possible without replicates! Set it manually
# dt2$common.dispersion <- 0.05
# plotBCV(dt2)

# Compute genewise exact tests for differences in the means 
# between two groups of negative-binomially distributed counts.
t1 <- edgeR::exactTest(dt2,
                       pair = c(1, 2))
# t1 <- edgeR::exactTest(dt2,
#                        pair = c(2, 3))
t1
topTags(t1, 
        n = 10)

sum1 <- decideTestsDGE(t1, 
                       adjust.method = "none",
                       p.value = 0.1)
summary(sum1)

plotSmear(t1)
abline(h = c(-0.2, 0.2), 
       col = "blue",
       lty = 2)

# GLM----
design.mat <- matrix(0, nrow = 8, ncol = 2)
design.mat[1:2, 1] <- design.mat[2:3, 2] <- 1
design.mat

m1 <- glmFit(dt2, design.mat)
summary(m1)

# Contrasts: HG - LG
lrt12 <- glmLRT(fit, 
                contrast = c(-1, 1))
lrt12
topTags(lrt12, 
        n = 10)
sum2 <- decideTestsDGE(lrt12, 
                       adjust.method = "none",
                       p.value = 0.05)
summary(sum2)


de2tags12 <- rownames(dt2)[as.logical(sum2)]
plotSmear(lrt12, 
          de.tags=de2tags12)
abline(h = c(-1, 1), 
       col = "blue",
       lty = 2)

# Log2 differences----
dt1$`diff(wj2,wj1)` <- log2(dt2$counts[, 2]/dt2$counts[, 1])

# Log2 means----
dt1$`mean(wj2,wj1)` <- log2(dt2$counts[, 2]*dt2$counts[, 1])/2

# Lot log2 differences vs. means----
plot(dt1$`diff(wj2,wj1)` ~ dt1$`mean(wj2,wj1)`)
plot(abs(dt1$`diff(wj2,wj1)`) ~ dt1$`mean(wj2,wj1)`)