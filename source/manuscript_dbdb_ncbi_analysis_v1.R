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
sampleNames(l1) <- c("DBDB 1",
                     "DBDB 2",
                     "DBDB 3",
                     "DBDB 4",
                     "DB+ 1",
                     "DB+ 2",
                     "DB+ 3",
                     "DB+ 4")
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


#dt1$diff.1.1 <- dt1$`DBDB 1` - dt1$`DB+ 1`
