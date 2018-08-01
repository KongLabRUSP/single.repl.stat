# |----------------------------------------------------------------------------------|
# | Project: Skin UVB SKH1 mouse model treated with UA/SFN                           |
# | Script: RNA-seq data preprocessing for N=1 manuscript                            |
# | Coordinator: Ran Yin, Renyi Wu                                                   |
# | Author: Davit Sargsyan                                                           |
# | Created: 07/23/2018                                                              |
# | Modified:                                                                        |
# |----------------------------------------------------------------------------------|
# sink(file = "tmp/log_data_preprocessing_uvb_rna_v1.txt")
date()

# Workflow: https://www.bioconductor.org/help/workflows/rnaseqGene/
# source("https://bioconductor.org/biocLite.R")
# biocLite("DESeq2")

require(data.table)
require(DESeq2)
require(readxl)
require(BiocParallel)
require(ggplot2)
require(knitr)

# NOTE on DESeq2 Output: 'baseMean' is the average of the normalized count values, 
# divided by the size factors, taken over all samples in the DESeqDataSet

# Load data----
dt1 <- fread("data/skin_uvb/featurescounts_uvb-skin_dedup_renyi_2-9-2018.csv",
             skip = 1)
dt1

# Remove unused columns----
dt1 <- dt1[, c(1, 7:ncol(dt1)), with = FALSE]
dt1

cnames <- colnames(dt1)[-1]
cnames <- gsub(x = cnames,
               pattern = ".dedup.bam",
               replacement = "")
cnames
colnames(dt1)[-1] <- cnames
dt1

# Select controls at weeks 2 and 25 only (max aging effect)----
tmp <- dt1[, c("Geneid",
               "02w_CON_0",
               "02w_CON_1",
               "25w_CON_0",
               "25w_CON_1")]

# Exclude genes with all zeros----
gene.keep <- tmp$Geneid[rowSums(tmp[, -1]) > 0]
tmp <- tmp[Geneid %in% gene.keep, ]
dt1 <- dt1[Geneid %in% gene.keep, ]

dt2 <- melt.data.table(data = tmp,
                       id.vars = 1,
                       measure.vars = 2:5,
                       variable.name = "Sample",
                       value.name = "GeneExp")
dt2$Grp <- factor(as.numeric(substr(dt2$Sample, 1, 2)))
dt2$Repl <- factor(as.numeric(substr(dt2$Sample, 9, 9)))
dt2$log2GeneExp <- log2(dt2$GeneExp + 1)
dt2

# Compute meas and sd for each gene----
dt3 <- aggregate(dt2$log2GeneExp,
                 by = list(dt2$Geneid),
                 FUN = "mean")
dt3$std <- aggregate(dt2$log2GeneExp,
                    by = list(dt2$Geneid),
                    FUN = "sd")$x
dt3 <- data.table(dt3)
colnames(dt3)[1:2] <- c("Geneid",
                        "mu")
dt3
plot(dt3$std ~ dt3$mu,
     pch = ".")

# Differences between the group means----
dt4 <- aggregate(dt2$log2GeneExp,
                 by = list(dt2$Geneid,
                           dt2$Grp),
                 FUN = "mean")
dt4 <- data.table(dt4)
colnames(dt4) <- c("Geneid",
                   "Grp",
                   "mu")
dt4 <- dcast.data.table(data = dt4,
                        Geneid ~ Grp,
                        value.var = "mu")
dt4$diff <- dt4$`25` - dt4$`2`
dt4

dt3 <- merge(dt3, 
             dt4,
             by = "Geneid")
dt3$absDiff <- abs(dt3$diff)
dt3

# Pairwise differences----
tmp <- dcast.data.table(data = dt2,
                        Geneid ~ Sample,
                        value.var = "log2GeneExp")
tmp$w25.0_w2.0 <- abs(tmp$`25w_CON_0` - tmp$`02w_CON_0`)
tmp$w25.0_w2.1 <- abs(tmp$`25w_CON_0` - tmp$`02w_CON_1`)
tmp$w25.1_w2.0 <- abs(tmp$`25w_CON_1` - tmp$`02w_CON_0`)
tmp$w25.1_w2.1 <- abs(tmp$`25w_CON_1` - tmp$`02w_CON_1`)
tmp

dt3 <- merge(dt3, 
             tmp[, c("Geneid",
                     "w25.0_w2.0", 
                     "w25.0_w2.1",
                     "w25.1_w2.0",
                     "w25.1_w2.1")],
             by = "Geneid")
dt3

# Plot----
plot(dt3$diff ~ dt3$mu,
     pch = ".")
abline(h = 0,
       col = "red",
       lty = 2)

# Smoothing with LOESS----
span <- 0.5
degree <- 1

# Weights: proportional to number of nearest neighbors
# i.e. the more genes have similar means, the lower the weight
dt3$wgt <- NA
for (i in 1:nrow(dt3)) {
  dt3$wgt[i] <- sum((dt3$mu - abs(dt3$mu[i])) < diff(range(dt3$mu))/100)
}
dt3$wgt <- dt3$wgt/sum(dt3$wgt)
sum(dt3$wgt)
hist(dt3$wgt)

# a. Predicted standard deviation, standard error and test stats for N=2----
dt3$std.loess <- loess(std ~ mu,
                       data = dt3,
                       span = span,
                       degree = degree,
                       family = "symmetric",
                       weights = wgt)$fitted
dt3$se.loess <- loess(std/sqrt(2) ~ mu,
                       data = dt3,
                       span = span,
                       degree = degree)$fitted
plot(dt3$se.loess ~ dt3$mu)
dt3$statN2 <- dt3$absDiff/dt3$std.loess
summary(dt3$statN2)
hist(dt3$statN2)
hist(log(dt3$statN2 + 1))

# b. Predicted standard deviation and test stats for N=1,w25.0 vs. w2.0 ----
dt3$std.w25.0_w2.0 <- loess(w25.0_w2.0 ~ mu,
                       data = dt3,
                       span = span,
                       degree = degree)$fitted
plot(dt3$std.w25.0_w2.0 ~ dt3$mu)
dt3$stat.w25.0_w2.0 <- dt3$w25.0_w2.0/dt3$std.w25.0_w2.0
summary(dt3$stat.w25.0_w2.0)
hist(dt3$stat.w25.0_w2.0)
hist(log(dt3$stat.w25.0_w2.0 + 1))

# c. Predicted standard deviation and test stats for N=1,w25.1 vs. w2.0 ----
dt3$std.w25.1_w2.0 <- loess(w25.1_w2.0 ~ mu,
                            data = dt3,
                            span = span,
                            degree = degree)$fitted
plot(dt3$std.w25.1_w2.0 ~ dt3$mu)
dt3$stat.w25.1_w2.0 <- dt3$w25.1_w2.0/dt3$std.w25.1_w2.0
summary(dt3$stat.w25.1_w2.0)
hist(dt3$stat.w25.1_w2.0)
hist(log(dt3$stat.w25.1_w2.0 + 1))

# d. Predicted standard deviation and test stats for N=1,w25.0 vs. w2.1 ----
dt3$std.w25.0_w2.1 <- loess(w25.0_w2.1 ~ mu,
                            data = dt3,
                            span = span,
                            degree = degree)$fitted
plot(dt3$std.w25.0_w2.1 ~ dt3$mu)
dt3$stat.w25.0_w2.1 <- dt3$w25.0_w2.1/dt3$std.w25.0_w2.1
summary(dt3$stat.w25.0_w2.1)
hist(dt3$stat.w25.0_w2.1)
hist(log(dt3$stat.w25.0_w2.1 + 1))

# e. Predicted standard deviation and test stats for N=1,w25.1 vs. w2.1 ----
dt3$std.w25.1_w2.1 <- loess(w25.1_w2.1 ~ mu,
                            data = dt3,
                            span = span,
                            degree = degree)$fitted
plot(dt3$std.w25.1_w2.1 ~ dt3$mu)
dt3$stat.w25.1_w2.1 <- dt3$w25.1_w2.1/dt3$std.w25.1_w2.1
summary(dt3$stat.w25.1_w2.1)
hist(dt3$stat.w25.1_w2.1)
hist(log(dt3$stat.w25.1_w2.1 + 1))

# f. Based on DESeq2----
mat <- data.frame(sample = c("02w_CON_0",
                             "02w_CON_1",
                             "25w_CON_0",
                             "25w_CON_1"),
                  time = factor(rep(c(2, 25), each = 2)),
                  repl = factor(rep(1:2, 2)))
mat

dtm <- as.matrix(dt1[, c("02w_CON_0",
                         "02w_CON_1",
                         "25w_CON_0",
                         "25w_CON_1")] + 1)
rownames(dtm) <- dt1$Geneid
dtm

dds <- DESeqDataSetFromMatrix(countData = dtm, 
                              colData = mat,
                              ~ time)
dds

# If all samples contain zeros, geometric means cannot be
# estimated. Change default 'type = "ratio"' to 'type = "poscounts"'.
# Type '?DESeq2::estimateSizeFactors' for more details.
dds <- estimateSizeFactors(dds,
                           type = "poscounts")
dds

# Run DESeq----
dds <- DESeq(dds,
             fitType = "local",
             parallel = TRUE)
res <- results(dds)
res
res <- data.table(Geneid = rownames(res),
                  std.deseq2 = sqrt(2)*res$lfcSE,
                  stat.deseq2 = res$stat,
                  deseq2.pval = res$pvalue)
res
resultsNames(dds)
colData(dds)

# Save DESeq2 p-values----
dt3 <- merge(dt3,
             res,
             by = "Geneid")
dt3

# Compare statistics----
plot(dt3$deseq2.pval ~ log(dt3$statN2 + 1), 
     pch = ".")
# NOTE: tails of the statistic's distributioon correspond to low p-values from DESeq2

plot(dt3$deseq2.pval ~ log(dt3$stat.w25.0_w2.0 + 1), 
     pch = ".")

plot(dt3$deseq2.pval ~ log(dt3$w25.0_w2.1 + 1), 
     pch = ".")

plot(dt3$deseq2.pval ~ log(dt3$w25.1_w2.0 + 1), 
     pch = ".")

plot(dt3$deseq2.pval ~ log(dt3$w25.1_w2.1 + 1), 
     pch = ".")

plot(log(dt3$statN2 + 1) ~ log(dt3$w25.0_w2.0 + 1), 
     pch = ".")
plot(log(dt3$statN2 + 1) ~ log(dt3$w25.0_w2.1 + 1), 
     pch = ".")
plot(log(dt3$statN2 + 1) ~ log(dt3$w25.1_w2.0 + 1), 
     pch = ".")
plot(log(dt3$statN2 + 1) ~ log(dt3$w25.1_w2.1 + 1), 
     pch = ".")

plot(dt3$std.deseq2 ~ dt3$mu, 
     pch = ".")
plot(dt3$std ~ dt3$mu, 
     pch = ".")
abline(a = 0, b = 2)

# Plot all together----
ggplot(dt3,
       aes(x = mu,
           y = absDiff)) +
  geom_point(aes(x = mu,
                 y = std),
             pch = ".",
             col = "black") +
  geom_line(aes(x = mu,
                y = std.loess,
                colour = "Gene Weighted SD"),
            # linetype = "dashed",
            size = 1.2) +
  geom_line(aes(x = mu,
                y = se.loess,
                colour = "Gene SE"),
            # linetype = "dashed",
            size = 1.2) +
  # geom_line(aes(x = mu,
  #               y = std.w25.0_w2.0),
  #           col = "blue",
  #           linetype = "dashed",
  #           size = 1.2) +
  # geom_line(aes(x = mu,
  #               y = std.w25.0_w2.1),
  #           col = "black",
  #           linetype = "dashed",
  #           size = 1.2) +
  # geom_line(aes(x = mu,
  #               y = std.w25.1_w2.0),
  #           col = "brown",
  #           linetype = "dashed",
  #           size = 1.2) +
  # geom_line(aes(x = mu,
  #               y = std.w25.1_w2.1),
  #           col = "orange",
  #           linetype = "dashed",
  #           size = 1.2) +
  geom_smooth(aes(x = mu,
                  y = std.deseq2,
                  colour = "DESeq2 SE")) +
  geom_abline(slope = c(0, 2),
              intercept = c(0, 0),
              linetype = "dashed") +
  scale_x_continuous("Mean") +
  scale_y_continuous("Absolute Log2 Differences") +
  scale_colour_manual("Estimates",
                      values = c(`Gene Weighted SD` = "green",
                                 `Gene SE` = "red",
                                 `DESeq2 SE` = "blue")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))
