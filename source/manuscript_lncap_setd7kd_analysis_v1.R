# |----------------------------------------------------------------------------------|
# | Project: Study of Diabetes in MES13 cells                                        |
# | Script: Analysis of microarray data from a LNCaP WT/KD experiment                |
# | Scientist: Chao Wang                                                             |
# | Author: Davit Sargsyan                                                           |
# | Created: 12/02/2017                                                              |
# |----------------------------------------------------------------------------------|
# NOTE: remove TNF samples
# Header----
# Save consol output to a log file
# sink(file = "tmp/log_WT_KDKD_analysis_v5.txt")

require(data.table)
require(ggplot2)

# Part I: Data----
# Preprocessed data was created by 'WT_KDKD_data_preprocessing_v1.R'
load("data/lncap_setd7_normilized_annotated.RData")
dt1

# Make all gene names unique----
setkey(dt1, SYMBOL)
dt1[, N := 1:.N,
    by = SYMBOL]
dt1$SYMBOL[dt1$N > 1] <- paste(dt1$SYMBOL[dt1$N > 1],
                               dt1$N[dt1$N > 1],
                               sep = "_")
dt1[, N := NULL]
dt1

# Remove TNF samples
dt1 <- droplevels(subset(dt1,
                         select = -c(6, 9)))
dt1

# Rename columns----
colnames(dt1)[5:8] <- c("WT_PEITC",
                        "WT",
                        "KD_PEITC",
                        "KD")

# Part II: long data----
dtl <- melt.data.table(dt1,
                       id.vars = "SYMBOL",
                       measure.vars = 5:8,
                       variable.name = "Group",
                       value.name = "Expression")
dtl$SYMBOL <- factor(dtl$SYMBOL)
dtl$Group <- factor(dtl$Group,
                    levels = c("WT",
                               "KD",
                               "WT_PEITC",
                               "KD_PEITC"))
dtl
summary(dtl)

# # Part III: Compare----
# Hitmap of differences----
# Differences----
dt1$`WT PEITC vs. WT` <- dt1$WT_PEITC - dt1$WT
dt1$`KD PEITC vs. KD` <- dt1$KD_PEITC - dt1$KD
dt1$`WT vs. KD` <- dt1$WT - dt1$KD

# Sums----
dt1$`Mean(WT PEITC, WT)` <- (dt1$WT_PEITC + dt1$WT)/2
dt1$`Mean(KD PEITC, KD)` <- (dt1$KD_PEITC + dt1$KD)/2
dt1$`Mean(WT, KD)` <- (dt1$WT + dt1$KD)/2

dt1

# # p-Values for single-replica samples----
# # Plot differences vs. means in WT PEITC vs. WT----
# plot(dt1$`WT PEITC vs. WT` ~ dt1$`Mean(WT PEITC, WT)`,
#      xlab = "Means",
#      ylab = "Differences",
#      main = "WT PEITC vs. WT")
# # Method 1: Regularize by the above within margin of epsilon
# epln <- 0.1
# dt1$sd <- NA
# 
# for (i in 1:nrow(dt1)) {
#   tmp <- subset(dt1,
#                 (`Mean(WT PEITC, WT)` <= dt1$`Mean(KD PEITC, KD)`[i] + epln) &
#                   (`Mean(WT PEITC, WT)` >= dt1$`Mean(KD PEITC, KD)`[i] - epln))
#   dt1$sd[i] <- sd(tmp$`WT PEITC vs. WT`)
# }
# 
# hist(dt1$`WT PEITC vs. WT`, 100)
# 
# m1 <- loess(dt1$sd ~ dt1$`Mean(WT PEITC, WT)`)
# dt1$`Fitted SD` <- predict(m1,
#                            newdata = data.frame(`Mean(WT PEITC, WT)` = dt1$`Mean(WT PEITC, WT)`))
# dt1
# dt1 <- dt1[order(dt1$`Mean(WT PEITC, WT)`), ]
# 
# plot(dt1$sd ~ dt1$`Mean(KD PEITC, KD)`,
#      xlab = "Means",
#      ylab = "SD",
#      main = "WT PEITC vs. WT")
# lines(dt1$`Fitted SD` ~ dt1$`Mean(WT PEITC, WT)`,
#       col = "red",
#       lw = 2)
# 
# # Assuming t follows normal distribution (it does not!)
# # dt1$t <- dt1$`WT PEITC vs. WT`/dt1$sd
# # qqnorm(dt1$t)
# # abline(0, 1)

# Method 2: spline over absolute values of differences
dt1$abs.wtpeitc.wt <- abs(dt1$`WT PEITC vs. WT`)


m2 <- smooth.spline(x = dt1$`Mean(WT PEITC, WT)`,
                    y = dt1$abs.wtpeitc.wt)
plot(m2)
m2

plot(dt1$abs.wtpeitc.wt ~ dt1$`Mean(WT PEITC, WT)`,
     main = "SD Estimation with Smoothing Spline",
     xlab = "Means",
     ylab = "Absolue Differences")
lines(m2, 
      col = "red",
      lw = 3)

res2 <- data.table(`Mean(WT PEITC, WT)` = m2$x,
                   `Spline SD` = m2$y)
res2

dt1 <- merge(dt1,
             res2,
             by = "Mean(WT PEITC, WT)")

# # How well the 2 methods agree?
# plot(dt1$`Fitted SD` ~ dt1$`Spline SD`,
#      xlim = c(0, 0.3),
#      ylim = c(0, 0.3),
#      main = "Comparison of SD From the Two Methods",
#      xlab = "Smoothing Spline of Absolute Differences",
#      ylab = "LOESS of Moving Window")
# abline(0, 1)

# # Use Method 2 to calculate stats using only the controls
# set.seed(2017)
# tmp <- data.table(s1 = sample(dt1$KD,
#                               nrow(dt1),
#                               replace = TRUE),
#                   s2 = sample(dt1$KD,
#                               nrow(dt1),
#                               replace = TRUE))
# plot(tmp)
# 
# tmp$diff <- tmp$s1 - tmp$s2
# tmp$mu <- (tmp$s1 + tmp$s2)/2
# plot(tmp$diff ~ tmp$mu,
#      main = "Random Samples with Replacement of WT",
#      xlab = "Means",
#      ylab = "Differences")
# 
# # WRONG! I cannot just sample at random from the vector 
# # as each gene has expression range!

dt1$t <- dt1$`WT PEITC vs. WT`/dt1$`Spline SD`
qqnorm(dt1$t)
qqline(dt1$t)

cutoffs <- quantile(dt1$t,
         probs = c(0.025, 
                   0.975))
hist(dt1$t, 100)
abline(v = cutoffs)

# CONTINUE HERE!!!



# Test statistic
dt1$t <- dt1$`WT PEITC vs. WT`/dt1$`Fitted SD`
qqnorm(dt1$t)
abline(0, 1)

dt1$p <- 2*pnorm(-abs(dt1$t))
hist(dt1$p)

tmp <- subset(dt1,
              p <= 0.05,
              select = c(1:6, 9, 12, 17:20))
tmp <- tmp[order(tmp$SYMBOL), ]
tmp
write.csv(tmp, file = "tmp/pvals.csv")

# sink()