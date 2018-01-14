# |----------------------------------------------------------------------------------|
# | Project: Study of Diabetes in MES13 cells                                        |
# | Script: Analysis of microarray data from the Golub's experiment                  |
# | Scientist: Chao Wang                                                             |
# | Author: Davit Sargsyan                                                           |
# | Created: 12/02/2017                                                              |
# |----------------------------------------------------------------------------------|
# NOTE: remove TNF samples
# Header----
# Save consol output to a log file
# sink(file = "tmp/log_WT_KDKD_analysis_v5.txt")

# source("https://bioconductor.org/biocLite.R")
# biocLite("multtest")

require(data.table)
require(ggplot2)
require(multtest)

# Part I: Data----
help(golub)
data(golub)
summary(golub)

tmp <- lapply(data.table(golub),
              function(a) {
                a <- a[a != min(a)]
                density(a)
              })

out <- list()
for (i in 1:length(tmp)) {
  out[[i]] <- data.table(smpl = i,
                         x = tmp[[i]]$x,
                         y = tmp[[i]]$y)
}
dt.dest <- rbindlist(out)
dt.dest

clr <- data.table(smpl = 1:length(golub.cl),
                  grp = golub.cl)

dt.dest <- merge(clr, dt.dest, by = "smpl")
dt.dest$grp <- factor(dt.dest$grp,
                      levels = c(0, 1),
                      labels = c("ALL", "AML"))
dt.dest

ggplot(dt.dest,
       aes(x = x,
           y = y,
           colour = grp)) +
  geom_line()

# Remove zeros----
dt2 <- lapply(data.table(golub),
              function(a) {
                a[a == min(a)] <- NA
                return(a)
              }) 
dt2 <- do.call("cbind", dt2)
dt2

plot(dt2[, 1] ~ dt2[, 38])
# First 27 columns are ALL, last 11 are AML
golub.cl

cor(dt2,use = "complete.obs")

lapply(tmp, plot, add = TRUE)
plot(tmp[[3]])
lapply()
hist(tmp, 100)
tmp <- tmp[tmp != min(tmp)]
hist(tmp, 100, freq = FALSE)
plot(density(tmp))

# Difference between 2 samples: ALL vs. AML
dt1 <- data.table(dt1, 
                  golub[, c(1, 38)])
dt1
hist(golub[, 1], 100)
plot(dt1$V1 ~ dt1$V2)

# Remove artificial data (most likely created for log of zeros by adding 1)
dt1 <- subset(dt1,
              V1 != min(V1) &
                V2 != min(V2))
plot(dt1)

# Standard deviation by gene
dt1 <- data.table(std.all = apply(golub[, golub.cl == 0],
                                  MARGIN = 1,
                                  FUN = sd),
                  std.aml = apply(golub[, golub.cl == 1],
                                  MARGIN = 1,
                                  FUN = sd),
                  mu.all = apply(golub[, golub.cl == 0],
                                 MARGIN = 1,
                                 FUN = mean),
                  mu.aml = apply(golub[, golub.cl == 1],
                                 MARGIN = 1,
                                 FUN = mean))
plot(dt1$std.all ~ dt1$mu.all)
plot(dt1$std.aml ~ dt1$mu.aml)

n.all <- sum(golub.cl == 0)
n.aml <- sum(golub.cl == 1)
dt1$std.pooled <- sqrt(((n.all - 1)*dt1$std.all^2 + (n.aml - 1)*dt1$std.aml^2)/(n.all + n.aml - 2))

# Part II: Compare----
# Means and differences
dt1$diff <- dt1$V1 - dt1$V2
dt1$avg.diff <- dt1$mu.all - dt1$mu.aml

dt1$mu <- (dt1$V1 + dt1$V2)/2
dt1$avg.mu <- (dt1$mu.all + dt1$mu.aml)/2
plot(dt1$mu ~ dt1$avg.mu)

dt1

# Plot differences vs. means vs. SD----
plot(dt1$avg.diff ~ dt1$avg.mu,
     xlab = "Means",
     ylab = "Differences",
     main = "ALL vs. AML")

plot(dt1$std.pooled ~ dt1$avg.mu,
     xlab = "Means",
     ylab = "Pooled SD",
     main = "ALL vs. AML")

plot(dt1$std.pooled ~ dt1$avg.diff,
     xlab = "Differences",
     ylab = "SD",
     main = "ALL vs. AML")

# Method 2: spline over absolute values of differences
dt1$abs.diff <- abs(dt1$avg.diff)

m2 <- smooth.spline(x = dt1$avg.mu,
                    y = dt1$abs.diff)
plot(m2)
m2

plot(dt1$abs.diff ~ dt1$avg.mu,
     main = "SD Estimation with Smoothing Spline",
     xlab = "Means",
     ylab = "Absolue Differences")
lines(m2, 
      col = "red",
      lw = 3)

res2 <- data.table(avg.mu = m2$x,
                   std.spline = m2$y)
res2

dt1 <- merge(dt1,
             res2,
             by = "avg.mu")
dt1

# How well the 2 methods agree?
plot(dt1$std.pooled ~ dt1$std.spline,
     xlim = c(0, 1.5),
     ylim = c(0, 1.5),
     main = "Comparison of SD From the Two Methods",
     xlab = "Pooled SD",
     ylab = "2-Sample Estimated SD")
abline(0, 1)

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