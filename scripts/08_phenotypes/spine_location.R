#!/usr/bin/env Rscript

# Step 08a: Analyze dorsal spine location and vertebral number in F2 hybrids.
# Original section: S* Analysis of dorsal spine location in the F2 cross
#
# Run with:
#   Rscript scripts/08_phenotypes/spine_location.R

rm(list = ls())

library(ggplot2)

setwd("C:/Users/bernerd/switchdrive/Institution/carlos/research/DSpheno_DB")

data <- read.table("d_dsPos.txt", h = TRUE)
iter <- 10000

bootstrap_mean_difference <- function(one, two, iter = 10000) {
  diffs <- numeric(iter)
  for (i in seq_len(iter)) {
    diffs[i] <- mean(sample(one, replace = TRUE)) - mean(sample(two, replace = TRUE))
  }
  diffs
}

# Histograms for 1-spined and 2-spined individuals.
one_spined <- data[which(data$Nspines == 1), ]
pdf("1DS.pdf")
hist(
  one_spined$ds1 + 0.001,
  xlim = c(3, 10),
  ylim = c(0, 30),
  breaks = seq(3, 10, by = 0.5),
  freq = TRUE,
  main = NA,
  ann = FALSE,
  las = 1
)
dev.off()

two_spined <- data[which(data$Nspine == 2), ]
pdf("2DS.pdf")
hist(
  c(two_spined$ds1 + 0.001, two_spined$ds2 + 0.001),
  xlim = c(3, 10),
  ylim = c(0, 30),
  breaks = seq(3, 10, by = 0.5),
  freq = TRUE,
  main = NA,
  ann = FALSE,
  las = 1
)
dev.off()

# Bootstrap CI for DS1 position difference.
one <- data[which(data$Nspine == 1), "ds1"]
two <- data[which(data$Nspine == 2), "ds1"]
difsDS1 <- bootstrap_mean_difference(one, two, iter)
print(mean(one) - mean(two))
print(quantile(difsDS1, prob = c(0.025, 0.975)))

# Bootstrap CI for DS3 position difference.
one <- data[which(data$Nspine == 1), "ds3"]
two <- data[which(data$Nspine == 2), "ds3"]
difsDS3 <- bootstrap_mean_difference(one, two, iter)
print(mean(one) - mean(two))
print(quantile(difsDS3, prob = c(0.025, 0.975)))

# DS3 position among 1-spined individuals with anterior vs posterior DS1.
oneAnt <- data[which((data$Nspine == 1) & (data$ds1 <= 5)), "ds3"]
onePost <- data[which((data$Nspine == 1) & (data$ds1 >= 6)), "ds3"]
difsDS3_ant_post <- bootstrap_mean_difference(oneAnt, onePost, iter)
print(mean(oneAnt) - mean(onePost))
print(quantile(difsDS3_ant_post, prob = c(0.025, 0.975)))

# Bootstrap CI for total vertebral number difference.
one <- data[which(data$Nspine == 1), "Nvert"]
two <- data[which(data$Nspine == 2), "Nvert"]
difsNvert <- bootstrap_mean_difference(one, two, iter)
print(mean(one) - mean(two))
print(quantile(difsNvert, prob = c(0.025, 0.975)))

# Plot point estimates and 95% CIs.
plot(0, 0, xlim = c(-0.38, 1.135), ylim = c(0, 2), type = "n")
segments(0, 0, 0, 2, col = "gray")
points(-0.077, 0.5, pch = 16, cex = 1.2)
segments(-0.350, 0.5, 0.192, 0.5, lwd = 2)
points(-0.027, 1, pch = 16, cex = 1.2)
segments(-0.188, 1, 0.136, 1, lwd = 2)
points(0.874, 1.5, pch = 16, cex = 1.2)
segments(0.654, 1.5, 1.105, 1.5, lwd = 2)

plot_data <- data.frame(
  name = c(rep("DS1", iter), rep("DS3", iter), rep("Nvert", iter)),
  value = c(difsDS1, difsDS3, difsNvert)
)

p <- ggplot(plot_data, aes(x = name, y = value, fill = name)) +
  geom_violin() +
  geom_segment(aes(x = 0.5, y = 0.654, xend = 0.5, yend = 1.105)) +
  geom_segment(aes(x = 1, y = -0.188, xend = 1, yend = 0.136)) +
  geom_segment(aes(x = 1.5, y = -0.350, xend = 1.5, yend = 0.192)) +
  geom_point(aes(x = 0.5, y = 0.874), color = "black") +
  geom_point(aes(x = 1, y = -0.027), color = "black") +
  geom_point(aes(x = 1.5, y = -0.077), color = "black")

plot(p)
