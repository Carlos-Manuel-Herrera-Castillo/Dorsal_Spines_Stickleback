#!/usr/bin/env Rscript

# Step 08c: Analyze first and third dorsal spine length in F2 hybrids.
# Original section: U* Analysis of first and third dorsal spine length
#
# Run with:
#   Rscript scripts/08_phenotypes/spine_length.R

rm(list = ls())

library(ggplot2)

setwd("C:/Users/bernerd/switchdrive/Institution/carlos/research/DSpheno_DB")

iter <- 10000

bootstrap_lm_effect <- function(data, response, size_var, iter = 10000) {
  formula <- as.formula(paste(response, "~", size_var, "+ phenotype"))
  red <- length(which(data$phenotype == "2"))
  compl <- length(which(data$phenotype == "3"))

  distr <- numeric(iter)
  for (i in seq_len(iter)) {
    boot_data <- data[
      c(
        sample(which(data$phenotype == "2"), size = red, replace = TRUE),
        sample(which(data$phenotype == "3"), size = compl, replace = TRUE)
      ),
    ]
    boot_model <- lm(formula, data = boot_data)
    distr[i] <- summary(boot_model)$coef[3, 1]
  }

  distr
}

standardize <- function(x) {
  (x - mean(x)) / sd(x)
}

# Publication version: size correction with PC1 from cube-root mass and standard length.
d <- read.table("d_dsLength.txt", h = TRUE)
d$phenotype <- as.character(d$phenotype)

size_pca <- prcomp(cbind(d$weight_g^(1 / 3), d$SL_mm), scale = TRUE)
pve <- size_pca$sdev^2 / sum(size_pca$sdev^2)
print(pve)

d$pc1 <- size_pca$x[, 1] * -1

# DS1 analysis.
d1 <- na.omit(d[, c("phenotype", "pc1", "DS1_mm")])
plot(
  d1[d1$phenotype == "2", "pc1"],
  d1[d1$phenotype == "2", "DS1_mm"],
  xlim = c(min(na.omit(d$pc1)), max(na.omit(d$pc1))),
  ylim = c(min(na.omit(d$DS1_mm)), max(na.omit(d$DS1_mm))),
  pch = 19,
  col = "blue"
)
points(d1[d1$phenotype == "3", "pc1"], d1[d1$phenotype == "3", "DS1_mm"], pch = 19, col = "red")

d1$DS1_mm <- standardize(d1$DS1_mm)

model_ds1_interaction <- lm(DS1_mm ~ pc1 * phenotype, data = d1)
print(summary(model_ds1_interaction))

model_ds1 <- lm(DS1_mm ~ pc1 + phenotype, data = d1)
print(summary(model_ds1))

ds1Distr <- bootstrap_lm_effect(d1, "DS1_mm", "pc1", iter)
q1 <- quantile(ds1Distr, probs = c(0.025, 0.975))
m1 <- median(ds1Distr)

# DS3 analysis.
d3 <- na.omit(d[, c("phenotype", "pc1", "DS3_mm")])
plot(
  d3[d3$phenotype == "2", "pc1"],
  d3[d3$phenotype == "2", "DS3_mm"],
  xlim = c(min(na.omit(d$pc1)), max(na.omit(d$pc1))),
  ylim = c(min(na.omit(d$DS3_mm)), max(na.omit(d$DS3_mm))),
  pch = 19,
  col = "blue"
)
points(d3[d3$phenotype == "3", "pc1"], d3[d3$phenotype == "3", "DS3_mm"], pch = 19, col = "red")

d3$DS3_mm <- standardize(d3$DS3_mm)

model_ds3_interaction <- lm(DS3_mm ~ pc1 * phenotype, data = d3)
print(summary(model_ds3_interaction))

model_ds3 <- lm(DS3_mm ~ pc1 + phenotype, data = d3)
print(summary(model_ds3))

ds3Distr <- bootstrap_lm_effect(d3, "DS3_mm", "pc1", iter)
q3 <- quantile(ds3Distr, probs = c(0.025, 0.975))
m3 <- median(ds3Distr)

plot_data <- data.frame(
  name = c(rep("ds1L", iter), rep("ds3L", iter)),
  value = c(ds1Distr, ds3Distr)
)

p <- ggplot(plot_data, aes(x = name, y = value, fill = name)) +
  geom_violin() +
  geom_segment(aes(x = 1, y = q1[1], xend = 1, yend = q1[2])) +
  geom_segment(aes(x = 2, y = q3[1], xend = 2, yend = q3[2])) +
  geom_point(aes(x = 1, y = m1), color = "black") +
  geom_point(aes(x = 2, y = m3), color = "black")

plot(p)
