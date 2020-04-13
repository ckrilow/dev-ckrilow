#!/usr/bin/env Rscript

library(ggplot2)

dat <- read.csv("../data/multiplet_rate_10x_3pt1.tsv", sep = "\t")
#dat$multiplet_pct <- dat$multiplet_pct / 100

print(summary(lm(multiplet_pct ~ no_of_cells_recovered, data = dat)))


pdf(file = paste0("multiplet_rate_10x_3pt1.pdf"), height = 5, width = 5)

plt <- ggplot2::ggplot(dat, ggplot2::aes(
    x = no_of_cells_recovered,
    y = multiplet_pct
))
plt <- plt + ggplot2::theme_bw(base_size = 12)
plt <- plt + ggplot2::geom_point()
plt <- plt + geom_smooth(method = "lm", se = T)
print(plt)

plt <- ggplot2::ggplot(dat, ggplot2::aes(
    x = no_of_cells_loaded,
    y = multiplet_pct
))
plt <- plt + ggplot2::theme_bw(base_size = 12)
plt <- plt + ggplot2::geom_point()
plt <- plt + geom_smooth(method = "lm", se = T)
print(plt)

dev.off()