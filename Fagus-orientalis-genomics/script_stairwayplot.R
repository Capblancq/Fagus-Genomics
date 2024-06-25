#################################################
#         Thibaut Capblancq January - 2023      #
#                                               #
#   STAIRWAY PLOT analysis visualization        #
#                                               #
#################################################

# Set working directory
setwd("~/Documents/Project-Fagus/Analyses/Study-Genetic-Diversity/Stairwayplot/")

# Required libraries
library(ggplot2)

#### For each species ####

# Stairway plot results
SPres <- lapply(list.files(recursive = T, path="./", pattern=".final.summary$", full.names = TRUE), function(x) read.table(x, header=TRUE))
names(SPres) <- unlist(lapply(strsplit(basename(list.files(recursive = T, path="./", pattern=".final.summary$", full.names = TRUE)), split = "_|\\."), function(x) x[2]))

# Confidence interval
ci.ne <- lapply(SPres, function(x) c(x$Ne_2.5., rev(x$Ne_97.5.)))
ci.ya <- lapply(SPres, function(x) c(x$year, rev(x$year)))

# Plot
pdf("./Stairwayplot_last4Myrs.pdf", width = 7, height = 4.5)
ggplot() +
  geom_path(data = SPres$hyrcanian, aes(x = year, y = Ne_median, colour = "Hyrcanian Forest"), lwd = 0.8) +
  geom_path(data = SPres$hyrcanian, aes(x = year, y = Ne_2.5., colour = "Hyrcanian Forest"), lty = "dotted", lwd = 0.4) +
  geom_path(data = SPres$hyrcanian, aes(x = year, y = Ne_97.5., colour = "Hyrcanian Forest"), lty = "dotted", lwd = 0.4) +
  geom_path(data = SPres$GCeast, aes(x = year, y = Ne_median, colour = "Greater Caucasus E"), lwd = 0.8) +
  geom_path(data = SPres$GCeast, aes(x = year, y = Ne_2.5., colour = "Greater Caucasus E"), lty = "dotted", lwd = 0.4) +
  geom_path(data = SPres$GCeast, aes(x = year, y = Ne_97.5., colour = "Greater Caucasus E"), lty = "dotted", lwd = 0.4) +
  geom_path(data = SPres$GCwest, aes(x = year, y = Ne_median, colour = "Greater Caucasus W"), lwd = 0.8) +
  geom_path(data = SPres$GCwest, aes(x = year, y = Ne_2.5., colour = "Greater Caucasus W"), lty = "dotted", lwd = 0.4) +
  geom_path(data = SPres$GCwest, aes(x = year, y = Ne_97.5., colour = "Greater Caucasus W"), lty = "dotted", lwd = 0.4) +
  geom_path(data = SPres$LC, aes(x = year, y = Ne_median, colour = "Lesser Caucasus"), lwd = 0.8) +
  geom_path(data = SPres$LC, aes(x = year, y = Ne_2.5., colour = "Lesser Caucasus"), lty = "dotted", lwd = 0.4) +
  geom_path(data = SPres$LC, aes(x = year, y = Ne_97.5., colour = "Lesser Caucasus"), lty = "dotted", lwd = 0.4) +
  scale_x_continuous(limits = c(800, 4000000), expand=c(0,0), trans = "log10", breaks = c(0, 1000, 20000, 200000, 2000000, 10000000, 40000000), labels = scales::comma) +
  #scale_y_continuous(expand=c(0.01,0), trans = "log10", breaks = c(0, 1000, 10000, 100000, 1000000, 10000000), labels = scales::comma) +
  scale_y_continuous(limits = c(0, 5000000), expand=c(0.01,0), breaks = c(0, 500000, 2000000, 4000000), labels = scales::comma) +
  labs(x = "Time (years ago)",
       y = "Effective population size (Ne)",
       colour = "") + 
  scale_colour_manual(values = c("Lesser Caucasus" = "#00798c", "Greater Caucasus W" = "#66a182", "Greater Caucasus E" = "#edae49", "Hyrcanian Forest" = "#d1495b")) +
  theme_bw(base_size = 12, base_family = "Times") +
  theme(legend.position = c(0.85,0.2),
  #theme(legend.position = c(0.15,0.88),
        legend.background = element_blank(),
    panel.grid = element_blank(),
    panel.background = element_blank(),
    axis.text.y = element_text(angle = 90, hjust = 0.5)
  )
dev.off()