#################################################
#         Thibaut Capblancq Jan - 2023          #
#                                               #
# Patterns of genetic diversity along           #
# F. orientalis genomes                         #
#                                               #
#################################################

# Set working directory
setwd("~/Documents/Project-Fagus/Analyses/Study-Genetic-Diversity/Genomic-Diversity/")

# Required libraries
library(reshape2)
library(ggplot2)
library(ggridges)
library(ggpubr)
library(purrr)

###############################################################
######################### DATA ################################

# Thetas tables from ANGSD
files <- list.files("./100kb_windows/", pattern = ".thetas.idx.pestPG$", full.names = T)
thetas <- lapply(files, read.table)
thetas <- lapply(thetas, function(x) {colnames(x) <- c("region", "Chr", "midPos", "tW", "tP", "tF", "tH", "tL", "Tajima", "fuf", "fud", "fayh", "zeng", "nSites"); return(x)})
names(thetas) <- unlist(lapply(strsplit(basename(files), split = "_|\\."), function(x) x[1]))

# Weight theats by the number of sites
thetas <- lapply(thetas, function(x) {x[,"tP"] <- x[,"tP"]/x[,"nSites"]; return(x)})
thetas <- lapply(thetas, function(x) {x[,"tW"] <- x[,"tW"]/x[,"nSites"]; return(x)})

# Remove small contigs
thetas <- lapply(thetas, function(x) x[grep("OU0", x[,"Chr"]),])

# Rename chromosomes
thetas <- lapply(thetas, function(x) {x[,"Chr"] <- substr(x[,"Chr"], 1, 8); return(x)})

# Add range information and filter out non used columns
thetas <- lapply(1:length(thetas), function(x) {thetas[[x]]$Range <- names(thetas)[x]; return(thetas[[x]][, c("Chr","midPos","tP","tW","Tajima","Range")])})
names(thetas) <- unlist(lapply(strsplit(basename(files), split = "_|\\."), function(x) x[1]))

# Fst tables
fst_GCW_GCE_H <- read.table("./100kb_windows/Fst_sliding_GCW_GCE_H", header = T)
colnames(fst_GCW_GCE_H) <- c("region","Chr","midPos","Nsites","Fst_GCW_GCE","Fst_GCW_H","Fst_GCE_H","PBS0","PBS1","PBS2")
fst_GCW_LC_H <- read.table("./100kb_windows/Fst_sliding_GCW_LC_H", header = T)
colnames(fst_GCW_LC_H) <- c("region","Chr","midPos","Nsites","Fst_GCW_LC","Fst_GCW_H","Fst_LC_H","PBS0","PBS1","PBS2")
fst <- merge(fst_GCW_GCE_H[,2:7], fst_GCW_LC_H[,c(2:5,7)], by = c("Chr","midPos"))

# Remove small contigs
fst <- fst[grep("OU0", fst[,"Chr"]),]

# Rename chromosomes
fst[,"Chr"] <- substr(fst[,"Chr"], 1, 8)

# Remove windows with 
fst <- fst[-which(fst[,"Nsites.x"]<2000 | fst[,"Nsites.y"]<2000),]

# LD per windows
files <- list.files("./", pattern = "_LDr2bp.txt$")
LD <- lapply(files, read.table)
LD <- lapply(LD, function(x) {colnames(x) <- c("Chr", "midPos", "MeanR2", "Nsites"); return(x)})
names(LD) <- unlist(lapply(strsplit(files, split = "\\.|_"), function(x) x[1]))
LD <- lapply(1:length(LD), function(x) {LD[[x]]$Range <- names(LD)[x]; return(LD[[x]][, c("Chr","midPos","MeanR2","Nsites","Range")])})
names(LD) <- unlist(lapply(strsplit(files, split = "\\.|_"), function(x) x[1]))
LD <- lapply(LD, function(x) {x[,"Chr"] <- substr(x[,"Chr"], 1, 8); return(x)})

###############################################################


###############################################################
##################      Statistics     ########################

# Pi
lapply(thetas, function(x) mean(x[,"tP"], na.rm = T))

# Waterson theta
lapply(thetas, function(x) mean(x[,"tW"], na.rm = T))

# Tajima's D
lapply(thetas, function(x) mean(x[,"Tajima"], na.rm = T))

# Fst
colMeans(fst[,-c(1:3,7)])

# LD
lapply(LD, function(x) mean(x[,"MeanR2"], na.rm = T))

###############################################################


###############################################################
########      Distribution across scaffolds       #############

# Create global table for diversity estimates
TAB.thetas <- do.call(rbind, thetas)
TAB.ld <- do.call(rbind, LD)
TAB.div <- merge(TAB.thetas, TAB.ld, by = c("Chr", "midPos", "Range"))
TAB.div <- melt(TAB.div, id.vars = c("Chr", "midPos", "Range", "Nsites"))

# Some formatting
TAB.div$Mean <- NA
TAB.div$Mean[TAB.div$variable=="tP"] <- mean(TAB.div$value[TAB.div$variable=="tP"], na.rm = T)
TAB.div$Mean[TAB.div$variable=="tW"] <- mean(TAB.div$value[TAB.div$variable=="tW"], na.rm = T)
TAB.div$Mean[TAB.div$variable=="Tajima"] <- mean(TAB.div$value[TAB.div$variable=="Tajima"], na.rm = T)
TAB.div$Mean[TAB.div$variable=="MeanR2"] <- mean(TAB.div$value[TAB.div$variable=="MeanR2"], na.rm = T)
TAB.div$variable <- as.character(TAB.div$variable)
TAB.div$species <- TAB.div$variable
TAB.div$variable[grep("tP", TAB.div$variable)] <- "π"
TAB.div$variable[grep("tW", TAB.div$variable)] <- "Theta W"
TAB.div$variable[grep("Tajima", TAB.div$variable)] <- "Tajima's D"
TAB.div$variable[grep("MeanR2", TAB.div$variable)] <- "LD mean R2"
TAB.div$variable <- factor(TAB.div$variable, levels = c("π", "Theta W", "Tajima's D", "LD mean R2"))
TAB.div$Range <- as.character(TAB.div$Range)
TAB.div$Range[TAB.div$Range=="LesserCaucasus"] <- "Lesser Caucasus"
TAB.div$Range[TAB.div$Range=="GreaterCaucasusWest"] <- "Greater Caucasus W"
TAB.div$Range[TAB.div$Range=="GreaterCaucasusEast"] <- "Greater Caucasus E"
TAB.div$Range <- factor(TAB.div$Range, levels = c("Lesser Caucasus", "Greater Caucasus W", "Greater Caucasus E", "Hyrcanian"))

# Global plot
pdf("Genomic_landscape.pdf", width = 12, height = 5)
ggplot(data = TAB.div, mapping = aes(x = midPos/1000000, y = value, color = Range)) +
  geom_hline(aes(yintercept = Mean), linewidth = 0.4, linetype="dotted") +
  geom_smooth(method = "loess", span = 0.1, linewidth = 0.7, se = FALSE) +
  #scale_color_manual(values = c(c("GreaterCaucasusWest" = "#BB3754FF", "GreaterCaucasusEast" = "#F98C0AFF", "LesserCaucasus" = "#FFF59D", "Hyrcanian" = "#56106EFF"))) +
  scale_color_manual(values = c(c("Greater Caucasus W" = "#66a182", "Greater Caucasus E" = "#edae49", "Lesser Caucasus" = "#00798c", "Hyrcanian" = "#d1495b"))) +
  facet_grid(variable ~ Chr, scales = "free", switch = "y", space = "free_x") +
  labs(x = "Nucleotide position (Mb)", 
       y = NULL) +
  theme(text = element_text(family = "Times"),
    strip.background = element_blank(),
    strip.placement = "outside",
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA),
    panel.background = element_blank(),
    legend.key = element_blank(),
    legend.title = element_blank(),
    legend.position = "bottom",
    axis.text.x = element_text()
  )
dev.off()

# Tajima's D only
TAB.div <- do.call(rbind, thetas)
TAB.div$Range <- factor(TAB.div$Range, levels = c("LesserCaucasus", "GreaterCaucasusWest", "GreaterCaucasusEast", "Hyrcanian"))
pdf("Genome_tajD.pdf", width = 20, height = 6)
ggplot(data = TAB.div, mapping = aes(x = midPos/1000000, y = Tajima, color = Range)) +
  geom_point(size = 0.3, alpha = 0.3) +
  scale_color_manual(values = c(c("GreaterCaucasusWest" = "#66a182", "GreaterCaucasusEast" = "#edae49", "LesserCaucasus" = "#00798c", "Hyrcanian" = "#d1495b"))) +
  facet_grid(Range ~ Chr, scales = "free_x", space = "free_x") +
  labs(x = "Nucleotide position (Mb)", 
       y = "Tajima's D") +
  theme( 
    text = element_text(family = "Times"),
    legend.position = "none",
    strip.background = element_blank(),
    strip.placement = "outside",
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA),
    panel.background = element_blank(),
    axis.text.x = element_text()
  )
dev.off()

# π only
TAB.div <- do.call(rbind, thetas)
TAB.div$Range <- factor(TAB.div$Range, levels = c("LesserCaucasus", "GreaterCaucasusWest", "GreaterCaucasusEast", "Hyrcanian"))
pdf("Genome_pi.pdf", width = 20, height = 6)
ggplot(data = TAB.div, mapping = aes(x = midPos/1000000, y = tP, color = Range)) +
  geom_point(size = 0.3, alpha = 0.3) +
  scale_color_manual(values = c(c("GreaterCaucasusWest" = "#66a182", "GreaterCaucasusEast" = "#edae49", "LesserCaucasus" = "#00798c", "Hyrcanian" = "#d1495b"))) +
  facet_grid(Range ~ Chr, scales = "free", space = "free_x") +
  labs(x = "Nucleotide position (Mb)", 
       y = "Nucleotide diversity π") +
  theme( 
    text = element_text(family = "Times"),
    legend.position = "none",
    strip.background = element_blank(),
    strip.placement = "outside",
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA),
    panel.background = element_blank(),
    axis.text.x = element_text()
  )
dev.off()

# LD
TAB.ld <- do.call(rbind, LD)
TAB.ld$Range <- factor(TAB.ld$Range, levels = c("LesserCaucasus", "GreaterCaucasusWest", "GreaterCaucasusEast", "Hyrcanian"))
TAB.ld <- TAB.ld[TAB.ld$Nsites>100,]

pdf("Genome_ld.pdf", width = 20, height = 6)
ggplot(data = TAB.ld, mapping = aes(x = midPos/1000000, y = HalfDecay, color = Range)) +
  geom_point(size = 0.5, alpha = 0.5) +
  scale_color_manual(values = c(c("GreaterCaucasusWest" = "#66a182", "GreaterCaucasusEast" = "#edae49", "LesserCaucasus" = "#00798c", "Hyrcanian" = "#d1495b"))) +
  facet_grid(Range ~ Chr, space = "free_x") +
  labs(x = "Nucleotide position (Mb)", 
       y = "LD half decay (bp)") +
  theme( 
    text = element_text(family = "Times"),
    legend.position = "none",
    strip.background = element_blank(),
    strip.placement = "outside",
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA),
    panel.background = element_blank(),
    axis.text.x = element_text()
  )
dev.off()

# Fst
TAB.fst <- melt(fst[,-c(3,7)], id.vars = c("Chr", "midPos"))
TAB.fst$variable <- as.character(TAB.fst$variable)
TAB.fst$variable[which(TAB.fst$variable == "Fst_GCW_GCE")] <- "GC West / GC East"
TAB.fst$variable[which(TAB.fst$variable == "Fst_GCW_H")] <- "GC West / Hyrcanian"
TAB.fst$variable[which(TAB.fst$variable == "Fst_GCW_LC")] <- "GC West / Lesser Caucasus"
TAB.fst$variable[which(TAB.fst$variable == "Fst_GCE_H")] <- "GC East / Hyrcanian"
TAB.fst$variable[which(TAB.fst$variable == "Fst_LC_H")] <- "Lesser Caucasus / Hyrcanian"
TAB.fst$variable <- factor(TAB.fst$variable, levels = c("GC West / Lesser Caucasus", "GC West / GC East", "GC East / Hyrcanian", "GC West / Hyrcanian", "Lesser Caucasus / Hyrcanian"))

# Plot
pdf("Genome_Fst.pdf", width = 12, height = 5)
ggplot(data = TAB.fst, mapping = aes(x = midPos/1000000, y = value, color = value)) +
  geom_hline(yintercept = 0.5, linewidth = 0.4, linetype="dotted") +
  geom_point(size = 0.8) +
  scale_color_viridis_c(option = "A", direction = -1) + 
  facet_grid(variable ~ Chr, scales = "free_x", space = "free_x") +
  labs(x = "Nucleotide position (Mb)", 
       y = "Fst",
       color = "Fst") +
  theme(
    text = element_text(family = "Times"),
    strip.background = element_blank(),
    strip.placement = "outside",
    panel.grid = element_blank(),
    panel.border = element_rect(colour = "black", fill=NA),
    panel.background = element_blank(),
    axis.text.x = element_text()
  )
dev.off()

###############################################################


#########################################
#########         PCA       #############

# Global table
tab.div <- Reduce(function(...) merge(..., by = c("Chr", "midPos")), thetas)
colnames(tab.div) <- c("Chr", "midPos", "tP_GCE", "tW_GCE", "TajD_GCE", "range_GCE", "tP_GCW", "tW_GCW", "TajD_GCW", "range_GCW", "tP_H", "tW_H", "TajD_H", "range_H", "tP_LC", "tW_LC", "TajD_LC", "range_LC")
tab.div <- tab.div[,-grep("range", colnames(tab.div))]
TAB.all <- merge(tab.div, fst, by = c("Chr", "midPos"))

# PCA
pca <- prcomp(TAB.all[,grep("TajD|Fst", colnames(TAB.all))], center = T, scale. = T)
#pca <- prcomp(TAB.all[,-which(colnames(TAB.all) %in% c("Chr","midPos","Nsites.x","Nsites.y"))], center = T, scale. = T)

# Formatting a table for ggplot
scores <- data.frame(TAB.all[,c("Chr", "midPos")], pca[["x"]][,c("PC1", "PC2")])
melt(scores, id.vars = c("Chr", "midPos", "PC1", "PC2"))

# Variable weights
variables <- data.frame(pca[["rotation"]][,c("PC1", "PC2")])
variables$param <- row.names(variables)

# Plot
pdf("PCA_windows.pdf", width = 12, height = 6)
ggplot() +
  geom_hline(yintercept=0, linetype="dotted", color = gray(.35), linewidth=0.5) +
  geom_vline(xintercept=0, linetype="dotted", color = gray(.35), linewidth=0.5) +
  geom_point(data = scores, aes(x = PC1, y = PC2), size = 1, alpha = 0.5) +
  geom_segment(data = variables, aes(x = 0, y = 0, xend = (PC1*50), yend = (PC2*50)), arrow = arrow(length = unit(1/2, "picas")), color = "black", size = 0.2) +
  geom_text(data = variables, aes(x = PC1*50, y = PC2*50, label = param)) +
  #facet_wrap(~variable) +
  #xlim(-5,8) +
  #ylim(-7,5) +
  theme_bw() +
  theme(panel.grid = element_blank(), panel.background = element_blank())
dev.off()
