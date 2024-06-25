#################################################
#        Thibaut Capblancq April - 2022         #
#                                               #
# Patterns of genetic diversity along           #
# Coenonympha genomes                           #
#                                               #
#################################################

# Set working directory
setwd("~/Documents/Project-Fagus/Analyses/Study-Genetic-Diversity/Site-Genetic-Diversity/")

# Required libraries
library(reshape2)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(raster)
library(ggpubr)

###############################################################
################# MEAN LOCALITY ESTIMATES #####################

# Metadata
meta <- read.table("../../F_orientalis_metadata.txt", header = T)
meta_loc <- aggregate(meta[,6:8], by = list(meta$Population), FUN = mean)
colnames(meta_loc)[1] <- "Population"

# Thetas tables from ANGSD
files <- list.files("./", pattern = ".pestPG$")
thetas <- lapply(files, function(x) read.table(x))
names <- unlist(lapply(strsplit(files, split = ".thetas"), function(x) x[1]))
names(thetas) <- names

# Averaging pi for each locality
pi <- unlist(lapply(thetas, function(x) mean(x[,5]/x[,14], na.rm = T)))
tajD <- unlist(lapply(thetas, function(x) mean(x[,9], na.rm = T)))

# Stats
min(pi)
mean(pi)
max(pi)
min(tajD)
mean(tajD)
max(tajD)

###############################################################


##############################
####     GENETIC LOAD     ####

# Files with freq major allele and annotation
list <- list.files(path = "./", pattern = "_FreqAnnotations.txt$", recursive=TRUE)
genes <- lapply(list, read.table)
names(genes) <- unlist(strsplit(list, split = "_FreqAnnotations.txt"))
genes <- lapply(genes, function(x) {colnames(x) <- c("Chrom_Pos", "Ploidy", "NbHaplotypes", "FreqRef", "Genotypes", "Annotation"); return(x)})

# Remove site with Na values
genes <- lapply(genes, function(x) {x <- x[!is.na(x$Genotypes),]; return(x)})
genes <- lapply(genes, function(x) {x <- x[!is.na(x$Annotation),]; return(x)})
genes <- lapply(genes, function(x) {x <- x[!is.na(x$FreqRef),]; return(x)})
genes <- lapply(genes, function(x) {x <- x[-which(x$NbHaplotypes<5),]; return(x)})

# Estimate Freq for ancestral homozygotes, heterozygotes and derived homozygotes
genes <- lapply(genes, function(x) {
  
  x$FreqAA <- as.integer(unlist(lapply(strsplit(x$Genotypes, split = "/"), function(x) x[1])))/(x$NbHaplotypes/2)
  x$FreqAD <- as.integer(unlist(lapply(strsplit(x$Genotypes, split = "/"), function(x) x[2])))/(x$NbHaplotypes/2)
  x$FreqDD <- as.integer(unlist(lapply(strsplit(x$Genotypes, split = "/"), function(x) x[3])))/(x$NbHaplotypes/2)
  x$FreqDerived <- x$FreqAD/2 + x$FreqDD
  
  return(x)})

# Remove sites that are monomorphic for the reference allele in the locality
genes <- lapply(genes, function(x) {x <- x[!(x$FreqDerived==0),]; return(x)})

# Number of synonymlous and non synonymous sites 
nbsyn <- lapply(genes, function(x) sum(x$Annotation=="synonymous_variant"))
nbnonsyn <- lapply(genes, function(x) sum(x$Annotation!="synonymous_variant"))

# Mean frequency of minor allele for synonymous and non-synonymous variants
meanfS <- lapply(genes, function(x) mean(x[which(x$Annotation=="synonymous_variant"), "FreqDerived"]))
meanfN <- lapply(genes, function(x) mean(x[which(x$Annotation!="synonymous_variant"), "FreqDerived"]))

# Genetic Load ratio
GenLoad <- NULL
GenLoad <- unlist(lapply(1:length(genes), function(x) { GenLoad[[names(genes)[x]]] <- (nbnonsyn[[x]]*meanfN[[x]])/(nbsyn[[x]]*meanfS[[x]]); return(GenLoad) }))

# Stats
min(GenLoad)
mean(GenLoad)
max(GenLoad)

# Realized load = nb heterozygote nonsynonymous genotypes / total nonsynonymous genotypes
GenLoad_fixed <- unlist(lapply(genes, function(x) mean(x$FreqDD)/(mean(x$FreqDD)+mean(x$FreqAD))))
GenLoad_het <- unlist(lapply(genes, function(x) mean(x$FreqAD)/(mean(x$FreqDD)+mean(x$FreqAD))))


##############################


##############################
####        FIGURE        ####

# Global table
TAB <- data.frame(meta_loc, Pi = pi, TajD = tajD, Load = GenLoad, Masked = GenLoad_het, Realized = GenLoad_fixed)

# Coorelation
summary(lm(TAB$Load ~ TAB$Longitude))
summary(lm(TAB$Pi ~ TAB$Longitude))
summary(lm(TAB$Load ~ TAB$Elevation))
summary(lm(TAB$Pi ~ TAB$Elevation))

# Mapping features
world <- map_data('world')
range <- crop(shapefile("../../Fagus_orientalis_range/Fagus_orientalis_EUFORGEN.shp"), extent(39, 52, 36, 46))
range <- fortify(range)

# Map
p1 <- ggplot(world, aes(long, lat)) + 
  geom_map(map=world, aes(map_id=region), fill=gray(.95), color="black") +
  geom_polygon(data = range, aes(long, lat, group = group), colour = NA, fill = "black", linewidth = 0.3, alpha = 0.2) + 
  coord_quickmap(xlim = c(41, 50), ylim = c(38.5, 43.5)) +
  geom_point(data=TAB, aes(x= Longitude, y= Latitude, size = Load, fill = Pi), shape = 21, alpha = 1) + 
  scale_size_continuous(range = c(1, 10)) +
  scale_fill_distiller(palette = "Spectral") +
  xlab("Longitude") + ylab("Latitude") +
  theme_bw(base_size = 12, base_family = "Times") +
  theme(panel.grid = element_blank(), legend.background = element_blank(), legend.box.background = element_blank(), panel.background = element_blank(), plot.background = element_blank(), strip.text = element_text(size=12))

# Map masked load for suppmat
pdf("Realized_Load.pdf", width = 8, height = 4.6)
ggplot(world, aes(long, lat)) + 
  geom_map(map=world, aes(map_id=region), fill=gray(.95), color="black") +
  geom_polygon(data = range, aes(long, lat, group = group), colour = NA, fill = "black", linewidth = 0.3, alpha = 0.2) + 
  coord_quickmap(xlim = c(41, 50), ylim = c(38.5, 43.5)) +
  geom_point(data=TAB, aes(x= Longitude, y= Latitude, size = GenLoad, fill = Realized), shape = 21, alpha = 1) + 
  scale_size_continuous(range = c(1, 10)) +
  scale_fill_distiller(palette = "Spectral") +
  xlab("Longitude") + ylab("Latitude") +
  theme_bw(base_size = 12, base_family = "Times") +
  theme(panel.grid = element_blank(), legend.background = element_blank(), legend.box.background = element_blank(), panel.background = element_blank(), plot.background = element_blank(), strip.text = element_text(size=12))
dev.off()

# Pi ~ elevation
TAB2 <- melt(TAB[,1:7], id.vars = c("Population","Latitude","Longitude","Elevation", "TajD"))
colnames(TAB2) <- c("Population","Latitude","Longitude","Elevation", "TajD", "Parameter", "Value_param")
TAB2 <- melt(TAB2, id.vars = c("Population","Latitude","TajD", "Parameter", "Value_param"))
colnames(TAB2) <- c("Population","Latitude","TajD", "Parameter", "Value_param", "Geography", "Value_geo")

p2 <- ggplot(data = TAB2, aes(y = Value_param, x = Value_geo)) + #, fill = TajD, size = TajD
  #geom_point(shape=19, size = 2) + 
  #geom_smooth(method="lm", colour = "grey20") + 
  geom_point(shape = 21, alpha = 1, size = 2.5, fill = "grey50") +
  #scale_fill_viridis_b() +
  scale_size_continuous(range = c(1, 7)) +
  facet_grid(Parameter~Geography, scales = "free", switch = "both") +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 12, base_family = "Times") +
  theme(strip.background = element_blank(),
        strip.placement = "outside",
        panel.grid = element_blank(), 
        legend.background = element_blank(), 
        legend.box.background = element_blank(), 
        panel.background = element_blank(), 
        plot.background = element_blank(), 
        strip.text = element_text(size=12))

# Figure
pdf("Diversity_Forientalis.pdf", width = 12, height = 4.6)
ggarrange(p1, p2, widths = c(1.5, 1), common.legend = F, legend = "right", labels = c("(a)", "(b)"), align = "h")
dev.off()


##############################
