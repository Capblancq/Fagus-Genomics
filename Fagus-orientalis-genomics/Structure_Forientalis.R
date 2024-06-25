#################################################
#     Thibaut Capblancq November - 2021         #
#                                               #
# PCA and NGSAdmix results for Fagus orientalis #
#                                               #
#################################################

# Required libraries
library(ade4)
library(ggplot2)
library(ggpubr)
library(wesanderson)
library(reshape2)
library(grid)
library(gridExtra)
library(plyr)
library(raster)
library(sp)
library(rgeos)
library(maps)
library(maptools)
library(mapplots)
library(scatterpie)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(ggspatial)
library(viridis)

# Set working directory
setwd("~/Documents/Project-Fagus/Analyses/Study-Genetic-Diversity/Genetic-Structure/")


#####################################################################
####################           PCA            #######################

## PCA from covariance matrix
TAB <- read.table("./Fagus_orientalis_poly.cov")
pca <- eigen(TAB)

## Info about individuals
names <- read.table("./FO_bam.list")
names <- unlist(strsplit(basename(as.character(names[,1])), split = ".final.bam"))
meta <- read.table("../../F_orientalis_metadata.txt", header=T)
meta <- meta[match(names, as.character(meta$Samples)),]

## Explained variance
var <- pca$values/sum(pca$values)
barplot(pca$values[1:20])
barplot(var[1:20])

## Table
TAB <- data.frame(PC1 = pca$vectors[,1], PC2 = pca$vectors[,2], Population = meta$Population, Range = meta$Area, Elevation = meta$Elevation)

## Removing an outlier individual AZ_10 862
TAB <- TAB[-which(meta$Samples == 862),]
TAB$Range <- factor(TAB$Range, levels = c("Lesser_Caucasus", "Likhi_Range", "Greater_Caucasus_W", "Greater_Caucasus_E", "Karabakh", "Hyrcanian_Forests"))
TAB <- TAB[order(TAB$Range),]

## Plot
p1 <- ggscatter(TAB, x = "PC1", y = "PC2",
                color = "Range",
                ellipse = FALSE, mean.point = FALSE,
                star.plot = TRUE, size = 2.5, star.plot.lwd = .2) +  
  #scale_color_manual(values = c("Greater_Caucasus" = "#BB3754FF", "Likhi_Range" = "#F98C0AFF", "Lesser_Caucasus" = "#FFF59D", "Hyrcanian_Forests" = "#56106EFF", "Karabakh" = "grey20")) +
  scale_color_manual(values = c("Greater_Caucasus_W" = "#66a182", "Greater_Caucasus_E" = "#edae49", "Lesser_Caucasus" = "#00798c", "Hyrcanian_Forests" = "#d1495b", "Karabakh" = "#2e4057", "Likhi_Range" = "black")) +
  geom_hline(yintercept=0, linetype="dotted", color = gray(.35), size=0.2) +
  geom_vline(xintercept=0, linetype="dotted", color = gray(.35), size=0.2) +
  labs(x = "PC1 (5.9%)", y = "PC2 (1%)") +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.background = element_blank(), legend.background = element_blank(), panel.grid = element_blank(), plot.background = element_blank(), legend.text=element_text(size=rel(.8)), strip.text = element_text(size=11), legend.position = c(0.82, 0.2))


#####################################################################


#####################################################################
##################           NGSadmix            ####################

## Import of the NGSadmix results
list <- list.files(pattern = ".qopt$")
proba<-lapply(list, function(x) read.table(x,h=F))

## Info about individuals
names <- read.table("./FO_bam.list")
names <- unlist(strsplit(basename(as.character(names[,1])), split = ".final.bam"))
meta <- read.table("../../F_orientalis_metadata.txt", header=T)
meta <- meta[match(names, as.character(meta$Samples)),]

## Removing an outlier individual AZ_10 862
proba<-lapply(proba, function(x) x[-which(meta$Samples == 862),])
meta <- meta[-which(meta$Samples == 862),]

## Sort by longitude and Area
long_area <- unlist(lapply(1:nrow(meta), function(x) min(meta$Longitude[meta$Area==meta$Area[x]])))
long_pop <- unlist(lapply(1:nrow(meta), function(x) min(meta$Longitude[meta$Population==meta$Population[x]])))

probas <- lapply(proba, function(x) x[order(long_area, long_pop, decreasing = F),])
groups <- meta$Area[order(long_area, long_pop, decreasing = F)]
pop <- meta$Population[order(long_area, long_pop, decreasing = F)]

## Barplot K from 2 to 4
seps = (as.integer(as.factor(groups)) - c(1, as.integer(as.factor(groups[-length(groups)]))))
cols=c('Proba1' = "#00798c", 'Proba2' = "#d1495b", 'Proba3' = "#66a182", 'Proba4' = "#edae49", 'Proba5' = "#6d597a", 'Proba6' = "#001219")

colnames(probas[[1]]) <- c('Proba1','Proba2')
colnames(probas[[2]]) <- c('Proba1','Proba3','Proba2')
colnames(probas[[3]]) <- c('Proba4','Proba2','Proba1','Proba3')
colnames(probas[[4]]) <- c('Proba1','Proba4','Proba5','Proba3', 'Proba2')
colnames(probas[[5]]) <- c('Proba4','Proba3','Proba1','Proba2', 'Proba5', 'Proba6')

pdf("./NGSadmix_K2-6.pdf", width = 10, height = 6)
opar <- par(no.readonly = T, xpd=FALSE)
layout(mat = matrix(c(rep(c(1:(5+1)),4)),byrow = F, ncol = 4), heights = c(rep(12,5), 6))
par(mar=c(0,3,1,4))
for (i in 1:5) {
  barplot(t(probas[[i]]), las=2, names = rep("", ncol(t(probas[[i]]))), axes = F, border = NA, space=0, col=cols[match(colnames(probas[[i]]),names(cols))])
  abline(v=(which(seps!=0)-1), lty=1, lwd=1)
  mtext(text = paste0('K = ', as.character(i+1)), side = 4, las=1)
  axis(side=2,pos=-2, las=2, at = pretty(c(0,1), n=5), labels=c('0%', rep("",4), '100%'))
}
coord <- 1:nrow(proba[[1]])

coord_site <- aggregate(coord, by = list(pop), FUN = mean)
axis(side=1, pos = -0.05, at = coord_site$x, labels=rep("", length(coord_site$Group.1)), cex.axis=0.5, padj = 0)
text(x = coord_site$x, 0.15, labels=as.character(coord_site$Group.1), cex=0.7, xpd = TRUE, srt = 90)

coord_region <- aggregate(coord, by = list(groups), FUN = mean)
axis(side=1, pos = -0.2, at = coord_region$x, labels=c("Greater Caucasus East", "Greater Caucasus West", "Hyrcanian Forests", "Karabakh", "Lesser Caucasus", "Likhi Range"), cex.axis=1, padj = 0)
dev.off()


## Likelihoods
lkh <- read.table("likelihoods.txt", header = T)
mean <- as.vector(aggregate(lkh$Likelihood, by = list(lkh$K), mean)[,2])
sd <- as.vector(aggregate(lkh$Likelihood, by = list(lkh$K), sd)[,2])
delta <- c()
for (i in 1:length(mean)){
  delta[i] <- abs(mean[i])/sd[i]
}

pdf("Likelihoods.pdf", height = 3.5, width = 3.5)
plot(2:6, delta[-1], xlab = 'K', ylab = 'Delta K', pch = 16, bty = "n", col = "gray25", axes = FALSE, cex = 1.5)
axis(1, col="gray25", col.ticks="gray25", col.axis="gray25")
axis(2, col="gray25", col.ticks="gray25", col.axis="gray25")
mtext("K", side=1, line=3, col="gray25")
mtext("Delta K", side=2, line=3, col="gray25")
lines(2:6, delta[-1], col = "gray25", lwd = 2)
dev.off()

## Maps and piecharts

# Mapping features
world <- map_data('world')
range <- crop(shapefile("../../Fagus_orientalis_range/Fagus_orientalis_EUFORGEN.shp"), extent(39, 52, 36, 46))
range <- fortify(range)

# Probabilities for K=2 and K=3
proba_2 <- proba[[1]]
rownames(proba_2) <- meta$ID
proba_3 <- proba[[2]]
rownames(proba_3) <- meta$ID
proba_4 <- proba[[3]]
rownames(proba_4) <- meta$ID

# Mean proba in each pop for K = 2 and K = 3
coord_pop <- aggregate(meta[,c("Longitude", "Latitude")], by = list(meta$Population), FUN = mean)
colnames(coord_pop) <- c("Population", "Longitude", "Latitude")
nb_ind <- aggregate(meta$Population, by = list(meta$Population), FUN = length)

prob_k2 <- aggregate(proba_2, by = list(meta$Population), FUN = mean)
colnames(prob_k2) <- c("Population", "Proba_1", "Proba_2")

prob_k3 <- aggregate(proba_3, by = list(meta$Population), FUN = mean)
colnames(prob_k3) <- c("Population", "Proba_2", "Proba_3", "Proba_1")

prob_k4 <- aggregate(proba_4, by = list(meta$Population), FUN = mean)
colnames(prob_k4) <- c("Population", "Proba_3", "Proba_1", "Proba_2", "Proba_4")

# Global table
TAB <- data.frame(Longitude = rep(coord_pop$Longitude, 3), Latitude = rep(coord_pop$Latitude, 3), Pop = rep(coord_pop$Population, 3), Proba_1 = c(prob_k2$Proba_1, prob_k3$Proba_1, prob_k4$Proba_1), Proba_2 = c(prob_k2$Proba_2, prob_k3$Proba_2, prob_k4$Proba_2), Proba_3 = c(rep(0, nrow(prob_k2)), prob_k3$Proba_3, prob_k4$Proba_3), Proba_4 = c(rep(0, 2*nrow(prob_k2)), prob_k4$Proba_4), Nb_ind = rep(nb_ind$x, 3), variable = rep(c("K = 2", "K = 3", "K = 4"), each = nrow(prob_k2)))

# Plotting
p2 <- ggplot(world, aes(long, lat)) + 
  geom_map(map=world, aes(map_id=region), fill=gray(.95), color="black") +
  geom_polygon(data = range, aes(long, lat, group = group), colour = NA, fill = "black", size = 0.3, alpha = 0.2) + 
  coord_quickmap(xlim = c(41, 50), ylim = c(38.5, 43.5)) +
  geom_scatterpie(data = TAB[which(TAB$variable=="K = 4"),], aes(x = Longitude, y = Latitude, group = Pop, r = log(Nb_ind)/8), cols = c("Proba_1", "Proba_2", "Proba_3", "Proba_4")) +
  scale_fill_manual(values = c(Proba_2 = "#00798c", Proba_1 = "#d1495b", Proba_3 = "#66a182", Proba_4 = "#edae49")) +
  #facet_wrap(~variable, ncol = 2) +
  xlab("Longitude") +
  ylab("Latitude") +
  guides(fill="none") +
  theme_bw(base_size = 11, base_family = "Times") +
  theme(panel.grid = element_blank(), strip.text = element_text(size=12))

# Figure
pdf("Structure_Forientalis.pdf", width = 11, height = 4.8)
ggarrange(p1, p2, widths = c(1, 1.4), common.legend = F, labels = c("(a)", "(b)"))
dev.off()


###################################################################################
