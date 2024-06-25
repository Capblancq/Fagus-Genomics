#################################################
#         Thibaut Capblancq Feb - 2023          #
#                                               #
# Sequencing coverage and polymorphic sites     #
# F. orientalis genomes                         #
#                                               #
#################################################

# Set working directory
setwd("~/Documents/Project-Fagus/Analyses/Study-Genetic-Diversity/Sequencing/")

# Required libraries


###############################################################
#############        MAPPING STATS       #####################

# Mean coverage per individuals after removing the PCR duplicates
cov_ind <- read.table("rmdup_mean_coverage.txt")
mean(cov_ind[,1])
min(cov_ind[,1])
max(cov_ind[,1])

# Mapping quality scores
mapping <- data.frame(read.table("rmdup.names.txt"), read.table("rmdup_reads_mapping_Qscores.txt"))
colnames(mapping) <- c("Sample", "Hits", "Total_mapped","Total_mapq10","Total_mapq20","Total_mapq30")
nodup <- 
mean(cov_ind[,1])
min(cov_ind[,1])
max(cov_ind[,1])
mean(mapping$Total_mapq20/mapping$Total_mapped)

# Number of reads
nbreads <- read.table("NBreads_1P.txt")
mean(nbreads[,1])
min(nbreads[,1])
max(nbreads[,1])

# PCR duplicates rates
dup <- data.frame(read.table("names.txt"), read.table("res.aln.reads.out"))
colnames(dup) <- c("Sample","Total", "Read1", "Read2","ProperlyPaired","ItselfMate","Singletons", "MateDiffChr", "MateDiffChr5")
dup <- dup[!duplicated(dup$Sample),]
dup <- dup[-which(dup$Sample==732),]
nodup <- data.frame(read.table("rmdup.names.txt"), read.table("rmdup.res.aln.reads.out"))
colnames(nodup) <- c("Sample","Total", "Read1", "Read2","ProperlyPaired","ItselfMate","Singletons", "MateDiffChr", "MateDiffChr5")
pcrdup <- merge(dup, nodup, by = c("Sample"))
pcrdup$ratio <- pcrdup$Total.y/pcrdup$Total.x
mean(1-pcrdup$ratio)


###############################################################

###############################################################
#############        MAPPING RESULTS      #####################

# Coverage tables from samtools bedcov
files <- list.files(".", pattern = "^coverage_")
cov <- lapply(files, read.table)
cov <- lapply(cov, function(x) {colnames(x) <- c("Chr", "Start", "Stop", "Reads", "Sites"); return(x)})
names(cov) <- unlist(lapply(strsplit(basename(files), split = "_|\\."), function(x) x[2]))

# Find mean coverage at covered sites
cov_all <- data.frame(Chr = cov[[1]][,1], midPos = cov[[1]][,2]+50000, Sites = cov[[1]][,5], Coverage = apply(do.call(cbind, lapply(cov, function(x) x$Reads/x$Sites)), 1, mean))

# Filter out windows with no site covered
cov_all$Coverage[cov_all$Sites==0] <- NA
cov_all$Sites[cov_all$Sites==0] <- NA

# Filtered out windows with coverage above the max used == 50
cov_all$Sites[cov_all$Coverage>50] <- NA
cov_all$Coverage[cov_all$Coverage>50] <- NA

# Nb of SNPs
NbSNPs <- read.table("NB_snps_windows.txt", header = T)
colnames(NbSNPs) <- c("Chr", "midPos", "NbSNPs")
NbSNPs$NbSNPs[NbSNPs$NbSNPs==0] <- NA

# Chromosome lengths
chrom_length <- read.table("Bagha_contig_lengths.txt")
colnames(chrom_length) <- c("Chrom", "Length")

# Coverage 
p.cov <- ggplot() +
  geom_rect(data = chrom_length, aes(xmin = 0, xmax = Length, ymin = (13-as.integer(as.factor(Chrom))-0.25), ymax = (13-as.integer(as.factor(Chrom))+0.25)), fill = "grey90") +
  geom_tile(data = cov_all, aes(x = midPos, y = (13-as.integer(as.factor(Chr))), fill = Coverage), width = 100000, height = 0.5) +
  scale_fill_viridis_b(option = "inferno", direction = -1, na.value = 'grey90', limits=c(0,50)) +
  xlab("Position along the chromosome (Mbp)") +
  ylab("") +
  labs(fill="Mean coverage") +
  theme_classic() +
  scale_y_continuous(breaks=seq(1,12,1), limits = c(0,12+1), labels = paste0("Chr",rev(seq(1,12,1)))) +
  scale_x_continuous(breaks = seq(0,80000000,10000000),
                     labels = seq(0,80,10)) +
  geom_segment(aes(x=0,xend=80000000,y=-Inf,yend=-Inf)) +
  theme(axis.line=element_blank(), axis.ticks.y = element_blank(), panel.grid = element_blank(), panel.background = element_blank(), legend.position = c(0.8, 0.25), legend.background = element_rect(fill = "transparent"))

# Covered sites
p.sites <- ggplot() +
  geom_rect(data = chrom_length, aes(xmin = 0, xmax = Length, ymin = (13-as.integer(as.factor(Chrom))-0.25), ymax = (13-as.integer(as.factor(Chrom))+0.25)), fill = "grey90") +
  geom_tile(data = cov_all, aes(x = midPos, y = (13-as.integer(as.factor(Chr))), fill = Sites/100000), width = 100000, height = 0.5) +
  scale_fill_viridis_b(option = "viridis", na.value = 'grey90', direction = -1, limits=c(0,1)) +
  xlab("Position along the chromosome (Mbp)") +
  ylab("") +
  labs(fill="Covered sites (%)") +
  theme_classic() +
  scale_y_continuous(breaks=seq(1,12,1), limits = c(0,12+1), labels = paste0("Chr",rev(seq(1,12,1)))) +
  scale_x_continuous(breaks = seq(0,80000000,10000000),
                     labels = seq(0,80,10)) +
  geom_segment(aes(x=0,xend=80000000,y=-Inf,yend=-Inf)) +
  theme(axis.line=element_blank(), axis.ticks.y = element_blank(), panel.grid = element_blank(), panel.background = element_blank(), legend.position = c(0.8, 0.25), legend.background = element_rect(fill = "transparent"))

# Number of SNPs
p.SNPs <- ggplot() +
  geom_rect(data = chrom_length, aes(xmin = 0, xmax = Length, ymin = (13-as.integer(as.factor(Chrom))-0.25), ymax = (13-as.integer(as.factor(Chrom))+0.25)), fill = "grey90") +
  geom_tile(data = NbSNPs, aes(x = midPos, y = (13-as.integer(as.factor(Chr))), fill = NbSNPs), width = 100000, height = 0.5) +
  scale_fill_viridis_b(option = "magma", na.value = 'grey90', direction = -1, limits=c(0,5000), values = c(0,0.1,0.2,0.5,1)) +
  xlab("Position along the chromosome (Mbp)") +
  ylab("") +
  labs(fill="Number of SNPs") +
  theme_classic() +
  scale_y_continuous(breaks=seq(1,12,1), limits = c(0,12+1), labels = paste0("Chr",rev(seq(1,12,1)))) +
  scale_x_continuous(breaks = seq(0,80000000,10000000),
                     labels = seq(0,80,10)) +
  geom_segment(aes(x=0,xend=80000000,y=-Inf,yend=-Inf)) +
  theme(axis.line=element_blank(), axis.ticks.y = element_blank(), panel.grid = element_blank(), panel.background = element_blank(), legend.position = c(0.8, 0.25), legend.background = element_rect(fill = "transparent"))

pdf("Sequencing_along_scaffolds.pdf", width = 14, height = 6)
ggarrange(p.cov, p.sites, p.SNPs, nrow = 1)
dev.off()

png("Sequencing_along_scaffolds.png", width = 14, height = 6, units = "in", res = 300)
ggarrange(p.cov, p.sites, p.SNPs, nrow = 1)
dev.off()

###############################################################
