# filter raw DArT data and generate PCA

# load required libraries
library(dartR)
library(adegenet)
library(melfuR)
library(vegan)
library(viridis)

# set working directory and location of data files
setwd("/path/to/your/working/directory")
DATA <- "/path/to/data"

# load data
load(paste0(DATA,"/KH_GP_raw.Rdata"))  
load(paste0(DATA,"/meta.Rdata"))


# filter by read depth
gl.report.rdepth(GP_dart)
rd <- gl.filter.rdepth(GP_dart, lower = 8, upper = 100, verbose = 3)

# filter loci with >20% missing data
gl.report.callrate(rd, method = "loc")
call80.loc <- gl.filter.callrate(rd, method = "loc", threshold = 0.8, mono.rm = TRUE, recursive = TRUE, recalc = TRUE)

# check missing data per individual
dartR.base::gl.report.callrate(call80.loc, method = "ind", ind.to.list = 50)

# filter for 99% repeatability and remove loci with NA in all individuals/populations
rep100 <- gl.filter.reproducibility(call80.loc, threshold = 0.99)
rep100 <- gl.filter.allna(rep100, by.pop = TRUE, verbose = 3)

# filter secondaries
gl.report.secondaries(rep100)
sec <- gl.filter.secondaries(rep100, method = "best")

# filter for Minor Allele Frequency (MAF) = 0.01
MAF01 <- gl.filter.maf(sec, threshold = 0.01)


# 3. PCA Analysis and Visualization
# impute missing data
geno <- gl2gi(MAF01)
imputed.geno <- impute_genotypes(geno, K = 5)

# generate matrix of allele counts
alleles <- imputed.geno@tab
snps <- alleles[, seq(1, ncol(alleles), 2)]
colnames(snps) <- locNames(imputed.geno)

# set plotting variables
cols <- viridis(n = 5)
meta$HighHe_color <- cols[as.factor(meta$basin_highHe)]
legend_colors_highHe <- unique(meta[c("basin_highHe", "HighHe_color")])

# run PCA
pc <- rda(snps)

# generate labels for axes
x.lab <- paste0("PC1 (", paste(round((pc$CA$eig[1] / pc$tot.chi * 100), 2)), "%)")
y.lab <- paste0("PC2 (", paste(round((pc$CA$eig[2] / pc$tot.chi * 100), 2)), "%)")

# generate plot
{pdf(file = "HighHe_pca.pdf", height = 6, width = 6)
plot(pc, choices = c(1, 2), type = "n", xlab = x.lab, ylab = y.lab, cex.lab = 1)
with(meta, points(pc, display = "sites", col = meta$HighHe_color, pch = 1, cex = 1.1))
legend("topleft", legend = legend_colors_highHe$basin_highHe, 
       col = legend_colors_highHe$HighHe_color, pch = 19, pt.cex = 1, cex = 1, 
       xpd = 1, box.lty = 0, bg = "transparent")
dev.off()}

