# hybrid index and triangle plot

# load required libraries
library(pophelper)
library(gghybrid)
library(tidyverse)
library(dartR.base)
library(introgress)
library(viridis)

# set working directory and location of data files
setwd("/path/to/your/working/directory")
DATA <- "/path/to/data"

# define custom version of introgress::triangle.plot modified for additional graphics control
triangle.plot.cb <- function (hi.index = NULL, int.het = NULL, pdf = TRUE, 
                              out.file = "tri_plot.pdf", cols = NULL, pch = 19, 
                              legend = NULL, leg.cols = NULL, leg.pch = NULL) {
  
  cols <- as.character(cols)
  
  pch <- as.numeric(as.character(pch))
  
  
  if (is.null(leg.pch)) {
    if (!is.null(legend)) {
      leg.pch <- rep(1, length(legend))
    } else {
      leg.pch <- 1  # Default value if legend is also NULL
    }
  } else {
    leg.pch <- as.numeric(as.character(leg.pch))
    if (!is.null(legend) && length(leg.pch) != length(legend)) {
      stop("leg.pch must be the same length as legend")
    }
  }
  
  if (is.null(hi.index) | is.null(int.het)) 
    stop("error, input data are not provided")
  if (is.data.frame(hi.index)) 
    hi.index <- hi.index[, 2]
  
  if (pdf) 
    pdf(file = paste(out.file))
  
  plot(hi.index, int.het, col = cols, pch = pch, xlab = "Hybrid index", ylab = "Interspecific heterozygosity", 
       xlim = c(0, 1), ylim = c(0, 1))
  
  lines(c(0, 0.5), c(0, 1))
  lines(c(0.5, 1), c(1, 0))
  
  legend("topleft", legend = legend, col = leg.cols, pch = leg.pch, pt.cex = 1.5, cex = 1, xpd = 1, box.lty = 0, bg = "transparent")
  
  if (pdf) 
    dev.off()
}




# load data
load(paste0(DATA,"/KH_GP_filter.Rdata"))  
load(paste0(DATA,"/meta.Rdata"))

# read admixture results
q5 <- readQ(paste0(DATA,"/GP_hybrids.5.Q"))
rownames(q5[[1]]) <- MAF01$ind.names
q5 <- alignK(q5)

# add q matrix to metadata
admix5 <- q5$GP_hybrids.5.Q
admix5$id <- rownames(admix5)
colnames(admix5) <- c("Bulloo", "LEB", "MDB2", "Fitzroy", "MDB1", "id")
admix5 <- merge(meta, admix5, by = "id")


# identify pure and hybrid individuals
MDB.inds <- admix5[(admix5$MDB1 > 0.99 | admix5$MDB2 > 0.99), ]
Fitzroy.inds <- admix5[admix5$Fitzroy > 0.99, ]
ids_to_exclude.FITZ <- unique(c(MDB.inds$id, Fitzroy.inds$id, admix5$id[admix5$LEB > 0.99], admix5$id[admix5$Bulloo > 0.99]))
# potential hybrids (i.e. anything not ref FITZ, LEB, BULLOO or MDB)
hybrid.inds <- meta$id[meta$basin == "MDB" & !meta$id %in% ids_to_exclude.FITZ]



# prepare genlight object
Fitzroy.MDB.hybrid.gl <- gl.keep.ind(MAF01, ind.list = c(MDB.inds$id, Fitzroy.inds$id, hybrid.inds), recalc = TRUE, mono.rm = TRUE)
hybrid.meta.FITZ <- meta[meta$id %in% Fitzroy.MDB.hybrid.gl@ind.names, ]

# get inds with high HE
high_het_inds <- data.frame(id=hybrid.meta.FITZ$id[hybrid.meta.FITZ$basin_highHe == "HighHe"])



# set population status
status_map <- c(rep("MDB", length(MDB.inds$id)),
                rep("Fitzroy", length(Fitzroy.inds$id)),
                rep("hybrid", length(hybrid.inds)))
names(status_map) <- c(MDB.inds$id, Fitzroy.inds$id, hybrid.inds)
hybrid.meta.FITZ$status <- status_map[hybrid.meta.FITZ$id]


# add population info to genlight object
strata(Fitzroy.MDB.hybrid.gl) <- hybrid.meta.FITZ[, 2:7]
setPop(Fitzroy.MDB.hybrid.gl) <- ~status

# prepare data for gghybrid
gl2structure(Fitzroy.MDB.hybrid.gl, add.columns = Fitzroy.MDB.hybrid.gl@pop,
             outpath = getwd(), outfile = "FitzroyMDBhybrid.str")
ghyb.fitzroy <- gghybrid::read.data("FitzroyMDBhybrid.str", nprecol = 2, precol.headers = 0, 
                                    MISSINGVAL = -9, NUMINDS = 558, NUMLOCI = 4364, ONEROW = 0, POPID = 1, MARKERNAME = TRUE)

gg.fitzroy <- data.prep(
  data = ghyb.fitzroy$data,
  loci = ghyb.fitzroy$loci,
  sourceAbsent = FALSE,
  max.S.MAF = 0.2, min.diff = 0.6,
  alleles = ghyb.fitzroy$alleles,
  S0 = "MDB",
  S1 = "Fitzroy",
  precols = ghyb.fitzroy$precols,
  return.genotype.table = TRUE,
  return.locus.table = TRUE
)

# estimate hybrid index (will take an hour or so depending on your computer)
hindlabel1.FITZ <- esth(data.prep.object = gg.fitzroy$data.prep,
                        read.data.precols = ghyb.fitzroy$precols,
                        include.Source = TRUE, 
                        nitt = 5000,
                        burnin = 1000,
                        test.subject = "INDLABEL", 
                        test.subject.compare = "POPID")



# get candidate hybrids
hybrid.inds.FITZ <- subset(gg.HI.FITZ, Source == "TEST" & h_cred_int_lower>S0.limit.FITZ & h_cred_int_upper<S1.limit.FITZ)

# Write the results to file
write.csv(hybrid.inds.FITZ, "gghybrid_res.csv")


# generate gghybrid plot
{pdf(file = "Fitzroy_gghybrid.pdf", width = 12, height = 6)
  plot_h(data = hindlabel1.FITZ$hi,
         test.subject = hindlabel1.FITZ$test.subject,
         mean.h.by = "POPID",
         sort.by = c("mean_h", "POPID", "h_posterior_mode"),
         col.group = "POPID",
         group.sep = "POPID",
         fill.source = TRUE,
         basic.lines = FALSE,
         source.col = c("blue", "red"),
         source.limits = c("blue", "red"),
         cex = 1, pch = 16,
         cex.lab = 1.5, cex.main = 1.5, ylim = c(0,1))
  dev.off()}


# prepare data for triangle plot
gg.HI.FITZ <- as.data.frame(hindlabel1.FITZ$hi)
gg.HI.tri.FITZ <- gg.HI.FITZ[match(Fitzroy.MDB.hybrid.gl$ind.names, gg.HI.FITZ$INDLABEL), ]
stopifnot(all(gg.HI.tri.FITZ$INDLABEL == Fitzroy.MDB.hybrid.gl$ind.names))
S0.limit.FITZ <- max(gg.HI.FITZ[gg.HI.FITZ$Source == "S0", ]$h_posterior_mode)
S1.limit.FITZ <- min(gg.HI.FITZ[gg.HI.FITZ$Source == "S1", ]$h_posterior_mode)


hybrid.meta.FITZ <- hybrid.meta.FITZ %>%
  mutate(He_status = case_when(
    het_test %in% "HighHe" ~ "HighHe",   
    TRUE ~ status
  ))



# assign various He, hybrid status and overlap for plotting
HighHe.hybrids <- intersect(hybrid.inds.FITZ$INDLABEL, high_het_inds$id)

hybrid.meta.FITZ$alt.status <- hybrid.meta.FITZ$He_status
hybrid.meta.FITZ <- hybrid.meta.FITZ %>%
  mutate(alt.status = ifelse(alt.status == "hybrid", "unknown", alt.status))

hybrid.meta.FITZ <- hybrid.meta.FITZ %>%
  mutate(alt.status = ifelse(alt.status == "unknown" & id %in% hybrid.inds.FITZ$INDLABEL, "Hindex", alt.status))


hybrid.meta.FITZ$overlap <- hybrid.meta.FITZ$alt.status
hybrid.meta.FITZ <- hybrid.meta.FITZ %>%
  mutate(overlap = ifelse(overlap == "HighHe" & id %in% HighHe.hybrids, "Hindex_HighHe", overlap))

# check meta and gl inds match perfectly
stopifnot(all(hybrid.meta.FITZ$id == Fitzroy.MDB.hybrid.gl$ind.names))


# estimate interspecific heterozygosity
strata(Fitzroy.MDB.hybrid.gl) <- hybrid.meta.FITZ[ ,2:10]
setPop(Fitzroy.MDB.hybrid.gl) <- ~basin
gg.loci.FITZ <- gg.fitzroy$locus.data$locus
gg.gl.FITZ <- gl.keep.loc(Fitzroy.MDB.hybrid.gl, loc.list=gg.loci.FITZ)
gg.gi.FITZ <- gl2gi(gg.gl.FITZ)

# get allele counts
alleles <- gg.gi.FITZ@tab
snps <- t(alleles[, seq(1, ncol(alleles), 2)])

# est. int het
int.het.FITZ <- calc.intersp.het(introgress.data = snps)
# Write the results to file
write.csv(data.frame(id=hybrid.meta.FITZ$id, IntHet= int.het.FITZ), "IndHet.csv")


# plot
setPop(Fitzroy.MDB.hybrid.gl) <- ~overlap
alt.species <- as.data.frame(Fitzroy.MDB.hybrid.gl@pop)
colnames(alt.species) <- "species"
unique_species <- as.data.frame(sort(unique(as.factor(alt.species$species))))
col_species <- factor(unique_species[,1])

alt.fitz_colors <- c("#440154FF", "darkorange", "dodgerblue3", "limegreen", "limegreen", "#FDE725FF")
names(alt.fitz_colors) = levels(col_species)

alt.fitz_pch <- c(19, 3, 3, 1, 4, 19)
names(alt.fitz_pch) = levels(col_species)



# map colours to points
alt.ind_cols <- plyr::mapvalues(alt.species$species, names(alt.fitz_colors), alt.fitz_colors)
alt.ind_pch <- plyr::mapvalues(alt.species$species, names(alt.fitz_pch), alt.fitz_pch)


alt.leg.labs <- c("MDB_reference", "Fitzroy_reference", "MDB", "HighHe", "Hindex", "HighHe_hybrid")
alt.leg.cols <- c("#440154FF", "#FDE725FF", "dodgerblue3", "darkorange", "limegreen", "limegreen")
alt.leg.pch <- c(19, 19, 3, 3, 1, 4)



{triangle.plot.cb(hi.index=gg.HI.tri.FITZ$h_posterior_mode, 
                  int.het=int.het.FITZ, 
                  pdf=TRUE, 
                  out.file="Fitzroy_tri_plot.pdf", 
                  cols = alt.ind_cols, pch=alt.ind_pch, 
                  legend=alt.leg.labs, 
                  leg.cols=alt.leg.cols, 
                  leg.pch=alt.leg.pch)
}

