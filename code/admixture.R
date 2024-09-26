# admixture analyses to reproduce figure S11 q4_He, q5_He, and q4_He_SL

# Load required libraries
library(adegenet)
library(pophelper)
library(viridis)
library(stringr)

# set working directory and location of data files
setwd("/path/to/your/working/directory")
DATA <- "/path/to/data"

# Load data
load(paste0(DATA,"/KH_GP_filter.Rdata"))  
load(paste0(DATA,"/meta.Rdata"))
load(paste0(DATA,"/plot_order.Rdata"))


# set and check metadata rows match the gl object sample order
meta <- meta[match(MAF01$ind.names, meta$id), ]
stopifnot(all(meta$id == MAF01$ind.names))

# set group labels and order for plotting
grp.labs <- as.data.frame(meta$het_test)
colnames(grp.labs) <- "grp"
grp.order <- as.character(plot_order[,1])

# function to run admixture
run_admixture <- function(gl, output.name, K, nproc, supervised = FALSE) {
  # Check that plink and admixture are installed and can be accessed in $PATH
  chk.exe <- function(exe) { if (Sys.which(exe) == "") stop(paste(exe, "not found in PATH")) }
  chk.exe("plink")
  chk.exe("admixture")
  
  # function to build ped and map files
  generate_ped_map_files <- function(genlight_obj, output.name) {
    # extract individual names
    individual_ids <- genlight_obj@ind.names
    
    # initialize ped file structure
    ped_data <- data.frame(
      FID = individual_ids,  # family id (can be the same as individual id)
      IID = individual_ids,  # individual id
      PID = 0,  # paternal id
      MID = 0,  # maternal id
      Sex = 0,  # sex (0 = unknown, 1 = male, 2 = female)
      Phenotype = 0  # phenotype (0 = missing, 1 = unaffected, 2 = affected)
    )
    
    # initialize a matrix to hold all genotype data
    num_individuals <- length(individual_ids)
    num_loci <- length(genlight_obj@loc.names)
    
    # function to convert genotypes to allele pairs
    convert_genotypes <- function(genotypes) {
      allele_pairs <- sapply(genotypes, function(x) {
        if (is.na(x)) c(0, 0)  # missing data
        else if (x == 0) c(1, 1)  # homozygous reference allele
        else if (x == 1) c(1, 2)  # heterozygous
        else if (x == 2) c(2, 2)  # homozygous alternate allele
        else c(0, 0)  # unexpected value, treat as missing
      })
      as.vector(allele_pairs)
    }
    
    # use lapply to process all individuals at once
    genotype_matrix <- do.call(rbind, lapply(genlight_obj@gen, function(x) convert_genotypes(as.integer(x))))
    
    # combine the ped data and the genotype matrix
    ped_data <- cbind(ped_data, genotype_matrix)
    
    # write ped file
    write.table(ped_data, file = paste0(output.name, ".ped"), sep = " ", quote = FALSE, row.names = FALSE, col.names = FALSE)
    
    # create the map file using @loc.names
    map_data <- data.frame(
      Chromosome = 0,  # chromosome number (placeholder)
      SNP = genlight_obj@loc.names,  # snp name
      GeneticDistance = 0,  # genetic distance (placeholder)
      Position = 0  # physical position (placeholder)
    )
    
    # write map file
    write.table(map_data, file = paste0(output.name, ".map"), sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    
    cat("ped and map files generated successfully.\n")
  }
  
  # Generate .ped and .map files
  generate_ped_map_files(gl, output.name)
  
  # Use plink to generate bed files
  system2("plink", args = c('--file', output.name, '--make-bed', '--out', output.name, "--allow-no-sex"))
  
  # Run ADMIXTURE
  if (supervised) {
    system2("admixture", args = c("--supervised", paste0(output.name, ".bed"), K, paste0("-j", nproc)), 
            stdout = paste0("K", sprintf("%02d", K), "_supervised.out"))
  } else {
    for (i in 1:K) {
      system2("admixture", args = c("--cv=10", paste0(output.name, ".bed"), i, paste0("-j", nproc)), 
              stdout = paste0("K", sprintf("%02d", i), ".out"))
    }
    
    # Get best K based on lowest CV error
    files <- list.files(pattern = "K[0-9]+.out")
    cv <- sapply(files, function(file) {
      as.numeric(str_extract(grep("CV error", readLines(file), value = TRUE), "\\d+\\.\\d+"))
    })
    write.table(cv, paste0(output.name, "_CV_error.txt"), col.names = FALSE, quote = FALSE)
    
    pdf(file = paste0(output.name, "_CV_error.pdf"), width = 6, height = 6)
    plot(seq_along(cv), cv, xlab = "K", ylab = "CV error")
    points(which.min(cv), min(cv), col = "red", lwd = 3, pch = 8)
    dev.off()
  }
}

# function to read Q files and align them
read_and_align_Q <- function(file_path) {
  q <- readQ(file_path)
  rownames(q[[1]]) <- MAF01$ind.names
  q <- alignK(q)
  return(q)
}

# function to plot Q matrix
plot_Q <- function(q, k, output_filename) {
  plotQ(q, 
        clustercol = viridis(k), 
        sharedindlab = FALSE, 
        showindlab = FALSE, 
        useindlab = TRUE, 
        indlabsize = 0.6, 
        indlabspacer = 0, 
        grplabsize = 1, 
        grplabangle = 60, 
        subsetgrp = grp.order,
        indlabangle = 60, 
        indlabvjust = 1, 
        grplab = grp.labs, 
        grplabpos = 0.7, 
        grplabheight = 0.4, 
        grplabspacer = 0, 
        ordergrp = TRUE, 
        showsp = FALSE, 
        linepos = 0.9, 
        splab = "", 
        height = 4,
        divsize = 0.1, 
        showyaxis = FALSE, 
        width = 40, 
        imgtype = "pdf", 
        outputfilename = output_filename, 
        exportpath = getwd())
}


# run unsupervised admixture for K2-20 (can skip and use results in the DATA dir)
# run_admixture(MAF01, "GP_hybrids", 20, 4)

# Plot q4_He and q5_He
# q4 <- read_and_align_Q("GP_hybrids.4.Q") # use this if you have run admixture yourself
q4 <- read_and_align_Q(paste0(DATA,"/GP_hybrids.4.Q"))
plot_Q(q4, 4, "q4_He")

# q5 <- read_and_align_Q("GP_hybrids.5.Q") # use this if you have run admixture yourself
q5 <- read_and_align_Q(paste0(DATA,"/GP_hybrids.5.Q"))
plot_Q(q5, 5, "q5_He")

# prepare data for supervised admixture
admix5 <- q5$GP_hybrids.5.Q
colnames(admix5) <- c("Bulloo", "LEB", "MDB2", "Fitzroy", "MDB1")
admix5 <- cbind(meta, admix5)

# define reference populations
ref_pops <- list(
  MDB_ref = admix5$id[admix5$MDB1 > 0.99 | admix5$MDB2 > 0.99],
  Fitzroy_ref = admix5$id[admix5$Fitzroy > 0.99],
  LEB_ref = admix5$id[admix5$LEB > 0.99],
  Bulloo_ref = admix5$id[admix5$Bulloo > 0.99]
)

# create reference population file
reference_pop <- data.frame(
  id = MAF01$ind.names,
  pop = rep("-", length(MAF01$ind.names))
)

for (pop in names(ref_pops)) {
  reference_pop$pop[reference_pop$id %in% ref_pops[[pop]]] <- pop
}

# write the population file
write.table(reference_pop$pop, file = "GP_hybrids.pop", row.names = FALSE, col.names = FALSE, quote = FALSE)

# run supervised admixture, (again can skip and use results in the DATA dir)
# run_admixture(MAF01, "GP_hybrids_ref", 4, 4, supervised = TRUE)

# Plot q4_He_SL
# q4_SL <- read_and_align_Q("GP_hybrids_ref.4.Q") # use this if you have run admixture yourself
q4_SL <- read_and_align_Q(paste0(DATA,"/GP_hybrids_ref.4.Q"))
plot_Q(q4_SL, 4, "q4_He_SL")


