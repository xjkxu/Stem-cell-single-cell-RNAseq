# Load required libraries
library(data.table)
library(readxl)
library(GenomicRanges)
library(parallel)
library(vioplot)
library(ggplot2)
library(tidyr)
library(dplyr)
home.dir <- "~/Desktop/MaizeGWASResults/"
GWAS.results.dir <- "GWASresults/GLM/"
h2.dir <- "h2_analysis/Genic/GeneSetNonStemCell/"
setwd(home.dir)

distance<-"2kb"
if(distance=="2kb"){
  regulatory_distance<-2000L
}
if(distance=="gene_only"){
  regulatory_distance<-0L
}
for(g in c(1:6)){
  geneset<-g
results.dir <- paste0(h2.dir,"GeneSet",geneset)
if (!dir.exists(results.dir)) { dir.create(results.dir) }

# Load gene coordinates
full_trait_names <- scan(text ="CobDiameter	CobWeight	EarDiameter	EarLength		
                         EarRankNumber EarRowNumber	EarWeight	KernelWeight20 SeedSetLength", 
                         what = character())
excel_file_path <- "Table_for_Brian_GWAS_update_other_cell_types_markers.xlsx"
multiplesheets <- function(fname) { 
  # getting info about all excel sheets 
  sheets <- readxl::excel_sheets(fname) 
  tibble <- lapply(sheets, function(x) readxl::read_excel(fname, sheet = x)) 
  data_frame <- lapply(tibble, as.data.frame) 
  # assigning names to data frames 
  names(data_frame) <- sheets 
  # print data frame 
  print(data_frame) 
} 
wb <- multiplesheets(excel_file_path)
if(geneset=="Combined"){
  print("combinding all gene sets")
  genes_cords <- do.call(rbind, wb[-1])
} else{genes_cords <- as.data.frame(wb[[geneset]])}
# Create a dataframe from cords
cords <- data.frame(Elements = genes_cords$`V5-position`) %>%
  separate(Elements, into = c("Chr", "Range"), sep = ":", convert = TRUE) %>%
  separate(Range, into = c("Start", "End"), sep = "-", convert = TRUE)
genes_cords<-cbind( genes_cords,cords)
df <- genes_cords[,c("V5-ID","Chr", "Start", "End")]
#df <- data.frame(matrix(nrow = nrow(genes_cords), ncol = 4))
colnames(df) <- c("Gene","Chromosome", "Start", "End")
# Split the values and assign to dataframe columns
#df$Gene <- genes_cords$V5.final
#df$Chromosome <- sapply(strsplit(genes_cords$V5.Position, "[:-]"), "[", 1)
#df$Start <- as.numeric(sapply(strsplit(genes_cords$V5.Position, "[:-]"), "[", 2))
##df$End <- as.numeric(sapply(strsplit(genes_cords$V5.Position, "[:-]"), "[", 3))
df$Stand <- "*"
# Create genomic ranges object
Genes_GR <- makeGRangesFromDataFrame(df,na.rm=TRUE,
                                     keep.extra.columns=T,
                                     ignore.strand=T,
                                     seqnames.field="Chromosome",
                                     start.field="Start",
                                     end.field="End",
                                     strand.field="strand")
# Load GWAS file paths
suffix <- ".assoc"  # we want to grab the association result files

file_names <- list.files(path =GWAS.results.dir,pattern = paste0("*", suffix))
# Print the file names
print(file_names)
# Function to remove first 7 and last 6 characters from a string
remove_chars <- function(x) {
  substr(x, 8, nchar(x) - 6)
}
# Apply the function to each element in the vector
traits <- sapply(file_names, remove_chars)

# Function to load GWAS summary statistics and find overlaps
run_overlap_analysis <- function(file){
  gwas_data <- fread(paste0(GWAS.results.dir,file),header=T)
  # Remove monomorphic markers and minor allele freq < 0.05
  gwas_data <- subset(gwas_data,MAF>=0.05)
  gwas_data$Strand="*"
  # Create genomic ranges object
  GWAS_GR <- makeGRangesFromDataFrame(gwas_data,
                                      keep.extra.columns=T,
                                      ignore.strand=T,
                                      seqnames.field="Chromosome",
                                      start.field="Basepair",
                                      end.field="Basepair",
                                      strand.field="Strand")
  
  # Find overlaps
  # Find overlaps between two GenomicRanges objects (Genes_GR and GWAS_GR)
  #regulatory_distance=2000L
  overlaps <- findOverlaps(Genes_GR, GWAS_GR,
                           maxgap=regulatory_distance, 
                           minoverlap=0L,
                           type=c("any"),
                           ignore.strand=TRUE)
  
  
  # Extract overlapping ranges from Genes_GR and GWAS_GR based on overlaps
  Genes_overlaps <- Genes_GR[queryHits(overlaps)]
  GWAS_overlaps <- GWAS_GR[subjectHits(overlaps)]
  
  # Access metadata columns for the overlapping ranges
  Genes_metadata_columns <- mcols(Genes_overlaps)
  GWAS_metadata_columns <- mcols(GWAS_overlaps)
  
  # Create a dataframe (overlap_df) containing information from both Genes_GR and GWAS_GR
  overlap_df <- data.frame(
    Gene = Genes_metadata_columns$Gene,        # Extract the 'Gene' metadata column from Genes_GR
    Chr = seqnames(Genes_overlaps),            # Extract chromosome information from Genes_overlaps
    Gene_Start = start(Genes_overlaps),         # Extract start position of genes from Genes_overlaps
    Gene_End = end(Genes_overlaps),             # Extract end position of genes from Genes_overlaps
    Marker = GWAS_metadata_columns$Predictor,  # Extract the 'Predictor' metadata column from GWAS_GR
    Position = start(GWAS_overlaps),            # Extract start position of markers from GWAS_overlaps
    Pvalue = GWAS_metadata_columns$Wald_P      # Extract the 'Wald_P' metadata column from GWAS_GR
  )
  # The resulting dataframe (overlap_df) now contains information about overlapping genes and markers
  return(overlap_df)
}

# Run function over all traits in parallel
overlap_results <- mclapply(file_names,
                            run_overlap_analysis,
                            mc.cores = detectCores()
)
names(overlap_results) <- traits
for(n in seq_along(overlap_results)){
  print(length(unique(overlap_results[[n]]$Gene)))
}
if (!dir.exists(paste0(results.dir,"/ColocalizedMarkers"))) { dir.create(paste0(results.dir,"/ColocalizedMarkers")) }
for (i in seq_along(overlap_results)) {
  filename <- paste0(results.dir,"/ColocalizedMarkers/",traits[i],"_ColocalizedMarkers.txt")
  write.table(unique(overlap_results[[i]]$Marker), file = filename,quote=F, col.names = F, row.names = F)
}

for(i in names(overlap_results)){
  print(paste(i))
  if(i==names(overlap_results)[1]){
  merged_overlap_results<-  overlap_results[[i]]
  colnames(merged_overlap_results)[which(colnames(merged_overlap_results)=="Pvalue") ]<-paste0("GLM3PC_",i,".assoc")
  merged_overlap_results[[(paste0("GLM3PC_",i,".assoc_FDR"))]] <- p.adjust(merged_overlap_results[, paste0("GLM3PC_",i,".assoc")], method = "BH")
   } else{
     merged_overlap_results_tmp<-  overlap_results[[i]]
     colnames(merged_overlap_results_tmp)[which(colnames(merged_overlap_results_tmp)=="Pvalue") ]<-paste0("GLM3PC_",i,".assoc")
     merged_overlap_results_tmp[[(paste0("GLM3PC_",i,".assoc_FDR"))]] <- p.adjust(merged_overlap_results_tmp[, paste0("GLM3PC_",i,".assoc")], method = "BH")
  
    merged_overlap_results <-merge(merged_overlap_results,
                              merged_overlap_results_tmp[,c("Marker",paste0("GLM3PC_",i,".assoc"),
                                                        paste0("GLM3PC_",i,".assoc_FDR"))],
                              by="Marker",all=T) 
   }
  print(paste(dim(merged_overlap_results)))
}
print(paste0("/Users/brice7/Desktop/MaizeGWASResults/",names(wb)[g],"_",distance,"_GWASResults.csv"))
write.table(merged_overlap_results,paste0("/Users/brice7/Desktop/MaizeGWASResults/",names(wb)[g],"_",distance,"_GWASResults.csv"),
            row.names = F,quote = F,sep = ",")

#Get h2 for these marker sets
#plink --bfile combined_data --extract "ColocalizedMarkers/CD_ColocalizedMarkers.txt" --make-bed --out "ColocalizedMarkers/data"
#./ldak5.2.linux --calc-kins-direct "ColocalizedMarkers/kin" --bfile "ColocalizedMarkers/data" --power 1
#./ldak5.2.linux --reml "ColocalizedMarkers/reml_CD" --pheno reseq.pheno.fam --mpheno 1 --grm "ColocalizedMarkers/kin"
#rm ColocalizedMarkers/kin**
#  rm ColocalizedMarkers/data**

###Run 1000 permutations with random gene sets ###

library(rtracklayer)
# Maize version 5 all genes
refv5 <- as.data.frame(fread("Zm-B73-REFERENCE-NAM-5.0_Zm00001e.1_ab.1_xref_gene_IDs.txt",header=T))
refv5 <- refv5[refv5$new_gene_model_chr %in% paste("chr",1:10,sep=""), ]
refv5$new_gene_model_chr <- as.numeric(gsub("chr", "", refv5$new_gene_model_chr))
refv5 <- refv5[,6:11]
colnames(refv5) <- c("Chromosome","start","end","id","gene_model_biotype_overlap_90%","strand")
refv5$start <- as.numeric(refv5$start)
refv5$end <- as.numeric(refv5$end)
refv5$Chromosome <- as.character(refv5$Chromosome)
# Create granges object 
genome_gr <-  makeGRangesFromDataFrame(refv5,
                                       keep.extra.columns=T,
                                       ignore.strand=F,
                                       seqnames.field="Chromosome",
                                       start.field="start",
                                       end.field="end",
                                       strand.field="strand")
# load all markers
if (!file.exists("All_GWAS_markers_map.txt")) { 
  print("Writing marker table to txt file")
  markers <- fread(file_names[1],header=T)
  markers <- markers[,c("Chromosome","Predictor","Basepair")]
  markers$Strand <- "*"
  write.table(markers$Predictor,file="All_GWAS_markers.txt",quote=F,row.names = F,col.names = F)
  write.table(markers,file="All_GWAS_markers_map.txt",quote=F,row.names = F,col.names = T)
} else {
  print("Reading in marker table")
  markers <- fread(file="All_GWAS_markers_map.txt",header=T)
}
markers_gr <- makeGRangesFromDataFrame(markers,
                                       keep.extra.columns=T,
                                       ignore.strand=T,
                                       seqnames.field="Chromosome",
                                       start.field="Basepair",
                                       end.field="Basepair",
                                       strand.field="Strand")# Select number genes equal to num_samples
# function to create random set of markers
random_set <- function(genome_gr){
  num_samples <- length(unique(overlap_results[[1]]$Gene))
  random_samples <- genome_gr[sample(length(genome_gr), num_samples)]

  overlaps <- findOverlaps(markers_gr,random_samples,
                           maxgap=regulatory_distance, 
                           minoverlap=0L,
                           type=c("any"),
                           ignore.strand=TRUE)
  # Find overlaps
  # Find overlaps between two GenomicRanges objects (Genes_GR and GWAS_GR)
  # Extract overlapping ranges from Genes_GR and GWAS_GR based on overlaps
  markers_overlaps <- markers_gr[queryHits(overlaps)]
  # Access metadata columns for the overlapping ranges
  markers_metadata_columns <- mcols(markers_overlaps)
  #ran_markers <- gsub('"', '', markers_metadata_columns$Predictor)
  return(unique(markers_metadata_columns$Predictor))
}


# Run function k times 
k <- 1000  
if (!dir.exists(paste0(results.dir,"/Permutations"))) { dir.create(paste0(results.dir,"/Permutations")) }

for(t in traits){
  if (!dir.exists(paste(results.dir,"/Permutations/Random_Marker_Sets_", t, sep=""))) { dir.create(paste(results.dir,"/Permutations/Random_Marker_Sets_", t, sep="")) }
  
  #parallel_results <- mclapply(1:k, function(i) {  # Set seed for reproducibility
  
  #}, mc.cores =detectCores())
  
  # create directory for random marker sets if it does not already exists
  # Write table for each random marker set 
  mclapply(1:k, function(i) {
    parallel_results<-random_set(genome_gr) 
    write.table(file=paste0(results.dir,"/Permutations/Random_Marker_Sets_",t,"/RandomSet_",i,".txt"),
                parallel_results,quote=F,row.names = F,col.names = F)
  }, mc.cores =detectCores())
  #confirm k files were written
  print(t)
  print(length(list.files(paste0(results.dir,"/Permutations/Random_Marker_Sets_",t)))==k)
}
}#end of g loop

#subset in plink and get h2 from LDAK
# Loop from 1 to 1000
#for i in {1..1000}; do
# Run PLINK command for each file
#plink --bfile combined_data --extract "Random_Marker_Sets/RandomSet_${i}.txt" --make-bed --out "data_${i}"
#./ldak5.2.linux --calc-kins-direct "kin_${i}" --bfile "data_${i}" --power 1
#./ldak5.2.linux --reml "reml${i}" --pheno reseq.pheno --mpheno 1 --grm "kin_${i}"
#done
#sbatch randomh2_CD.sh
#sbatch randomh2_CW.sh
#sbatch randomh2_ED.sh
#sbatch randomh2_EL.sh
#sbatch randomh2_EW.sh
#sbatch randomh2_ERkN.sh
#sbatch randomh2_ERN.sh
#sbatch randomh2_KW20.sh
#sbatch randomh2_SSL.sh
#read in reml output for each trait

read_and_process_file <- function(r,results_dir) {
  #tmp <- fread(paste0(results_dir, "/ColocalizedMarkers/reml_", r, ".reml"), fill = TRUE)[13:20, ]
  #return(data.frame(trait = r, h2 = -2*as.numeric(tmp[2, 2])))
  
  tmp <- fread(paste0(results_dir, "/ColocalizedMarkers/reml_", r, ".reml"), fill = TRUE)
  
  return(data.frame(trait = r, h2 = as.numeric(tmp[18, 2])))
}

home.dir <- "~/Desktop/MaizeGWASResults/"
GWAS.results.dir <- "GWASresults/GLM/"
setwd(home.dir)
suffix <- ".assoc"  # we want to grab the association result files
file_names <- list.files(path =GWAS.results.dir,pattern = paste0("*", suffix))
# Print the file names
print(file_names)
# Function to remove first 7 and last 6 characters from a string
remove_chars <- function(x) {
  substr(x, 8, nchar(x) - 6)
}
# Apply the function to each element in the vector
traits <- sapply(file_names, remove_chars)
full_trait_names <- scan(text ="CobDiameter	CobWeight	EarDiameter	EarLength		
                         EarRankNumber EarRowNumber	EarWeight	KernelWeight20 SeedSetLength", 
                         what = character())
traits <-data.frame(traits=traits,full_trait_names=full_trait_names)

excel_file_path <- "Table_for_Brian_GWAS_update_other_cell_types_markers.xlsx"
multiplesheets <- function(fname) { 
  # getting info about all excel sheets 
  sheets <- readxl::excel_sheets(fname) 
  tibble <- lapply(sheets, function(x) readxl::read_excel(fname, sheet = x)) 
  data_frame <- lapply(tibble, as.data.frame) 
  # assigning names to data frames 
  names(data_frame) <- sheets 
  # print data frame 
  print(data_frame) 
} 
wb <- multiplesheets(excel_file_path)
# Define the list of genesets
genesets <- c( 1,2, 3, 4,5,6)
# Iterate over each geneset
for(d in c("Genic","2kb")){
  h2.dir <- paste0("h2_analysis/",d,"/GeneSetNonStemCell/")
  h2_dataframes <- list()
  #setwd(paste(h2.dir))
  for (geneset in genesets) {
    h2_dataframes[["trait"]] <- as.factor(traits$traits)
    results.dir.tmp <- paste0(h2.dir,"GeneSet",geneset)
    # Process h2 for the current geneset
    h2_list <- mclapply(traits$traits, read_and_process_file, results_dir = results.dir.tmp, mc.cores = detectCores())
    h2 <- do.call(rbind, h2_list)
    # Store h2 dataframe for the current geneset in the list
    h2_dataframes[[as.character(geneset)]] <- as.numeric(h2$h2)
  }
  # Combine the results into a single data frame
  if(d == "Genic"){
    combined_h2 <- as.data.frame(do.call(cbind, h2_dataframes))
    combined_h2$trait <- unlist(traits$traits)
    combined_h2$regions <- d
  }else{
    combined_h2_tmp <- as.data.frame(do.call(cbind, h2_dataframes))
    combined_h2_tmp$trait <- unlist(traits$traits)
    combined_h2_tmp$regions <- d
    combined_h2 <- rbind(combined_h2,combined_h2_tmp)
  }
}
combined_h2$control <- "NonStemCellGenes"

write.table(combined_h2)
# #read in reml output for each trait and permutation
process_permutations <- function(results_dir, t, num_permutations) {
  # Define a function to process each permutation
  process_permutation <- function(r, results_dir, t) {
    if (file.exists(paste0(results_dir, "/Permutations/", t, "/reml", r, ".reml"))) {
      tmp <- fread(paste0(results_dir, "/Permutations/", t, "/reml", r, ".reml"), fill = TRUE)[13:20, ]
      #h2_value <- -2*as.numeric(tmp[2, 2])
      #tmp <- fread(paste0(results_dir, "/Permutations/", t, "/reml", r, ".vars"), fill = TRUE)
      h2_value <-as.numeric(tmp[6, 2])
    } else {
      h2_value <- NA
    }
    return(h2_value)
  }
  # Use mclapply to parallelize the processing of permutations
  h2_perm <- mclapply(1:num_permutations, function(r) process_permutation(r, results_dir, t))
  # Convert the list of results into a data frame
  h2_perm_df <- data.frame(h2 = unlist(h2_perm))
  return(h2_perm_df)
}

traits
genesets <- c(1, 2, 3, 4,5,6)
permutations_list <- list()
gene_space <- c("Genic","2kb")
for(gs in gene_space){
  h2.dir <- paste0("h2_analysis/",gs,"/GeneSetNonStemCell/")
  for(g in genesets){
    print(paste(gs,g))
    results.dir <- paste0(h2.dir,"GeneSet",g)
    h2_perm_list <- lapply(traits$traits, function(t) {
      process_permutations(results.dir, t, 1000)
    })
    h2_perm_df <- do.call(cbind, h2_perm_list)
    colnames(h2_perm_df) <- traits$traits
    
    #wb$`3_Stem cell DEGs-IM-V5only`#Plot
    melt_permutation <- melt(setDT(h2_perm_df))
    melt_permutation$geneset <- g
    melt_permutation$genespace <- gs
    permutations_list[[g]] <- melt_permutation
  }
  if(gs == gene_space[1]){
    permutations <- do.call(rbind,permutations_list)
  } else{
    permutations_tmp <- do.call(rbind,permutations_list)
    permutations <- rbind(permutations,permutations_tmp)
  }
}
permutations.tmp <- permutations
permutations$region<-as.factor(permutations$genespace)
str(permutations)
permutations$region <- factor(permutations$region, levels = rev(levels(permutations$region)))
permutations <- merge(permutations,traits,by.x="variable",by.y="traits",all.x=T)
combined_h2_melt <- melt(as.data.table(combined_h2))
combined_h2_melt$regions <- factor(combined_h2_melt$regions)

permutations$ID <- paste(permutations$variable,permutations$genespace,permutations$geneset,sep="_")
combined_h2_melt$ID <-paste(combined_h2_melt$trait,combined_h2_melt$regions,combined_h2_melt$variable,sep="_")
permutations <- merge(permutations, combined_h2_melt,by="ID")

#get percentiles
ID <- unique(combined_h2_melt$ID)
ecdf_function <- ecdf(permutations$value.x[which(permutations$ID==ID)])
ecdf_function(combined_h2_melt$value[which(combined_h2_melt$ID==ID)]) * 100


percentile_df <- data.frame(ID = character(),h2 = numeric(), Percentile_95 = numeric(), 
                            Percentile_90 = numeric(), stringsAsFactors = FALSE)
unique_IDs <- unique(combined_h2_melt$ID)
for (ID in unique_IDs) {
  print(paste(ID))
  # Subset the values for the current ID from permutations and combined_h2_melt
  permutations_subset <- permutations$value.x[which(permutations$ID == ID)]
  if(length(which((is.na(permutations_subset )==T)))>0){  # If all values were NA
    print(paste("Skipping ID:", ID, "due to NA values"))
    next
  }  
  mean_value <- mean(permutations$value.x[permutations$ID == ID], na.rm = TRUE)
  
  if (is.na(mean_value) || is.nan(mean_value)) {  # If all values were NA
    print(paste("Skipping ID:", ID, "due to NA values"))
    next
  }  
  combined_h2_melt_subset <- combined_h2_melt$value[which(combined_h2_melt$ID == ID)]
  # Calculate the ECDF function for the current ID
  ecdf_function <- ecdf(permutations_subset)
  # Calculate the percentile of each value in combined_h2_melt for the current ID
  percentiles <- ecdf_function(combined_h2_melt_subset) * 100
  # Calculate the 95th and 90th percentiles
  percentile_95 <- quantile(permutations_subset, 0.95)
  percentile_90 <- quantile(permutations_subset, 0.90)
  # Create a dataframe with ID and Percentile
  ID_percentile_df <- data.frame(ID = rep(ID, length(combined_h2_melt_subset)),
                                 h2 = combined_h2_melt_subset,
                                 Percentile = percentiles, 
                                 Percentile_95 = percentile_95, 
                                 Percentile_90 = percentile_90)
  # Append the results to the main dataframe
  percentile_df <- rbind(percentile_df, ID_percentile_df)
}
#write.csv(percentile_df, "percentile_2log_data.csv", row.names = FALSE)

# Plot violin plot per gene set 
permutations <- merge(permutations,percentile_df[,c("ID","Percentile_95")])

combined_h2_melt <-  merge(combined_h2_melt ,percentile_df[,c("ID","Percentile_95")],by="ID")
tmp<-combined_h2_melt[combined_h2_melt$trait=="CD",]

perm_plot <- ggplot(permutations, aes(x = value.x, y = as.factor(geneset), fill = genespace)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  facet_grid(variable.x ~ genespace, scales = "free_y") +  # Row: variable.x, Column: genespace
  theme_minimal() +
  theme(
    panel.spacing = unit(1, "cm"),  # Space between facets
    strip.text = element_text(size = 7, face = "bold"),  # Improve readability of facet labels
    axis.text.y = element_text(size = 5, face = "bold")  # Adjust y-axis text size
    
  ) +
  labs(
    title = "",
    x = "h2",
    y = "Geneset"
  )


save(permutations,file="NonStemCellGeneSeth2Permuations.Rdata")



for(i in genesets){
  permutations_tmp <- permutations[which(permutations$geneset==i),]
  #permutations_tmp$region <- factor(permutations_tmp$region, levels = rev(levels(permutations_tmp$region)))
  permutations_tmp <-  permutations_tmp[order(permutations_tmp$trait, permutations_tmp$region), ]
  unique_IDs_tmp <-unique(permutations_tmp$ID)
  # Create a new column with index positions of matching elements
  permutations_tmp$index <- match(permutations_tmp$ID, unique_IDs_tmp)
  # Reorder the dataframe based on the order of the matching vector
  permutations_tmp <-  permutations_tmp[order( permutations_tmp$index), ]
  # Reorder the dataframe based on the order of the matching vector
  if(i != "Combined"){
    title <- names(wb)[i]
  } else if (i == "Combined"){
    title <- "All_Gene_Sets_Combined"
  }
  # Calculate unique trait names and regions
  unique_traits <- unique(permutations_tmp$trait)
  unique_regions <- unique(permutations_tmp$region)
  # Determine the number of unique traits and regions
  num_traits <- length(unique_traits)
  num_regions <- length(unique_regions)
  # Create positions for violins
  positions <- seq(0.5, num_traits * num_regions - 0.5, length.out = num_traits * num_regions)
  # Create a vector from 1 to 18
  vector <- 1:c(num_traits*num_regions)
  # Create a vector of colors alternating for each value
  colors <- ifelse(vector %% 2 == 1, "#66c2a4", "#fc8e62")
  # Create violin plot
  jpeg(paste0(title,"_varG_results.jpeg"), width = 8, height = 3.5, units = "in",res=300)
  par(xpd=TRUE)
  vioplot(
    permutations_tmp$value.x ~ permutations_tmp$index),
    col = colors,
    names = rep(unique_traits, each = length(unique_regions)),
    at = positions,
    horizontal = FALSE,
    main = title,
    xaxt = "n",
    xlab = "Trait",
    ylab = "sigma_g",font.axis = 2,
    drawRect = FALSE  # Remove the white dot at the mean
  )
  
  # Customize axis labels
  axis(side = 1, at = seq(1, length(unique_traits) * length(unique_regions), by = 2), labels = unique_traits,font.axis = 2 )
  positions_h2 <- seq(0.5, length(unique_traits) * length(unique_regions) - 0.5, length.out = length(unique_traits) * length(unique_regions))
  
  combined_h2_tmp <- combined_h2_melt[combined_h2_melt$ID %in% unique_IDs_tmp, ]
  # Create a new column with index positions of matching elements
  combined_h2_tmp$index <- match(combined_h2_tmp$ID, unique_IDs_tmp)
  # Reorder the dataframe based on the order of the matching vector
  combined_h2_tmp <- combined_h2_tmp[order(combined_h2_tmp$index), ]
  points(positions_h2, combined_h2_tmp$value, col = "black", pch = 19, cex = 0.6)
  points(positions_h2, combined_h2_tmp$Percentile_95, col = "red", pch = "-", cex = 1.5)
  legend("bottomright", 
         legend = c("Genic", "2kb"), 
         fill = c("#66c2a4", "#fc8e62"),
         title = "",
         horiz = TRUE,bg = "transparent",  # Remove legend background
         box.lwd = 0)
  dev.off()
  gc()
}

#########plot each trait with its own axis
library(gridExtra)
unique_traits <- unique(permutations_tmp$traits)
plot_list <- list()  # Empty list to store each trait-specific plot
for(i in genesets){
  for(t in unique(permutations$variable.x)){
    # Filter data for the current geneset
    permutations_tmp <- permutations[permutations$geneset == i, ]
    permutations_tmp <- permutations_tmp[permutations_tmp$variable.x == t, ]
    
    permutations_tmp <- permutations_tmp[order(permutations_tmp$traits, permutations_tmp$region), ]
    
    combined_h2_tmp <- combined_h2_melt[combined_h2_melt$ID %in% unique(permutations_tmp$ID), ]
    # Create a new column with index positions of matching elements
    combined_h2_tmp$index <- match(combined_h2_tmp$ID, unique(permutations_tmp$ID))
    
    unique_IDs_tmp <-unique_IDs_tmp <- unique(permutations_tmp$ID)
    # Create a new column with index positions of matching elements
    permutations_tmp$index <- match(permutations_tmp$ID, unique_IDs_tmp)
    # Reorder the dataframe based on the order of the matching vector
    permutations_tmp <-  permutations_tmp[order( permutations_tmp$index), ]
    
    permutations_tmp$fill_factor <- factor(permutations_tmp$index %% 2)
    
    # Create the violin plot for each trait
    p <- ggplot(permutations_tmp, aes(x = factor(index), y = value.x, fill = factor(index %% 2))) +
      geom_violin(trim = FALSE, show.legend = FALSE) +
      scale_fill_manual(values = colors) +
      geom_point(data=combined_h2_tmp,show.legend = FALSE,
                 aes(x = factor(index), y = value), color = "black", size = 3) +
      geom_point(data=combined_h2_tmp,show.legend = FALSE,
                 aes(x = factor(index), y = Percentile_95),  col = "red", pch = "-", cex = 8) +
      labs(title = paste(t), x = "Geneset", y = "sigma_2_g") +
      scale_x_discrete(labels = c("1" = "Genic", "2" = "2kb")) +  # Custom labels
      
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "bold"),
            plot.title = element_text(size = 14, face = "bold"))
    #print(p)
    
    # Save each individual plot to the list
    plot_list[[t]] <- p
  }
  # Arrange all trait-specific plots for the current geneset in a grid layout and save as a JPEG
  jpeg(paste0(i, "_sigma_2_g.jpeg"), width = 10, height = 8, units = "in", res = 300)
  grid.arrange(grobs = plot_list, ncol = 3)  # Adjust `ncol` for layout (e.g., 2 columns)
  dev.off()
}





