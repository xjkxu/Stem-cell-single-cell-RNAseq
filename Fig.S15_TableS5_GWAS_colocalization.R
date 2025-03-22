
library(ggplot2)
library(data.table)
library(RColorBrewer)
library(dplyr)
library(parallel)
library(GenomicRanges)
library(ggrepel)

home.dir <- "~/Desktop/MaizeGWASResults/"
setwd(home.dir)
GWAS.results.dir <- "GWASresults/GLM/"
if (!dir.exists(paste0(GWAS.results.dir,"ManhattanPlots"))) { dir.create(paste0(GWAS.results.dir,"ManhattanPlots")) }

# Load your GWAS results data 
# Load GWAS file paths
suffix <- ".assoc"  # we want to grab the association result files
file_names <- list.files(path =GWAS.results.dir,pattern = paste0("*", suffix))
remove_chars <- function(x) {
  substr(x, 8, nchar(x) - 6)
}
# Function to read a GWAS result file and extract predictor names and p-values
read_gwas_file <- function(file_path) {
  gwas_data <- fread(paste0(GWAS.results.dir,file_path), header = TRUE)  # Adjust read function based on file format
  #remove maf 0.05 pvalues (set to NA)
  rows_to_modify <- gwas_data$MAF < 0.05
  gwas_data$Wald_P[rows_to_modify] <- NA
  return(data.frame(Predictor = gwas_data$Predictor, Pvalue = gwas_data$Wald_P))
}

# Read and combine predictor names and p-values from all files
combined_data <- mclapply(file_names, read_gwas_file)

# Combine into a single dataframe and confirm they are combinding the same markers
combined_df <- do.call(cbind, combined_data)
predictors<- combined_df[, seq(1, ncol(combined_df), by = 2)]
are_equal <- all(apply(predictors, 2, function(x) all(x == predictors[, 1])))
# Print result
if (are_equal) {
  cat("All odd-numbered columns are equal vectors.\n")
} else {
  cat("Not all odd-numbered columns are equal vectors.\n")
}
rm(predictors)

pvalues <- combined_df[, c(1,seq(2, ncol(combined_df), by = 2))]
# Add a column with file names as identifier
colnames(pvalues)[2:ncol(pvalues)] <- file_names

#Add in gene info
geneset="Combined"
results.dir <- paste0("GeneSet",geneset)
if (!dir.exists(results.dir)) { dir.create(results.dir) }

# Load gene coordinates
excel_file_path <- "Table_for_Brian_GWAS_update_v5.xlsx"
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
wb$`Read Me`<-NULL
excel_file_path <- "Table_for_Brian_GWAS_update_other_cell_types_markers.xlsx"
wb2 <- multiplesheets(excel_file_path)
wb <- c(wb,wb2)
rm(wb2)

##########loop through all genesets
regulatory_distance<-2000L

for(geneset in c(1:length(wb),"Combined")){
 if(geneset=="Combined"){
  print("combinding all gene sets")
  genes_cords <- do.call(rbind, wb[1:4])
} else{genes_cords <- as.data.frame(wb[[geneset]])} 

# Create a dataframe
df <- data.frame(matrix(nrow = nrow(genes_cords), ncol = 4))
colnames(df) <- c("Gene","Chromosome", "Start", "End")
# Split the values and assign to dataframe columns
df$Gene <- genes_cords$V5.final
df$Chromosome <- sapply(strsplit(genes_cords$V5.Position, "[:-]"), "[", 1)
df$Start <- as.numeric(sapply(strsplit(genes_cords$V5.Position, "[:-]"), "[", 2))
df$End <- as.numeric(sapply(strsplit(genes_cords$V5.Position, "[:-]"), "[", 3))
df$Stand <- "*"
df <- df[!duplicated(df$Gene),]

# Create genomic ranges object
Genes_GR <- makeGRangesFromDataFrame(df,na.rm=TRUE,
                                     keep.extra.columns=T,
                                     ignore.strand=T,
                                     seqnames.field="Chromosome",
                                     start.field="Start",
                                     end.field="End",
                                     strand.field="strand")

#load GWAS summary statistics and find overlaps
if (file.exists("All_GWAS_markers_map.txt")) { 
  print("Reading in marker table")
  markers <- fread(file="All_GWAS_markers_map.txt",header=T)
} else {
  print("Cannot find marker table")
}

# Create genomic ranges object
GWAS_GR <- makeGRangesFromDataFrame(markers,
                                    keep.extra.columns=T,
                                    ignore.strand=T,
                                    seqnames.field="Chromosome",
                                    start.field="Basepair",
                                    end.field="Basepair",
                                    strand.field="Strand")

# Find overlaps between two GenomicRanges objects (Genes_GR and GWAS_GR)
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
  Position = start(GWAS_overlaps)           # Extract start position of markers from GWAS_overlaps
)
#confirm no duplicate gene marker combos
tmp <- overlap_df 
tmp$tmp <- paste0(tmp$Gene,tmp$Marker)
print(paste("Length of unique marker gene combo is equal to nrow of dataframe:",length(unique(tmp$tmp))==nrow(tmp)))
rm(tmp)


#find markers in multiple genes
duplicates <- duplicated(overlap_df$Marker,na.rm = TRUE) | duplicated(overlap_df$Marker, fromLast = TRUE)
overlap_df$tmp <- NA

# Sort and collapse the 'Gene' column for each 'Marker'
overlap_df <- overlap_df %>%
  group_by(Marker) %>%
  mutate(tmp = paste0(sort(Gene), collapse = ";"))

# Remove duplicate rows and keep unique combinations of 'Marker' and 'tmp'
overlap_df_duplicates_markers_removed <- unique(select(overlap_df, Marker, tmp))

#remove duplicates
overlap_df_duplicates_markers_removed <- overlap_df_duplicates_markers_removed[!duplicated(overlap_df_duplicates_markers_removed$Marker),]
overlap_df <- merge(overlap_df_duplicates_markers_removed,pvalues,by.x="Marker",by.y="Predictor",all.y=T)

rm(markers,pvalues,GWAS_GR,GWAS_metadata_columns,
   GWAS_overlaps,df,combined_df,combined_data,wb,overlaps,genes_cords,
   Genes_GR,Genes_overlaps,Genes_metadata_columns,overlap_df_duplicates_markers_removed)


#save(overlap_df,file=paste0(GWAS.results.dir,"CandidateGenes_GLM_GWASoutput_overlap_2kbregulatory.Rdata"))
#load(paste0(GWAS.results.dir,"CandidateGenes_GLM_GWASoutput_overlap_2kbregulatory.Rdata"))
#colnames(overlap_df)[2] <- "Gene"

#add in columns for which genes pass FDR and bonferoni thresholds
process_file <- function(overlap_df, col_name) {  
  col_index <-which(colnames(overlap_df)==col_name)
  #get FDR
  tmp <- overlap_df[, c(1:2, col_index)]
  tmp <- tmp[complete.cases(tmp[col_name]), ]
  tmp[[paste0(col_name,"_FDR")]] <- p.adjust(tmp[, col_name], method = "BH")
  tmp <- tmp[complete.cases(tmp["tmp"]), ]
  colnames(tmp)
  return(tmp)
  #bonferoni threshold
}
FDR_results <- mclapply(file_names, function(col_name) {
  process_file(overlap_df, col_name)
}) 

for(i in 1:length(FDR_results)){
  print(i)
  if(i==1){
    FDR_results_full <- FDR_results[[1]]
  } else {
    FDR_results_tmp <- FDR_results[[i]]
    FDR_results_full <- merge(FDR_results_full,FDR_results_tmp[,c(1,3,4)],by="Marker",all=T)
  }
}

#rm(FDR_results_tmp,FDR_results)
#write.table(FDR_results_full,file="~/Google Drive/My Drive/MaizeData/EarTraits_Xu_SingleCell/Genes_GWAS_overlap_with_FDR_values_2kb.csv",
#            row.names = F,col.names = T,quote = F, sep = "\t")
#merge with reference genome genes
#run script maizev5genes_WGS_overlap.R
load("maizev5genes_WGS_overlap_genic_2kbregulatory.Rdata")
overlap_df_ref <- overlap_df_ref[!duplicated(overlap_df_ref$Marker),]

overlap_df <- merge(overlap_df_ref[,c("Gene","Marker")],overlap_df ,by.y="Marker",by.x="Marker",all.y=T)
dim(overlap_df)
length(unique(overlap_df$Marker))
duplicates <- duplicated(overlap_df$Marker,na.rm = TRUE) | duplicated(overlap_df$Marker, fromLast = TRUE)
length(which(duplicates==T))

#chi-square to find out if number of genes with hits is higher than expected
# Create a contingency table

# Count the number of items in column1 less than the constant
result_df <-  data.frame(
  Trait = character(0),
  UniqueGeneCount  = numeric(0),
  RefGeneCount = numeric(0)
)

#### get gene count
  #geneset total
  #geneset passing fdr
  #fullgenome total
  #fullgenome passing fdr

gene_colocalization_count <- function(col_index) {
  constant_value <- 0.05
  col_name <- names(overlap_df)[col_index]
  tmp <- overlap_df[, c(1:3, col_index)]
  tmp <- tmp[complete.cases(tmp[col_name]), ]
  tmp$FDR <- p.adjust(tmp[, col_name], method = "BH")
  tmp <- tmp[which(tmp[,"FDR"] <= constant_value), ]
  
  unique_gene_count <- length(unique(na.omit(tmp$Gene)))
  unique_gene_count_ref <- length(unique(na.omit(tmp$Gene)))
  
  result_df <- data.frame(Trait = col_name,
                          UniqueGeneCount = unique_gene_count,
                          RefGeneCount = unique_gene_count_ref)
  return(result_df)
}
#gene_colocalization_count(col_index=8)
gene_colocalization_list <- mclapply(4:(ncol(overlap_df)),
                                     gene_colocalization_count,
                                     mc.cores =detectCores())
gene_colocalization_df <- do.call(rbind, gene_colocalization_list)
target_genes <- na.omit(unique(overlap_df$Gene.y))
target_genes<- target_genes[-c(grep(";", target_genes))]
gene_colocalization_df$TotalTargetGenes <- length(target_genes)
total_genes <- na.omit(unique(overlap_df$Gene.x))
gene_colocalization_df$TotalRefGenes <- length(total_genes)
gene_colocalization_df <- gene_colocalization_df[,c("Trait","UniqueGeneCount",
                                                    "TotalTargetGenes","RefGeneCount",
                                                    "TotalRefGenes")]
gene_colocalization_df$ChiSquare <- NA
for(i in 1:nrow(gene_colocalization_df)){
  df1 <- gene_colocalization_df[i,2:3]
  df2 <- gene_colocalization_df[i,4:5]
  colnames(df1) <- c("observed","total")
  colnames(df2) <- c("observed","total")
  observed_data <- rbind(df1,df2) 
  test <- chisq.test(observed_data)
  gene_colocalization_df$ChiSquare[i] <- test$p.value
}
gene_colocalization_df$UniqueGeneCountRatio <-gene_colocalization_df$UniqueGeneCount/gene_colocalization_df$TotalTargetGenes
gene_colocalization_df$RefGeneCountRatio <-gene_colocalization_df$RefGeneCount/gene_colocalization_df$TotalRefGenes

write.table(gene_colocalization_df,file="~/Google Drive/My Drive/MaizeData/EarTraits_Xu_SingleCell/gene_colocalization_count_chisq_2kb.csv",row.names = F,col.names = T,quote=F)
#downsize for testing
#results <- results[sample(1:nrow(results),5000),]
#file <- file_names[1]

# Define colors



create_manhattan_plot <- function(file_path) {
  # Read data
  print(paste("Reading in GWAS results for file",file_path)) 
  
  results <- fread(paste0(GWAS.results.dir,file_path), header = TRUE)
  results <- subset(results, MAF >= 0.05)
  results$log_10 <- -log10(results$Wald_P)
  colnames(results)[which(colnames(results) == "Chromosome")] <- "Chr"
  results$Chr <- as.factor(results$Chr)
  results$FDR <- p.adjust(results$Wald_P, method = "BH")
  if(threshold_type=="FDR"){
    threshold_line <- results[which(results$FDR<= threshold),]
    threshold_line <- -log10(threshold_line$Wald_P[which.max(threshold_line$FDR)])
  }else if (threshold_type=="Bon"){
    threshold_line <- -log10(threshold/nrow(results))
  } else {
    return("No Threshold Given")
  }
  #tmp <- overlap_df[complete.cases(overlap_df["Gene"]),]
  #tmp$group <- "Candidate"
  #tmp <- tmp[,c("Marker","group")]
  #results <- merge(results,tmp,by.x="Predictor",by.y="Marker",all.x=T)
  #results$group[which(is.na(results$group)==T)]="NonCandidate"
  #results$alpha <- ifelse(results$group == "Candidate", 1, 0.2)
  #results$alpha <- ifelse(results$group == "NonCandidate", 0.2,1)
  #results <- results[which(results$group=="Candidate"),]
  overlap_df_test<-overlap_df[complete.cases(overlap_df["Gene2"]),]
  overlap_df_test<-overlap_df_test[complete.cases(overlap_df_test[file_path]),]
  overlap_df_test$log10 <- -log10(overlap_df_test[file_path])
  hits <- overlap_df_test[overlap_df_test$log10 >= threshold_line,]
  hits <- merge(hits[,c("Marker","Gene2" )],results[,c("Chr","Predictor","Basepair","log_10")],by.x="Marker",by.y="Predictor")
  # Sort dataframe descending order
  hits <- hits[order(-hits$log_10), ]
  # Remove duplicates in col1 keeping the first occurrence (highest value in col2)
  hits <- hits[!duplicated(hits$Gene), ]
  
  # Create Manhattan plot
  manplot <- ggplot(data = results, aes(x = Basepair / 10^6, y = log_10, color = Chr)) +
    #geom_point(aes(alpha = factor(group))) +
    geom_point() +
    #scale_alpha_manual(values = c("Candidate" = 1, "NonCandidate" = 0.2), guide = FALSE) + # Manually set alpha values
    scale_colour_manual(values = setNames(colors, paste0("", seq_along(colors)))) +
    facet_wrap(~Chr, nrow = 1, scales = "free_x", strip.position = "bottom", labeller = label_parsed) +
    geom_hline(yintercept = threshold_line, lty = 2, linewidth = 0.75, color = "red", alpha = 0.75) +
    theme_classic(base_size = 16) +
    theme(axis.text = element_text(face = "bold"),
          strip.text = element_text(size = 15, face = "bold"),
          axis.title = element_text(face = "bold"),
          axis.text.x = element_blank(), axis.ticks.x = element_blank(), strip.placement = "outside",
          strip.background = element_blank(),
          panel.border = element_blank(), panel.background = element_blank(), panel.grid = element_blank(),
          panel.spacing.x = unit(0, "line"),
          legend.position = "none") +
    labs(title = paste0(file_path), x = "Position (Mbp)", y = "-log10(pvalue)") +
    scale_y_continuous(expand = c(0, 0), limits = c(0, max(results$log_10) + 1))
  
  
  manplot <- manplot +
    geom_point(data = hits, aes(x = Basepair / 10^6, y = log_10), color = "black", size = 2.9)
  
  # Conditionally add geom_text_repel only if the logical statement is True
  if (gene_naming) {
    manplot <- manplot +
      geom_text_repel(data = hits, aes(x = Basepair / 10^6, y = log_10, label = Gene2), 
                      vjust = 0, hjust = 0, size = 3, nudge_y = 0.35, min.segment.length = 0.5,
                      fontface = "bold", color = "black", angle = 0)
  }
  #geom_text(data = results %>% filter(Predictor %in% hits$Marker), aes(x = Basepair / 10^6, y = log_10, label =), vjust = -0.5, hjust = 0.5, size = 3)
  
  # Save as PDF
  #print(paste("Saving PDF",file_path))
  # pdf_file <- paste0("ManhattanPlots/",file_path, ".pdf")
  # ggsave(plot = manplot, filename = pdf_file, units = "in", height = 5, width = 20, dpi = 320)
  
  # Save as JPEG
  print(paste("Saving JPEG",file_path))
  if (gene_naming) {
    jpeg_file <- paste0(GWAS.results.dir,"ManhattanPlots/",file_path, "_",region,"_FDR_NamedGenes.jpeg")
  } else {
    jpeg_file <- paste0(GWAS.results.dir,"ManhattanPlots/",file_path, "_",region,"_FDR.jpeg")
  }
  ggsave(plot = manplot, filename = jpeg_file, units = "in", height = 5.5, width = 20, dpi = 320, quality = 10)
}
gc()
#create_manhattan_plot(file_path = file_names)
home.dir <- "~/Desktop/MaizeGWASResults/"
setwd(home.dir)
GWAS.results.dir <- "GWASresults/GLM/"
if (!dir.exists(paste0(GWAS.results.dir,"ManhattanPlots"))) { dir.create(paste0(GWAS.results.dir,"ManhattanPlots")) }

region <- "Genic"

load(paste0(GWAS.results.dir,"CandidateGenes_GLM_GWASoutput_overlap_",region,".Rdata"))
colnames(overlap_df)[2] <- "Gene"
#overlap_df$Gene<-gsub("00001eb", "_", overlap_df$Gene)
# Load your GWAS results data 
# Load GWAS file paths
suffix <- ".assoc"  # we want to grab the association result files
file_names <- list.files(path =GWAS.results.dir,pattern = paste0("*", suffix))

colors <- brewer.pal(n = 10, name = "Set1")
colors[6] <- "#DAA520"
colors[10] <- "#4EB3D3"

interesting_genes <- fread("ERN_interesting_genes.csv",header=T)
interesting_genes <- fread("SSL_interesting_genes.csv",header=T)
interesting_genes <- interesting_genes[interesting_genes$Region == region, ]

overlap_df$Gene2 <- NA
for(i in 1:nrow(interesting_genes)){
  overlap_df$Gene2[which(overlap_df$Gene==interesting_genes$ZmID[i])] <- interesting_genes$Name[i]
}

threshold <- 0.05
threshold_type <- "FDR"
gene_naming <-F
# Run the function in parallel using mcapply
mclapply(file_names[6], create_manhattan_plot, mc.cores = detectCores())
# Print the results (file paths of saved plots)

