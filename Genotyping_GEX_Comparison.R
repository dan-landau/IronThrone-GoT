#Read in options from Rscript command to set filtered barcode text file location, molecule info h5 file, IronThrone output table, Gene Name, output file for metadata, threshold
options <- commandArgs(trailingOnly = TRUE)
bc_loc <- options[1]
h5_file <- options[2]
got_df_loc <- options[3]
target_gene <- options[4]
output_dir <- options[5]
quant_thresh <- options[6]


#UMI to Binary and Binary to UMI conversion functions
umi_bin_to_seq <- function(umi_decimal, umi_len = 12){
  umi_bin <- R.utils::intToBin(umi_decimal)
  umi <- NA
  if (nchar(umi_bin) < (umi_len*2)){
    difference <- (umi_len*2)-nchar(umi_bin)
    umi_bin <- paste0(paste0(rep(0, difference), collapse = ""), umi_bin)
  }
  if (nchar(umi_bin) == (umi_len*2)){
    umi_bits <- substring(umi_bin, seq(1, nchar(umi_bin), 2), seq(2, nchar(umi_bin), 2))
    umi_char <- plyr::mapvalues(umi_bits, from = c("00", "01", "10", "11"), to = c("A", "C", "G", "T"))
    umi <- paste(umi_char, sep = "", collapse = "")
  }
  return(umi)
}


umi_seq_to_bin_decimal <- function(umi_seq){
  umi_char <- unlist(strsplit(umi_seq, ""))
  umi_bits <- plyr::mapvalues(umi_char, to = c("00", "01", "10", "11"), from = c("A", "C", "G", "T"))
  umi_bin <- paste0(umi_bits, collapse = "")
  umi_decimal <- strtoi(umi_bin, base = 2)
  return(umi_decimal)
}


library(parallel)
library(tidyverse)
library(rhdf5)
library(stringdist)

#Load barcodes that will be used to generate seurat object
seurat_bcs <- scan(file = bc_loc, what = "character")
seurat_bcs <- gsub("-.*", "", seurat_bcs)

#Retrieve GEX molecule information
gex_molecules <- list()
gex_molecules[["barcode_idx"]] <- h5read(h5_file, name = "barcode_idx") + 1
gex_molecules[["barcodes"]] <- h5read(h5_file, name = "barcodes")
gex_molecules[["umi_counts"]] <- h5read(h5_file, name = "count")
gex_molecules[["feature_idx"]] <- h5read(h5_file, name = "feature_idx") + 1
gex_molecules[["feature_name"]] <- h5read(h5_file, name = "features/name")
gex_molecules[["feature_id"]] <- h5read(h5_file, name = "features/id")
gex_molecules[["umi"]] <- h5read(h5_file, name = "umi")

#Create the genotyping metadata data frame
md <- data.frame(BC = seurat_bcs)

#Read in IronThrone genotyping library results
got_df <- read.delim(got_df_loc, stringsAsFactors = FALSE)
got_df <- got_df[got_df$UMI != "",]

#Unfiltered Genotyping using complete IronThrone output
temp_metadata <- merge(x = got_df, by.x = c("BC"), y = md, by.y = c("BC"), all.y = TRUE)

md$unfilt.WT.calls <- as.numeric(temp_metadata$WT.calls)
md$unfilt.MUT.calls <- as.numeric(temp_metadata$MUT.calls)
md$unfilt.Total.calls <- as.numeric(md$unfilt.WT.calls) + as.numeric(md$unfilt.MUT.calls)

#Any MUT calls is MUT, 0 MUT and any WT calls is WT
md$unfilt.Genotype <- ifelse(is.na(md$unfilt.WT.calls), "No Data",
                                  ifelse(md$unfilt.MUT.calls>0, "MUT",
                                         ifelse(md$unfilt.WT.calls>=1, "WT", "NA")))


#Split got df into a per-UMI data frame
split_got_df <- data.frame(matrix(nrow = length(unlist(strsplit(got_df[,"UMI"],";"))), ncol = 0))
for (i in colnames(got_df)){
  per_umi <- length(grep(";", got_df[,i])) > 0
  if (per_umi){
    split_got_df[,i] <- unlist(strsplit(got_df[,i],";"))
  } else {
    split_got_df[,i] <- rep(got_df[,i], times = got_df$WT.calls + got_df$MUT.calls + got_df$amb.calls)
  }
}
split_got_df$num.WT.in.dups <- as.numeric(split_got_df$num.WT.in.dups)
split_got_df$num.MUT.in.dups <- as.numeric(split_got_df$num.MUT.in.dups)
split_got_df$num.amb.in.dups <- as.numeric(split_got_df$num.amb.in.dups)
split_got_df$WT.calls <- as.numeric(split_got_df$WT.calls)
split_got_df$MUT.calls <- as.numeric(split_got_df$MUT.calls)
split_got_df$amb.calls <- as.numeric(split_got_df$amb.calls)

split_got_df$BC_UMI <- paste0(split_got_df$BC, "_", split_got_df$UMI)
split_got_df$total_dups <- split_got_df$num.WT.in.dups + split_got_df$num.MUT.in.dups + split_got_df$num.amb.in.dups
split_got_df$total_dups_wt_mut <- split_got_df$num.WT.in.dups + split_got_df$num.MUT.in.dups
split_got_df$UMI_bin <- unlist(mclapply(split_got_df$UMI, umi_seq_to_bin_decimal, mc.cores = detectCores()))
split_got_df$BC_UMI_bin <- paste0(split_got_df$BC, "_", split_got_df$UMI_bin)

#Add info about BC/UMI pairs as found in scRNA molecule info for target gene
target_gene_idx <- which(gex_molecules$feature_name == target_gene)
target_gene_id <- gex_molecules$feature_id[target_gene_idx]
target_gene_entries <- which(gex_molecules$feature_idx == target_gene_idx)
df_10x <- data.frame("BC" = gex_molecules$barcodes[gex_molecules$barcode_idx[target_gene_entries]])

umi_filt <- gex_molecules$umi[target_gene_entries]
df_10x$UMI <- unlist(mclapply(umi_filt, mc.cores = detectCores(), FUN = function(x){
  umi_bin_to_seq(x, umi_len = 12)
}))
df_10x$counts <- gex_molecules$umi_counts[target_gene_entries]
df_10x$BC_UMI <- paste0(df_10x$BC, "_", df_10x$UMI)

#Is there an exact match to a DNMT3A molecule in GEX?
split_got_df$Exact_Match <- split_got_df$BC_UMI %in% df_10x$BC_UMI

#Is there an approximate match to a DNMT3A molecule in GEX?
split_got_df$Approx_Match <- unlist(mclapply(split_got_df$BC_UMI, mc.cores = detectCores(), FUN = function(x){
  ain(x, df_10x$BC_UMI, method = "lv", maxDist = 2)
}))


#Molecule info for all genes in 10X molecule info
target_bc_idx <- which(gex_molecules$barcodes %in% got_df$BC)
molecules <- gex_molecules$barcode_idx %in% target_bc_idx
df_all_gene <- data.frame("BC_IDX" = gex_molecules$barcode_idx[molecules])
df_all_gene$BC <- gex_molecules$barcodes[df_all_gene$BC_IDX]
df_all_gene$UMI_bin <- gex_molecules$umi[molecules]
df_all_gene$BC_UMI_bin <- paste0(df_all_gene$BC, "_", df_all_gene$UMI_bin)
df_all_gene$gene_idx <- gex_molecules$feature_idx[molecules]
df_all_gene$gene <- gex_molecules$feature_name[df_all_gene$gene_idx]
df_all_gene$count <- gex_molecules$umi_counts[molecules]
to_keep <- df_all_gene$BC_UMI_bin %in% (split_got_df %>% pull(BC_UMI_bin))

df_all_gene_to_keep <- df_all_gene[to_keep,]
df_all_gene_to_keep$UMI <- unlist(mclapply(df_all_gene_to_keep$UMI_bin, mc.cores = detectCores(), FUN = function(x){
  umi_bin_to_seq(x, umi_len = 12)
}))
df_all_gene_to_keep$BC_UMI <- paste0(df_all_gene_to_keep$BC, "_", df_all_gene_to_keep$UMI)

df_all_gene_collapse <- df_all_gene_to_keep
sort(table(df_all_gene_to_keep$BC_UMI), decreasing = TRUE)

for (k in unique(df_all_gene_collapse$BC_UMI)){
  target_rows <- which(df_all_gene_collapse$BC_UMI == k)
  if (length(target_rows) > 1){
    sub_df <- df_all_gene_collapse[target_rows,]
    if (target_gene %in% sub_df$gene){
      if(length(grep("_CITE", sub_df$gene)) > 0){
        sub_df$gene <- paste0("Multiple_", target_gene, "_CITE")
      } else{
        sub_df$gene <- paste0("Multiple_", target_gene)
      }
    } else {
      sub_df$gene <- "Multiple"
    }
    min_row <- min(target_rows)
    target_rows <- target_rows[target_rows != min_row]
    df_all_gene_collapse[min_row,] <- sub_df[1,]
    df_all_gene_collapse <- df_all_gene_collapse[-target_rows,]
  }
}

split_got_df_gene <- (merge(split_got_df, df_all_gene_collapse[,c("BC_UMI", "gene")], by = "BC_UMI", all.x = TRUE, all.y = FALSE, sort = FALSE))

split_got_df_gene$In_GEX <- !is.na(split_got_df_gene$gene)

split_got_df_gene$Gene_Group <- ifelse(split_got_df_gene$Exact_Match,
                                       "Exact",
                                       ifelse(split_got_df_gene$Approx_Match,
                                              "Approx",
                                              ifelse(split_got_df_gene$In_GEX,
                                                     "Other Gene",
                                                     "No Gene")))

#Approx Matches + UMIs that don't match any non DNMT3A gene in the 10X library
concatenate_got <- function(BC, split_df){
  single_bc_mat <- split_df[split_df[,"BC"] == BC,]
  single_bc_vec <- apply(single_bc_mat, MARGIN = 2, FUN = function(x) paste0(x, collapse = ";"))
  single_bc_vec["BC"] <- BC
  single_bc_vec["WT.calls"] <- sum(single_bc_mat[,"call.in.dups"] == "WT")
  single_bc_vec["MUT.calls"] <- sum(single_bc_mat[,"call.in.dups"] == "MUT")
  single_bc_vec["amb.calls"] <- sum(single_bc_mat[,"call.in.dups"] == "AMB")
  single_bc_df <- t(as.data.frame(single_bc_vec, stringsAsFactors = FALSE))
  rownames(single_bc_df) <- NULL
  return(single_bc_df)
}

unique_bc_approx_no_gene <- unique(split_got_df_gene %>% filter(Gene_Group != "Other Gene") %>% pull(BC))
split_got_df_approx_no_gene <- split_got_df_gene %>% filter(Gene_Group != "Other Gene")
concat_got_df_approx_no_gene <- as.data.frame(Reduce(rbind, mclapply(unique_bc_approx_no_gene, FUN = function(x) (concatenate_got(BC = x, split_df = split_got_df_approx_no_gene)), mc.cores = detectCores())), stringsAsFactors = FALSE)
concat_got_df_approx_no_gene$Genotype <- ifelse(is.na(concat_got_df_approx_no_gene$WT.calls), "No Data",
                                                ifelse(concat_got_df_approx_no_gene$MUT.calls>0, "MUT",
                                                       ifelse(concat_got_df_approx_no_gene$WT.calls>=1, "WT", "NA")))


temp_metadata <- merge(md, concat_got_df_approx_no_gene[, c("BC", "Genotype", "WT.calls", "MUT.calls")], by = "BC", all.x = TRUE, all.y = FALSE)
rownames(temp_metadata) <- temp_metadata$BC
md$filt.Genotype <- temp_metadata$Genotype
md$WT.calls.filt <- as.numeric(temp_metadata$WT.calls)
md$MUT.calls.filt <- as.numeric(temp_metadata$MUT.calls)
md$Total.calls.filt <- md$WT.calls.filt + md$MUT.calls.filt

#Approx Matches + UMIs that don't match any non DNMT3A gene in the 10X library, threshold applied
other_gene_counts <- split_got_df_gene %>% filter(Gene_Group == "Other Gene") %>% pull(total_dups_wt_mut)
threshold <- quantile(other_gene_counts, probs = quant_thresh)

thresh_plot <- ggplot(split_got_df_gene, aes(y = log10(total_dups), x = Gene_Group, fill = Gene_Group))+ geom_violin(position = position_dodge(0.9), trim = FALSE) + 
  geom_boxplot(width=0.1, position = position_dodge(0.9), alpha = 0.5) + 
  geom_hline(yintercept = log10(threshold)) + 
  theme_bw() + 
  labs(y = "log10(Supporting Read Counts per UMI)", x = "Amplicon Match to GEX Library")

if(file.exists(paste0(output_dir, "/threshold_plot.pdf"))){
  print("Warning: File threshold_plot.pdf already exists in output directory")
} else{
  ggsave(filename = paste0(output_dir, "/threshold_plot.pdf"), plot = thresh_plot, device = "pdf")
}


split_got_df_gene$Keep <- ifelse(split_got_df_gene$Gene_Group %in% c("Exact", "Approx"), TRUE,
                                 ifelse(split_got_df_gene$Gene_Group == "Other Gene", FALSE,
                                        ifelse(split_got_df_gene$total_dups_wt_mut > threshold, TRUE, FALSE)))
split_got_df_gene_thresh <- split_got_df_gene

unique_bc_approx_no_gene_thresh <- unique(split_got_df_gene_thresh %>% filter(Keep) %>% pull(BC))
split_got_df_approx_no_gene_thresh <- split_got_df_gene_thresh %>% filter(Keep)
concat_got_df_approx_no_gene_thresh <- as.data.frame(Reduce(rbind, mclapply(unique_bc_approx_no_gene_thresh, FUN = function(x) (concatenate_got(BC = x, split_df = split_got_df_approx_no_gene_thresh)), mc.cores = 16)), stringsAsFactors = FALSE)
concat_got_df_approx_no_gene_thresh$Genotype <- ifelse(is.na(concat_got_df_approx_no_gene_thresh$WT.calls), "No Data",
                                                       ifelse(concat_got_df_approx_no_gene_thresh$MUT.calls>0, "MUT",
                                                              ifelse(concat_got_df_approx_no_gene_thresh$WT.calls>=1, "WT", "NA")))


temp_metadata <- merge(md, concat_got_df_approx_no_gene_thresh[, c("BC", "Genotype", "WT.calls", "MUT.calls")], by = "BC", all.x = TRUE, all.y = FALSE)
rownames(temp_metadata) <- md$BC
md$thresh.filt.Genotype <- temp_metadata$Genotype
md$WT.calls.thresh.filt <- as.numeric(temp_metadata$WT.calls)
md$MUT.calls.thresh.filt <- as.numeric(temp_metadata$MUT.calls)
md$Total.calls.thresh.filt <- md$WT.calls.thresh.filt + md$MUT.calls.thresh.filt

if(file.exists(paste0(output_dir, "/metadata.Rdata"))){
  print("Warning: File metadata.Rdata already exists in output directory")
} else{
  save(md, file = paste0(output_dir, "/metadata.Rdata"))
}


