#Read in options from Rscript command to set working directory, pcr ratio threshold, levenshtein distance for UMI collapse, and pcr duplicate threshold ####
options <- commandArgs(trailingOnly = TRUE)
wd <- options[1]
pcr_thresh <- as.numeric(options[2])
ld <- as.numeric(options[3])
dupcut <- as.numeric(options[4])
threads <- as.numeric(options[5])

setwd(wd)
library(parallel)

#Concatenate output files for parallel IronThrone runs #####
split_got <- data.frame()
for (i in list.files()){
    split_got <- rbind(split_got, read.delim(paste0(i,"/myGoT.summTable.txt"), stringsAsFactors = FALSE))
}

#Remove final MUT/WT/Amb UMI tallies, this will be recalculated after collapsing UMIs
split_got2 <- as.data.frame(split_got[,1:(ncol(split_got)-3)], stringsAsFactors = FALSE)

#Melt the data frame so each BC/UMI pair is a single row entry
split_got_df <- data.frame(matrix(nrow = length(unlist(strsplit(split_got2[,"UMI"],";"))), ncol = 0))
for (i in colnames(split_got2)){
  split_got_df[,i] <- unlist(strsplit(split_got2[,i],";"))
}

split_got_df_num <- split_got_df[,c("BC", "UMI", "num.WT.in.dups", "num.MUT.in.dups", "num.amb.in.dups")]
split_got_df_num[,"num.WT.in.dups"] <- as.numeric(split_got_df_num[,"num.WT.in.dups"])
split_got_df_num[,"num.MUT.in.dups"] <- as.numeric(split_got_df_num[,"num.MUT.in.dups"])
split_got_df_num[,"num.amb.in.dups"] <- as.numeric(split_got_df_num[,"num.amb.in.dups"])

#Function that takes the melted data frame from above and collapses down all UMIs associated with a single barcode (including duplicates) across all parallel runs of IronTHrone
concatenate_got <- function(BC, split_df){
  single_bc_mat <- split_df[split_df[,"BC"] == BC,]
  single_bc_mat_collapse <- data.frame(matrix(nrow = length(unique(single_bc_mat[,"UMI"]))))
  single_bc_mat_collapse[,1] <- NULL
  single_bc_mat_collapse[,"BC"] <- BC
  single_bc_mat[,"BC"] <- NULL
  single_bc_mat_collapse <- cbind(single_bc_mat_collapse, aggregate(single_bc_mat[,-1], by = list(single_bc_mat[,"UMI"]), sum))
  colnames(single_bc_mat_collapse) <- c("BC", "UMI", "num.WT.in.dups", "num.MUT.in.dups", "num.amb.in.dups")
  
  wt_frac <- single_bc_mat_collapse[,"num.WT.in.dups"]/(single_bc_mat_collapse[,"num.WT.in.dups"] + single_bc_mat_collapse[,"num.MUT.in.dups"])
  mut_frac <- single_bc_mat_collapse[,"num.MUT.in.dups"]/(single_bc_mat_collapse[,"num.WT.in.dups"] + single_bc_mat_collapse[,"num.MUT.in.dups"])
  single_bc_mat_collapse[,"call.in.dups"] <- ifelse(single_bc_mat_collapse[,"num.WT.in.dups"] + single_bc_mat_collapse[,"num.MUT.in.dups"] == 0, "AMB", 
                                                    ifelse(wt_frac > pcr_thresh, "WT", 
                                                           ifelse(mut_frac > pcr_thresh, "MUT", "AMB")))
  
  single_bc_mat_collapse[,"num.WT.in.dups"] <- as.character(single_bc_mat_collapse[,"num.WT.in.dups"])
  single_bc_mat_collapse[,"num.MUT.in.dups"] <- as.character(single_bc_mat_collapse[,"num.MUT.in.dups"])
  single_bc_mat_collapse[,"num.amb.in.dups"] <- as.character(single_bc_mat_collapse[,"num.amb.in.dups"])
  single_bc_vec <- apply(single_bc_mat_collapse, MARGIN = 2, FUN = function(x) paste0(x, collapse = ";"))
  single_bc_vec["BC"] <- BC
  single_bc_df <- t(as.data.frame(single_bc_vec, stringsAsFactors = FALSE))
  rownames(single_bc_df) <- NULL
  return(single_bc_df)
}

#Identify all unique barcodes in the data frame and run the concatenating function 
unique_bc <- unique(split_got_df_num[,"BC"])
concat_got_df <- as.data.frame(Reduce(rbind, mclapply(unique_bc, FUN = function(x) (concatenate_got(BC = x, split_df = split_got_df_num)), mc.cores = threads)), stringsAsFactors = FALSE)
write.table(concat_got_df, file = "../myGoT.summTable.concat.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)


#Collapse UMIs #####
raw_GoT_table <- concat_got_df
pcr_ratio_thresh <- pcr_thresh


#Function to collapse UMIs for a single cell barcode
list_collapse <- function(single_got_row){
  UMIs <- unlist(strsplit(single_got_row[,"UMI"], ";")) #Separate all UMIs
  num.WT.in.dups <- as.numeric(unlist(strsplit(single_got_row[,"num.WT.in.dups"], ";"))) #Separate number of WT PCR duplicates for each UMI
  num.MUT.in.dups <- as.numeric(unlist(strsplit(single_got_row[,"num.MUT.in.dups"], ";"))) #Separate number of MUT PCR duplicates for each UMI
  num.amb.in.dups <- as.numeric(unlist(strsplit(single_got_row[,"num.amb.in.dups"], ";"))) #Separate number of ambiguous PCR duplicates for each UMI
  call.in.dups <- unlist(strsplit(single_got_row[,"call.in.dups"], ";")) #Separate WT/MUT/Amb call for each UMI
  
  
  match_list <- lapply(UMIs, FUN =  function(x) agrep(x, UMIs, max.distance = ld)) #For each UMI, use the agrep function to identify all similar UMIs within the Levenshtein distance threshold
  num_of_matches <- unlist(lapply(match_list, FUN = function(x) length(x))) #Count the number of similar UMIs for each UMI
  
  #Continue to collapse while there are similar UMIs remaining in the set
  while (sum(num_of_matches) > length(num_of_matches)){
    to_collapse <- which(num_of_matches == max(num_of_matches))[1] #Identify the UMI with the most number of similar UMIs in the set
    matches_t0 <- numeric()
    matches_t1 <- match_list[[to_collapse]] #Identify the group of UMIs to be collapsed
    while (length(matches_t1) > length(matches_t0)){
      to_add <- unique(unlist(match_list[c(matches_t1)]))
      matches_t0 <- matches_t1
      matches_t1 <- to_add
    }
    
    WT_dups <- num.WT.in.dups[matches_t1] #Identify the number of WT PCR duplicates for each UMI to be collapsed
    MUT_dups <- num.MUT.in.dups[matches_t1] #Identify the number of MUT PCR duplicates for each UMI to be collapsed
    AMB_dups <- num.amb.in.dups[matches_t1] #Identify the number of Amb PCR duplicates for each UMI to be collapsed
    num.WT.in.dups[to_collapse] <- sum(WT_dups) #Add the WT PCR duplicates for all UMIs to be collapsed
    num.MUT.in.dups[to_collapse] <- sum(MUT_dups) #Add the MUT PCR duplicates for all UMIs to be collapsed
    num.amb.in.dups[to_collapse] <- sum(AMB_dups) #Add the Amb PCR duplicates for all UMIs to be collapsed
    if (sum(WT_dups) + sum(MUT_dups) == 0){
      call.in.dups[to_collapse] <- "AMB" #If no WT or MUT supporting PCR duplicates, the UMI is ambiguous
    } else {
      pcr_ratio <- (max(c(sum(WT_dups), sum(MUT_dups)))[1])/(sum(WT_dups) + sum(MUT_dups)) #Calculate the ratio of supporting PCR duplicates in favor of the majority genotype (WT or MUT)
      if(pcr_ratio > pcr_ratio_thresh){
        call.in.dups[to_collapse] <- ifelse(sum(WT_dups) > sum(MUT_dups), "WT", "MUT") #If the ratio passes the specified threshold, call the UMI in favor of the majority genotype
      } else {
        call.in.dups[to_collapse] <- "AMB" #If the ratio is below the specified threshold, the UMI is ambiguous
      }
    }
    
    #Remove the UMIs and corresponding data for all but the UMI to which the group was collapsed
    matches_rm <- matches_t1[matches_t1 != to_collapse] 
    UMIs <- UMIs[-matches_rm]
    num.WT.in.dups <- num.WT.in.dups[-matches_rm]
    num.MUT.in.dups <- num.MUT.in.dups[-matches_rm]
    num.amb.in.dups <- num.amb.in.dups[-matches_rm]
    call.in.dups <- call.in.dups[-matches_rm]
    
    #Recalculate the number of matches within the UMI set for the single cell barcode
    match_list <- lapply(UMIs, FUN =  function(x) agrep(x, UMIs, max.distance = ld))
    num_of_matches <- unlist(lapply(match_list, FUN = function(x) length(x)))
  }
  #Reformat the data for remaining UMIs following collapsing back into the format of the original IronThone data frame
  sort_val <- order(num.WT.in.dups + num.MUT.in.dups + num.amb.in.dups, decreasing = TRUE)
  single_got_row[,"UMI"] <- paste0(UMIs[sort_val], collapse = ";")
  single_got_row[,"num.WT.in.dups"] <- paste0(num.WT.in.dups[sort_val], collapse = ";")
  single_got_row[,"num.MUT.in.dups"] <- paste0(num.MUT.in.dups[sort_val], collapse = ";")
  single_got_row[,"num.amb.in.dups"] <- paste0(num.amb.in.dups[sort_val], collapse = ";")
  single_got_row[,"call.in.dups"] <- paste0(call.in.dups[sort_val], collapse = ";")
  single_got_row[,"WT.calls"] <- sum(call.in.dups[sort_val] == "WT")
  single_got_row[,"MUT.calls"] <- sum(call.in.dups[sort_val] == "MUT")
  single_got_row[,"amb.calls"] <- sum(call.in.dups[sort_val] == "AMB")
  
  #Keep only those UMIs with a total number of supporting PCR duplicates above the pre-specified threshold
  dup_thresh <- dupcut
  
  UMIs <- unlist(strsplit(single_got_row[,"UMI"], ";"))
  num.WT.in.dups <- as.numeric(unlist(strsplit(single_got_row[,"num.WT.in.dups"], ";")))
  num.MUT.in.dups <- as.numeric(unlist(strsplit(single_got_row[,"num.MUT.in.dups"], ";")))
  num.amb.in.dups <- as.numeric(unlist(strsplit(single_got_row[,"num.amb.in.dups"], ";")))
  call.in.dups <- unlist(strsplit(single_got_row[,"call.in.dups"], ";"))
  sum_reads <- num.WT.in.dups + num.MUT.in.dups
  threshold_filter <- sum_reads >= dup_thresh
  
  single_got_row[,"UMI"] <- paste0(UMIs[threshold_filter], collapse = ";")
  single_got_row[,"num.WT.in.dups"] <- paste0(num.WT.in.dups[threshold_filter], collapse = ";")
  single_got_row[,"num.MUT.in.dups"] <- paste0(num.MUT.in.dups[threshold_filter] , collapse = ";")
  single_got_row[,"num.amb.in.dups"] <- paste0(num.amb.in.dups[threshold_filter], collapse = ";")
  single_got_row[,"call.in.dups"] <- paste0(call.in.dups[threshold_filter], collapse = ";")
  single_got_row[,"WT.calls"] <- sum(call.in.dups[threshold_filter] == "WT")
  single_got_row[,"MUT.calls"] <- sum(call.in.dups[threshold_filter] == "MUT")
  single_got_row[,"amb.calls"] <- sum(call.in.dups[threshold_filter] == "AMB")
  
  return(single_got_row)
}

#Run the single-barcode UMI collapsing function on the entire dataset
UMI_collapse <- function(GoT_table){
  GoT_table_to_collapse <- data.frame(BC = GoT_table[,"BC"],
                                      UMI = GoT_table[,"UMI"],
                                      num.WT.in.dups = GoT_table[,"num.WT.in.dups"],
                                      num.MUT.in.dups = GoT_table[,"num.MUT.in.dups"],
                                      num.amb.in.dups = GoT_table[,"num.amb.in.dups"],
                                      call.in.dups = GoT_table[,"call.in.dups"],
                                      WT.calls = "",
                                      MUT.calls = "",
                                      amb.calls = "",
                                      stringsAsFactors = FALSE)
  #Convert the GoT data frame into a list of single-row data frames where each entry is the data for a single cell barcode
  GoT_list <- split(GoT_table_to_collapse, seq(nrow(GoT_table_to_collapse)))
  collapsed_GoT_list <- mclapply(GoT_list, FUN = list_collapse, mc.cores = threads)
  #Convert collapsed list back to a data frame
  collapsed_GoT_table <- do.call("rbind", collapsed_GoT_list)
  return(collapsed_GoT_table)
}

max_collapsed_GoT_table <- UMI_collapse(raw_GoT_table)

#Save output
write.table(max_collapsed_GoT_table, file = "../myGoT.summTable.concat.umi_collapsed.txt", sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)