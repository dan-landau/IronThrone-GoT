options <- commandArgs(trailingOnly = TRUE)
wd <- options[1]
pcr_thresh <- as.numeric(options[2])
ld <- as.numeric(options[3])
dupcut <- as.numeric(options[4])

setwd(wd)

library(parallel)

#Concatenate
for (i in list.files()){
  if (as.numeric(i)==0){
    split_got <- read.delim(paste0(i,"/myGoT.summTable.txt"), stringsAsFactors = FALSE)
  } else {
    split_got <- rbind(split_got, read.delim(paste0(i,"/myGoT.summTable.txt"), stringsAsFactors = FALSE))
  }
}

split_got2 <- as.data.frame(split_got[,1:(ncol(split_got)-3)], stringsAsFactors = FALSE)
split_got_df <- data.frame(matrix(nrow = length(unlist(strsplit(split_got2[,"UMI"],";")))))
for (i in colnames(split_got2)){
  split_got_df[,i] <- unlist(strsplit(split_got2[,i],";"))
}
split_got_df <- split_got_df[,2:ncol(split_got_df)]


concatenate_got <- function(BC, split_df){
  single_bc_mat <- split_df[split_df[,"BC"] == BC,]
  single_bc_vec <- apply(single_bc_mat, MARGIN = 2, FUN = function(x) paste0(x, collapse = ";"))
  single_bc_vec["BC"] <- BC
  single_bc_df <- t(as.data.frame(single_bc_vec, stringsAsFactors = FALSE))
  rownames(single_bc_df) <- NULL
  return(single_bc_df)
}


unique_bc <- unique(split_got_df[,"BC"])


concat_got_df <- as.data.frame(Reduce(rbind, mclapply(unique_bc, FUN = function(x) (concatenate_got(BC = x, split_df = split_got_df)))), stringsAsFactors = FALSE)

write.table(concat_got_df, file = "../myGoT.summTable.concat.txt", sep = "\t", row.names = FALSE, col.names = TRUE)


#Collapse UMIs
raw_GoT_table <- concat_got_df

pcr_ratio_thresh <- pcr_thresh

list_collapse <- function(single_got_row){
  UMIs <- unlist(strsplit(single_got_row[,"UMI"], ";"))
  num.WT.in.dups <- as.numeric(unlist(strsplit(single_got_row[,"num.WT.in.dups"], ";")))
  num.MUT.in.dups <- as.numeric(unlist(strsplit(single_got_row[,"num.MUT.in.dups"], ";")))
  num.amb.in.dups <- as.numeric(unlist(strsplit(single_got_row[,"num.amb.in.dups"], ";")))
  call.in.dups <- unlist(strsplit(single_got_row[,"call.in.dups"], ";"))
  
  
  match_list <- lapply(UMIs, FUN =  function(x) agrep(x, UMIs, max.distance = ld))
  num_of_matches <- unlist(lapply(match_list, FUN = function(x) length(x)))
  
  while (sum(num_of_matches) > length(num_of_matches)){
    
    
    to_collapse <- which(num_of_matches == max(num_of_matches))[1]
    matches_t0 <- numeric()
    matches_t1 <- match_list[[to_collapse]]
    while (length(matches_t1) > length(matches_t0)){
      to_add <- unique(unlist(match_list[c(matches_t1)]))
      matches_t0 <- matches_t1
      matches_t1 <- to_add
    }
    
    WT_dups <- num.WT.in.dups[matches_t1]
    MUT_dups <- num.MUT.in.dups[matches_t1]
    AMB_dups <- num.amb.in.dups[matches_t1]
    num.WT.in.dups[to_collapse] <- sum(WT_dups)
    num.MUT.in.dups[to_collapse] <- sum(MUT_dups)
    num.amb.in.dups[to_collapse] <- sum(AMB_dups)
    if (sum(WT_dups) + sum(MUT_dups) == 0){
      call.in.dups[to_collapse] <- "AMB"
    } else {
      pcr_ratio <- (max(c(sum(WT_dups), sum(MUT_dups)))[1])/(sum(WT_dups) + sum(MUT_dups))
      if(pcr_ratio > pcr_ratio_thresh){
        call.in.dups[to_collapse] <- ifelse(sum(WT_dups) > sum(MUT_dups), "WT", "MUT")
      } else {
        call.in.dups[to_collapse] <- "AMB"
      }
    }
    
    matches_rm <- matches_t1[matches_t1 != to_collapse]
    UMIs <- UMIs[-matches_rm]
    num.WT.in.dups <- num.WT.in.dups[-matches_rm]
    num.MUT.in.dups <- num.MUT.in.dups[-matches_rm]
    num.amb.in.dups <- num.amb.in.dups[-matches_rm]
    call.in.dups <- call.in.dups[-matches_rm]
    
    match_list <- lapply(UMIs, FUN =  function(x) agrep(x, UMIs, max.distance = ld))
    num_of_matches <- unlist(lapply(match_list, FUN = function(x) length(x)))
  }
  single_got_row[,"UMI"] <- paste0(UMIs, collapse = ";")
  single_got_row[,"num.WT.in.dups"] <- paste0(num.WT.in.dups, collapse = ";")
  single_got_row[,"num.MUT.in.dups"] <- paste0(num.MUT.in.dups, collapse = ";")
  single_got_row[,"num.amb.in.dups"] <- paste0(num.amb.in.dups, collapse = ";")
  single_got_row[,"call.in.dups"] <- paste0(call.in.dups, collapse = ";")
  single_got_row[,"WT.calls"] <- sum(call.in.dups == "WT")
  single_got_row[,"MUT.calls"] <- sum(call.in.dups == "MUT")
  single_got_row[,"amb.calls"] <- sum(call.in.dups == "AMB")
  
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
  GoT_list <- split(GoT_table_to_collapse, seq(nrow(GoT_table_to_collapse)))
  #rename_vector <- colnames(GoT_table_to_collapse)
  collapsed_GoT_list <- mclapply(GoT_list, FUN = list_collapse)
  collapsed_GoT_table <- do.call("rbind", collapsed_GoT_list)
  return(collapsed_GoT_table)
}

max_collapsed_GoT_table_higher2 <- UMI_collapse(raw_GoT_table)

write.table(max_collapsed_GoT_table_higher2, file = "../myGoT.summTable.concat.umi_collapsed.txt", sep = "\t", row.names = FALSE, col.names = TRUE)