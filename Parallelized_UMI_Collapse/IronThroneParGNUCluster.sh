#!/bin/bash

cd $(dirname $0)

fastqR1=$(pwd)'/'$(ls *R1*)
fastqR2=$(pwd)'/'$(ls *R2*)
target_lines=500000
config=$(pwd)'/'$(ls *.config)
whitelist=$(pwd)'/'$(ls *[b,B]arcode*)
umilen=12
bclen=16
run=linear
mmtch=0.2
postP=0.99
dupcut=2
sample=myGoT
outdir=./
log=myGoT.log
keepouts=0
verbose=0
skip_shuf=0
pcr_read_threshold=0.5
skip_iron_throne=0
levenshtein_distance=0.1


#Set Up Command Line Options
while [ "$1" != "" ]; do
	case $1 in
		-f1 | --fastqR1 )		shift
					fastqR1=$1
					;;
		-f2 | --fastqR2 )		shift
					fastqR2=$1
					;;
		-tl | --target_lines )	shift
					target_lines=$1
					;;
		-c | --config )		shift
					config=$1
					;;
		-w | --whitelist )	shift
					whitelist=$1
					;;
		-u | --umilen )		shift
					umilen=$1
					;;
		-b | --bclen )		shift
					bclen=$1
					;;
		-r | --run )		shift
					run=$1
					;;
		-m | --mmtch )		shift
					mmtch=$1
					;;
		-p | --postP )		shift
					postP=$1
					;;
		-d | --dupcut )		shift
					dupcut=$1
					;;
		-s | --sample )		shift
					sample=$1
					;;
		-o | --outdir )		shift
					outdir=$1
					;;
		-l | --log )		shift
					log=$1
					;;
		-k | --keepouts )	shift
					keepouts=$1
					;;
		-v | --verbose )	shift
					verbose=$1
					;;
		-z | --skip_shuf )	shift
					skip_shuf=$1
					;;
		-pcr | --pcr_read_threshold )	shift
					pcr_read_threshold=$1
					;;
		-x | --skip_iron_throne )	shift
					skip_iron_throne=$1
					;;
		-ld | --levenshtein_distance )	shift
					levenshtein_distance=$1
					;;
	esac
	shift
done

target_reads=$((target_lines / 4))

if (( $(echo "$config" | wc -w) != 1 ))
then
        echo Single config file needs to be specified or spaces need to be removed from filename
        exit 1
fi

if (( $(echo "$whitelist" | wc -w) != 1 ))
then
        echo Single whitelist file needs to be specified or spaces need to be removed from filename. If filename contains the string "barcode" and is located in the same directory as this script, it will be called automatically.
        exit 1
fi

config=$(readlink -f $config)
whitelist=$(readlink -f $whitelist)

#If skip_shuf option is not passed, will do the following section to shuffle and split the reads
if ((skip_shuf != 1))
then
	if (( $(echo "$fastqR1" | wc -w) != 1 ))
	then
		echo Single R1 file needs to be specified or spaces need to be removed from filename
		exit 1
	fi

	if (( $(echo "$fastqR2" | wc -w) != 1 ))
	then
	        echo Single R2 file needs to be specified or spaces need to be removed from filename
	        exit 1
	fi
	#Join R1 and R2 into a single file with each line containing tab-separated corresponding lines of R1 and R2
	paste $fastqR1 $fastqR2 > combined.fastq
	echo fastq files joined


	#Randomly sort lines of combined R1/R2 file
	if (($(grep ";" combined.fastq | wc -l) == 0))
	then
		awk '{printf("%s%s",$0,(NR%4==0)?"\n":";")}' combined.fastq | sort -R | tr ";" "\n" > combined_shuffled.fastq
		echo fastq files shuffled
	elif (($(grep "|" combined.fastq | wc -l) == 0))
	then
		awk '{printf("%s%s",$0,(NR%4==0)?"\n":"|")}' combined.fastq | sort -R | tr "|" "\n" > combined_shuffled.fastq
		echo fastq files shuffled
	else
		echo "New awk-line character needed"
		exit 1
	fi

	#Separate shuffled file back into R1 and R2
	cut -f1 -d$'\t' combined_shuffled.fastq > shuffled.R1.fastq
	cut -f2 -d$'\t' combined_shuffled.fastq > shuffled.R2.fastq
	echo shuffled fastq files cut back into R1 and R2

	#Split files into pieces of $target_line lines
	total_lines=$(wc -l < combined_shuffled.fastq)
	total_reads=$((total_lines / 4))
	reads_mod=$((total_reads % target_reads))
	#Adjust target_lines value to ensure number of reads in final file is > 0.9 of all other files
	if ((reads_mod*10 < target_reads*9 && reads_mod != 0))
	then
		remainder_reads=$((total_reads % target_reads))
		total_files=$((total_reads / target_reads))
		to_add=$(( (remainder_reads / total_files) + 1))
		target_reads=$((target_reads + to_add))
		target_lines=$((target_reads * 4))
	fi


	split -d -a 3 -l $target_lines shuffled.R1.fastq shuffled.R1
	split -d -a 3 -l $target_lines shuffled.R2.fastq shuffled.R2

	total_files=$(ls shuffled.R1[0-9]* | wc -l)


	for file in $(ls | grep '.*R[0-9][0-9][0-9][0-9]'); do mv "$file" "$file.fastq"; done
	mkdir shuffled_split
	for file in $(ls | grep '.*R[0-9][0-9][0-9][0-9]'); do mv "$file" "./shuffled_split/"; done
	echo R1 and R2 split into $total_files pieces

	mkdir preprocessing_fastqs
	mv combined.fastq preprocessing_fastqs/
	mv combined_shuffled.fastq preprocessing_fastqs/
	mv shuffled.R1.fastq preprocessing_fastqs/
	mv shuffled.R2.fastq preprocessing_fastqs/
fi




#Create main output folder
mkdir -p Output
cd Output/
main_output_folder=$(pwd)

cd ..
cd shuffled_split/

if ((skip_iron_throne != 1))
then
	echo Begin job parallelization
fi

touch ../Parallel_Command_List.txt
>../Parallel_Command_List.txt

#Loop through split R1 and R2 files
total_files=0
for i in $(ls | grep '.*R[0-9][0-9][0-9][0-9]' | sed 's/\.fastq//g' | sed 's/.*R[0-9]//g' | sort | uniq);
do
R1=$(pwd)'/'$(ls | grep "R1${i}");
R2=$(pwd)'/'$(ls | grep "R2${i}");
output=${main_output_folder}'/'${i}
mkdir -p ${output};

echo "module load got/0.1; \
IronThrone-GoT \
-r ${run} \
-f1 ${R1} \
-f2 ${R2} \
-c ${config} \
-w ${whitelist} \
-u ${umilen} \
-b ${bclen} \
-o ${output} \
-m ${mmtch} \
-p ${postP} \
-d 1 \
-s ${sample} \
-l ${log} \
-k ${keepouts} \
-v ${verbose}" >> ../Parallel_Command_List.txt


done

#Back to main level folder
cd ..

if ((skip_iron_throne != 1))
then
	parallel :::: Parallel_Command_List.txt

	echo All instances of IronThrone complete
fi

cat <<EOF > Combine_IronThrone_Parallel_Output.R

setwd("$main_output_folder")

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

pcr_ratio_thresh <- ${pcr_read_threshold}

list_collapse <- function(single_got_row){
  UMIs <- unlist(strsplit(single_got_row[,"UMI"], ";"))
  num.WT.in.dups <- as.numeric(unlist(strsplit(single_got_row[,"num.WT.in.dups"], ";")))
  num.MUT.in.dups <- as.numeric(unlist(strsplit(single_got_row[,"num.MUT.in.dups"], ";")))
  num.amb.in.dups <- as.numeric(unlist(strsplit(single_got_row[,"num.amb.in.dups"], ";")))
  call.in.dups <- unlist(strsplit(single_got_row[,"call.in.dups"], ";"))


  match_list <- lapply(UMIs, FUN =  function(x) agrep(x, UMIs, max.distance = ${levenshtein_distance}))
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

    match_list <- lapply(UMIs, FUN =  function(x) agrep(x, UMIs, max.distance = ${levenshtein_distance}))
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

  dup_thresh <- ${dupcut}

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


EOF

Rscript Combine_IronThrone_Parallel_Output.R

echo All IronThrone outputs concatenated into myGoT.summTable.concat.txt

outdir=$(readlink -f $outdir)

if [ ! -f $outdir'/myGoT.summTable.concat.txt' ]
then
	mv myGoT.summTable.concat.txt $outdir
	mv myGoT.summTable.concat.umi_collapsed.txt $outdir
fi

rm Combine_IronThrone_Parallel_Output.R
rm Parallel_Command_List.txt
