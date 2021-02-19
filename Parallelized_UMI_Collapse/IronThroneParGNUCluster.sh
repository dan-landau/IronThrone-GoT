#!/bin/bash

#Navigate to directory out of which parallelization script is being run
cd $(dirname $0)


#Set Up Command Line Option Defaults
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
low_mem=0


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
		-lm | --low_mem )		shift
					low_mem=$1
					;;
	esac
	shift
done


#Convert desired line number of split fastq into number of reads
target_reads=$((target_lines / 4))


#Ensure config and barcode whitelist file are specified or present in the working directory
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


#Set precise paths for config and whitelist files
config=$(readlink -f $config)
whitelist=$(readlink -f $whitelist)

#If skip_shuf option is not passed, will do the following section to shuffle and split the reads
if ((skip_shuf != 1))
then
	#Check for R1 and R2 fastq files
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
	if ($low_mem == 1)
	then
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
	else
		if (($(grep ";" combined.fastq | wc -l) == 0))
		then
			awk '{printf("%s%s",$0,(NR%4==0)?"\n":";")}' combined.fastq | shuf | tr ";" "\n" > combined_shuffled.fastq
			echo fastq files shuffled
		elif (($(grep "|" combined.fastq | wc -l) == 0))
		then
			awk '{printf("%s%s",$0,(NR%4==0)?"\n":"|")}' combined.fastq | shuf | tr "|" "\n" > combined_shuffled.fastq
			echo fastq files shuffled
		else
			echo "New awk-line character needed"
			exit 1
		fi
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

	#Move split fastq files into shuffled_split directory
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


#Create text file of commands for GNU Parallel to execute
touch ../Parallel_Command_List.txt
>../Parallel_Command_List.txt

#Loop through split R1 and R2 files, creating directories for each split's individual IronThrone run and adding a command to the parallel command list with the corresponding R1 and R2 filenames
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


#Run list of IronThrone commands on split fastqs using GNU Parallel
if ((skip_iron_throne != 1))
then
	parallel :::: Parallel_Command_List.txt

	echo All instances of IronThrone complete
fi


#Call R script to concatenate and collapse output
Rscript Combine_IronThrone_Parallel_Output.R $main_output_folder ${pcr_read_threshold} ${levenshtein_distance} ${dupcut}

echo All IronThrone outputs concatenated into myGoT.summTable.concat.txt

outdir=$(readlink -f $outdir)

if [ ! -f $outdir'/myGoT.summTable.concat.txt' ]
then
	mv myGoT.summTable.concat.txt $outdir
	mv myGoT.summTable.concat.umi_collapsed.txt $outdir
fi


rm Parallel_Command_List.txt
