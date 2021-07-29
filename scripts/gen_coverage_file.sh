#!/usr/bin/env bash

##############################################################################################################################################################
#
# This script is meant to be run on the contigs and the sequencing reads to generate coverage information
# Ideally it should take in the assembly file of all of your samples, followed by the reads of all the samples that went into the assembly.

# Some of the programs this pipeline uses are from the binning.sh file in MetaWRAP.
##############################################################################################################################################################


help_message () {
	echo ""
	echo "Usage: bash gen_coverage_file.sh [options] -a assembly.fa -o output_dir readsA_1.fastq readsA_2.fastq ... [readsX_1.fastq readsX_2.fastq]"
	echo "Note1: Make sure to provide all your separately replicate read files, not the joined file."
	echo "Note2: You may provide single end or interleaved reads as well with the use of the correct option"
	echo "Note3: If the output already has the .bam alignments files from previous runs, the module will skip re-aligning the reads"
	echo ""
	echo "Options:"
	echo ""
	echo "	-a STR    metagenomic assembly file"
	echo "	-o STR    output directory"
	echo "	-t INT    number of threads (default=1)"
	echo "	-m INT		amount of RAM available (default=4)"
	echo "	-l INT		minimum contig length to bin (default=1000bp)."
	echo "	--single-end	non-paired reads mode (provide *.fastq files)"
	echo "	--interleaved	the input read files contain interleaved paired-end reads"
	echo "";}



########################################################################################################
########################               LOADING IN THE PARAMETERS                ########################
########################################################################################################

# setting scripts and databases from config file (should be in same folder as main script)
#config_file=$(which config-metawrap)
#source $config_file

SOFT="$( cd "$( dirname "$0"  )" && pwd  )"


#/home/wangzy/tools/test_tool/metabinner/scripts

comm () { ${SOFT}/print_comment.py "$1" "-"; }
error () { ${SOFT}/print_comment.py "$1" "*"; exit 1; }
warning () { ${SOFT}/print_comment.py "$1" "*"; }
announcement () { ${SOFT}/print_comment.py "$1" "#"; }



# default params
threads=1; mem=4; len=1000; out=false; ASSEMBLY=false
# long options defaults
read_type=paired


# load in params

# loop through input params
while true; do
        case "$1" in
                -t) threads=$2; shift 2;;
		-m) mem=$2; shift 2;;
                -o) out=$2; shift 2;;
                -a) ASSEMBLY=$2; shift 2;;
		-l) len=$2; shift 2;;
                -h | --help) help_message; exit 1; shift 1;;
		--single-end) read_type=single; shift 1;;
		--interleaved) read_type=interleaved; shift 1;;
                --) help_message; exit 1; shift; break ;;
                *) break;;
        esac
done

########################################################################################################
########################           MAKING SURE EVERYTHING IS SET UP             ########################
########################################################################################################
# Make sure at least one binning method was chosen

# check if all parameters are entered
if [ $out = false ] || [ $ASSEMBLY = false ] ; then 
	comm "Non-optional parameters -a and/or -o were not entered"
	help_message; exit 1
fi

#check if the assembly file exists
if [ ! -s $ASSEMBLY ]; then error "$ASSEMBLY does not exist. Exiting..."; fi

comm "Entered read type: $read_type"

if [ $read_type = paired ]; then
	# check for at least one pair of read fastq files:
	F="no"; R="no"
	for num in "$@"; do
		if [[ $num == *"_1.fastq" ]]; then F="yes"; fi
		if [[ $num == *"_2.fastq" ]]; then R="yes"; fi
	done
	if [ $F = "no" ] || [ $R = "no" ]; then
		comm "Unable to find proper fastq read pair in the format *_1.fastq and *_2.fastq"
		help_message; exit 1
	fi
else
	# check for at least one fastq read
	F="no"
	for num in "$@"; do
		if [[ $num == *".fastq" ]]; then F="yes"; fi
	done
	if [ $F = "no" ]; then
		comm "Unable to find read files in format *.fastq (for single-end or interleaved reads)"
		help_message; exit 1
	fi
fi

if [ $read_type = paired ]; then
	#determine number of fastq read files provided:
	num_of_F_read_files=$(for I in "$@"; do echo $I | grep _1.fastq; done | wc -l)
	num_of_R_read_files=$(for I in "$@"; do echo $I | grep _2.fastq; done | wc -l)

	comm "$num_of_F_read_files forward and $num_of_R_read_files reverse read files detected"
	if [ ! $num_of_F_read_files == $num_of_R_read_files ]; then error "Number of F and R reads must be the same!"; fi
fi

#
## Checks for correctly configures meta-scripts folder
#if [ ! -s $SOFT/sort_contigs.py ]; then
#	error "The folder $SOFT doesnt exist. Please make sure config.sh is in the same filder as the mains scripts and all the paths in the config.sh file are correct"
#fi


########################################################################################################
########################                    BEGIN PIPELINE!                     ########################
########################################################################################################


########################################################################################################
########################         ALIGNING READS TO MAKE COVERAGE FILES          ########################
########################################################################################################
announcement "ALIGNING READS TO MAKE COVERAGE FILES"

# setting up the output folder
if [ ! -d $out ]; then mkdir $out;
else 
	echo "Warning: $out already exists."
fi

if [ ! -d ${out}/work_files ]; then mkdir ${out}/work_files; fi

if [ -f ${out}/work_files/assembly.fa ]; then
	comm "Looks like the assembly file is already coppied, but will re-transfer just in case to avoid truncation problems."
	cp $ASSEMBLY ${out}/work_files/assembly.fa
else
	comm "making copy of assembly file $ASSEMBLY"
	cp $ASSEMBLY ${out}/work_files/assembly.fa
fi

tmp=${ASSEMBLY##*/}
sample=${tmp%.*}

# Index the assembly
if [ -f ${out}/work_files/assembly.fa.bwt ]; then
	comm "Looks like there is a index of the assembly already. Skipping..."
else
	comm "Indexing assembly file"
	bwa index ${out}/work_files/assembly.fa
	if [[ $? -ne 0 ]] ; then error "Something went wrong with indexing the assembly. Exiting."; fi
fi

# If there are several pairs of reads passed, they are processed sepperately
for num in "$@"; do
	# paired end reads
	if [ $read_type = paired ]; then
		if [[ $num == *"_1.fastq"* ]]; then 
			reads_1=$num
			reads_2=${num%_*}_2.fastq
			if [ ! -s $reads_1 ]; then error "$reads_1 does not exist. Exiting..."; fi
			if [ ! -s $reads_2 ]; then error "$reads_2 does not exist. Exiting..."; fi
	
			tmp=${reads_1##*/}
			sample=${tmp%_*}
			
			if [[ ! -f ${out}/work_files/${sample}.bam ]]; then
				comm "Aligning $reads_1 and $reads_2 back to assembly"
				bwa mem -v 1 -t $threads ${out}/work_files/assembly.fa $reads_1 $reads_2 > ${out}/work_files/${sample}.sam
				if [[ $? -ne 0 ]]; then error "Something went wrong with aligning $reads_1 and $reads_2 reads to the assembly. Exiting"; fi

				comm "Sorting the $sample alignment file"
				samtools sort -T ${out}/work_files/tmp-samtools -@ $threads -O BAM -o ${out}/work_files/${sample}.bam ${out}/work_files/${sample}.sam
				if [[ $? -ne 0 ]]; then error "Something went wrong with sorting the alignments. Exiging..."; fi
				rm ${out}/work_files/${sample}.sam
			else
				comm "skipping aligning $sample reads to assembly because ${out}/work_files/${sample}.bam already exists."
			fi
		fi

	# single end or interleaved reads
	else
		if [[ $num == *".fastq"* ]]; then
			reads=$num
			if [ ! -s $reads ]; then error "$reads does not exist. Exiting..."; fi
			tmp=${reads##*/}
			sample=${tmp%.*}
			if [[ ! -f ${out}/work_files/${sample}.bam ]]; then
				comm "Aligning $reads back to assembly, and sorting the alignment"
				if [ $read_type = single ]; then
					bwa mem -t $threads ${out}/work_files/assembly.fa $reads > ${out}/work_files/${sample}.sam
					if [[ $? -ne 0 ]]; then error "Something went wrong with aligning the reads to the assembly!"; fi
				elif [ $read_type = interleaved ]; then
					bwa mem -v 1 -p -t $threads ${out}/work_files/assembly.fa $reads > ${out}/work_files/${sample}.sam
					if [[ $? -ne 0 ]]; then error "Something went wrong with aligning the reads to the assembly!"; fi
				else
					error "something from with the read_type (=$read_type)"
				fi
				
				comm "Sorting the $sample alignment file"
				samtools sort -T ${out}/work_files/tmp-samtools -@ $threads -O BAM -o ${out}/work_files/${sample}.bam ${out}/work_files/${sample}.sam
				if [[ $? -ne 0 ]]; then error "Something went wrong with sorting the alignments. Exiging..."; fi
				rm ${out}/work_files/${sample}.sam
			else
				comm "skipping aligning $sample reads to assembly because ${out}/work_files/${sample}.bam already exists."
			fi
		fi
	fi
done




#if [ $maxbin2 = true ]; then
#        ########################################################################################################
#        ########################                   Making contig depth                    ########################
#        ########################################################################################################
#        announcement "Making contig depth "

comm "making contig depth file..."
      ${SOFT}/jgi_summarize_bam_contig_depths --outputDepth ${out}/work_files/mb2_master_depth.txt --noIntraDepthVariance ${out}/work_files/*.bam
      if [[ $? -ne 0 ]]; then error "Something went wrong with making contig depth file. Exiting."; fi

#	#calculate total numper of columns
#	A=($(head -n 1 ${out}/work_files/mb2_master_depth.txt))
#	N=${#A[*]}
#
#	# split the contig depth file into multiple files
#	comm "split master contig depth file into individual files for maxbin2 input"
#	if [ -f ${out}/work_files/mb2_abund_list.txt ]; then rm ${out}/work_files/mb2_abund_list.txt; fi
#	for i in $(seq 4 $N); do
#		sample=$(head -n 1 ${out}/work_files/mb2_master_depth.txt | cut -f $i)
#		echo "processing $sample depth file..."
#		grep -v totalAvgDepth ${out}/work_files/mb2_master_depth.txt | cut -f 1,$i > ${out}/work_files/mb2_${sample%.*}.txt
#		if [[ $out == /* ]]; then
#			echo ${out}/work_files/mb2_${sample%.*}.txt >> ${out}/work_files/mb2_abund_list.txt
#		else
#			echo $(pwd)/${out}/work_files/mb2_${sample%.*}.txt >> ${out}/work_files/mb2_abund_list.txt
#		fi
#	done

#comm "cleaning up *.bam to save space..."
#rm ${out}/work_files/*bam
cat ${out}/work_files/mb2_master_depth.txt | cut -f -1,4- > ${out}/coverage_profile.tsv
cat ${out}/work_files/mb2_master_depth.txt | awk '{if ($2>1000) print $0 }' | cut -f -1,4- > ${out}/coverage_profile_f1k.tsv
########################################################################################################
########################      BINNING PIPELINE SUCCESSFULLY FINISHED!!!         ########################
########################################################################################################
announcement "The process of generating coverage pipeline finished!!!"

