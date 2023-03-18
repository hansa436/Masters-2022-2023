#!/bin/bash

# Split Reads for Concatination Run

echo setting paramters
 


#parameters
#FQ1=Fastq_output run1

FQ1=run1
INPUT=fastq_file_location
adaptor_sequence=Adaptor

#make directory for Fastq location
mkdir ${FQ1}
mkdir ${FQ2}

# set mismatch allowance
mm=4

#OF1=Output File location for run 1
#OF2=Output File location for run 2
#mm= mismatches
OF1=${FQ1}/${mm}mm


#make directory for Output files
mkdir ${OF1}


#RD1=Run_id for run1
#RD2=Run_id for run2
RD1=R1_${mm}mm



echo Starting ${RD1} analysis

echo Producing Fasta/Fastq files

# make fastq/fasta files
## only need to do thos of it hasn't previously been done

# run 1

#only required if not made perviously/*
 ## make a fastq file with all the reads
#{
zcat ${INPUT}/* > ${FQ1}/${RD1}.fastq
	
## Turn the fastq into fasta file for easy parsing of seqkit subseq

	### cat reds the file
	### paste - - - - means create 4 coloumns
	### sed removed the @ infront of the fastq Id seq and replaces it with > to turn it into the fasta formt
	### cut cuts feilds 1 and 2
	### tr changes tab (|t) separation to new line (\n) to make it compatible with fasta format

cat ${FQ1}/${RD1}.fastq | paste - - - - | sed 's/^@/>/g'| cut -f1-2 | tr '\t' '\n' > ${FQ1}/${RD1}.fasta
#}

#####################################################################
# Finding adaptor

echo Finding Adaptors

## find the strand and location of the adaptor for each read
./seqkit locate -m ${mm} -f adaptor_seq.fasta ${FQ1}/${RD1}.fastq | cut -f 1,5,6 > ${OF1}/location_of_adapters_${RD1}.txt

## find how many times each stand has an adapter
cut -f1 ${OF1}/location_of_adapters_${RD1}.txt | sort | uniq -c | sort -bgr > ${OF1}/number_of_adaptors_per_read_${RD1}.txt


# Split the reads run 

## create a bed file using seqkit locate to find adaptors

./seqkit locate --bed -m ${mm} -f adaptor_seq.fasta ${FQ1}/${RD1}.fastq > ${OF1}/location_of_adapters_${RD1}.bed

echo Finding reads

## take the complement of the bed file to find everything that is not an adaptor

## to get this to weork we need to create a list of the lengths of each reads and the ids for the genome file

## get list of length
## for every sequence print the lengths
cat ${FQ1}/${RD1}.fastq | awk '{if(NR%4==2) print length($1)}' > ${OF1}/${RD1}_length.txt

##get read IDs only
## Find lines that start wit @ as these signify the read_ids 
#use grep to extract he first 37 positions which corresponds to the read ids
grep ^@ ${FQ1}/${RD1}.fastq | grep -Po "^....................................." | sed 's/@//g' > ${OF1}/${RD1}_read_ids.txt

## Form one file with two columns
paste ${OF1}/${RD1}_read_ids.txt ${OF1}/${RD1}_length.txt > ${OF1}/sudo_genome_${RD1}.txt

## sort each file for this to work
bedtools sort -i ${OF1}/location_of_adapters_${RD1}.bed > ${OF1}/location_of_adapters_${RD1}_sorted.bed
sort ${OF1}/sudo_genome_${RD1}.txt > ${OF1}/sudo_genome_${RD1}_sorted.txt

echo splitting reads 

bedtools complement -i ${OF1}/location_of_adapters_${RD1}_sorted.bed -g ${OF1}/sudo_genome_${RD1}_sorted.txt > ${OF1}/read_section_without_adaptor_${RD1}.txt


## split using seqkit subseq
	#input is bed and fasta

./seqkit subseq --bed ${OF1}/read_section_without_adaptor_${RD1}.txt ${FQ1}/${RD1}.fasta > ${OF1}/reads_split_based_on_adaptor_${RD1}.fasta



# Map to human genome

echo mapping to human genome

## convert fasta to fastq
	# -F '#' adds a fake quality score for fastq conversion
seqtk seq -F '#' ${OF1}/reads_split_based_on_adaptor_${RD1}.fasta > ${OF1}/reads_split_based_on_adaptor_${RD1}.fastq

## Map using Minimap
./minimap2/minimap2 -I8G -ax map-ont ../genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa ${OF1}/reads_split_based_on_adaptor_${RD1}.fastq > ${OF1}/mapped_split_reads_${RD1}.bam

echo extracting coverage

## look at the coverage 
samtools sort ${OF1}/mapped_split_reads_${RD1}.bam > ${OF1}/mapped_split_reads_${RD1}_sorted.bam
samtools coverage ${OF1}/mapped_split_reads_${RD1}_sorted.bam > ${OF1}/coverage_mapped_reads_${RD1}.txt

echo filtering by 50bp min

# filter by length
mkdir ${OF1}/50bp_filtered
## remove everything lss than 50bp
cat ${OF1}/reads_split_based_on_adaptor_${RD1}.fastq | paste - - - - | awk 'length($2)  >= 50' | sed 's/\t/\n/g' > ${OF1}/50bp_filtered/length_filtered_gDNA_reads_${RD1}.fastq

echo mapping filtered reads genome

## map to the genoome

./minimap2/minimap2 -I8G -ax map-ont ../umi/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa ${OF1}/50bp_filtered/length_filtered_gDNA_reads_${RD1}.fastq > ${OF1}/50bp_filtered/mapped_length_filtered_reads_${RD1}.bam

echo extracting filtered reads coverage

## look at the coverage 
samtools sort ${OF1}/50bp_filtered/mapped_length_filtered_reads_${RD1}.bam > ${OF1}/50bp_filtered/mapped_length_filtered_reads_${RD1}_sorted.bam
samtools coverage ${OF1}/50bp_filtered/mapped_length_filtered_reads_${RD1}_sorted.bam > ${OF1}/50bp_filtered/coverage_mapped_filtered_reads_{$RD1}.txt

echo removing all temp/intermediate files

# Remove all tempory/intermediate files
rm ${OF1}/mapped_split_reads_${RD1}.bam
rm ${OF1}/sudo_genome_${RD1}.txt 
rm ${OF1}/location_of_adapters_${RD1}.bed
rm ${OF1}/50bp_filtered/mapped_length_filtered_reads_${RD1}.bam

echo ${RD1} analysis complete


