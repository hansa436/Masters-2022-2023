#!/bin/bash

ID=YourID

Input=filename.fastq

Adaptor_seq=AAAA
Complement_adaptor_seq=TTTT

mkdir ${ID}

######################################

echo "Extraction of ctDNA commensing" >>${ID}/Pipeline_log.txt

#this find the positions of the adaptor sequnces in the read and outputs it into a file
# extract the read ids (awk) for subsetting to get sense and antisense reads

./seqkit locate -d -i -P -p ${Adaptor_seq} ${Input} |  awk '{print $1}' OFS='\t' ${ID}/linker_position_sense.tsv> ${ID}/ctDNA_UMI_read_IDs_sense.txt

./seqkit locate -d -i -P -p ${Complement_adaptor_seq} ${Input}  | awk '{print $1}' OFS='\t' ${ID}/linker_position_antisense.tsv> ${ID}/ctDNA_UMI_read_IDs_antisense.txt

# subsets the fastq file with all reads to produce a fastq file with only reads that contains the ctDNA

seqtk subseq ${Input} ${ID}/ctDNA_UMI_read_IDs_sense.txt > ${ID}/reads_with_ctDNA_sense.fastq

seqtk subseq ${Input} ${ID}/ctDNA_UMI_read_IDs_antisense.txt > ${ID}/reads_with_ctDNA_antisense.fastq

# reverse complements the antisense reads to sense orientation (-r reverses the fastq sequences)

./seqkit seq -r -p ${ID}/reads_with_ctDNA_antisense.fastq > ${ID}/reads_with_ctDNA_reverse_complement.fastq

# merge complemented file into file with sense orientated reads

cat ${ID}/reads_with_ctDNA_reverse_complement.fastq ${ID}/reads_with_ctDNA_sense.fastq > ${ID}/reads_with_ctDNA.fastq

# re-extract the insert position now all reads are in the same orientation and amke a psotion file based on that

./seqkit locate -d -i -P -p ${Adaptor_seq} ${ID}/reads_with_ctDNA.fastq  > ${ID}/linker_position.tsv

# extract only relevant info into a a file with just the readID start and end position

awk '{print $1, $5,$6}' OFS='\t' ${ID}/linker_position.tsv> ${ID}/start_end_pos.bed

#this adds 10 to before the linker (UMI) and 25bp to after  the linker (ctDNA) to extract the UMI and ctDNA insert sequneces

awk '{$2=$2-11}1' OFS='\t' ${ID}/start_end_pos.bed | awk '{$3=$3+26}1' OFS='\t' > ${ID}/alt_start_end_pos.bed

# this removes the header (makes downstream analysis easier)

sed -i '1d' ${ID}/alt_start_end_pos.bed

# Removes any reads that are less than 1 (can occur when adding 10bp to start postion ) and creates a bed file 

awk '(NR==1)||($2 > 1)' ${ID}/alt_start_end_pos.bed > ${ID}/ctDNA_UMI_location.bed

echo "Extraction of reads with ctDNA completed" >> ${ID}/Pipeline_log.txt


#######################################

# this section is optional, only required if the data needs demultiplexing based on a barcode sequence

echo "Extraction of reads with barcode commencing" >> insert/${ID1}/Pipeline_log.txt

./seqkit grep -m 3 -s -i -p ${Barcode} ${ID}/reads_with_ctDNA.fastq > ${ID}/reads_with_ctDNA_barcode.fastq

echo "Extraction of reads with barcode completed" >> insert/${ID1}/Pipeline_log.txt


#######################

echo "Mapping of reads to the genome commencing" >> ${ID}/Pipeline_log.txt

# Map the genes to the genome to get the biologically relevant genes (Minimap2 for ONT fgneerated data BWA from Illumina generated data)
# the two mappigns are the same, it is just dependednt on if a barcode was used or not, just # out the one not required.

#./minimap2/minimap2 -I8G -ax map-ont ./genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa ${ID}/reads_with_ctDNA.fastq > ${ID}/amplicons.bam

./minimap2/minimap2 -I8G -ax map-ont ./genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa ${ID}/reads_with_ctDNA_barcode.fastq > ${ID}/amplicons.bam

# sort bam file

samtools sort ${ID}/amplicons.bam > ${ID}/mapped_amplicons_sorted.bam

samtools index ${ID}/mapped_amplicons_sorted.bam

#to get get the bed coverage

samtools bedcov hg38_QIAseq_DNA_panel_no_chr.roi.bed ${ID}/mapped_amplicons_sorted.bam > ${ID}/hg38_coverage_report.bed

#removes the blank extra coloumn in there:

awk -F "\t" 'OFS="\t" {print $1,$2,$3,$4,$6}' ${ID}/hg38_coverage_report.bed >${ID}/hg38_bedcov_report.bed

# remove the original output:

rm ${ID}/hg38_coverage_report.bed

echo "Mapping of reads to the genome completed" >> ${ID}/Pipeline_log.txt

###########################

echo "Extraction of biologically relevant reads commencing" >> ${ID}/Pipeline_log.txt

# Extract the biolgically relevant reads

mkdir genes/${ID}

samtools view ${ID}/mapped_amplicons_sorted.bam 1:26696363-26780766| cut -d$'\t' -f1 |sort | uniq> genes/${ID}/ARID1A_IDs.txt

samtools view ${ID}/mapped_amplicons_sorted.bam 2:140233175-1142130739| cut -d$'\t' -f1 |sort | uniq> genes/${ID}/LRP1B_IDs.txt

samtools view ${ID}/mapped_amplicons_sorted.bam 2:211383604-212538540 | cut -d$'\t' -f1 |sort | uniq > genes/${ID}/ERBB4_IDs.txt

samtools view ${ID}/mapped_amplicons_sorted.bam 3:41224058-41239352 | cut -d$'\t' -f1 |sort | uniq > genes/${ID}/CTNNB1_IDs.txt

samtools view ${ID}/mapped_amplicons_sorted.bam 3:49360198-49375599| cut -d$'\t' -f1 |sort | uniq> genes/${ID}/RHOA_IDs.txt

samtools view ${ID}/mapped_amplicons_sorted.bam 3:179198815-179234374| cut -d$'\t' -f1 |sort | uniq> genes/${ID}/PIK3CA_IDs.txt

samtools view ${ID}/mapped_amplicons_sorted.bam 4:152322870-152411813| cut -d$'\t' -f1 |sort | uniq> genes/${ID}/FBXWF_IDs.txt

samtools view ${ID}/mapped_amplicons_sorted.bam 5:112707707-112844136| cut -d$'\t' -f1 |sort | uniq > genes/${ID}/APC_IDs.txt

samtools view ${ID}/mapped_amplicons_sorted.bam 7:55019267-55205627| cut -d$'\t' -f1  |sort | uniq> genes/${ID}/EGFR_IDs.txt

samtools view ${ID}/mapped_amplicons_sorted.bam 7:140716017-140924718| cut -d$'\t' -f1 |sort | uniq > genes/${ID}/BRAF_IDs.txt

samtools view ${ID}/mapped_amplicons_sorted.bam 9:21968218-21994464| cut -d$'\t' -f1 |sort | uniq > genes/${ID}/CDKN2A_IDs.txt

samtools view ${ID}/mapped_amplicons_sorted.bam 10:87863886-87965537| cut -d$'\t' -f1 |sort | uniq > genes/${ID}/PTEN_IDs.txt

samtools view ${ID}/mapped_amplicons_sorted.bam 11:108227609-108365518| cut -d$'\t' -f1 |sort | uniq > genes/${ID}/ATM_IDs.txt

samtools view ${ID}/mapped_amplicons_sorted.bam 12:25209784-25245394| cut -d$'\t' -f1 |sort | uniq > genes/${ID}/KRAS_IDs.txt

samtools view ${ID}/mapped_amplicons_sorted.bam 12:56080290-56102065| cut -d$'\t' -f1 |sort | uniq > genes/${ID}/ERBB3_IDs.txt

samtools view ${ID}/mapped_amplicons_sorted.bam 13:32316450-32398780| cut -d$'\t' -f1 |sort | uniq > genes/${ID}/BRCA2_IDs.txt

samtools view ${ID}/mapped_amplicons_sorted.bam 16:68737405-68833509| cut -d$'\t' -f1 |sort | uniq > genes/${ID}/CHD1_IDs.txt

samtools view ${ID}/mapped_amplicons_sorted.bam 17:7668197-7676604| cut -d$'\t' -f1 |sort | uniq > genes/${ID}/TP53_IDs.txt

samtools view ${ID}/mapped_amplicons_sorted.bam 17:39699549-39728054| cut -d$'\t' -f1 |sort | uniq > genes/${ID}/ERBB2_IDs.txt

samtools view ${ID}/mapped_amplicons_sorted.bam 17:58354932-58415587| cut -d$'\t' -f1 |sort | uniq > genes/${ID}/RNF43_IDs.txt

samtools view ${ID}/mapped_amplicons_sorted.bam 18:51047036-51078477| cut -d$'\t' -f1 |sort | uniq > genes/${ID}/SMAD4_IDs.txt

cat genes/${ID}/*.txt >${ID}/all_targeted_reads.txt

seqtk subseq ${ID}/reads_with_ctDNA.fastq ${ID}/all_targeted_reads.txt > ${ID}/bio_relevant_reads.fastq

echo "Extraction of biologically relevant reads completed" >> ${ID}/Pipeline_log.txt

#################################

echo "Extraction of UMI/CTDNA sequences commencing" >> ${ID}/Pipeline_log.txt

# UMI/ctDNA sequences

# only extract sequence of insterest
# First line outputs the read ID and the sequence generating a fatsa file
# second line subsets it by the umi-ctDNA locations in the bed files to generate a fasta with only the UMI-linker-ctDNA sequence 

cat ${ID}/bio_relevant_reads.fastq | awk '{if(NR%4==1) {printf(">%s\n",substr($0,2));} else if(NR%4==2) print;}' > ${ID}/reads_with_ID_linker.fasta
bedtools getfasta -fi ${ID}/reads_with_ID_linker.fasta -bed ${ID}/ctDNA_UMI_location.bed > ${ID}/UMI_linker_seq.fasta

# This extacts only the sequence removing read names
awk 'NR % 2 == 0' ${ID}/UMI_linker_seq.fasta> ${ID}/UMI_insert_seq.txt

# extract the Umi sequences and the insert sequences then put them together without the linker for comparison of thow many times the same UMI-linker sequence occurs
# The UMI sequnce is from 1-11bp, the ctDNA is from 22 till end of read, so by extracting those regions using cut the Linker region is removed from the analysis
# to see how many times this occured in the dataset the sequences are sorted by uniq -c which counts every time the same sequence is seen, and sorting by -bgr lists the reads in deceding oorder of most copies to least copies 
cut -c1-11,22- ${ID}/UMI_insert_seq.txt > ${ID}/UMI_insert_sequences.txt
sort ${ID}/UMI_insert_sequences.txt | uniq -c |sort -bgr > ${ID}/UMI_Insert_sequences_Unique.txt

echo "Extraction of UMI/ctDNA sequences Completed" >> ${ID}/Pipeline_log.txt

echo "Script is fully complete!" >> ${ID}/Pipeline_log.txt
