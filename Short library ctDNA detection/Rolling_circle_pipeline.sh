#!/bin/bash


ID=mapped_R2C2_consesus
in=R2C2_Consensus.fasta
genome= hg38.fa
genes=panel_genes

mkdir ${ID}


# Map the genes to the genome to get the biologically relevant genes

echo "Mapping intiated"
../../minimap2/minimap2 -I8G -ax map-ont --secondary=no ${genome} ${in} > ${ID}/${ID}_mapped_amplicons.bam

echo "Mapping complete"


# sort bam file

echo "Sorting intiated"

samtools sort ${ID}/${ID}_mapped_amplicons.bam > ${ID}/${ID}_mapped_amplicons_sorted.bam

samtools index ${ID}/${ID}_mapped_amplicons_sorted.bam

echo "Sorting complete"

# filter bam file

# the other way to do it is by sbestting the mapping file by the sam-flags( good if you're trying to map to a specific coordinate

# Samtools -F works as an inverse search so it will only keep things that don't have that parameter
# (optional -b) will output in a bam format
# sam flags ( the second coloumn of the bam file)
        ## 0= mapped
        ## 4= unmapped
        ## 16= read on the reverse strand
        ## 256= not primary aligment
        ## 272= not primary match reverse strand
        ## 2048= suplementary alignment
        # 2064= Supplmentary alignment on the reverse strand
# if you ever want to look them up you can use this website and type in the flag: https:/broadinstitute.github.io/picard/explain-flags.html
# have also had issues doing it a on one line so its easier to do it and make mulitple output files
# so remove the forward flag i.e 256 or 2048 also removes the complmenatry flag i.e 272 or 2064


echo "Sam flag filtering intiated"

samtools view -F4, -F256, -F2048 ${ID}_mapped_amplicons_sorted.bam > ${ID}/flag_filters_primary_only.bam


# sort bam file

echo "Sorting filtered intiated"


samtools sort ${ID}/flag_filters_primary_only.bam > ${ID}/flag_filters_primary_only_sorted.bam

samtools index ${ID}/flag_filters_primary_only_sorted.bam

echo "Sorting complete"



rm ${ID}/flag_filters_primary_only.bam

echo "Sam flag filtering completed"


# to get the bed coverage
# using the panel of /filtered_genes I'm interested in I can get the coverage


echo "Coverage intiated"
samtools bedcov ../../panel/hg38_QIAseq_DNA_panel_no_chr.roi.bed ${ID}/${ID}_mapped_amplicons_sorted.bam> ${ID}/hg38_coverage_report.bed

# To remove the weird extra coloumn in there:

awk -F "\t" 'OFS="\t" {print $1,$2,$3,$4,$6}' ${ID}/hg38_coverage_report.bed >${ID}/hg38_bedcov_report_panel.bed

# remove the orignial output:
rm ${ID}/hg38_coverage_report.bed

echo "Coverage complete"

# Extract the biolgically relevant reads

echo "Extraction of biologically relevant reads intiated"

echo "Extraction of panel reads intiated"
# samtools view: this opens up the maped bam file
# Coordinates: the coorinates or your region of interest
# cut -d sets delinit=atr which will be tab
# cut -f1 will extract the read ids only
# sort | uniq will provide only one of the same read ID to revent double mapping of the same read

mkdir ${ID}/${genes}

samtools view ${ID}/${ID}_mapped_amplicons_sorted.bam 1:26696363-26780766| cut -d$'\t' -f1 |sort | uniq > ${ID}/${genes}/ARID1A_IDs.txt

samtools view ${ID}/${ID}_mapped_amplicons_sorted.bam 2:140233175-142130739| cut -d$'\t' -f1 |sort | uniq > ${ID}/${genes}/LRP1B_IDs.txt

samtools view ${ID}/${ID}_mapped_amplicons_sorted.bam 2:211383604-212538540 | cut -d$'\t' -f1 |sort | uniq > ${ID}/${genes}/ERBB4_IDs.txt

samtools view ${ID}/${ID}_mapped_amplicons_sorted.bam 3:41224058-41239352 | cut -d$'\t' -f1 |sort | uniq > ${ID}/${genes}/CTNNB1_IDs.txt

samtools view ${ID}/${ID}_mapped_amplicons_sorted.bam 3:49360198-49375599| cut -d$'\t' -f1 |sort | uniq > ${ID}/${genes}/RHOA_IDs.txt

samtools view ${ID}/${ID}_mapped_amplicons_sorted.bam 3:179198815-179234374| cut -d$'\t' -f1 |sort | uniq > ${ID}/${genes}/PIK3CA_IDs.txt

samtools view ${ID}/${ID}_mapped_amplicons_sorted.bam 4:152322870-152411813| cut -d$'\t' -f1 |sort | uniq > ${ID}/${genes}/FBXWF_IDs.txt

samtools view ${ID}/${ID}_mapped_amplicons_sorted.bam 5:112707707-112844136| cut -d$'\t' -f1 |sort | uniq > ${ID}/${genes}/APC_IDs.txt

samtools view ${ID}/${ID}_mapped_amplicons_sorted.bam 7:55019267-55205627| cut -d$'\t' -f1  |sort | uniq > ${ID}/${genes}/EGFR_IDs.txt

samtools view ${ID}/${ID}_mapped_amplicons_sorted.bam 7:140716017-140924718| cut -d$'\t' -f1 |sort | uniq > ${ID}/${genes}/BRAF_IDs.txt

samtools view ${ID}/${ID}_mapped_amplicons_sorted.bam 9:21968218-21994464| cut -d$'\t' -f1 |sort | uniq > ${ID}/${genes}/CDKN2A_IDs.txt

samtools view ${ID}/${ID}_mapped_amplicons_sorted.bam 10:87863886-87965537| cut -d$'\t' -f1 |sort | uniq > ${ID}/${genes}/PTEN_IDs.txt

samtools view ${ID}/${ID}_mapped_amplicons_sorted.bam 11:108227609-108365518| cut -d$'\t' -f1 |sort | uniq > ${ID}/${genes}/ATM_IDs.txt

samtools view ${ID}/${ID}_mapped_amplicons_sorted.bam 12:25209784-25245394| cut -d$'\t' -f1 |sort | uniq > ${ID}/${genes}/KRAS_IDs.txt

samtools view ${ID}/${ID}_mapped_amplicons_sorted.bam 12:56080290-56102065| cut -d$'\t' -f1 |sort | uniq > ${ID}/${genes}/ERBB3_IDs.txt

samtools view ${ID}/${ID}_mapped_amplicons_sorted.bam 13:32316450-32398780| cut -d$'\t' -f1 |sort | uniq > ${ID}/${genes}/BRCA2_IDs.txt

samtools view ${ID}/${ID}_mapped_amplicons_sorted.bam 16:68737405-68833509| cut -d$'\t' -f1 |sort | uniq > ${ID}/${genes}/CHD1_IDs.txt

samtools view ${ID}/${ID}_mapped_amplicons_sorted.bam 17:7668197-7676604| cut -d$'\t' -f1 |sort | uniq > ${ID}/${genes}/TP53_IDs.txt

samtools view ${ID}/${ID}_mapped_amplicons_sorted.bam 17:39699549-39728054| cut -d$'\t' -f1 |sort | uniq > ${ID}/${genes}/ERBB2_IDs.txt

samtools view ${ID}/${ID}_mapped_amplicons_sorted.bam 17:58354932-58415587| cut -d$'\t' -f1 |sort | uniq > ${ID}/${genes}/RNF43_IDs.txt

samtools view ${ID}/${ID}_mapped_amplicons_sorted.bam 18:51047036-51078477| cut -d$'\t' -f1 |sort | uniq > ${ID}/${genes}/SMAD4_IDs.txt

cat ${ID}/${genes}/*.txt > ${ID}/${genes}_mapped_reads_filtered.txt

echo "Extraction of biologically relevant reads complete"

echo "Subsetting intial file by reads in panel intiated"

echo "Gastric filtered"

seqtk subseq ${in} ${ID}/mapped_reads_flag_filtered.txt> ${ID}/Reads_in_gastic_panel.fasta

echo "Subset complete"


echo "Filter bam by panel commensing/getting mapping quality"

#echo "Gastric filtered"

#samtools view ${ID}/${ID}_mapped_amplicons_sorted.bam | fgrep -f ${ID}/Reads_in_gastic_panel_fasta.fasta > ${ID}/filtered_read_panel.sam

#cut -f 1,5 filtered_read_panel.sam > mapping_scores_filtered.txt

echo "Variant calling commensing"

# look for variants

bcftools mpileup -Ou -f ${genome} ${ID}/flag_filters_primary_only_sorted.bam| bcftools call -vmO z -o mapped_R2C2_consesus/hg38_consensus_reads.vcf.gz
echo "variant calling complete"



echo "hg38_consensus_reads.vcf.gz completed"

#### Variant subset by poziton

echo "Starting gastric panel extraction per gene"

tabix -h hg38_consensus_reads.vcf.gz 1:26696363-26780766 > mutation_info/consensus/ARID1A.vcf
tabix -h hg38_consensus_reads.vcf.gz 2:140233175-142130739 > mutation_info/consensus/LRP1B.vcf
tabix -h hg38_consensus_reads.vcf.gz 2:211383604-212538540 > mutation_info/consensus/ERBB4.vcf
tabix -h hg38_consensus_reads.vcf.gz 3:41224058-41239352 > mutation_info/consensus/CTNNB1.vcf
tabix -h hg38_consensus_reads.vcf.gz 3:49360198-49375599 > mutation_info/consensus/RHOA.vcf
tabix -h hg38_consensus_reads.vcf.gz 3:179198815-179234374 > mutation_info/consensus/PIK3CA.vcf
tabix -h hg38_consensus_reads.vcf.gz 4:152322870-152411813 > mutation_info/consensus/FBXWF.vcf
tabix -h hg38_consensus_reads.vcf.gz 5:112707707-112844136> mutation_info/consensus/APC.vcf
tabix -h hg38_consensus_reads.vcf.gz 7:55019267-55205627 > mutation_info/consensus/EGFR.vcf
tabix -h hg38_consensus_reads.vcf.gz 7:140716017-140924718 > mutation_info/consensus/BRAF.vcf
tabix -h hg38_consensus_reads.vcf.gz 9:21968218-21994464 > mutation_info/consensus/CDKN2A.vcf
tabix -h hg38_consensus_reads.vcf.gz 10:87863886-87965537 > mutation_info/consensus/PTEN.vcf
tabix -h hg38_consensus_reads.vcf.gz 11:108227609-108365518 > mutation_info/consensus/ATM.vcf
tabix -h hg38_consensus_reads.vcf.gz 12:25209784-25245394 > mutation_info/consensus/KRAS.vcf
tabix -h hg38_consensus_reads.vcf.gz 12:56080290-56102065 > mutation_info/consensus/ERBB3.vcf
tabix -h hg38_consensus_reads.vcf.gz 13:32316450-32398780 > mutation_info/consensus/BRCA2.vcf
tabix -h hg38_consensus_reads.vcf.gz 16:68737405-68833509 > mutation_info/consensus/CHD1.vcf
tabix -h hg38_consensus_reads.vcf.gz 17:7668197-7676604 > mutation_info/consensus/TP53.vcf
tabix -h hg38_consensus_reads.vcf.gz 17:39699549-39728054 > mutation_info/consensus/ERBB2.vcf
tabix -h hg38_consensus_reads.vcf.gz 17:58354932-58415587 > mutation_info/consensus/RNF43.vcf
tabix -h hg38_consensus_reads.vcf.gz 18:51047036-51078477 > mutation_info/consensus/SMAD4.vcf

echo "Completed gastric panel extraction per gene"
