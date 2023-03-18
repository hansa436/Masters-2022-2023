
#!/bin/bash

# For maping files to put into IGV to see where things load

hairpin=hairpin_name
fastq=fastq_file.fastq

#to Loaad into IGV the fasta (reference) files need a .fai index
samtools faidx ${hairpin}_long.fa
samtools faidx ${hairpin}_short.fa

#for bwa mapping these need to be index
bwa index ${hairpin}_long.fa
bwa index ${hairpin}_short.fa


# The layout of each thing is
# 1. map the fstq using the related fasta file
# 2. sort it
# 3. sort by maped genes (-F 4 flag)
# 4. Sort by primary alignments only (-F 256 flag)
# 5. make an index file for the mapped primary lignments
# 6. remove unwanted intermediate files
# 7. Profit (Jk load into IGV for analyis)


# this maps the read to the fatsa file, then outputs that to a bam file
#the bam file is then sorted and indexed so it can go into IGV
#by using the primary mapping ones we can exclude anything that doesn't map and only take things that map

bwa mem ${hairpin)_short.fa mixexorw_S1_L001_R1_001.fastq.gz > ${hairpin)_short.bam
samtools sort ${hairpin)_short.bam> ${hairpin)_short_sorted.bam
samtools view -b -F 4 ${hairpin)_short_sorted.bam > ${hairpin)_mapped.bam
samtools view -b -F 256 ${hairpin)_mapped.bam > ${hairpin)_short_primary_mapped.bam
samtools index ${hairpin)_short_primary_mapped.bam
rm ${hairpin)_short.bam ${hairpin)_mapped.bam

bwa mem ${hairpin)_long.fa mixexorw_S1_L001_R1_001.fastq.gz > rwbasic_long.bam
samtools sort rwbasic_long.bam> ${hairpin)_long_sorted.bam
samtools view -b -F 4 ${hairpin)_long_sorted.bam > ${hairpin)_mapped.bam
samtools view -b -F 256 ${hairpin)_mapped.bam > ${hairpin)_long_primary_mapped.bam
samtools index ${hairpin)_long_primary_mapped.bam
rm {hairpin)_long.bam ${hairpin)_mapped.bam





