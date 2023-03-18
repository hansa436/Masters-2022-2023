#!/bin/bash


Hairpin=Hairpin_seq
Anitsense_RNF= Antisense RNF sequence
Anitsense_RNF_10bp_before= 10bp sequence of antisense RNF before of the hairpin
Anitsense_RNF_10bp_after= 10bp sequence of antisense RNF after of the hairpin
Sense_RNF= Sense RNF sequence
Sense_RNF_10bp_before= 10bp sequnce of sense RNF before of the hairpin
Sense_RNF_10bp_after= 10bp sequence of sense RNF after of the hairpin
Fastq=fastq_file.fastq
Output_full_reads=full_reads_output_file.txt
Output_components=components_output_file.txt

echo "Full hairpin reads" >${Output_full_reads}

echo "RWbasic_short_S_H_AS">>${Output_full_reads} 
./seqkit grep -s -m12 -i -p ${Sense_RNF}${Hairpin}${Antisense_RNF} ${Fastq} |  grep -c @M >>${Output_full_reads}

echo "RWbasic_short_AS_H_S">>${Output_full_reads}
./seqkit grep -s -m12 -i -p ${Antisense_RNF}${Hairpin}${Sense_RNF} ${Fastq} |  grep -c @M  >>${Output_full_reads}

echo "RWbasic_long_S_H_AS">>${Output_full_reads}
./seqkit grep -s -m27 -i -p ${Sense_RNF}${Hairpin}${Antisense_RNF} ${Fastq} |  grep -c @M >>${Output_full_reads}

echo "RWbasic_long_AS_H_S">>${Output_full_reads}
./seqkit grep -s -m27 -i -p ${Antisense_RNF}${Hairpin}${Sense_RNF} ${Fastq} |  grep -c @M >>${Output_full_reads}



echo "hpREADs componets" > ${Output_components}

echo " " >> ${Output_components}

echo "RWbasics long" >> ${Output_components}

echo "RWbasics full S_H" >> ${Output_components}

../seqkit grep -s -m15 -i -p ${Sense_RNF}${Hairpin} ${Fastq} |  grep -c @M  >> ${Output_components}

echo "RWbasics full AS_H" >> ${Output_components}

../seqkit grep -s -m15 -i -p ${Antisense_RNF}${Hairpin} ${Fastq} |  grep -c @M >> ${Output_components}

echo "RWbasics full H_S" >> ${Output_components}

../seqkit grep -s -m15 -i -p ${Hairpin}${Sense_RNF} ${Fastq} |  grep -c @M >> ${Output_components}

echo "RWbasics full H_AS" >> ${Output_components}

../seqkit grep -s -m15 -i -p ${Hairpin}${Antisense_RNF} ${Fastq} |  grep -c @M >> ${Output_components}

echo "RWbasics 10bp S_H_AS" >> ${Output_components}

../seqkit grep -s -m5 -i -p ${Sense_RNF_10bp_before}${Hairpin}${Antisense_RNF_10bp_after} ${Fastq} |  grep -c @M >> ${Output_components}

echo "RWbasics 10bp AS_H_S" >> ${Output_components}

../seqkit grep -s -m5 -i -p ${Antisense_RNF_10bp_before}${Hairpin}${Sense_RNF_10bp_after} ${Fastq} |  grep -c @M >> ${Output_components}

echo "RWbasics 10bp H_S" >> ${Output_components}

../seqkit grep -s -m4 -i -p ${Hairpin}${Sense_RNF_10bp_after} ${Fastq} |  grep -c @M >> ${Output_components}

echo "RWbasics 10bp H_AS" >> ${Output_components}

../seqkit grep -s -m4 -i -p ${Hairpin}${Antisense_RNF_10bp_after} ${Fastq} |  grep -c @M >> ${Output_components}

echo "RWbasics 10bp S_H" >> ${Output_components}

../seqkit grep -s -m4 -i -p ${Sense_RNF_10bp_before}${Hairpin} ${Fastq} |  grep -c @M >> ${Output_components}

echo "RWbasics 10bp AS_H" >> ${Output_components}

../seqkit grep -s -m4 -i -p ${Antisense_RNF_10bp_before}${Hairpin} ${Fastq} |  grep -c @M >> ${Output_components}

echo "RWbasics full S RNF" >> ${Output_components}

../seqkit grep -s -m12 -i -P -p ${Sense_RNF} ${Fastq} |  grep -c @M >> ${Output_components}

echo "RWbasics full AS RNF" >> ${Output_components}

../seqkit grep -s -m12 -i -P -p ${Antisense_RNF} ${Fastq} |  grep -c @M >> ${Output_components}

echo "RWbasics full Hairpin" >> ${Output_components}

../seqkit grep -s -m3 -i -p ${Hairpin} ${Fastq} |  grep -c @M >> ${Output_components}



echo " " >> ${Output_components}

echo "RWbasics  short" >> ${Output_components}

echo "RWbasics full S_H" >> ${Output_components}

../seqkit grep -s -m8 -i -p ${Sense_RNF}${Hairpin} ${Fastq} |  grep -c @M >> ${Output_components}

echo "RWbasics full AS_H" >> ${Output_components}

../seqkit grep -s -m8 -i -p ${Antisense_RNF}${Hairpin} ${Fastq} |  grep -c @M >> ${Output_components}

echo "RWbasics full H_S" >> ${Output_components}

../seqkit grep -s -m8 -i -p ${Hairpin}${Sense_RNF} ${Fastq} |  grep -c @M >> ${Output_components}

echo "RWbasics full H_AS" >> ${Output_components}

../seqkit grep -s -m8 -i -p ${Hairpin}${Antisense_RNF} ${Fastq} |  grep -c @M >> ${Output_components}

echo "RWbasics 10bp S_H_AS" >> ${Output_components}

../seqkit grep -s -m5 -i -p ${Sense_RNF_10bp_before}${Hairpin}${Antisense_RNF_10bp_after} ${Fastq} |  grep -c @M >> ${Output_components}

echo "RWbasics 10bp AS_H_S" >> ${Output_components}

../seqkit grep -s -m5 -i -p ${Antisense_RNF_10bp_before}${Hairpin}${Sense_RNF_10bp_after} ${Fastq} |  grep -c @M >> ${Output_components}

echo "RWbasics 10bp H_S" >> ${Output_components}

../seqkit grep -s -m4 -i -p ${Hairpin}${Sense_RNF_10bp_after} ${Fastq} |  grep -c @M >> ${Output_components}

echo "RWbasics 10bp H_AS" >> ${Output_components}

../seqkit grep -s -m4 -i -p ${Hairpin}${Antisense_RNF_10bp_after} ${Fastq} |  grep -c @M >> ${Output_components}

echo "RWbasics 10bp S_H" >> ${Output_components}

../seqkit grep -s -m4 -i -p ${Sense_RNF_10bp_before}${Hairpin} ${Fastq} |  grep -c @M >> ${Output_components}

echo "RWbasics 10bp AS_H" >> ${Output_components}

../seqkit grep -s -m4 -i -p ${Antisense_RNF_10bp_before}${Hairpin} ${Fastq} |  grep -c @M >> ${Output_components}

echo "RWbasics full S RNF" >> ${Output_components}

../seqkit grep -s -m4 -i -P -p ${Sense_RNF} ${Fastq} |  grep -c @M >> ${Output_components}

echo "RWbasics full AS RNF" >> ${Output_components}

../seqkit grep -s -m4 -i -P -p ${Antisense_RNF} ${Fastq} |  grep -c @M >> ${Output_components}

echo "RWbasics full Hairpin" >> ${Output_components}

../seqkit grep -s -m3 -i -p ${Hairpin} ${Fastq} |  grep -c @M >> ${Output_components}