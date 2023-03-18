#!/bin/bash

# Megldodon script
# 	Extracted from:https://github.com/nanoporetech/megalodon
# Example command to output basecalls, mappings, and CpG 5mC and 5hmC methylation in both per-read (``mod_mappings``) and aggregated (``mods``) formats
#   Compute settings: GPU devices 0 and 1 with 20 CPU cores


raw_fast5s=27_5mc_only.fast5
reference=../genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa
# output=./megalodon/27ac
#guppy_path=../../ont-guppy-cpu/bin/

megalodon \
    ${raw_fast5s} \
 	--guppy-server-path ont-guppy-cpu/bin/guppy_basecall_server \
	--guppy-config dna_r9.4.1_450bps_fast.cfg \
   	--remora-modified-bases dna_r9.4.1_e8 fast 0.0.0 5hmc_5mc CG 0 \
	--outputs basecalls mappings mod_mappings mods \
    	--reference ${reference} \
    	--devices 0 1 \
	--guppy-server-path /data/ont-guppy/bin/guppy_basecall_server

