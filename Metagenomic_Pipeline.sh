#Metagenomics pipeline
#Initiated: October 30, 2019 by Rachael Storo
#Last Edit: October 30, 2019 by Rachael Storo
#Purpose: Create a script for the processing of environmental metagenomic samples
#Disclaimer: No pipeline or script should be followed blindly. Ensure that this script and the programs used make sense for the samples you have.
#You may or may not utilize all steps of this pipeline- that depends entirely on your question. Understanding the pipeline and your question will save you time by not running steps you don't need.
#Note: If you are performing this on the cluster, all programs are installed. Otherwise, ensure you install all programs and dependencies.
#This is written for paired-end data.

#Step One- Demultiplex
# Not all data will require this step- in fact, many metagenomes are sequenced individually or are given to you already demultiplexed.

usearch -fastx_demux R1.fq -reverse R2.fq -index I1.fq -barcodes bar.fa \
  -fastqout fwd_demux.fq -output2 rev_demux.fq

#Step Two- Quality trimming and filtering

java -jar trimmomatic-0.39.jar PE input_forward.fq.gz input_reverse.fq.gz output_forward_paired.fq.gz output_forward_unpaired.fq.gz output_reverse_paired.fq.gz output_reverse_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:36

#Step Three- No Assembly

metaphlan2.py SRS014476-Supragingival_plaque.fasta.gz  --input_type fasta > SRS014476-Supragingival_plaque_profile.txt

#Step Four- Assembly

megahit -1 pe_1.fq -2 pe_2.fq -o out

#Step Five- Gene Calling and Taxonomic Profiling

prokka contigs.fa

#Step Six- Read Based Analysis

phyloseq

#Step Seven- Mapping
bowtie2

#Step Eight- Recovering Metagenome Assembled Genomes

anvio

#Step Nine- Phylogenomics

anvio

#Step Ten- Comparative Genomics

anvio

#Step Eleven- Pangenomics

anvio
