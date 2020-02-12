#Oil Fluff Metagenomics Pipeline
#Initiated: October 30, 2019 by Rachael Storo
#Last Edit: February 11, 2020 by Rachael Storo

#Purpose: Create a script for the processing of oil fluff metagenomic samples

#Read this paper before use of this pipeline, it will provide a visual context to what the majority of these steps are doing;
#Eren, A. & Esen, Özcan & Quince, Christopher & Vineis, Joseph & Morrison, Hilary & Sogin, Mitchell & Delmont, Tom. (2015). Anvi'o: An advanced analysis and visualization platformfor 'omics data. PeerJ. 3. e1319. 10.7717/peerj.1319.

#Also be sure to familiarize yourself with Anvi'o tutorials- Meren and co have some really useful ones.


#Notes for use: This is written in steps, hopefully to be streamlined at some point. Each step includes a header
#For the Sapelo2 cluster at UGA. If you are not using the cluster, do not include this code.
#If you are performing this on the cluster, all programs are installed. Otherwise, ensure you install all programs and dependencies.
#This is written for paired-end data.


#Disclaimer: No pipeline or script should be followed blindly. Ensure that this script and the programs used make sense for the samples you have.
#You may or may not utilize all steps of this pipeline- that depends entirely on your question. Understanding the pipeline and your question will save you time by not running steps you don't need.

#Note: This script was adapted from the Oil_Fluff.sh


##########################################################################################################################


#Preprocessing:
#Most of this protocol utilizes FASTQ files, but some need FASTA. If you only have FASTQ files,
#you can convert them to FASTA as follows when executed in the directory with the FASTQ files:

for i in $(ls *.fastq);
do
cat $i | paste - - - - | sed 's/^@/>/'| cut -f1-2 | tr '\t' '\n' > ${i%.fastq}.fa
done


##########################################################################################################################



#Step One- Quality filtering
#You first need to generate a TAB-delimited samples.txt file to point out where are your raw R1 and R2 files for each sample:
#This step requires a file called samples.txt which is outlined here: http://merenlab.org/tutorials/assembly-based-metagenomics/ in the quality filtering subheader.



iu-gen-configs samples.txt -o 01_QC
iu-filter-quality-minoche 01_QC/2013.ini




##########################################################################################################################

#Step Two- Assembly
#You can choose to assemble each metagenome individually.
#Then, you should reformat them to simplify the names for anvio


megahit -1 {R1_File} -2 {R2_File}--min-contig-len 1000 -m 0.8 --k-min 27 --presets meta-large -o 02_ASSEMBLY/ -t 8
mkdir 03_CONTIGS

anvi-script-reformat-fasta 02_ASSEMBLY/final.contigs.fa -o 03_CONTIGS/contigs.fa --min-len 2500 --simplify-names --report name_conversions.txt



##########################################################################################################################
