#File Created: 1/12/2020 by Rachael Storo
#Last Edited: 1/16/2020 by Rachael Storo

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

#Note: This script was adapted from the metagenomic_pipeline.sh and this is the actual code run in the oil fluff analyses





#Step one- Preprocessing:
#Most of this protocol utilizes FASTQ files, but some need FASTA. If you only have FASTQ files, you can convert them to FASTA as follows when executed in the directory with FASTQ files:
for i in $(ls *.fastq);
do
cat $i | paste - - - - | sed 's/^@/>/'| cut -f1-2 | tr '\t' '\n' > ${i%.fastq}.fa
done


#Step Two- Simplify Names

#PBS -S /bin/bash
#PBS -N simplify
#PBS -q joye_q
#PBS -l nodes=1:ppn=4
#PBS -l walltime=20:00:00
#PBS -l mem=200gb


BASEDIR=/scratch/rck80079/Oil_Fluff_Data/Oil_Fluff
cd $BASEDIR

singularity exec /usr/local/singularity-images/anvio-6.1.simg anvi-script-reformat-fasta JW860_ONT_EN527_44.fa \
                           --simplify-names \
                           -o JW860_ONT_EN527_44_fixy.fa

mv JW860_ONT_EN527_44_fixy.fa JW860_ONT_EN527_44.fa


#Step three- build database

#PBS -S /bin/bash
#PBS -N db
#PBS -q joye_q
#PBS -l nodes=1:ppn=10
#PBS -l walltime=100:00:00
#PBS -l mem=200gb


BASEDIR=/scratch/rck80079/Oil_Fluff_Data/Oil_Fluff
cd $BASEDIR

singularity exec /usr/local/singularity-images/anvio-6.1.simg anvi-gen-contigs-database -f JW860_ONT_EN527_44.fa -o JW860_ONT_EN527_44.db
singularity exec /usr/local/singularity-images/anvio-6.1.simg anvi-run-hmms -c JW860_ONT_EN527_44.db --num-threads 10
singularity exec /usr/local/singularity-images/anvio-6.1.simg anvi-display-contigs-stats JW860_ONT_EN527_44.db >minion_stats.txt



#Step 4- Get out Ribosomal rRNA seqs
#You can take and blast the output fasta file from this step against the nr database to see what you get.

#PBS -S /bin/bash
#PBS -N rRNA
#PBS -q joye_q
#PBS -l nodes=1:ppn=4
#PBS -l walltime=20:00:00
#PBS -l mem=200gb


BASEDIR=/scratch/rck80079/Oil_Fluff_Data/Oil_Fluff
cd $BASEDIR

singularity exec /usr/local/singularity-images/anvio-6.1.simg anvi-get-sequences-for-hmm-hits -c JW860_ONT_EN527_44.db \
                                --hmm-source Ribosomal_RNAs \
                                --gene-name Bacterial_16S_rRNA,Bacterial_23S_rRNA,Archaeal_16S_rRNA,Archaeal_23S_rRNA \
                                -o JW860_ONT_EN527_44_rRNA.fa
