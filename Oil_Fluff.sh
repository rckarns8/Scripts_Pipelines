#Oil Fluff Metagenomics Pipeline
#Initiated: October 30, 2019 by Rachael Storo
#Last Edit: December 13, 2019 by Rachael Storo

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





#Preprocessing:
#Most of this protocol utilizes FASTQ files, but some need FASTA. If you only have FASTQ files, you can convert them to FASTA as follows when executed in the directory with FASTQ files:
for i in $(ls *.fastq);
do
cat $i | paste - - - - | sed 's/^@/>/'| cut -f1-2 | tr '\t' '\n' > ${i%.fastq}.fa
done





#Step Two- Quality filtering and trimming
#Quality filtering will remove adapters and primers to give you clean fasta files.
#You should change the file TruSeq3-PE.fa to match whatever primers you are using.
#This step requires a file called samples.txt which is outlined here: http://merenlab.org/tutorials/assembly-based-metagenomics/ in the quality filtering subheader.

#PBS -S /bin/bash
#PBS -N Quality_filtering
#PBS -q joye_q
#PBS -l nodes=1:ppn=4
#PBS -l walltime=20:00:00
#PBS -l mem=200gb


BASEDIR=/scratch/rck80079/OC26_Leftover
cd $BASEDIR


singularity exec /usr/local/singularity-images/anvio-6.1.simg iu-gen-configs samples.txt -o 01_QC
singularity exec /usr/local/singularity-images/anvio-6.1.simg iu-filter-quality-minoche 01_QC/EN600_74.ini
singularity exec /usr/local/singularity-images/anvio-6.1.simg iu-filter-quality-minoche 01_QC/AT26-13_43.ini
singularity exec /usr/local/singularity-images/anvio-6.1.simg iu-filter-quality-minoche 01_QC/EN510_39.ini
singularity exec /usr/local/singularity-images/anvio-6.1.simg iu-filter-quality-minoche 01_QC/AT18-02_MUC19-8_0-3.ini
singularity exec /usr/local/singularity-images/anvio-6.1.simg iu-filter-quality-minoche 01_QC/EN586_65.ini
singularity exec /usr/local/singularity-images/anvio-6.1.simg iu-filter-quality-minoche 01_QC/EN559_168.ini
singularity exec /usr/local/singularity-images/anvio-6.1.simg iu-filter-quality-minoche 01_QC/AT26-13_43_JW826.ini
singularity exec /usr/local/singularity-images/anvio-6.1.simg iu-filter-quality-minoche 01_QC/EN600_79.ini
singularity exec /usr/local/singularity-images/anvio-6.1.simg iu-filter-quality-minoche 01_QC/EN600_79.ini
singularity exec /usr/local/singularity-images/anvio-6.1.simg iu-filter-quality-minoche 01_QC/EN559_174.ini
singularity exec /usr/local/singularity-images/anvio-6.1.simg iu-filter-quality-minoche 01_QC/EN586_70.ini
singularity exec /usr/local/singularity-images/anvio-6.1.simg iu-filter-quality-minoche 01_QC/AT26-13_48.ini
singularity exec /usr/local/singularity-images/anvio-6.1.simg iu-filter-quality-minoche 01_QC/EN527_50.ini
singularity exec /usr/local/singularity-images/anvio-6.1.simg iu-filter-quality-minoche 01_QC/EN527_45.ini
singularity exec /usr/local/singularity-images/anvio-6.1.simg iu-filter-quality-minoche 01_QC/OC468-2_C86_STE25.ini
singularity exec /usr/local/singularity-images/anvio-6.1.simg iu-filter-quality-minoche 01_QC/OC468-2_C83_STE25.ini
singularity exec /usr/local/singularity-images/anvio-6.1.simg iu-filter-quality-minoche 01_QC/OC468-2_C82_STE25.ini
singularity exec /usr/local/singularity-images/anvio-6.1.simg iu-filter-quality-minoche 01_QC/PE1031_96.ini
singularity exec /usr/local/singularity-images/anvio-6.1.simg iu-filter-quality-minoche 01_QC/OC468-2_C93.ini
singularity exec /usr/local/singularity-images/anvio-6.1.simg iu-filter-quality-minoche 01_QC/PE1031_54.ini
singularity exec /usr/local/singularity-images/anvio-6.1.simg iu-filter-quality-minoche 01_QC/OC468-2_C89.ini
singularity exec /usr/local/singularity-images/anvio-6.1.simg iu-filter-quality-minoche 01_QC/EN510_34.ini
singularity exec /usr/local/singularity-images/anvio-6.1.simg iu-filter-quality-minoche 01_QC/EN496_122.ini
singularity exec /usr/local/singularity-images/anvio-6.1.simg iu-filter-quality-minoche 01_QC/AT18-02_MUC19-8_16-21.ini









#Step Three- No Assembly
#You can choose to do a read-based analysis depending on your question. This would be a great
#way to take a quick look at the data before doing an assembly-based analysis, and sometimes is sufficient
#for answering your scientific question. The read-based analysis includes steps three and four.
#you can run MetaPhlAn2 with paired end data, but run each fasta file separately.



#Step Four- (co)Assembly
#You can choose to assemble each metagenome individually or to co-assemble them all together,
#and then map individual reads to that co-assembly. Read up on the theory before deciding which is right for your data.
#R1s.txt and R2s.txt are files which simply list the file location of the R1 and R2 reads you are co assembling.


#PBS _S /bin/bash
#PBS -N 2010 assembly
#PBS -q highmem_q
#PBS -l nodes=1:ppn=4
#PBS -l walltime=50:00:00
#PBS -l mem=800gb

BASEDIR=/scratch/rck80079/Oil_Fluff_Data/Oil_Fluff
cd $BASEDIR

module load MEGAHIT/1.1.3-foss-2016b-Python-2.7.14

megahit -1 /scratch/rck80079/Oil_Fluff_Data/Oil_Fluff/01_QC/2010/2010-QUALITY_PASSED_R1.fa -2 /scratch/rck80079/Oil_Fluff_Data/Oil_Fluff/01_QC/2010/2010-QUALITY_PASSED_R2.fa --min-contig-len 1000 -m 0.8 --k-min 27 --presets meta-large -o 02_ASSEMBLY/2010/ -t 4

singularity exec /usr/local/singularity-images/anvio-5.4.simg anvi-script-reformat-fasta 02_ASSEMBLY/2010/final.contigs.fa -o 03_CONTIGS/PE1031_Sed_53_OC26.fa --min-len 2500 --simplify-names --report name_conversions.txt




#PBS _S /bin/bash
#PBS -N 2010-2 assembly
#PBS -q highmem_q
#PBS -l nodes=1:ppn=4
#PBS -l walltime=50:00:00
#PBS -l mem=800gb

BASEDIR=/scratch/rck80079/Oil_Fluff_Data/Oil_Fluff
cd $BASEDIR

module load MEGAHIT/1.1.3-foss-2016b-Python-2.7.14

megahit -1 /scratch/rck80079/Oil_Fluff_Data/Oil_Fluff/01_QC/2010_2/2010_2-QUALITY_PASSED_R1.fa -2 /scratch/rck80079/Oil_Fluff_Data/Oil_Fluff/01_QC/2010_2/2010_2-QUALITY_PASSED_R2.fa --min-contig-len 1000 -m 0.8 --k-min 27 --presets meta-large -o 02_ASSEMBLY/2010_2 -t 4

singularity exec /usr/local/singularity-images/anvio-5.4.simg anvi-script-reformat-fasta 02_ASSEMBLY/2010_2/final.contigs.fa -o 03_CONTIGS/PE1031_Sed_95_OC26.fa --min-len 2500 --simplify-names --report name_conversions.txt


#PBS _S /bin/bash
#PBS -N 2011 assembly
#PBS -q highmem_q
#PBS -l nodes=1:ppn=4
#PBS -l walltime=50:00:00
#PBS -l mem=800gb

BASEDIR=/scratch/rck80079/Oil_Fluff_Data/Oil_Fluff
cd $BASEDIR

module load MEGAHIT/1.1.3-foss-2016b-Python-2.7.14

megahit -1 /scratch/rck80079/Oil_Fluff_Data/Oil_Fluff/01_QC/2011/2011-QUALITY_PASSED_R1.fa -2 /scratch/rck80079/Oil_Fluff_Data/Oil_Fluff/01_QC/2011/2011-QUALITY_PASSED_R2.fa --min-contig-len 1000 -m 0.8 --k-min 27 --presets meta-large -o 02_ASSEMBLY/2011 -t 4


singularity exec /usr/local/singularity-images/anvio-5.4.simg anvi-script-reformat-fasta 02_ASSEMBLY/2011/final.contigs.fa -o 03_CONTIGS/EN496_Sed_117_OC26.fa --min-len 2500 --simplify-names --report name_conversions.txt



#PBS _S /bin/bash
#PBS -N 2012 assembly
#PBS -q highmem_q
#PBS -l nodes=1:ppn=4
#PBS -l walltime=50:00:00
#PBS -l mem=800gb

BASEDIR=/scratch/rck80079/Oil_Fluff_Data/Oil_Fluff
cd $BASEDIR

module load MEGAHIT/1.1.3-foss-2016b-Python-2.7.14

megahit -1 /scratch/rck80079/Oil_Fluff_Data/Oil_Fluff/01_QC/2012/2012-QUALITY_PASSED_R1.fa -2 /scratch/rck80079/Oil_Fluff_Data/Oil_Fluff/01_QC/2012/2012-QUALITY_PASSED_R2.fa --min-contig-len 1000 -m 0.8 --k-min 27 --presets meta-large -o 02_ASSEMBLY/2012/ -t 4

singularity exec /usr/local/singularity-images/anvio-5.4.simg anvi-script-reformat-fasta 02_ASSEMBLY/2012/final.contigs.fa -o 03_CONTIGS/EN510_Sed_33_OC26.fa --min-len 2500 --simplify-names --report name_conversions.txt

#PBS _S /bin/bash
#PBS -N 2013 assembly
#PBS -q highmem_q
#PBS -l nodes=1:ppn=4
#PBS -l walltime=50:00:00
#PBS -l mem=800gb

BASEDIR=/scratch/rck80079/Oil_Fluff_Data/Oil_Fluff
cd $BASEDIR

module load MEGAHIT/1.1.3-foss-2016b-Python-2.7.14

megahit -1 /scratch/rck80079/Oil_Fluff_Data/Oil_Fluff/01_QC/2013/2013-QUALITY_PASSED_R1.fa -2 /scratch/rck80079/Oil_Fluff_Data/Oil_Fluff/01_QC/2013/2013-QUALITY_PASSED_R2.fa --min-contig-len 1000 -m 0.8 --k-min 27 --presets meta-large -o 02_ASSEMBLY/2013/ -t 4


singularity exec /usr/local/singularity-images/anvio-5.4.simg anvi-script-reformat-fasta 02_ASSEMBLY/2013/final.contigs.fa -o 03_CONTIGS/EN527_Sed_44_OC26.fa --min-len 2500 --simplify-names --report name_conversions.txt

#PBS _S /bin/bash
#PBS -N 2014 assembly
#PBS -q highmem_q
#PBS -l nodes=1:ppn=4
#PBS -l walltime=50:00:00
#PBS -l mem=800gb

BASEDIR=/scratch/rck80079/Oil_Fluff_Data/Oil_Fluff
cd $BASEDIR

module load MEGAHIT/1.1.3-foss-2016b-Python-2.7.14

megahit -1 /scratch/rck80079/Oil_Fluff_Data/Oil_Fluff/01_QC/2014/2014-QUALITY_PASSED_R1.fa -2 /scratch/rck80079/Oil_Fluff_Data/Oil_Fluff/01_QC/2014/2014-QUALITY_PASSED_R2.fa --min-contig-len 1000 -m 0.8 --k-min 27 --presets meta-large -o 02_ASSEMBLY/2014/ -t 4

singularity exec /usr/local/singularity-images/anvio-5.4.simg anvi-script-reformat-fasta 02_ASSEMBLY/2014/final.contigs.fa -o 03_CONTIGS/AT26-13_Sed_42_OC26.fa --min-len 2500 --simplify-names --report name_conversions.txt


#PBS _S /bin/bash
#PBS -N 2015 assembly
#PBS -q highmem_q
#PBS -l nodes=1:ppn=4
#PBS -l walltime=50:00:00
#PBS -l mem=800gb

BASEDIR=/scratch/rck80079/Oil_Fluff_Data/Oil_Fluff
cd $BASEDIR

module load MEGAHIT/1.1.3-foss-2016b-Python-2.7.14

megahit -1 /scratch/rck80079/Oil_Fluff_Data/Oil_Fluff/01_QC/2015/2015-QUALITY_PASSED_R1.fa -2 /scratch/rck80079/Oil_Fluff_Data/Oil_Fluff/01_QC/2015/2015-QUALITY_PASSED_R2.fa --min-contig-len 1000 -m 0.8 --k-min 27 --presets meta-large -o 02_ASSEMBLY/2015/ -t 4

singularity exec /usr/local/singularity-images/anvio-5.4.simg anvi-script-reformat-fasta 02_ASSEMBLY/2015/final.contigs.fa -o 03_CONTIGS/EN559_Sed_167_OC26.fa --min-len 2500 --simplify-names --report name_conversions.txt



#PBS _S /bin/bash
#PBS -N 2016 assembly
#PBS -q highmem_q
#PBS -l nodes=1:ppn=4
#PBS -l walltime=50:00:00
#PBS -l mem=800gb

BASEDIR=/scratch/rck80079/Oil_Fluff_Data/Oil_Fluff
cd $BASEDIR

module load MEGAHIT/1.1.3-foss-2016b-Python-2.7.14

megahit -1 /scratch/rck80079/Oil_Fluff_Data/Oil_Fluff/01_QC/2016/2016-QUALITY_PASSED_R1.fa -2 /scratch/rck80079/Oil_Fluff_Data/Oil_Fluff/01_QC/2016/2016-QUALITY_PASSED_R2.fa --min-contig-len 1000 -m 0.8 --k-min 27 --presets meta-large -o 02_ASSEMBLY/2016/ -t 4

singularity exec /usr/local/singularity-images/anvio-5.4.simg anvi-script-reformat-fasta 02_ASSEMBLY/2016/final.contigs.fa -o 03_CONTIGS/EN586_Sed_64_OC26.fa --min-len 2500 --simplify-names --report name_conversions.txt


#PBS _S /bin/bash
#PBS -N 2017 assembly
#PBS -q highmem_q
#PBS -l nodes=1:ppn=4
#PBS -l walltime=50:00:00
#PBS -l mem=800gb

BASEDIR=/scratch/rck80079/Oil_Fluff_Data/Oil_Fluff
cd $BASEDIR

module load MEGAHIT/1.1.3-foss-2016b-Python-2.7.14

megahit -1 /scratch/rck80079/Oil_Fluff_Data/Oil_Fluff/01_QC/2017/2017-QUALITY_PASSED_R1.fa -2 /scratch/rck80079/Oil_Fluff_Data/Oil_Fluff/01_QC/2017/2017-QUALITY_PASSED_R2.fa --min-contig-len 1000 -m 0.8 --k-min 27 --presets meta-large -o 02_ASSEMBLY/2017/ -t 4

singularity exec /usr/local/singularity-images/anvio-5.4.simg anvi-script-reformat-fasta 02_ASSEMBLY/2017/final.contigs.fa -o 03_CONTIGS/EN600_Sed_73_OC26.fa --min-len 2500 --simplify-names --report name_conversions.txt

#################################################################################
#Other SAMPLES


#PBS _S /bin/bash
#PBS -N EN600_74
#PBS -q highmem_q
#PBS -l nodes=1:ppn=4
#PBS -l walltime=50:00:00
#PBS -l mem=800gb

BASEDIR=/scratch/rck80079/OC26_Leftover
cd $BASEDIR

module load MEGAHIT/1.1.3-foss-2016b-Python-2.7.14

megahit -1 /scratch/rck80079/OC26_Leftover/01_QC/EN600_74-QUALITY_PASSED_R1.fa -2 /scratch/rck80079/OC26_Leftover/01_QC/EN600_74-QUALITY_PASSED_R2.fa --min-contig-len 1000 -m 0.8 --k-min 27 --presets meta-large -o 02_ASSEMBLY/EN600_74/ -t 4

singularity exec /usr/local/singularity-images/anvio-5.4.simg anvi-script-reformat-fasta 02_ASSEMBLY/EN600_74/final.contigs.fa -o 03_CONTIGS/EN600_74.fa --min-len 2500 --simplify-names --report name_conversions.txt



#PBS _S /bin/bash
#PBS -N AT26-13_43
#PBS -q highmem_q
#PBS -l nodes=1:ppn=4
#PBS -l walltime=50:00:00
#PBS -l mem=800gb

BASEDIR=/scratch/rck80079/OC26_Leftover
cd $BASEDIR

module load MEGAHIT/1.1.3-foss-2016b-Python-2.7.14

megahit -1 /scratch/rck80079/OC26_Leftover/01_QC/AT26-13_43-QUALITY_PASSED_R1.fa -2 /scratch/rck80079/OC26_Leftover/01_QC/AT26-13_43-QUALITY_PASSED_R2.fa --min-contig-len 1000 -m 0.8 --k-min 27 --presets meta-large -o 02_ASSEMBLY/AT26-13_43/ -t 4

singularity exec /usr/local/singularity-images/anvio-5.4.simg anvi-script-reformat-fasta 02_ASSEMBLY/AT26-13_43/final.contigs.fa -o 03_CONTIGS/AT26-13_43.fa --min-len 2500 --simplify-names --report name_conversions.txt


#PBS _S /bin/bash
#PBS -N EN510_39
#PBS -q highmem_q
#PBS -l nodes=1:ppn=4
#PBS -l walltime=50:00:00
#PBS -l mem=800gb

BASEDIR=/scratch/rck80079/OC26_Leftover
cd $BASEDIR

module load MEGAHIT/1.1.3-foss-2016b-Python-2.7.14

megahit -1 /scratch/rck80079/OC26_Leftover/01_QC/EN510_39/EN510_39-QUALITY_PASSED_R1.fa -2 /scratch/rck80079/OC26_Leftover/01_QC/EN510_39/EN510_39-QUALITY_PASSED_R2.fa --min-contig-len 1000 -m 0.8 --k-min 27 --presets meta-large -o 02_ASSEMBLY/EN510_39/ -t 4

singularity exec /usr/local/singularity-images/anvio-5.4.simg anvi-script-reformat-fasta 02_ASSEMBLY/EN510_39/final.contigs.fa -o 03_CONTIGS/EN510_39.fa --min-len 2500 --simplify-names --report name_conversions.txt



#PBS _S /bin/bash
#PBS -N AT18-02_MUC19-8_0-3
#PBS -q highmem_q
#PBS -l nodes=1:ppn=4
#PBS -l walltime=50:00:00
#PBS -l mem=800gb

BASEDIR=/scratch/rck80079/OC26_Leftover
cd $BASEDIR

module load MEGAHIT/1.1.3-foss-2016b-Python-2.7.14

megahit -1 /scratch/rck80079/OC26_Leftover/01_QC/AT18-02_MUC19-8_0-3/AT18-02_MUC19-8_0-3-QUALITY_PASSED_R1.fa -2 /scratch/rck80079/OC26_Leftover/01_QC/AT18-02_MUC19-8_0-3/AT18-02_MUC19-8_0-3-QUALITY_PASSED_R2.fa --min-contig-len 1000 -m 0.8 --k-min 27 --presets meta-large -o 02_ASSEMBLY/AT18-02_MUC19-8_0-3/ -t 4

singularity exec /usr/local/singularity-images/anvio-5.4.simg anvi-script-reformat-fasta 02_ASSEMBLY/AT18-02_MUC19-8_0-3/final.contigs.fa -o 03_CONTIGS/AT18-02_MUC19-8_0-3.fa --min-len 2500 --simplify-names --report name_conversions.txt


#PBS _S /bin/bash
#PBS -N EN586_65
#PBS -q highmem_q
#PBS -l nodes=1:ppn=4
#PBS -l walltime=50:00:00
#PBS -l mem=800gb

BASEDIR=/scratch/rck80079/OC26_Leftover
cd $BASEDIR

module load MEGAHIT/1.1.3-foss-2016b-Python-2.7.14

megahit -1 /scratch/rck80079/OC26_Leftover/01_QC/EN586_65/EN586_65-QUALITY_PASSED_R1.fa -2 /scratch/rck80079/OC26_Leftover/01_QC/EN586_65/EN586_65-QUALITY_PASSED_R2.fa --min-contig-len 1000 -m 0.8 --k-min 27 --presets meta-large -o 02_ASSEMBLY/EN586_65/ -t 4

singularity exec /usr/local/singularity-images/anvio-5.4.simg anvi-script-reformat-fasta 02_ASSEMBLY/EN586_65/final.contigs.fa -o 03_CONTIGS/EN586_65.fa --min-len 2500 --simplify-names --report name_conversions.txt



#PBS _S /bin/bash
#PBS -N EN559_168
#PBS -q highmem_q
#PBS -l nodes=1:ppn=4
#PBS -l walltime=50:00:00
#PBS -l mem=800gb

BASEDIR=/scratch/rck80079/OC26_Leftover
cd $BASEDIR

module load MEGAHIT/1.1.3-foss-2016b-Python-2.7.14

megahit -1 /scratch/rck80079/OC26_Leftover/01_QC/EN559_168/EN559_168-QUALITY_PASSED_R1.fa -2 /scratch/rck80079/OC26_Leftover/01_QC/EN559_168/EN559_168-QUALITY_PASSED_R2.fa --min-contig-len 1000 -m 0.8 --k-min 27 --presets meta-large -o 02_ASSEMBLY/EN559_168/ -t 4

singularity exec /usr/local/singularity-images/anvio-5.4.simg anvi-script-reformat-fasta 02_ASSEMBLY/EN559_168/final.contigs.fa -o 03_CONTIGS/EN559_168.fa --min-len 2500 --simplify-names --report name_conversions.txt


#PBS _S /bin/bash
#PBS -N AT26-13_43_JW826
#PBS -q highmem_q
#PBS -l nodes=1:ppn=4
#PBS -l walltime=50:00:00
#PBS -l mem=800gb

BASEDIR=/scratch/rck80079/OC26_Leftover
cd $BASEDIR

module load MEGAHIT/1.1.3-foss-2016b-Python-2.7.14

megahit -1 /scratch/rck80079/OC26_Leftover/01_QC/AT26-13_43_JW826/AT26-13_43_JW826-QUALITY_PASSED_R1.fa -2 /scratch/rck80079/OC26_Leftover/01_QC/AT26-13_43_JW826/AT26-13_43_JW826-QUALITY_PASSED_R2.fa --min-contig-len 1000 -m 0.8 --k-min 27 --presets meta-large -o 02_ASSEMBLY/AT26-13_43_JW826/ -t 4

singularity exec /usr/local/singularity-images/anvio-5.4.simg anvi-script-reformat-fasta 02_ASSEMBLY/AT26-13_43_JW826/final.contigs.fa -o 03_CONTIGS/AT26-13_43_JW826.fa --min-len 2500 --simplify-names --report name_conversions.txt


#PBS _S /bin/bash
#PBS -N EN600_79
#PBS -q highmem_q
#PBS -l nodes=1:ppn=4
#PBS -l walltime=50:00:00
#PBS -l mem=800gb

BASEDIR=/scratch/rck80079/OC26_Leftover
cd $BASEDIR

module load MEGAHIT/1.1.3-foss-2016b-Python-2.7.14

megahit -1 /scratch/rck80079/OC26_Leftover/01_QC/EN600_79/EN600_79-QUALITY_PASSED_R1.fa -2 /scratch/rck80079/OC26_Leftover/01_QC/EN600_79/EN600_79-QUALITY_PASSED_R2.fa --min-contig-len 1000 -m 0.8 --k-min 27 --presets meta-large -o 02_ASSEMBLY/EN600_79/ -t 4

singularity exec /usr/local/singularity-images/anvio-5.4.simg anvi-script-reformat-fasta 02_ASSEMBLY/EN600_79/final.contigs.fa -o 03_CONTIGS/EN600_79.fa --min-len 2500 --simplify-names --report name_conversions.txt


#PBS _S /bin/bash
#PBS -N EN559_174
#PBS -q highmem_q
#PBS -l nodes=1:ppn=4
#PBS -l walltime=50:00:00
#PBS -l mem=800gb

BASEDIR=/scratch/rck80079/OC26_Leftover
cd $BASEDIR

module load MEGAHIT/1.1.3-foss-2016b-Python-2.7.14

megahit -1 /scratch/rck80079/OC26_Leftover/01_QC/EN559_174/EN559_174-QUALITY_PASSED_R1.fa -2 /scratch/rck80079/OC26_Leftover/01_QC/EN559_174/EN559_174-QUALITY_PASSED_R2.fa --min-contig-len 1000 -m 0.8 --k-min 27 --presets meta-large -o 02_ASSEMBLY/EN559_174/ -t 4

singularity exec /usr/local/singularity-images/anvio-5.4.simg anvi-script-reformat-fasta 02_ASSEMBLY/EN559_174/final.contigs.fa -o 03_CONTIGS/EN559_174.fa --min-len 2500 --simplify-names --report name_conversions.txt


#PBS _S /bin/bash
#PBS -N EN586_70
#PBS -q highmem_q
#PBS -l nodes=1:ppn=4
#PBS -l walltime=50:00:00
#PBS -l mem=800gb

BASEDIR=/scratch/rck80079/OC26_Leftover
cd $BASEDIR

module load MEGAHIT/1.1.3-foss-2016b-Python-2.7.14

megahit -1 /scratch/rck80079/OC26_Leftover/01_QC/EN586_70/EN586_70-QUALITY_PASSED_R1.fa -2 /scratch/rck80079/OC26_Leftover/01_QC/EN586_70/EN586_70-QUALITY_PASSED_R2.fa --min-contig-len 1000 -m 0.8 --k-min 27 --presets meta-large -o 02_ASSEMBLY/EN586_70/ -t 4

singularity exec /usr/local/singularity-images/anvio-5.4.simg anvi-script-reformat-fasta 02_ASSEMBLY/EN586_70/final.contigs.fa -o 03_CONTIGS/EN586_70.fa --min-len 2500 --simplify-names --report name_conversions.txt


#PBS _S /bin/bash
#PBS -N AT26-13_48
#PBS -q highmem_q
#PBS -l nodes=1:ppn=4
#PBS -l walltime=50:00:00
#PBS -l mem=800gb

BASEDIR=/scratch/rck80079/OC26_Leftover
cd $BASEDIR

module load MEGAHIT/1.1.3-foss-2016b-Python-2.7.14

megahit -1 /scratch/rck80079/OC26_Leftover/01_QC/AT26-13_48/AT26-13_48-QUALITY_PASSED_R1.fa -2 /scratch/rck80079/OC26_Leftover/01_QC/AT26-13_48/AT26-13_48-QUALITY_PASSED_R2.fa --min-contig-len 1000 -m 0.8 --k-min 27 --presets meta-large -o 02_ASSEMBLY/AT26-13_48/ -t 4

singularity exec /usr/local/singularity-images/anvio-5.4.simg anvi-script-reformat-fasta 02_ASSEMBLY/AT26-13_48/final.contigs.fa -o 03_CONTIGS/AT26-13_48.fa --min-len 2500 --simplify-names --report name_conversions.txt


#PBS _S /bin/bash
#PBS -N EN527_50
#PBS -q highmem_q
#PBS -l nodes=1:ppn=4
#PBS -l walltime=50:00:00
#PBS -l mem=800gb

BASEDIR=/scratch/rck80079/OC26_Leftover
cd $BASEDIR

module load MEGAHIT/1.1.3-foss-2016b-Python-2.7.14

megahit -1 /scratch/rck80079/OC26_Leftover/01_QC/EN527_50/EN527_50-QUALITY_PASSED_R1.fa -2 /scratch/rck80079/OC26_Leftover/01_QC/EN527_50/EN527_50-QUALITY_PASSED_R2.fa --min-contig-len 1000 -m 0.8 --k-min 27 --presets meta-large -o 02_ASSEMBLY/EN527_50/ -t 4

singularity exec /usr/local/singularity-images/anvio-5.4.simg anvi-script-reformat-fasta 02_ASSEMBLY/EN527_50/final.contigs.fa -o 03_CONTIGS/EN527_50.fa --min-len 2500 --simplify-names --report name_conversions.txt


#PBS _S /bin/bash
#PBS -N EN527_45
#PBS -q highmem_q
#PBS -l nodes=1:ppn=4
#PBS -l walltime=50:00:00
#PBS -l mem=800gb

BASEDIR=/scratch/rck80079/OC26_Leftover
cd $BASEDIR

module load MEGAHIT/1.1.3-foss-2016b-Python-2.7.14

megahit -1 /scratch/rck80079/OC26_Leftover/01_QC/EN527_45/EN527_45-QUALITY_PASSED_R1.fa -2 /scratch/rck80079/OC26_Leftover/01_QC/EN527_45/EN527_45-QUALITY_PASSED_R2.fa --min-contig-len 1000 -m 0.8 --k-min 27 --presets meta-large -o 02_ASSEMBLY/EN527_45/ -t 4

singularity exec /usr/local/singularity-images/anvio-5.4.simg anvi-script-reformat-fasta 02_ASSEMBLY/EN527_45/final.contigs.fa -o 03_CONTIGS/EN527_45.fa --min-len 2500 --simplify-names --report name_conversions.txt


#PBS _S /bin/bash
#PBS -N OC468-2_C86_STE25
#PBS -q highmem_q
#PBS -l nodes=1:ppn=4
#PBS -l walltime=50:00:00
#PBS -l mem=800gb

BASEDIR=/scratch/rck80079/OC26_Leftover
cd $BASEDIR

module load MEGAHIT/1.1.3-foss-2016b-Python-2.7.14

megahit -1 /scratch/rck80079/OC26_Leftover/01_QC/OC468-2_C86_STE25/OC468-2_C86_STE25-QUALITY_PASSED_R1.fa -2 /scratch/rck80079/OC26_Leftover/01_QC/OC468-2_C86_STE25/OC468-2_C86_STE25-QUALITY_PASSED_R2.fa --min-contig-len 1000 -m 0.8 --k-min 27 --presets meta-large -o 02_ASSEMBLY/OC468-2_C86_STE25/ -t 4

singularity exec /usr/local/singularity-images/anvio-5.4.simg anvi-script-reformat-fasta 02_ASSEMBLY/OC468-2_C86_STE25/final.contigs.fa -o 03_CONTIGS/OC468-2_C86_STE25.fa --min-len 2500 --simplify-names --report name_conversions.txt


#PBS _S /bin/bash
#PBS -N OC468-2_C83_STE25
#PBS -q highmem_q
#PBS -l nodes=1:ppn=4
#PBS -l walltime=50:00:00
#PBS -l mem=800gb

BASEDIR=/scratch/rck80079/OC26_Leftover
cd $BASEDIR

module load MEGAHIT/1.1.3-foss-2016b-Python-2.7.14

megahit -1 /scratch/rck80079/OC26_Leftover/01_QC/OC468-2_C83_STE25/OC468-2_C83_STE25-QUALITY_PASSED_R1.fa -2 /scratch/rck80079/OC26_Leftover/01_QC/OC468-2_C83_STE25/OC468-2_C83_STE25-QUALITY_PASSED_R2.fa --min-contig-len 1000 -m 0.8 --k-min 27 --presets meta-large -o 02_ASSEMBLY/OC468-2_C83_STE25/ -t 4

singularity exec /usr/local/singularity-images/anvio-5.4.simg anvi-script-reformat-fasta 02_ASSEMBLY/OC468-2_C83_STE25/final.contigs.fa -o 03_CONTIGS/OC468-2_C83_STE25.fa --min-len 2500 --simplify-names --report name_conversions.txt


#PBS _S /bin/bash
#PBS -N OC468-2_C82_STE25
#PBS -q highmem_q
#PBS -l nodes=1:ppn=4
#PBS -l walltime=50:00:00
#PBS -l mem=800gb

BASEDIR=/scratch/rck80079/OC26_Leftover
cd $BASEDIR

module load MEGAHIT/1.1.3-foss-2016b-Python-2.7.14

megahit -1 /scratch/rck80079/OC26_Leftover/01_QC/OC468-2_C82_STE25/OC468-2_C82_STE25-QUALITY_PASSED_R1.fa -2 /scratch/rck80079/OC26_Leftover/01_QC/OC468-2_C82_STE25/OC468-2_C82_STE25-QUALITY_PASSED_R2.fa --min-contig-len 1000 -m 0.8 --k-min 27 --presets meta-large -o 02_ASSEMBLY/OC468-2_C82_STE25/ -t 4

singularity exec /usr/local/singularity-images/anvio-5.4.simg anvi-script-reformat-fasta 02_ASSEMBLY/OC468-2_C82_STE25/final.contigs.fa -o 03_CONTIGS/OC468-2_C82_STE25.fa --min-len 2500 --simplify-names --report name_conversions.txt


#PBS _S /bin/bash
#PBS -N PE1031_96
#PBS -q highmem_q
#PBS -l nodes=1:ppn=4
#PBS -l walltime=50:00:00
#PBS -l mem=800gb

BASEDIR=/scratch/rck80079/OC26_Leftover
cd $BASEDIR

module load MEGAHIT/1.1.3-foss-2016b-Python-2.7.14

megahit -1 /scratch/rck80079/OC26_Leftover/01_QC/PE1031_96/PE1031_96-QUALITY_PASSED_R1.fa -2 /scratch/rck80079/OC26_Leftover/01_QC/PE1031_96/PE1031_96-QUALITY_PASSED_R2.fa --min-contig-len 1000 -m 0.8 --k-min 27 --presets meta-large -o 02_ASSEMBLY/PE1031_96/ -t 4

singularity exec /usr/local/singularity-images/anvio-5.4.simg anvi-script-reformat-fasta 02_ASSEMBLY/PE1031_96/final.contigs.fa -o 03_CONTIGS/PE1031_96.fa --min-len 2500 --simplify-names --report name_conversions.txt


#PBS _S /bin/bash
#PBS -N OC468-2_C93
#PBS -q highmem_q
#PBS -l nodes=1:ppn=4
#PBS -l walltime=50:00:00
#PBS -l mem=800gb

BASEDIR=/scratch/rck80079/OC26_Leftover
cd $BASEDIR

module load MEGAHIT/1.1.3-foss-2016b-Python-2.7.14

megahit -1 /scratch/rck80079/OC26_Leftover/01_QC/OC468-2_C93/OC468-2_C93-QUALITY_PASSED_R1.fa -2 /scratch/rck80079/OC26_Leftover/01_QC/OC468-2_C93/OC468-2_C93-QUALITY_PASSED_R2.fa --min-contig-len 1000 -m 0.8 --k-min 27 --presets meta-large -o 02_ASSEMBLY/OC468-2_C93/ -t 4

singularity exec /usr/local/singularity-images/anvio-5.4.simg anvi-script-reformat-fasta 02_ASSEMBLY/OC468-2_C93/final.contigs.fa -o 03_CONTIGS/OC468-2_C93.fa --min-len 2500 --simplify-names --report name_conversions.txt


#PBS _S /bin/bash
#PBS -N PE1031_54
#PBS -q highmem_q
#PBS -l nodes=1:ppn=4
#PBS -l walltime=50:00:00
#PBS -l mem=800gb

BASEDIR=/scratch/rck80079/OC26_Leftover
cd $BASEDIR

module load MEGAHIT/1.1.3-foss-2016b-Python-2.7.14

megahit -1 /scratch/rck80079/OC26_Leftover/01_QC/PE1031_54/PE1031_54-QUALITY_PASSED_R1.fa -2 /scratch/rck80079/OC26_Leftover/01_QC/PE1031_54/PE1031_54-QUALITY_PASSED_R2.fa --min-contig-len 1000 -m 0.8 --k-min 27 --presets meta-large -o 02_ASSEMBLY/PE1031_54/ -t 4

singularity exec /usr/local/singularity-images/anvio-5.4.simg anvi-script-reformat-fasta 02_ASSEMBLY/PE1031_54/final.contigs.fa -o 03_CONTIGS/PE1031_54.fa --min-len 2500 --simplify-names --report name_conversions.txt


#PBS _S /bin/bash
#PBS -N OC468-2_C89
#PBS -q highmem_q
#PBS -l nodes=1:ppn=4
#PBS -l walltime=50:00:00
#PBS -l mem=800gb

BASEDIR=/scratch/rck80079/OC26_Leftover
cd $BASEDIR

module load MEGAHIT/1.1.3-foss-2016b-Python-2.7.14

megahit -1 /scratch/rck80079/OC26_Leftover/01_QC/OC468-2_C89/OC468-2_C89-QUALITY_PASSED_R1.fa -2 /scratch/rck80079/OC26_Leftover/01_QC/OC468-2_C89/OC468-2_C89-QUALITY_PASSED_R2.fa --min-contig-len 1000 -m 0.8 --k-min 27 --presets meta-large -o 02_ASSEMBLY/OC468-2_C89/ -t 4

singularity exec /usr/local/singularity-images/anvio-5.4.simg anvi-script-reformat-fasta 02_ASSEMBLY/OC468-2_C89/final.contigs.fa -o 03_CONTIGS/OC468-2_C89.fa --min-len 2500 --simplify-names --report name_conversions.txt



#PBS _S /bin/bash
#PBS -N EN510_34
#PBS -q highmem_q
#PBS -l nodes=1:ppn=4
#PBS -l walltime=50:00:00
#PBS -l mem=800gb

BASEDIR=/scratch/rck80079/OC26_Leftover
cd $BASEDIR

module load MEGAHIT/1.1.3-foss-2016b-Python-2.7.14

megahit -1 /scratch/rck80079/OC26_Leftover/01_QC/EN510_34/EN510_34-QUALITY_PASSED_R1.fa -2 /scratch/rck80079/OC26_Leftover/01_QC/EN510_34/EN510_34-QUALITY_PASSED_R2.fa --min-contig-len 1000 -m 0.8 --k-min 27 --presets meta-large -o 02_ASSEMBLY/EN510_34/ -t 4

singularity exec /usr/local/singularity-images/anvio-5.4.simg anvi-script-reformat-fasta 02_ASSEMBLY/EN510_34/final.contigs.fa -o 03_CONTIGS/EN510_34.fa --min-len 2500 --simplify-names --report name_conversions.txt



#PBS _S /bin/bash
#PBS -N EN496_122
#PBS -q highmem_q
#PBS -l nodes=1:ppn=4
#PBS -l walltime=50:00:00
#PBS -l mem=800gb

BASEDIR=/scratch/rck80079/OC26_Leftover
cd $BASEDIR

module load MEGAHIT/1.1.3-foss-2016b-Python-2.7.14

megahit -1 /scratch/rck80079/OC26_Leftover/01_QC/EN496_122/EN496_122-QUALITY_PASSED_R1.fa -2 /scratch/rck80079/OC26_Leftover/01_QC/EN496_122/EN496_122-QUALITY_PASSED_R2.fa --min-contig-len 1000 -m 0.8 --k-min 27 --presets meta-large -o 02_ASSEMBLY/EN496_122/ -t 4

singularity exec /usr/local/singularity-images/anvio-5.4.simg anvi-script-reformat-fasta 02_ASSEMBLY/EN496_122/final.contigs.fa -o 03_CONTIGS/EN496_122.fa --min-len 2500 --simplify-names --report name_conversions.txt


#PBS _S /bin/bash
#PBS -N AT18-02_MUC19-8_16-21
#PBS -q highmem_q
#PBS -l nodes=1:ppn=4
#PBS -l walltime=50:00:00
#PBS -l mem=800gb

BASEDIR=/scratch/rck80079/OC26_Leftover
cd $BASEDIR

module load MEGAHIT/1.1.3-foss-2016b-Python-2.7.14

megahit -1 /scratch/rck80079/OC26_Leftover/01_QC/AT18-02_MUC19-8_16-21/AT18-02_MUC19-8_16-21-QUALITY_PASSED_R1.fa -2 /scratch/rck80079/OC26_Leftover/01_QC/AT18-02_MUC19-8_16-21/AT18-02_MUC19-8_16-21-QUALITY_PASSED_R2.fa --min-contig-len 1000 -m 0.8 --k-min 27 --presets meta-large -o 02_ASSEMBLY/AT18-02_MUC19-8_16-21/ -t 4

singularity exec /usr/local/singularity-images/anvio-5.4.simg anvi-script-reformat-fasta 02_ASSEMBLY/AT18-02_MUC19-8_16-21/final.contigs.fa -o 03_CONTIGS/AT18-02_MUC19-8_16-21.fa --min-len 2500 --simplify-names --report name_conversions.txt


#Step Five- Mapping
#Map the sample reads to the co-assembly to generate covereage information.

#PBS -S /bin/bash
#PBS -q joye_q
#PBS -N Mapping
#PBS -l nodes=1:ppn=4
#PBS -l walltime=20:00:00
#PBS -l mem=200gb

BASEDIR=/scratch/rck80079/Oil_Fluff_Data/Oil_Fluff/
cd $BASEDIR
mkdir 04_MAPPING
time bowtie2-build 03_CONTIGS/contigs.fa 04_MAPPING/contigs
module load Bowtie2/2.3.4.1-foss-2016b
module load SAMtools/1.6-foss-2016b

time bowtie2-build 03_CONTIGS/contigs.fa 04_MAPPING/contigs

# A simple loop to serially map all samples.
# referenced from within http://merenlab.org/tutorials/assembly_and_mapping/

# how many threads should each mapping task use?
NUM_THREADS=4

for sample in `awk '{print $1}' samples.txt`
do
    if [ "$sample" == "sample" ]; then continue; fi

    # you need to make sure you "ls 01_QC/*QUALITY_PASSED_R1*" returns R1 files for all your samples in samples.txt
    R1s=`ls 01_QC/01_QC/$sample*QUALITY_PASSED_R1* | python -c 'import sys; print ",".join([x.strip() for x in sys.stdin.readlines()])'`
    R2s=`ls 01_QC/01_QC/$sample*QUALITY_PASSED_R2* | python -c 'import sys; print ",".join([x.strip() for x in sys.stdin.readlines()])'`

    time bowtie2 --threads $NUM_THREADS -x 04_MAPPING/contigs -1 $R1s -2 $R2s --no-unal -S 04_MAPPING/$sample.sam
    time samtools view -F 4 -bS 04_MAPPING/$sample.sam > 04_MAPPING/$sample-RAW.bam
    singularity exec /usr/local/singularity-images/anvio-6.1.simg anvi-init-bam 04_MAPPING/$sample-RAW.bam -o 04_MAPPING/$sample.bam
    rm 04_MAPPING/$sample.sam 04_MAPPING/$sample-RAW.bam
done

# mapping is done, and we no longer need bowtie2-build files
rm 04_MAPPING/*.bt2








#Step Six- Assign taxonomy

#PBS -S /bin/bash
#PBS -q highmem_q
#PBS -N kraken2
#PBS -l nodes=1:ppn=20
#PBS -l walltime=10:00:00
#PBS -l mem=300gb
BASEDIR=/scratch/rck80079/Oil_Fluff_Data/Oil_Fluff/01_QC/2010

cd /scratch/rck80079/Oil_Fluff_Data/Oil_Fluff/01_QC/
module load Kraken2/2.0.7-beta-foss-2018a-Perl-5.26.1

#You only need to build this once, and you can use it over and over.
kraken2-build --standard --threads 20 --db . --download-taxonomy wget ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz
wget ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_wgs.accession2taxid.gz
touch accmap.dlflag

mv nucl_gb.accession2taxid.gz /scratch/rck80079/Oil_Fluff_Data/Oil_Fluff/01_QC/kraken_database/taxonomy
mv nucl_wgs.accession2taxid.gz /scratch/rck80079/Oil_Fluff_Data/Oil_Fluff/01_QC/kraken_database/taxonomy
cd $BASEDIR
kraken2 --db kraken_database --quick --threads 20 --paired 2010-QUALITY_PASSED_R1.fa 2010-QUALITY_PASSED_R2.fa




#PBS -S /bin/bash
#PBS -q highmem_q
#PBS -N assign taxa
#PBS -l nodes=1:ppn=20
#PBS -l walltime=10:00:00
#PBS -l mem=300gb
BASEDIR=/scratch/rck80079/Oil_Fluff_Data/Oil_Fluff/01_QC/2010
cd $BASEDIR


module load Bracken/2.2-foss-2016b-Python-2.7.14
bracken -d kraken_database -i




#Step Seven- Make Anvio Database for visualization
#A note on MAGs- the accepted quality of MAGs has been put forward in the Woyke et al 2018 paper
#High-quality MAGs have greater than 90% completion, and <5% contamination. These are not to be
#considered the same as isolate genomes, but are a pretty reasonable representative genome of closely
#related genomic lineages.

#PBS -S /bin/bash
#PBS -q joye_q
#PBS -N makeDB
#PBS -l nodes=1:ppn=20
#PBS -l walltime=200:00:00
#PBS -l mem=250gb

BASEDIR=/scratch/rck80079/Oil_Fluff_Data/Oil_Fluff/
cd $BASEDIR



singularity exec /usr/local/singularity-images/anvio-6.1.simg anvi-script-reformat-fasta 03_CONTIGS/contigs.fa -o contigs-fixed.fa -l 0 --simplify-names
singularity exec /usr/local/singularity-images/anvio-6.1.simg anvi-gen-contigs-database -f 03_CONTIGS/contigs.fa -o contigs.db -n 'Oil Fluff contigs datbase'
singularity exec /usr/local/singularity-images/anvio-6.1.simg anvi-run-hmms -c contigs.db
singularity exec /usr/local/singularity-images/anvio-6.1.simg anvi-display-contigs-stats contigs.db
singularity exec /usr/local/singularity-images/anvio-6.1.simg anvi-run-ncbi-cogs -c CONTIGS.db --num-threads 20
singularity exec /usr/local/singularity-images/anvio-6.1.simg anvi-get-sequences-for-gene-calls -c CONTIGS.db -o gene_calls.fa
singularity exec /usr/local/singularity-images/anvio-6.1.simg makeDB.sh -e -t 20



#PBS -S /bin/bash
#PBS -q joye_q
#PBS -N profile
#PBS -l nodes=1:ppn=1
#PBS -l walltime=50:00:00
#PBS -l mem=50gb

BASEDIR=/scratch/rck80079/Oil_Fluff_Data/Oil_Fluff/
cd $BASEDIR


singularity exec /usr/local/singularity-images/anvio-6.1.simg anvi-profile -i 04_MAPPING/2010.bam -c contigs.db
singularity exec /usr/local/singularity-images/anvio-6.1.simg anvi-merge 2010/PROFILE.db 2010_2/PROFILE.db 2011/PROFILE.db 2012/PROFILE.db 2013/PROFILE.db 2014/PROFILE.db 2015/PROFILE.db 2016/PROFILE.db 2017/PROFILE.db -o OIL-SAMPLES-MERGED -c contigs.db
singularity exec /usr/local/singularity-images/anvio-6.1.simganvi-interactive -p OIL-SAMPLES-MERGED/PROFILE.db -c contigs.db
singularity exec /usr/local/singularity-images/anvio-6.1.simganvi-summarize -p OIL-SAMPLES-MERGED/PROFILE.db -c contigs.db -o SAMPLES-SUMMARY -C CONCOCT

#singularity exec /usr/local/singularity-images/anvio-6.1.simg anvi-refine -p OIL-SAMPLES-MERGED/PROFILE.db -c contigs.db -C CONCOCT -b Group_6




#Step Seven- Gene Calling
#You can do this step after Recovering MAGs from your metagenomes or on each metagenome.

#PBS -S /bin/bash
#PBS -N Step_4
#PBS -q joye_q
#PBS -l nodes=1:ppn=1
#PBS -l walltime=100:00:00
#PBS -l mem=50gb



BASEDIR=/work/sbjlab/rck/Sed_Assemblies/AT26-13-87/
cd $BASEDIR

ml prokka/1.13-foss-2016b-BioPerl-1.7.1
time prokka contigs.fa

BASEDIR=/work/sbjlab/rck/Sed_Assemblies/AT26-13-89/
cd $BASEDIR

ml prokka/1.13-foss-2016b-BioPerl-1.7.1
time prokka contigs.fa

BASEDIR=/work/sbjlab/rck/Sed_Assemblies/AT26-13-91/
cd $BASEDIR

ml prokka/1.13-foss-2016b-BioPerl-1.7.1
time prokka contigs.fa




#Step Eight- Phylogenomics
#This step is written to be done on a merged database and profile of all metagenomes.

#PBS -S /bin/bash
#PBS -q joye_q
#PBS -N jobname
#PBS -l nodes=1:ppn=1
#PBS -l walltime=100:00:00
#PBS -l mem=50gb

cd $PBS_O_WORKDIR

singularity exec /usr/local/singularity-images/anvio-5.4.simg anvi-get-sequences-for-hmm-hits -c CONTIGS.db \
      -p PROFILE.db \
      -o seqs-for-phylogenomics.fa \
      --list-hmm-sources

singularity exec /usr/local/singularity-images/anvio-5.4.simg anvi-get-sequences-for-hmm-hits -c CONTIGS.db \
      -p PROFILE.db \
      -o seqs-for-phylogenomics.fa \
      --hmm-source Bacteria_71 \
      -C default \
      --gene-names Ribosomal_L1,Ribosomal_L2,Ribosomal_L3,Ribosomal_L4,Ribosomal_L5,Ribosomal_L6
singularity exec /usr/local/singularity-images/anvio-5.4.simg anvi-get-sequences-for-hmm-hits -c CONTIGS.db \
      -p PROFILE.db \
      -o seqs-for-phylogenomics.fa \
      --hmm-source Bacteria_71 \
      -C default \
      --gene-names Ribosomal_L1,Ribosomal_L2,Ribosomal_L3,Ribosomal_L4,Ribosomal_L5,Ribosomal_L6 \
      --concatenate-genes \
      --return-best-hit \
      --get-aa-sequences

singularity exec /usr/local/singularity-images/anvio-5.4.simg anvi-gen-phylogenomic-tree -f seqs-for-phylogenomics.fa \
      -o phylogenomic-tree.txt

anvi-interactive --tree phylogenomic-tree.txt \
      -p temp-profile.db \
      --title "Pylogenomics of Bins" \
      --manual

anvi-interactive -p PROFILE.db \
      -c CONTIGS.db \
      -C default \
      --tree phylogenomic-tree.txt








#Step Nine- Comparative (Pan) Genomics

#PBS -S /bin/bash
#PBS -q joye_q
#PBS -N jobname
#PBS -l nodes=1:ppn=1
#PBS -l walltime=100:00:00
#PBS -l mem=50gb

cd $PBS_O_WORKDIR
singularity exec /usr/local/singularity-images/anvio-5.4.simg anvi-import-collection additional-files/collections/sample.txt \
    --bins-info additional-files/collections/info.txt \
    -p PROFILE.db \
    -c CONTIGS.db \
    -C AT26-13-Co-Assembly

singularity exec /usr/local/singularity-images/anvio-5.4.simg anvi-gen-genomes-storage -i additional-files/pangenomics/internal-genomes.txt \
    -e additional-files/pangenomics/external-genomes.txt \
    -o AT26-13-Co-Assembly-GENOMES.db

anvi-pan-genome -g AT26-13-Co-Assembly-GENOMES.db \
    -n AT26-13-Co-Assembly \
    -o PAN \
    --num-threads 10

anvi-display-pan -g AT26-13-Co-Assembly-GENOMES.db \
    -p PAN/AT26-13-Co-Assembly-PAN.db \
    --title "AT26-13-Co-Assembly Pan"


anvi-import-state -p PAN/AT26-13-Co-Assembly-PAN.db \
    --state additional-files/state-files/state-pan.json \
    --name default

anvi-display-pan -g AT26-13-Co-Assembly-GENOMES.db \
    -p PAN/AT26-13-Co-Assembly-PAN.db \
    --title "AT26-13-Co-Assembly Pan"


anvi-compute-genome-similarity -e additional-files/pangenomics/external-genomes.txt \
    -i additional-files/pangenomics/internal-genomes.txt \
    --program pyANI \
    -o ANI \
    -T 6 \
    --pan-db PAN/AT26-13-Co-Assembly-PAN.db

anvi-display-pan -g AT26-13-Co-Assembly-GENOMES.db \
    -p PAN/AT26-13-Co-Assembly-PAN.db \
    --title "AT26-13-Co-Assembly Pan"

anvi-import-misc-data -p PAN/AT26-13-Co-Assembly-PAN.db \
   --target-data-table layers \
   additional-files/pangenomics/additional-layers-data.txt


anvi-display-pan -g AT26-13-Co-Assembly-GENOMES.db \
   -p PAN/AT26-13-Co-Assembly-PAN.db \
   --title "AT26-13-Co-Assembly Pan"

anvi-import-collection additional-files/pangenomics/pan-collection.txt \
    --bins-info additional-files/pangenomics/pan-collection-info.txt \
    -p PAN/AT26-13-Co-Assembly-PAN.db \
    -C default

anvi-display-pan -g Enterococcus-GENOMES.db \
    -p PAN/AT26-13-Co-Assembly-PAN.db \
    --title "AT26-13-Co-Assembly Pan"

anvi-summarize -p PAN/AT26-13-Co-Assembly-PAN.db \
    -g AT26-13-Co-Assembly-GENOMES.db \
    -C default \
    -o PAN_SUMMARY

open PAN_SUMMARY/index.html
gzip -d PAN_SUMMARY/AT26-13-Co-Assembly_protein_clusters_summary.txt.gz








#Step Ten- SNV analysis
#Profiling SNVs allows for you to examine the microbial population genetics in your metagenomes

#PBS -S /bin/bash
#PBS -q joye_q
#PBS -N jobname
#PBS -l nodes=1:ppn=1
#PBS -l walltime=100:00:00
#PBS -l mem=50gb

cd $PBS_O_WORKDIR
singularity exec /usr/local/singularity-images/anvio-5.4.simg anvi-import-collection additional-files/collections/merens.txt \
      --bins-info additional-files/collections/merens-info.txt \
      -p PROFILE.db \
      -c CONTIGS.db \
      -C default

# importing taxonomy for gene calls
singularity exec /usr/local/singularity-images/anvio-5.4.simg anvi-import-taxonomy-for-genes -c CONTIGS.db \
      -i additional-files/centrifuge-files/centrifuge_report.tsv \
      additional-files/centrifuge-files/centrifuge_hits.tsv \
      -p centrifuge

# importing the state file so things look pretty
singularity exec /usr/local/singularity-images/anvio-5.4.simg anvi-import-state --state additional-files/state-files/state-merged.json \
      --name default \
      -p PROFILE.db

singularity exec /usr/local/singularity-images/anvio-5.4.simg anvi-split -p PROFILE.db \
      -c CONTIGS.db \
      -C default \
      -b AT26-13-Co-Assembly \
      -o MAGs

singularity exec /usr/local/singularity-images/anvio-5.4.simg anvi-estimate-genome-completeness -p MAGs/E_facealis/PROFILE.db \
      -c MAGs/AT26-13-Co-Assembly/CONTIGS.db \
      -C DEFAULT

anvi-interactive -p MAGs/AT26-13-Co-Assembly/PROFILE.db \
      -c MAGs/AT26-13-Co-Assembly/CONTIGS.db

singularity exec /usr/local/singularity-images/anvio-5.4.simg anvi-gen-variability-profile -c MAGs/AT26-13-Co-Assembly/CONTIGS.db \
      -p MAGs/AT26-13-Co-Assembly/PROFILE.db \
      -C DEFAULT \
      -b ALL_SPLITS \
      --samples-of-interest additional-files/samples.txt \
      --min-coverage-in-each-sample 20 \
      --min-occurrence 3 \
      --include-split-names \
      --quince-mode \
      -o AT26-13-Co-Assembly-SNVs.txt

cat additional-files/samples.txt
singularity exec /usr/local/singularity-images/anvio-5.4.simg anvi-script-snvs-to-interactive AT26-13-Co-Assembly-SNVs.txt -o AT26-13-Co-Assembly_snvs
anvi-interactive --profile AT26-13-Co-Assembly_snvs/profile.db \
                 --tree AT26-13-Co-Assembly_snvs/tree.txt \
                 --view-data AT26-13-Co-Assembly_snvs/view.txt \
                 --title "SNV Profile for the AT26-13-Co-Assembly bin" \
                 --manual

singularity exec /usr/local/singularity-images/anvio-5.4.simg anvi-import-state -p AT26-13-Co-Assembly_snvs/profile.db \
                 --state additional-files/state-files/state-snvs.json \
                 --name default

anvi-interactive -d AT26-13-Co-Assembly_snvs/view.txt \
                 -t AT26-13-Co-Assembly_snvs/tree.txt \
                 -p AT26-13-Co-Assembly_snvs/profile.db \
                 --title "SNV Profile for the AT26-13-Co-Assembly bin" \
                 --manual

singularity exec /usr/local/singularity-images/anvio-5.4.simg anvi-gen-fixation-index-matrix --variability-profile AT26-13-Co-Assembly-SNVs.txt \
               --output-file FST_AT26-13-Co-Assembly.txt

singularity exec /usr/local/singularity-images/anvio-5.4.simg anvi-matrix-to-newick FST_AT26-13-Co-Assembly.txt \
               --output-file FST_AT26-13-Co-Assembly.newick

anvi-interactive -t FST_AT26-13-Co-Assembly.newick \
                 -p FST_AT26-13-Co-Assembly.db \
                 --manual
