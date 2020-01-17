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


BASEDIR=/scratch/rck80079/Oil_Fluff_Data/Oil_Fluff
cd $BASEDIR


singularity exec /usr/local/singularity-images/anvio-6.1.simg iu-gen-configs samples.txt -o 01_QC
singularity exec /usr/local/singularity-images/anvio-6.1.simg iu-filter-quality-minoche 01_QC/2010.ini
singularity exec /usr/local/singularity-images/anvio-6.1.simg iu-filter-quality-minoche 01_QC/2010_2.ini
singularity exec /usr/local/singularity-images/anvio-6.1.simg iu-filter-quality-minoche 01_QC/2011.ini
singularity exec /usr/local/singularity-images/anvio-6.1.simg iu-filter-quality-minoche 01_QC/2012.ini
singularity exec /usr/local/singularity-images/anvio-6.1.simg iu-filter-quality-minoche 01_QC/2013.ini
singularity exec /usr/local/singularity-images/anvio-6.1.simg iu-filter-quality-minoche 01_QC/2014.ini
singularity exec /usr/local/singularity-images/anvio-6.1.simg iu-filter-quality-minoche 01_QC/2015.ini
singularity exec /usr/local/singularity-images/anvio-6.1.simg iu-filter-quality-minoche 01_QC/2016.ini
singularity exec /usr/local/singularity-images/anvio-6.1.simg iu-filter-quality-minoche 01_QC/2017.ini




#Step Three- No Assembly
#You can choose to do a read-based analysis depending on your question. This would be a great
#way to take a quick look at the data before doing an assembly-based analysis, and sometimes is sufficient
#for answering your scientific question. The read-based analysis includes steps three and four.
#you can run MetaPhlAn2 with paired end data, but run each fasta file separately. MetaPhlAn2 does not rely on mapping information.
#You can merge the outputs into one file if you do this.

#PBS -S /bin/bash
#PBS -q joye_q
#PBS -N Step_3
#PBS -l nodes=1:ppn=1
#PBS -l walltime=240:00:00
#PBS -l mem=100g
#PBS -j oe


BASEDIR=/scratch/rck80079/Oil_Fluff_Data/Oil_Fluff/01_QC/2010
cd $BASEDIR

module load MetaPhlAn2/2.7.8-foss-2016b-Python-2.7.14
metaphlan2.py 2010-QUALITY_PASSED_R1.fa  --input_type fasta > 2010-QUALITY_PASSED_R1_profile.fa
metaphlan2.py 2010-QUALITY_PASSED_R2.fa  --input_type fasta > 2010-QUALITY_PASSED_R2_profile.fa
cat 2010-QUALITY_PASSED_R1_profile.fa 2010-QUALITY_PASSED_R2_profile.fa > 2010_Metaphlan_profile.txt



#Take only unique taxa, combine duplicates
awk '!h[$4,$NF]++ { print $4, $NF }' FS="|" 2010_Metaphlan_profile.txt > 2010_temp.txt
awk '!h[$1,$NF]++ { print $1, $NF }' FS=" " 2010_temp.txt  > 2010_clean_profile.txt
ls -1 | grep "^o" 2010_clean_profile.txt > 2010_order_profile.txt
sed 's/ \+/,/g' 2010_order_profile.txt > 2010_order_profiles.csv
awk -F, '!seen[$1]++' 2010_order_profiles.csv > 2010_order_profile.csv
rm 2010_temp.txt 2010_clean_profile.txt  2010_order_profiles.csv


#PBS -S /bin/bash
#PBS -q highmem_q
#PBS -N Step_3_2010_2
#PBS -l nodes=1:ppn=1
#PBS -l walltime=240:00:00
#PBS -l mem=100g
#PBS -j oe


BASEDIR=/scratch/rck80079/Oil_Fluff_Data/Oil_Fluff/01_QC/2010_2
cd $BASEDIR

module load MetaPhlAn2/2.7.8-foss-2016b-Python-2.7.14
metaphlan2.py 2010_2-QUALITY_PASSED_R1.fa  --input_type fasta > 2010_2-QUALITY_PASSED_R1_profile.fa
metaphlan2.py 2010_2-QUALITY_PASSED_R2.fa  --input_type fasta > 2010_2-QUALITY_PASSED_R2_profile.fa
cat 2010_2-QUALITY_PASSED_R1_profile.fa 2010_2-QUALITY_PASSED_R2_profile.fa > 2010_2_Metaphlan_profile.txt



#Take only unique taxa, combine duplicates
awk '!h[$4,$NF]++ { print $4, $NF }' FS="|" 2010_2_Metaphlan_profile.txt > 2010_2_temp.txt
awk '!h[$1,$NF]++ { print $1, $NF }' FS=" " 2010_2_temp.txt  > 2010_2_clean_profile.txt
ls -1 | grep "^o" 2010_2_clean_profile.txt > 2010_2_order_profile.txt
sed 's/ \+/,/g' 2010_2_order_profile.txt > 2010_2_order_profiles.csv
awk -F, '!seen[$1]++' 2010_2_order_profiles.csv > 2010_2_order_profile.csv
rm 2010_2_temp.txt 2010_2_clean_profile.txt  2010_2_order_profiles.csv

#PBS -S /bin/bash
#PBS -q highmem_q
#PBS -N Step_3_2011
#PBS -l nodes=1:ppn=1
#PBS -l walltime=240:00:00
#PBS -l mem=100g
#PBS -j oe


BASEDIR=/scratch/rck80079/Oil_Fluff_Data/Oil_Fluff/01_QC/2011
cd $BASEDIR

module load MetaPhlAn2/2.7.8-foss-2016b-Python-2.7.14
metaphlan2.py 2011-QUALITY_PASSED_R1.fa  --input_type fasta > 2011-QUALITY_PASSED_R1_profile.fa
metaphlan2.py 2011-QUALITY_PASSED_R2.fa  --input_type fasta > 2011-QUALITY_PASSED_R2_profile.fa
cat 2011-QUALITY_PASSED_R1_profile.fa 2011-QUALITY_PASSED_R2_profile.fa > 2011_Metaphlan_profile.txt



#Take only unique taxa, combine duplicates
awk '!h[$4,$NF]++ { print $4, $NF }' FS="|" 2011_Metaphlan_profile.txt > 2011_temp.txt
awk '!h[$1,$NF]++ { print $1, $NF }' FS=" " 2011_temp.txt  > 2011_clean_profile.txt
ls -1 | grep "^o" 2011_clean_profile.txt > 2011_order_profile.txt
sed 's/ \+/,/g' 2011_order_profile.txt > 2011_order_profiles.csv
awk -F, '!seen[$1]++' 2011_order_profiles.csv > 2011_order_profile.csv
rm 2010_temp.txt 2011_clean_profile.txt  2011_order_profiles.csv

#PBS -S /bin/bash
#PBS -q highmem_q
#PBS -N Step_3_2012
#PBS -l nodes=1:ppn=1
#PBS -l walltime=240:00:00
#PBS -l mem=100g
#PBS -j oe


BASEDIR=/scratch/rck80079/Oil_Fluff_Data/Oil_Fluff/01_QC/2012
cd $BASEDIR

module load MetaPhlAn2/2.7.8-foss-2016b-Python-2.7.14
metaphlan2.py 2012-QUALITY_PASSED_R1.fa  --input_type fasta > 2012-QUALITY_PASSED_R1_profile.fa
metaphlan2.py 2012-QUALITY_PASSED_R2.fa  --input_type fasta > 2012-QUALITY_PASSED_R2_profile.fa
cat 2012-QUALITY_PASSED_R1_profile.fa 2012-QUALITY_PASSED_R2_profile.fa > 2012_Metaphlan_profile.txt



#Take only unique taxa, combine duplicates
awk '!h[$4,$NF]++ { print $4, $NF }' FS="|" 2012_Metaphlan_profile.txt > 2012_temp.txt
awk '!h[$1,$NF]++ { print $1, $NF }' FS=" " 2012_temp.txt  > 2012_clean_profile.txt
ls -1 | grep "^o" 2012_clean_profile.txt > 2012_order_profile.txt
sed 's/ \+/,/g' 2012_order_profile.txt > 2012_order_profiles.csv
awk -F, '!seen[$1]++' 2012_order_profiles.csv > 2012_order_profile.csv
rm 2012_temp.txt 2012_clean_profile.txt  2012_order_profiles.csv

#PBS -S /bin/bash
#PBS -q highmem_q
#PBS -N Step_3_2013
#PBS -l nodes=1:ppn=1
#PBS -l walltime=240:00:00
#PBS -l mem=100g
#PBS -j oe


BASEDIR=/scratch/rck80079/Oil_Fluff_Data/Oil_Fluff/01_QC/2013
cd $BASEDIR

module load MetaPhlAn2/2.7.8-foss-2016b-Python-2.7.14
metaphlan2.py 2013-QUALITY_PASSED_R1.fa  --input_type fasta > 2013-QUALITY_PASSED_R1_profile.fa
metaphlan2.py 2013-QUALITY_PASSED_R2.fa  --input_type fasta > 2013-QUALITY_PASSED_R2_profile.fa
cat 2013-QUALITY_PASSED_R1_profile.fa 2013-QUALITY_PASSED_R2_profile.fa > 2013_Metaphlan_profile.txt



#Take only unique taxa, combine duplicates
awk '!h[$4,$NF]++ { print $4, $NF }' FS="|" 2013_Metaphlan_profile.txt > 2013_temp.txt
awk '!h[$1,$NF]++ { print $1, $NF }' FS=" " 2013_temp.txt  > 2013_clean_profile.txt
ls -1 | grep "^o" 2013_clean_profile.txt > 2013_order_profile.txt
sed 's/ \+/,/g' 2013_order_profile.txt > 2013_order_profiles.csv
awk -F, '!seen[$1]++' 2013_order_profiles.csv > 2013_order_profile.csv
rm 2013_temp.txt 2013_clean_profile.txt  2013_order_profiles.csv

#PBS -S /bin/bash
#PBS -q highmem_q
#PBS -N Step_3_2014
#PBS -l nodes=1:ppn=1
#PBS -l walltime=240:00:00
#PBS -l mem=100g
#PBS -j oe


BASEDIR=/scratch/rck80079/Oil_Fluff_Data/Oil_Fluff/01_QC/2014
cd $BASEDIR

module load MetaPhlAn2/2.7.8-foss-2016b-Python-2.7.14
metaphlan2.py 2014-QUALITY_PASSED_R1.fa  --input_type fasta > 2014-QUALITY_PASSED_R1_profile.fa
metaphlan2.py 2014-QUALITY_PASSED_R2.fa  --input_type fasta > 2014-QUALITY_PASSED_R2_profile.fa
cat 2014-QUALITY_PASSED_R1_profile.fa 2014-QUALITY_PASSED_R2_profile.fa > 2014_Metaphlan_profile.txt



#Take only unique taxa, combine duplicates
awk '!h[$4,$NF]++ { print $4, $NF }' FS="|" 2014_Metaphlan_profile.txt > 2014_temp.txt
awk '!h[$1,$NF]++ { print $1, $NF }' FS=" " 2014_temp.txt  > 2014_clean_profile.txt
ls -1 | grep "^o" 2014_clean_profile.txt > 2014_order_profile.txt
sed 's/ \+/,/g' 2014_order_profile.txt > 2014_order_profiles.csv
awk -F, '!seen[$1]++' 2014_order_profiles.csv > 2014_order_profile.csv
rm 2014_temp.txt 2014_clean_profile.txt  2014_order_profiles.csv

#PBS -S /bin/bash
#PBS -q highmem_q
#PBS -N Step_3_2015
#PBS -l nodes=1:ppn=1
#PBS -l walltime=240:00:00
#PBS -l mem=100g
#PBS -j oe


BASEDIR=/scratch/rck80079/Oil_Fluff_Data/Oil_Fluff/01_QC/2015
cd $BASEDIR

module load MetaPhlAn2/2.7.8-foss-2016b-Python-2.7.14
metaphlan2.py 2015-QUALITY_PASSED_R1.fa  --input_type fasta > 2015-QUALITY_PASSED_R1_profile.fa
metaphlan2.py 2015-QUALITY_PASSED_R2.fa  --input_type fasta > 2015-QUALITY_PASSED_R2_profile.fa
cat 2015-QUALITY_PASSED_R1_profile.fa 2015-QUALITY_PASSED_R2_profile.fa > 2015_Metaphlan_profile.txt



#Take only unique taxa, combine duplicates
awk '!h[$4,$NF]++ { print $4, $NF }' FS="|" 2015_Metaphlan_profile.txt > 2015_temp.txt
awk '!h[$1,$NF]++ { print $1, $NF }' FS=" " 2015_temp.txt  > 2015_clean_profile.txt
ls -1 | grep "^o" 2015_clean_profile.txt > 2015_order_profile.txt
sed 's/ \+/,/g' 2015_order_profile.txt > 2015_order_profiles.csv
awk -F, '!seen[$1]++' 2015_order_profiles.csv > 2015_order_profile.csv
rm 2015_temp.txt 2015_clean_profile.txt  2015_order_profiles.csv

#PBS -S /bin/bash
#PBS -q highmem_q
#PBS -N Step_3_2016
#PBS -l nodes=1:ppn=1
#PBS -l walltime=240:00:00
#PBS -l mem=100g
#PBS -j oe


BASEDIR=/scratch/rck80079/Oil_Fluff_Data/Oil_Fluff/01_QC/2016
cd $BASEDIR

module load MetaPhlAn2/2.7.8-foss-2016b-Python-2.7.14
metaphlan2.py 2016-QUALITY_PASSED_R1.fa  --input_type fasta > 2016-QUALITY_PASSED_R1_profile.fa
metaphlan2.py 2016-QUALITY_PASSED_R2.fa  --input_type fasta > 2016-QUALITY_PASSED_R2_profile.fa
cat 2016-QUALITY_PASSED_R1_profile.fa 2016-QUALITY_PASSED_R2_profile.fa > 2016_Metaphlan_profile.txt



#Take only unique taxa, combine duplicates
awk '!h[$4,$NF]++ { print $4, $NF }' FS="|" 2016_Metaphlan_profile.txt > 2016_temp.txt
awk '!h[$1,$NF]++ { print $1, $NF }' FS=" " 2016_temp.txt  > 2016_clean_profile.txt
ls -1 | grep "^o" 2016_clean_profile.txt > 2016_order_profile.txt
sed 's/ \+/,/g' 2016_order_profile.txt > 2016_order_profiles.csv
awk -F, '!seen[$1]++' 2016_order_profiles.csv > 2016_order_profile.csv
rm 2016_temp.txt 2016_clean_profile.txt  2016_order_profiles.csv

#PBS -S /bin/bash
#PBS -q highmem_q
#PBS -N Step_3_2017
#PBS -l nodes=1:ppn=1
#PBS -l walltime=240:00:00
#PBS -l mem=100g
#PBS -j oe


BASEDIR=/scratch/rck80079/Oil_Fluff_Data/Oil_Fluff/01_QC/2017
cd $BASEDIR

module load MetaPhlAn2/2.7.8-foss-2016b-Python-2.7.14
metaphlan2.py 2017-QUALITY_PASSED_R1.fa  --input_type fasta > 2017-QUALITY_PASSED_R1_profile.fa
metaphlan2.py 2017-QUALITY_PASSED_R2.fa  --input_type fasta > 2017-QUALITY_PASSED_R2_profile.fa
cat 2017-QUALITY_PASSED_R1_profile.fa 2017-QUALITY_PASSED_R2_profile.fa > 2017_Metaphlan_profile.txt



#Take only unique taxa, combine duplicates
awk '!h[$4,$NF]++ { print $4, $NF }' FS="|" 2017_Metaphlan_profile.txt > 2017_temp.txt
awk '!h[$1,$NF]++ { print $1, $NF }' FS=" " 2017_temp.txt  > 2017_clean_profile.txt
ls -1 | grep "^o" 2017_clean_profile.txt > 2017_order_profile.txt
sed 's/ \+/,/g' 2017_order_profile.txt > 2017_order_profiles.csv
awk -F, '!seen[$1]++' 2017_order_profiles.csv > 2017_order_profile.csv
rm 2017_temp.txt 2017_clean_profile.txt  2017_order_profiles.csv


#This part of step 3 is to combine outputs for all samples for visualization in Rstudio interactive on the cluster.
#This output can be viewed on your local machine.

#PBS -S /bin/bash
#PBS -q joye_q
#PBS -N Step_3
#PBS -l nodes=1:ppn=1
#PBS -l walltime=20:00:00
#PBS -l mem=10g
#PBS -j oe
BASEDIR=/scratch/rck80079/Oil_Fluff_Data/Oil_Fluff/01_QC
cd $BASEDIR

cat /scratch/rck80079/Oil_Fluff_Data/Oil_Fluff/01_QC/2010_2/2010_2_order_profile.csv /scratch/rck80079/Oil_Fluff_Data/Oil_Fluff/01_QC/2010/2010_order_profile.csv /scratch/rck80079/Oil_Fluff_Data/Oil_Fluff/01_QC/2011/2011_order_profile.csv /scratch/rck80079/Oil_Fluff_Data/Oil_Fluff/01_QC/2012/2012_order_profile.csv /scratch/rck80079/Oil_Fluff_Data/Oil_Fluff/01_QC/2013/2013_order_profile.csv /scratch/rck80079/Oil_Fluff_Data/Oil_Fluff/01_QC/2014/2014_order_profile.csv /scratch/rck80079/Oil_Fluff_Data/Oil_Fluff/01_QC/2015/2015_order_profile.csv /scratch/rck80079/Oil_Fluff_Data/Oil_Fluff/01_QC/2016/2016_order_profile.csv /scratch/rck80079/Oil_Fluff_Data/Oil_Fluff/01_QC/2017/2017_order_profile.csv >combined_Oil_fluff_Metaphlan_profile.csv

qlogin

module load R/3.4.4-foss-2016b-X11-20160819-GACRC
R
# at this point, be sure to install the packages d3heatmap and htmlwidgets. Instructions for install are here: https://blog.rstudio.com/2015/06/24/d3heatmap/


#install.packages("d3heatmap")
#install.packages("htmlwidgets")

R --no-save < /scratch/rck80079/Oil_Fluff_Data/Oil_Fluff/01_QC/heatmap.r
scp *.html /scratch/rck80079/Oil_Fluff_Data/Oil_Fluff/01_QC



#Step Four- (co)Assembly
#You can choose to assemble each metagenome individually or to co-assemble them all together,
#and then map individual reads to that co-assembly. Read up on the theory before deciding which is right for your data.
#R1s.txt and R2s.txt are files which simply list the file location of the R1 and R2 reads you are co assembling.

#PBS _S /bin/bash
#PBS -N Step_5_co-assembly
#PBS -q highmem_q
#PBS -l nodes=1:ppn=32
#PBS -l walltime=200:00:00
#PBS -l mem=800gb

BASEDIR=/scratch/rck80079/Oil_Fluff_Data/Oil_Fluff
cd $BASEDIR

module load MEGAHIT/1.1.3-foss-2016b-Python-2.7.14

megahit -1 /scratch/rck80079/Oil_Fluff_Data/Oil_Fluff/01_QC/2010/2010-QUALITY_PASSED_R1.fa,/scratch/rck80079/Oil_Fluff_Data/Oil_Fluff/01_QC/2010_2/2010_2-QUALITY_PASSED_R1.fa,/scratch/rck80079/Oil_Fluff_Data/Oil_Fluff/01_QC/2011/2011-QUALITY_PASSED_R1.fa,/scratch/rck80079/Oil_Fluff_Data/Oil_Fluff/01_QC/2012/2012-QUALITY_PASSED_R1.fa,/scratch/rck80079/Oil_Fluff_Data/Oil_Fluff/01_QC/2013/2013-QUALITY_PASSED_R1.fa,/scratch/rck80079/Oil_Fluff_Data/Oil_Fluff/01_QC/2014/2014-QUALITY_PASSED_R1.fa,/scratch/rck80079/Oil_Fluff_Data/Oil_Fluff/01_QC/2015/2015-QUALITY_PASSED_R1.fa,/scratch/rck80079/Oil_Fluff_Data/Oil_Fluff/01_QC/2016/2016-QUALITY_PASSED_R1.fa,/scratch/rck80079/Oil_Fluff_Data/Oil_Fluff/01_QC/2017/2017-QUALITY_PASSED_R1.fa -2 /scratch/rck80079/Oil_Fluff_Data/Oil_Fluff/01_QC/2010/2010-QUALITY_PASSED_R2.fa,/scratch/rck80079/Oil_Fluff_Data/Oil_Fluff/01_QC/2010_2/2010_2-QUALITY_PASSED_R2.fa,/scratch/rck80079/Oil_Fluff_Data/Oil_Fluff/01_QC/2011/2011-QUALITY_PASSED_R2.fa,/scratch/rck80079/Oil_Fluff_Data/Oil_Fluff/01_QC/2012/2012-QUALITY_PASSED_R2.fa,/scratch/rck80079/Oil_Fluff_Data/Oil_Fluff/01_QC/2013/2013-QUALITY_PASSED_R2.fa,/scratch/rck80079/Oil_Fluff_Data/Oil_Fluff/01_QC/2014/2014-QUALITY_PASSED_R2.fa,/scratch/rck80079/Oil_Fluff_Data/Oil_Fluff/01_QC/2015/2015-QUALITY_PASSED_R2.fa,/scratch/rck80079/Oil_Fluff_Data/Oil_Fluff/01_QC/2016/2016-QUALITY_PASSED_R2.fa,/scratch/rck80079/Oil_Fluff_Data/Oil_Fluff/01_QC/2017/2017-QUALITY_PASSED_R2.fa --min-contig-len 1000 -m 0.8 --k-min 27 --presets meta-large -o 02_ASSEMBLY/ -t 32
mkdir 03_CONTIGS

singularity exec /usr/local/singularity-images/anvio-5.4.simg anvi-script-reformat-fasta 02_ASSEMBLY/final.contigs.fa -o 03_CONTIGS/contigs.fa --min-len 2500 --simplify-names --report name_conversions.txt








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










#Step Six- Recovering Metagenome Assembled Genomes
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
#PBS -q highmem_q
#PBS -N kaiju
#PBS -l nodes=1:ppn=10
#PBS -l walltime=200:00:00
#PBS -l mem=400gb

BASEDIR=/scratch/rck80079/Oil_Fluff_Data/Oil_Fluff/
cd $BASEDIR

module load kaiju/1.6.2

#Only make this DB once, you can reuse it.
mkdir kaijudb
cd kaijudb
kaiju-makedb -t 10 -s refseq

kaiju -t kaijudbnodes.dmp \
      -f kaijudb/kaiju_db.fmi \
      -i gene_calls.fa \
      -o gene_calls_nr.out \
      -z 16 \
      -v \
      1>job.out 2>job.err

singularity exec /usr/local/singularity-images/anvio-6.1.simg addTaxonNames -t kaijudb/nodes.dmp \
      -n kaijudb/names.dmp \
      -i gene_calls_nr.out \
      -o gene_calls_nr.names \
      -r superkingdom,phylum,order,class,family,genus,species

singularity exec /usr/local/singularity-images/anvio-6.1.simg anvi-import-taxonomy-for-genes -i gene_calls_nr.names \
      -c contigs.db \
      --just-do-it \
      -p kaiju






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
