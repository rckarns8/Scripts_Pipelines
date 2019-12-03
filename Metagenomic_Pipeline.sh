#Metagenomics pipeline
#Initiated: October 30, 2019 by Rachael Storo
#Last Edit: December 2, 2019 by Rachael Storo

#Purpose: Create a script for the processing of environmental metagenomic samples

#Read this paper before use of this pipeline, it will provide a visual context to what the majority of these steps are doing;
#Eren, A. & Esen, Özcan & Quince, Christopher & Vineis, Joseph & Morrison, Hilary & Sogin, Mitchell & Delmont, Tom. (2015). Anvi'o: An advanced analysis and visualization platformfor 'omics data. PeerJ. 3. e1319. 10.7717/peerj.1319.

#Also be sure to familiarize yourself with Anvi'o tutorials- Meren and co have some really useful ones.


#Notes for use: This is written in steps, hopefully to be streamlined at some point. Each step includes a header
#For the Sapelo2 cluster at UGA. If you are not using the cluster, do not include this code.
#If you are performing this on the cluster, all programs are installed. Otherwise, ensure you install all programs and dependencies.
#This is written for paired-end data.


#Disclaimer: No pipeline or script should be followed blindly. Ensure that this script and the programs used make sense for the samples you have.
#You may or may not utilize all steps of this pipeline- that depends entirely on your question. Understanding the pipeline and your question will save you time by not running steps you don't need.

#I will test this pipeline with some metagenomic data that my lab has on hand, which are called
#AT26-13_87, AT26-13_89, and AT26-13_91. These are samples from a sediment core, at different sediment depths.
#from a site in the Gulf of Mexico called GC600, a natural oil seep site, sampled on 06/04/2014. Sample 87 is the closest to the
#sediment/water interface at 0-3 cm sediment depth, sample 89 is 6-9 cm sediment depth, and 91 is 12-15 cm sediment depth.



#NOTE: anywhere with '<text>' requires you to change the file names. Be sure to remove <> before running code.





#Preprocessing:
#Most of this protocol utilizes FASTQ files, but some need FASTA. If you only have FASTQ files, you can convert them to FASTA as follows when executed in the directory with FASTQ files:
for i in $(ls *.fastq);
do
cat $i | paste - - - - | sed 's/^@/>/'| cut -f1-2 | tr '\t' '\n' > ${i%.fastq}.fa
done







#Step One- Demultiplex
# Not all data will require this step- in fact, many metagenomes are sequenced individually or are given to you already demultiplexed.

#PBS -S /bin/bash
#PBS -q batch
#PBS -N Step_1
#PBS -l nodes=1:ppn=1:AMD
#PBS -l walltime=20:00:00
#PBS -l mem=50gb

BASEDIR=</work/sbjlab/rck/Sed_Assemblies/AT26-13-87/>
cd $BASEDIR

module load USEARCH/10.0.240-i86linux32
usearch -fastx_demux <R1.fq> -reverse <R2.fq> -index <I1.fq> -barcodes <bar.fa> \
        -fastqout fwd_demux.fq -output2 rev_demux.fq







#Step Two- Quality trimming and filtering
#Quality filtering will remove adapters and primers to give you clean fasta files.
#You should change the file TruSeq3-PE.fa to match whatever primers you are using.
#This step requires a file called samples.txt which is outlined here: http://merenlab.org/tutorials/assembly-based-metagenomics/ in the quality filtering subheader.

#PBS -S /bin/bash
#PBS -N Step_2
#PBS -q batch
#PBS -l nodes=1:ppn=4:AMD
#PBS -l walltime=20:00:00
#PBS -l mem=50gb


BASEDIR=</work/sbjlab/rck/Sed_Assemblies/AT26-13-87/>
cd $BASEDIR

module load Trimmomatic/0.36-Java-1.8.0_144

time java -jar /usr/local/apps/eb/Trimmomatic/0.36-Java-1.8.0_144/trimmomatic-0.36.jar PE -threads 4 <AT26-13_87_Sed_S49_L003_R1_001.fastq> <AT26-13_87_Sed_S49_L003_R2_001.fastq> <AT26-13_87_forward_paired.fq.gz> <AT26-13_87_forward_unpaired.fq.gz> <AT26-13_87_reverse_paired.fq.gz> <AT26-13_87_reverse_unpaired.fq.gz> ILLUMINACLIP:<TruSeq3-PE.fa>:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:36

singularity exec /usr/local/singularity-images/anvio-5.4.simg iu-gen-configs <samples.txt> -o 01_QC
singularity exec /usr/local/singularity-images/anvio-5.4.simg iu-filter-quality-minoche 01_QC/<AT26-13-87.ini>






#Step Three- No Assembly
#You can choose to do a read-based analysis depending on your question. This would be a great
#way to take a quick look at the data before doing an assembly-based analysis, and sometimes is sufficient
#for answering your scientific question. The read-based analysis includes steps three and four.
#you can run MetaPhlAn2 with paired end data, but run each fasta file separately. MetaPhlAn2 does not rely on mapping information.
#You can merge the outputs into one file if you do this.

#PBS -S /bin/bash
#PBS -q highmem_q
#PBS -N Step_3
#PBS -l nodes=1:ppn=1
#PBS -l walltime=240:00:00
#PBS -l mem=100g
#PBS -j oe


BASEDIR=</work/sbjlab/rck/Sed_Assemblies/AT26-13-87/>
cd $BASEDIR

module load MetaPhlAn2/2.7.8-foss-2016b-Python-2.7.14
metaphlan2.py <AT26-13_87_Sed_S49_L003_R1_001.fa>  --input_type fasta > <AT26-13-87_R1_profile.txt>
metaphlan2.py <AT26-13_87_Sed_S49_L003_R2_001.fa>  --input_type fasta > <AT26-13-87_R2_profile.txt>
cat <AT26-13-87_R1_profile.txt> <AT26-13-87_R2_profile.txt> > <87_Metaphlan_profile.txt>


awk '!h[$4,$NF]++ { print $4, $NF }' FS="|" <87_Metaphlan_profile.txt> > <87_temp.txt>
awk '!h[$1,$NF]++ { print $1, $NF }' FS=" " <87_temp.txt>  > <87_clean_profile.txt>
ls -1 | grep "^o" <87_clean_profile.txt> > <87_order_profile.txt>
sed 's/ \+/,/g' <87_order_profile.txt> > <87_order_profiles.csv>
awk -F, '!seen[$1]++' <87_order_profiles.csv> > <87_order_profile.csv>
rm <87_temp.txt> <87_clean_profile.txt>  <87_order_profiles.csv>



#PBS -S /bin/bash
#PBS -q highmem_q
#PBS -N Step_3
#PBS -l nodes=1:ppn=1
#PBS -l walltime=240:00:00
#PBS -l mem=100g
#PBS -j oe


BASEDIR=</work/sbjlab/rck/Sed_Assemblies/AT26-13-89/>
cd $BASEDIR

module load MetaPhlAn2/2.7.8-foss-2016b-Python-2.7.14
metaphlan2.py <AT26-13_89_S50_L003_R1_001.fa>  --input_type fasta > <AT26-13-89_R1_profile.txt>
metaphlan2.py <AT26-13_89_S50_L003_R2_001.fa>  --input_type fasta > <AT26-13-89_R2_profile.txt>
cat <AT26-13-89_R1_profile.txt> <AT26-13-89_R2_profile.txt> > <89_Metaphlan_profile.txt>

awk '!h[$4,$NF]++ { print $4, $NF }' FS="|" 89_Metaphlan_profile.txt >89_temp.txt
awk '!h[$1,$NF]++ { print $1, $NF }' FS=" " 89_temp.txt  >89_clean_profile.txt
ls -1 | grep "^o" 89_clean_profile.txt >89_order_profile.txt
sed 's/ \+/,/g' 89_order_profile.txt > 89_order_profiles.csv
awk -F, '!seen[$1]++' 89_order_profiles.csv > 89_order_profile.csv
rm 89_temp.txt 89_clean_profile.txt  89_order_profiles.csv



#PBS -S /bin/bash
#PBS -q highmem_q
#PBS -N Step_3
#PBS -l nodes=1:ppn=1
#PBS -l walltime=240:00:00
#PBS -l mem=100g
#PBS -j oe


BASEDIR=</work/sbjlab/rck/Sed_Assemblies/AT26-13-91/>
cd $BASEDIR

module load MetaPhlAn2/2.7.8-foss-2016b-Python-2.7.14
metaphlan2.py <AT26-13_91_Sed_S51_L003_R1_001.fa>  --input_type fasta > <AT26-13-91_R1_profile.txt>
metaphlan2.py <AT26-13_91_Sed_S51_L003_R2_001.fa>  --input_type fasta > <AT26-13-91_R2_profile.txt>
cat <AT26-13-91_R1_profile.txt> <AT26-13-91_R2_profile.txt> > <91_Metaphlan_profile.txt>



awk '!h[$4,$NF]++ { print $4, $NF }' FS="|" <91_Metaphlan_profile.txt> > <91_temp.txt>
awk '!h[$1,$NF]++ { print $1, $NF }' FS=" " <91_temp.txt>  > <91_clean_profile.txt>
ls -1 | grep "^o" <91_clean_profile.txt> > <91_order_profile.txt>
sed 's/ \+/,/g' <91_order_profile.txt> > <91_order_profiles.csv>
awk -F, '!seen[$1]++' <91_order_profiles.csv> > <91_order_profile.csv>
rm <91_temp.txt> <91_clean_profile.txt>  <91_order_profiles.csv>



#This part of step 3 is to combine outputs for all samples for visualization in Rstudio interactive on the cluster.
#This output can be viewed on your local machine.

#PBS -S /bin/bash
#PBS -q batch
#PBS -N Step_3
#PBS -l nodes=1:ppn=1
#PBS -l walltime=20:00:00
#PBS -l mem=10g
#PBS -j oe
BASEDIR=/work/sbjlab/rck/Sed_Assemblies
cd $BASEDIR

cat /work/sbjlab/rck/Sed_Assemblies/AT26-13-87/87_order_profile.csv /work/sbjlab/rck/Sed_Assemblies/AT26-13-89/89_order_profile.csv /work/sbjlab/rck/Sed_Assemblies/AT26-13-91/91_order_profile.csv >combined_87-89-91_Metaphlan_profile.csv

qlogin

module load R/3.4.4-foss-2016b-X11-20160819-GACRC
R
# at this point, be sure to install the packages d3heatmap and htmlwidgets. Instructions for install are here: https://blog.rstudio.com/2015/06/24/d3heatmap/
install.packages("d3heatmap")
install.packages("htmlwidgets")
R --no-save < /work/sbjlab/rck/Scripts_Pipelines/heatmap.r
scp *.html /work/sbjlab/rck/Scripts_Pipelines/




#Step Four- (co)Assembly
#You can choose to assemble each metagenome individually or to co-assemble them all together,
#and then map individual reads to that co-assembly. Read up on the theory before deciding which is right for your data.
#R1s.txt and R2s.txt are files which simply list the file location of the R1 and R2 reads you are co assembling.

#PBS _S /bin/bash
#PBS -N Step_5
#PBS -q highmem_q
#PBS -l nodes=1:ppn=32
#PBS -l walltime=400:00:00
#PBS -l mem=100gb

BASEDIR=</home/rck80079/Sed_Assemblies/>
cd $BASEDIR

module load MEGAHIT/1.1.3-foss-2016b-Python-2.7.14

megahit -1 <R1s.txt> -2 <R2s.txt> --min-contig-len 1000 -m 0.85 -o 02_ASSEMBLY/ -t 32
mkdir 03_CONTIGS

singularity exec /usr/local/singularity-images/anvio-5.4.simg anvi-script-reformat-fasta 02_ASSEMBLY/final.contigs.fa -o 03_CONTIGS/contigs.fa --min-len 2500 --simplify-names --report name_conversions.txt








#Step Five- Mapping
#Map the sample reads to the co-assembly to generate covereage information.

#PBS -S /bin/bash
#PBS -q batch
#PBS -N jobname
#PBS -l nodes=1:ppn=1
#PBS -l walltime=10:00:00
#PBS -l mem=50gb

BASEDIR=</home/rck80079/Sed_Assemblies/>
cd $BASEDIR

module load Bowtie2/2.3.4.1-foss-2016b
module load SAMtools/1.6-foss-2016b


mkdir 04_MAPPING
time bowtie2-build 03_CONTIGS/contigs.fa 04_MAPPING/contigs
time bowtie2 --threads $NUM_THREADS -x 04_MAPPING/contigs -1 <01_QC/AT26-13-87_R1.fastq> -2 <01_QC/AT26-13-87_R2.fastq> -S <04_MAPPING/AT26-13-87.sam>
time samtools view -F 4 -bS 04_MAPPING/<AT26-13-87.sam> > 04_MAPPING/<AT26-13-87-RAW.bam>
singularity exec /usr/local/singularity-images/anvio-5.4.simg anvi-init-bam <04_MAPPING/AT26-13-87-RAW.bam> -o <04_MAPPING/AT26-13-87.bam>
rm 04_MAPPING/<AT26-13-87.sam> 04_MAPPING/<AT26-13-87-RAW.bam>

#PBS -S /bin/bash
#PBS -q batch
#PBS -N jobname
#PBS -l nodes=1:ppn=1
#PBS -l walltime=10:00:00
#PBS -l mem=50gb

BASEDIR=</home/rck80079/Sed_Assemblies/>
cd $BASEDIR

module load Bowtie2/2.3.4.1-foss-2016b
module load SAMtools/1.6-foss-2016b


mkdir 04_MAPPING
time bowtie2-build 03_CONTIGS/contigs.fa 04_MAPPING/contigs
time bowtie2 --threads $NUM_THREADS -x 04_MAPPING/contigs -1 <01_QC/AT26-13-89_R1.fastq> -2 <01_QC/AT26-13-89_R2.fastq> -S 04_MAPPING/<AT26-13-89.sam>
time samtools view -F 4 -bS 04_MAPPING/<AT26-13-89.sam> > 04_MAPPING/<AT26-13-89-RAW.bam>
singularity exec /usr/local/singularity-images/anvio-5.4.simg anvi-init-bam 04_MAPPING/<AT26-13-89-RAW.bam> -o 04_MAPPING/<AT26-13-89.bam>
rm 04_MAPPING/<AT26-13-89.sam> 04_MAPPING/<AT26-13-89-RAW.bam>

#PBS -S /bin/bash
#PBS -q batch
#PBS -N jobname
#PBS -l nodes=1:ppn=1
#PBS -l walltime=10:00:00
#PBS -l mem=50gb

BASEDIR=</home/rck80079/Sed_Assemblies/>
cd $BASEDIR

module load Bowtie2/2.3.4.1-foss-2016b
module load SAMtools/1.6-foss-2016b


mkdir 04_MAPPING
time bowtie2-build 03_CONTIGS/contigs.fa 04_MAPPING/contigs
time bowtie2 --threads $NUM_THREADS -x 04_MAPPING/contigs -1 <01_QC/AT26-13-91_R1.fastq> -2 <01_QC/AT26-13-91_R2.fastq> -S 04_MAPPING/<AT26-13-91.sam>
time samtools view -F 4 -bS 04_MAPPING/<AT26-13-91.sam> > 04_MAPPING/<AT26-13-91-RAW.bam>
singularity exec /usr/local/singularity-images/anvio-5.4.simg anvi-init-bam 04_MAPPING/<AT26-13-91-RAW.bam> -o 04_MAPPING/<AT26-13-91.bam>
rm 04_MAPPING/<AT26-13-91.sam> 04_MAPPING/<AT26-13-91-RAW.bam>







#Step Six- Recovering Metagenome Assembled Genomes
#A note on MAGs- the accepted quality of MAGs has been put forward in the Woyke et al 2018 paper
#High-quality MAGs have greater than 90% completion, and <5% contamination. These are not to be
#considered the same as isolate genomes, but are a pretty reasonable representative genome of closely
#related genomic lineages.

#PBS -S /bin/bash
#PBS -q highmem_q
#PBS -N jobname
#PBS -l nodes=1:ppn=1
#PBS -l walltime=50:00:00
#PBS -l mem=50gb

BASEDIR=</home/rck80079/Sed_Assemblies/>
cd $BASEDIR


module load kaiju/1.6.2

singularity exec /usr/local/singularity-images/anvio-5.4.simg anvi-script-reformat-fasta 03_CONTIGS/contigs.fa -o contigs-fixed.fa -l 0 --simplify-names
singularity exec /usr/local/singularity-images/anvio-5.4.simg anvi-gen-contigs-database -f 03_CONTIGS/contigs.fa -o contigs.db -n 'AT26-13-87-89-91 contigs datbase'
singularity exec /usr/local/singularity-images/anvio-5.4.simg anvi-run-hmms -c contigs.db
singularity exec /usr/local/singularity-images/anvio-5.4.simg anvi-display-contigs-stats contigs.db
singularity exec /usr/local/singularity-images/anvio-5.4.simg anvi-run-ncbi-cogs -c CONTIGS.db --num-threads 20
singularity exec /usr/local/singularity-images/anvio-5.4.simg anvi-get-sequences-for-gene-calls -c CONTIGS.db -o gene_calls.fa
singularity exec /usr/local/singularity-images/anvio-5.4.simg makeDB.sh -e -t 20
kaiju -t /path/to/nodes.dmp \
      -f /path/to/kaiju_db.fmi \
      -i gene_calls.fa \
      -o gene_calls_nr.out \
      -z 16 \
      -v \
      1>job.out 2>job.err

singularity exec /usr/local/singularity-images/anvio-5.4.simg addTaxonNames -t /path/to/nodes.dmp \
      -n /path/to/names.dmp \
      -i gene_calls_nr.out \
      -o gene_calls_nr.names \
      -r superkingdom,phylum,order,class,family,genus,species

singularity exec /usr/local/singularity-images/anvio-5.4.simg anvi-import-taxonomy-for-genes -i gene_calls_nr.names \
      -c contigs.db \
      --just-do-it \
      -p kaiju

singularity exec /usr/local/singularity-images/anvio-5.4.simg anvi-init-bam <AT26-13-87-RAW.bam> -o <AT26-13-87.bam>
singularity exec /usr/local/singularity-images/anvio-5.4.simg anvi-init-bam <AT26-13-89-RAW.bam> -o <AT26-13-89.bam>
singularity exec /usr/local/singularity-images/anvio-5.4.simg anvi-init-bam <AT26-13-91-RAW.bam> -o <AT26-13-91.bam>
singularity exec /usr/local/singularity-images/anvio-5.4.simg anvi-profile -i <AT26-13-87.bam> -c contigs.db
singularity exec /usr/local/singularity-images/anvio-5.4.simg anvi-merge <AT26-13-87/PROFILE.db> <AT26-13-89/PROFILE.db> <AT26-13-91/PROFILE.db> -o SAMPLES-MERGED -c contigs.db
singularity exec /usr/local/singularity-images/anvio-5.4.simganvi-interactive -p SAMPLES-MERGED/PROFILE.db -c contigs.db
singularity exec /usr/local/singularity-images/anvio-5.4.simganvi-summarize -p SAMPLES-MERGED/PROFILE.db -c contigs.db -o SAMPLES-SUMMARY -C CONCOCT
singularity exec /usr/local/singularity-images/anvio-5.4.simg anvi-refine -p MERGED_PROFILE/PROFILE.db -c contigs.db -C CONCOCT -b Group_6




#Step Seven- Gene Calling
#You can do this step after Recovering MAGs from your metagenomes or on each metagenome.

#PBS -S /bin/bash
#PBS -N Step_4
#PBS -q highmem_q
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
#PBS -q highmem_q
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
#PBS -q highmem_q
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
#PBS -q highmem_q
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
