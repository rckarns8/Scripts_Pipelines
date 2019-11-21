#Metagenomics pipeline
#Initiated: October 30, 2019 by Rachael Storo
#Last Edit: November 21, 2019 by Rachael Storo
#Purpose: Create a script for the processing of environmental metagenomic samples
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



#Step One- Demultiplex
# Not all data will require this step- in fact, many metagenomes are sequenced individually or are given to you already demultiplexed.
#For the purpose of future use, I will include this step, but data I am currently working with do not require it.

##PBS -S /bin/bash
##PBS -q batch
##PBS -N j_usearch
##PBS -l nodes=1:ppn=1:AMD
##PBS -l walltime=100:00:00
##PBS -l mem=100gb

#cd $PBS_O_WORKDIR

#module load USEARCH/10.0.240-i86linux32
#usearch -fastx_demux R1.fq -reverse R2.fq -index I1.fq -barcodes bar.fa \
        #-fastqout fwd_demux.fq -output2 rev_demux.fq







#Step Two- Quality trimming and filtering
#Quality filtering will remove adapters and primers to give you clean fasta files.

#PBS -S /bin/bash
#PBS -N j_trim
#PBS -q batch
#PBS -l nodes=1:ppn=4:AMD
#PBS -l walltime=100:00:00
#PBS -l mem=10gb

cd $PBS_O_WORKDIR

module load Trimmomatic/0.36-Java-1.8.0_144

time java -jar /usr/local/apps/eb/Trimmomatic/0.36-Java-1.8.0_144/trimmomatic-0.36.jar PE -threads 4 input_forward.fq.gz input_reverse.fq.gz output_forward_paired.fq.gz output_forward_unpaired.fq.gz output_reverse_paired.fq.gz output_reverse_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:36








#Step Three- No Assembly
#You can choose to do a read-based analysis depending on your question. This would be a great
#way to take a quick look at the data before doing an assembly-based analysis, and sometimes is sufficient
#for answering your scientific question. The read-based analysis includes steps three through five.

#PBS -S /bin/bash
#PBS -q batch
#PBS -N jobname
#PBS -l nodes=1:ppn=1
#PBS -l walltime=240:00:00
#PBS -l mem=100g
#PBS -j oe

#PBS -M username@uga.edu
#PBS -m ae

cd $PBS_O_WORKDIR

module load MetaPhlAn2/2.7.8-foss-2016b-Python-2.7.14
metaphlan2.py Sample.fasta.gz  --input_type fasta > Sample_profile.txt








#Step Four- Gene Calling and Taxonomic Profiling
#In addition to being part of the read based analysis, you can also do this step after Recovering
#MAGs from your metagenomes (Step eight below).

#PBS -S /bin/bash
#PBS -N j_prokka
#PBS -q batch
#PBS -l nodes=1:ppn=1
#PBS -l walltime=120:00:00
#PBS -l mem=20gb

cd $PBS_O_WORKDIR
ml prokka/1.13-foss-2016b-BioPerl-1.7.1
time prokka contigs.fa








#Step Five- (co)Assembly
#You can choose to assemble each metagenome individually or to co-assemble them all together,
#and then map individual reads to that co-assembly. Read up on the theory before deciding which is right for your data.

#PBS _S /bin/bash
#PBS -N j_megahit
#PBS -q gpu_q
#PBS -l nodes=1:ppn=12
#PBS -l walltime=480:00:00
#PBS -l mem=10gb

cd $PBS_O_WORKDIR

module load MEGAHIT/1.1.3-foss-2016b-Python-2.7.14

R1s=`ls 01_QC/*QUALITY_PASSED_R1* | python -c 'import sys; print ",".join([x.strip() for x in sys.stdin.readlines()])'`
R2s=`ls 01_QC/*QUALITY_PASSED_R2* | python -c 'import sys; print ",".join([x.strip() for x in sys.stdin.readlines()])'`
megahit -1 $R1s -2 $R2s --min-contig-len 1000 -m 0.85 -o 02_ASSEMBLY/ -t 40
mkdir 03_CONTIGS

singularity exec /usr/local/singularity-images/anvio-5.1.simg anvi-script-reformat-fasta 02_ASSEMBLY/final.contigs.fa -o 03_CONTIGS/contigs.fa --min-len 2500 --simplify-names --report name_conversions.txt








#Step Six- Mapping
#Map the individual sample reads to the co-assembly to generate covereage information.

#PBS -S /bin/bash
#PBS -q batch
#PBS -N jobname
#PBS -l nodes=1:ppn=1
#PBS -l walltime=4:00:00
#PBS -l mem=2gb

cd $PBS_O_WORKDIR

module load Bowtie2/2.3.4.1-foss-2016b
module load SAMtools/1.6-foss-2016b


mkdir 04_MAPPING
time bowtie2-build 03_CONTIGS/contigs.fa 04_MAPPING/contigs
time bowtie2 --threads $NUM_THREADS -x 04_MAPPING/contigs -1 01_QC/Sample_01-QUALITY_PASSED_R1.fastq -2 01_QC/Sample_01-QUALITY_PASSED_R2.fastq -S 04_MAPPING/Sample_01.sam
time samtools view -F 4 -bS 04_MAPPING/Sample_01.sam > 04_MAPPING/Sample_01-RAW.bam
singularity exec /usr/local/singularity-images/anvio-5.1.simg anvi-init-bam 04_MAPPING/Sample_01-RAW.bam -o 04_MAPPING/Sample_01.bam
rm 04_MAPPING/Sample_01.sam 04_MAPPING/Sample_01-RAW.bam








#Step Seven- Recovering Metagenome Assembled Genomes
#A note on MAGs- the accepted quality of MAGs has been put forward in the Woyke et al 2018 paper
#High-quality MAGs have greater than 90% completion, and <5% contamination. These are not to be
#considered the same as isolate genomes, but are a pretty reasonable representative genome of closely
#related genomic lineages.

#PBS -S /bin/bash
#PBS -q batch
#PBS -N jobname
#PBS -l nodes=1:ppn=1
#PBS -l walltime=4:00:00
#PBS -l mem=2gb

cd $PBS_O_WORKDIR
module load kaiju/1.6.2

singularity exec /usr/local/singularity-images/anvio-5.1.simg anvi-script-reformat-fasta contigs.fa -o contigs-fixed.fa -l 0 --simplify-names
singularity exec /usr/local/singularity-images/anvio-5.1.simg anvi-gen-contigs-database -f contigs.fa -o contigs.db -n 'An example contigs datbase'
singularity exec /usr/local/singularity-images/anvio-5.1.simg anvi-run-hmms -c contigs.db
singularity exec /usr/local/singularity-images/anvio-5.1.simg anvi-display-contigs-stats contigs.db
singularity exec /usr/local/singularity-images/anvio-5.1.simg anvi-run-ncbi-cogs -c CONTIGS.db --num-threads 20
singularity exec /usr/local/singularity-images/anvio-5.1.simg anvi-get-sequences-for-gene-calls -c CONTIGS.db -o gene_calls.fa
singularity exec /usr/local/singularity-images/anvio-5.1.simg makeDB.sh -e -t 20
kaiju -t /path/to/nodes.dmp \
      -f /path/to/kaiju_db.fmi \
      -i gene_calls.fa \
      -o gene_calls_nr.out \
      -z 16 \
      -v \
      1>job.out 2>job.err

singularity exec /usr/local/singularity-images/anvio-5.1.simg addTaxonNames -t /path/to/nodes.dmp \
      -n /path/to/names.dmp \
      -i gene_calls_nr.out \
      -o gene_calls_nr.names \
      -r superkingdom,phylum,order,class,family,genus,species

singularity exec /usr/local/singularity-images/anvio-5.1.simg anvi-import-taxonomy-for-genes -i gene_calls_nr.names \
      -c contigs.db \
      --just-do-it \
      -p kaiju

singularity exec /usr/local/singularity-images/anvio-5.1.simg anvi-init-bam SAMPLE-01-RAW.bam -o SAMPLE-01.bam
for sample in `cat SAMPLE_IDs`; do anvi-init-bam $sample-RAW.bam -o $sample.bam; done
singularity exec /usr/local/singularity-images/anvio-5.1.simg anvi-profile -i SAMPLE-01.bam -c contigs.db
singularity exec /usr/local/singularity-images/anvio-5.1.simg anvi-merge SAMPLE-01/PROFILE.db SAMPLE-02/PROFILE.db SAMPLE-03/PROFILE.db -o SAMPLES-MERGED -c contigs.db
singularity exec /usr/local/singularity-images/anvio-5.1.simganvi-interactive -p SAMPLES-MERGED/PROFILE.db -c contigs.db
singularity exec /usr/local/singularity-images/anvio-5.1.simganvi-summarize -p SAMPLES-MERGED/PROFILE.db -c contigs.db -o SAMPLES-SUMMARY -C CONCOCT
singularity exec /usr/local/singularity-images/anvio-5.1.simg anvi-refine -p MERGED_PROFILE/PROFILE.db -c contigs.db -C CONCOCT -b Group_6








#Step Eight- Phylogenomics

#PBS -S /bin/bash
#PBS -q batch
#PBS -N jobname
#PBS -l nodes=1:ppn=1
#PBS -l walltime=4:00:00
#PBS -l mem=2gb

cd $PBS_O_WORKDIR

singularity exec /usr/local/singularity-images/anvio-5.1.simg anvi-get-sequences-for-hmm-hits -c CONTIGS.db \
      -p PROFILE.db \
      -o seqs-for-phylogenomics.fa \
      --list-hmm-sources

singularity exec /usr/local/singularity-images/anvio-5.1.simg anvi-get-sequences-for-hmm-hits -c CONTIGS.db \
      -p PROFILE.db \
      -o seqs-for-phylogenomics.fa \
      --hmm-source Bacteria_71 \
      -C default \
      --gene-names Ribosomal_L1,Ribosomal_L2,Ribosomal_L3,Ribosomal_L4,Ribosomal_L5,Ribosomal_L6
singularity exec /usr/local/singularity-images/anvio-5.1.simg anvi-get-sequences-for-hmm-hits -c CONTIGS.db \
      -p PROFILE.db \
      -o seqs-for-phylogenomics.fa \
      --hmm-source Bacteria_71 \
      -C default \
      --gene-names Ribosomal_L1,Ribosomal_L2,Ribosomal_L3,Ribosomal_L4,Ribosomal_L5,Ribosomal_L6 \
      --concatenate-genes \
      --return-best-hit \
      --get-aa-sequences

singularity exec /usr/local/singularity-images/anvio-5.1.simg anvi-gen-phylogenomic-tree -f seqs-for-phylogenomics.fa \
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
#PBS -q batch
#PBS -N jobname
#PBS -l nodes=1:ppn=1
#PBS -l walltime=4:00:00
#PBS -l mem=2gb

cd $PBS_O_WORKDIR
singularity exec /usr/local/singularity-images/anvio-5.1.simg anvi-import-collection additional-files/collections/sample.txt \
    --bins-info additional-files/collections/info.txt \
    -p PROFILE.db \
    -c CONTIGS.db \
    -C Sample_name

singularity exec /usr/local/singularity-images/anvio-5.1.simg anvi-gen-genomes-storage -i additional-files/pangenomics/internal-genomes.txt \
    -e additional-files/pangenomics/external-genomes.txt \
    -o Sample_name-GENOMES.db

anvi-pan-genome -g Sample_name-GENOMES.db \
    -n Sample_name \
    -o PAN \
    --num-threads 10

anvi-display-pan -g Sample_name-GENOMES.db \
    -p PAN/Sample_name-PAN.db \
    --title "Sample_name Pan"


anvi-import-state -p PAN/Sample_name-PAN.db \
    --state additional-files/state-files/state-pan.json \
    --name default

anvi-display-pan -g Sample_name-GENOMES.db \
    -p PAN/Sample_name-PAN.db \
    --title "Sample_name Pan"


anvi-compute-genome-similarity -e additional-files/pangenomics/external-genomes.txt \
    -i additional-files/pangenomics/internal-genomes.txt \
    --program pyANI \
    -o ANI \
    -T 6 \
    --pan-db PAN/Sample_name-PAN.db

anvi-display-pan -g Sample_name-GENOMES.db \
    -p PAN/Sample_name-PAN.db \
    --title "Sample_name Pan"

anvi-import-misc-data -p PAN/Sample_name-PAN.db \
   --target-data-table layers \
   additional-files/pangenomics/additional-layers-data.txt


anvi-display-pan -g Sample_name-GENOMES.db \
   -p PAN/Sample_name-PAN.db \
   --title "Sample_name Pan"

anvi-import-collection additional-files/pangenomics/pan-collection.txt \
    --bins-info additional-files/pangenomics/pan-collection-info.txt \
    -p PAN/Sample_name-PAN.db \
    -C default

anvi-display-pan -g Enterococcus-GENOMES.db \
    -p PAN/Sample_name-PAN.db \
    --title "Sample_name Pan"

anvi-summarize -p PAN/Sample_name-PAN.db \
    -g Sample_name-GENOMES.db \
    -C default \
    -o PAN_SUMMARY

open PAN_SUMMARY/index.html
gzip -d PAN_SUMMARY/Sample_name_protein_clusters_summary.txt.gz








#Step Ten- SNV analysis
#Profiling SNVs allows for you to examine the microbial population genetics in your metagenomes

#PBS -S /bin/bash
#PBS -q batch
#PBS -N jobname
#PBS -l nodes=1:ppn=1
#PBS -l walltime=4:00:00
#PBS -l mem=2gb

cd $PBS_O_WORKDIR
singularity exec /usr/local/singularity-images/anvio-5.1.simg anvi-import-collection additional-files/collections/merens.txt \
      --bins-info additional-files/collections/merens-info.txt \
      -p PROFILE.db \
      -c CONTIGS.db \
      -C default

# importing taxonomy for gene calls
singularity exec /usr/local/singularity-images/anvio-5.1.simg anvi-import-taxonomy-for-genes -c CONTIGS.db \
      -i additional-files/centrifuge-files/centrifuge_report.tsv \
      additional-files/centrifuge-files/centrifuge_hits.tsv \
      -p centrifuge

# importing the state file so things look pretty
singularity exec /usr/local/singularity-images/anvio-5.1.simg anvi-import-state --state additional-files/state-files/state-merged.json \
      --name default \
      -p PROFILE.db

singularity exec /usr/local/singularity-images/anvio-5.1.simg anvi-split -p PROFILE.db \
      -c CONTIGS.db \
      -C default \
      -b Sample_name \
      -o MAGs

singularity exec /usr/local/singularity-images/anvio-5.1.simg anvi-estimate-genome-completeness -p MAGs/E_facealis/PROFILE.db \
      -c MAGs/Sample_name/CONTIGS.db \
      -C DEFAULT

anvi-interactive -p MAGs/Sample_name/PROFILE.db \
      -c MAGs/Sample_name/CONTIGS.db

singularity exec /usr/local/singularity-images/anvio-5.1.simg anvi-gen-variability-profile -c MAGs/Sample_name/CONTIGS.db \
      -p MAGs/Sample_name/PROFILE.db \
      -C DEFAULT \
      -b ALL_SPLITS \
      --samples-of-interest additional-files/samples.txt \
      --min-coverage-in-each-sample 20 \
      --min-occurrence 3 \
      --include-split-names \
      --quince-mode \
      -o Sample_name-SNVs.txt

cat additional-files/samples.txt
singularity exec /usr/local/singularity-images/anvio-5.1.simg anvi-script-snvs-to-interactive Sample_name-SNVs.txt -o Sample_name_snvs
anvi-interactive --profile Sample_name_snvs/profile.db \
                 --tree Sample_name_snvs/tree.txt \
                 --view-data Sample_name_snvs/view.txt \
                 --title "SNV Profile for the Sample_name bin" \
                 --manual

singularity exec /usr/local/singularity-images/anvio-5.1.simg anvi-import-state -p Sample_name_snvs/profile.db \
                 --state additional-files/state-files/state-snvs.json \
                 --name default

anvi-interactive -d Sample_name_snvs/view.txt \
                 -t Sample_name_snvs/tree.txt \
                 -p Sample_name_snvs/profile.db \
                 --title "SNV Profile for the Sample_name bin" \
                 --manual

singularity exec /usr/local/singularity-images/anvio-5.1.simg anvi-gen-fixation-index-matrix --variability-profile Sample_name-SNVs.txt \
               --output-file FST_Sample_name.txt

singularity exec /usr/local/singularity-images/anvio-5.1.simg anvi-matrix-to-newick FST_Sample_name.txt \
               --output-file FST_Sample_name.newick

anvi-interactive -t FST_Sample_name.newick \
                 -p FST_Sample_name.db \
                 --manual
