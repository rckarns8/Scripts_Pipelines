#Metagenomics pipeline
#Initiated: October 30, 2019 by Rachael Storo
#Last Edit: October 31, 2019 by Rachael Storo
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

#Step Four- (co)Assembly
R1s=`ls 01_QC/*QUALITY_PASSED_R1* | python -c 'import sys; print ",".join([x.strip() for x in sys.stdin.readlines()])'`
R2s=`ls 01_QC/*QUALITY_PASSED_R2* | python -c 'import sys; print ",".join([x.strip() for x in sys.stdin.readlines()])'`
megahit -1 $R1s -2 $R2s --min-contig-len 1000 -m 0.85 -o 02_ASSEMBLY/ -t 40
mkdir 03_CONTIGS
anvi-script-reformat-fasta 02_ASSEMBLY/final.contigs.fa -o 03_CONTIGS/contigs.fa --min-len 2500 --simplify-names --report name_conversions.txt


#Step Five- Gene Calling and Taxonomic Profiling

prokka contigs.fa

#Step Six- Read Based Analysis

phyloseq??

#Step Seven- Mapping
#Map the individual sample reads to the co-assembly
mkdir 04_MAPPING
bowtie2-build 03_CONTIGS/contigs.fa 04_MAPPING/contigs
bowtie2 --threads $NUM_THREADS -x 04_MAPPING/contigs -1 01_QC/Sample_01-QUALITY_PASSED_R1.fastq -2 01_QC/Sample_01-QUALITY_PASSED_R2.fastq -S 04_MAPPING/Sample_01.sam
samtools view -F 4 -bS 04_MAPPING/Sample_01.sam > 04_MAPPING/Sample_01-RAW.bam
anvi-init-bam 04_MAPPING/Sample_01-RAW.bam -o 04_MAPPING/Sample_01.bam
rm 04_MAPPING/Sample_01.sam 04_MAPPING/Sample_01-RAW.bam

#Step Eight- Recovering Metagenome Assembled Genomes

anvi-script-reformat-fasta contigs.fa -o contigs-fixed.fa -l 0 --simplify-names
anvi-gen-contigs-database -f contigs.fa -o contigs.db -n 'An example contigs datbase'
anvi-run-hmms -c contigs.db
anvi-display-contigs-stats contigs.db
anvi-run-ncbi-cogs -c CONTIGS.db --num-threads 20
anvi-get-sequences-for-gene-calls -c CONTIGS.db -o gene_calls.fa
makeDB.sh -e -t 20
kaiju -t /path/to/nodes.dmp \
      -f /path/to/kaiju_db.fmi \
      -i gene_calls.fa \
      -o gene_calls_nr.out \
      -z 16 \
      -v

addTaxonNames -t /path/to/nodes.dmp \
      -n /path/to/names.dmp \
      -i gene_calls_nr.out \
      -o gene_calls_nr.names \
      -r superkingdom,phylum,order,class,family,genus,species

anvi-import-taxonomy-for-genes -i gene_calls_nr.names \
      -c contigs.db \
      --just-do-it \
      -p kaiju

anvi-init-bam SAMPLE-01-RAW.bam -o SAMPLE-01.bam
for sample in `cat SAMPLE_IDs`; do anvi-init-bam $sample-RAW.bam -o $sample.bam; done
anvi-profile -i SAMPLE-01.bam -c contigs.db
anvi-merge SAMPLE-01/PROFILE.db SAMPLE-02/PROFILE.db SAMPLE-03/PROFILE.db -o SAMPLES-MERGED -c contigs.db
anvi-interactive -p SAMPLES-MERGED/PROFILE.db -c contigs.db
anvi-summarize -p SAMPLES-MERGED/PROFILE.db -c contigs.db -o SAMPLES-SUMMARY -C CONCOCT



anvi-refine -p MERGED_PROFILE/PROFILE.db -c contigs.db -C CONCOCT -b Group_6
#Step Nine- Phylogenomics
anvi-get-sequences-for-hmm-hits -c CONTIGS.db \
      -p PROFILE.db \
      -o seqs-for-phylogenomics.fa \
      --list-hmm-sources


anvi-get-sequences-for-hmm-hits -c CONTIGS.db \
      -p PROFILE.db \
      -o seqs-for-phylogenomics.fa \
      --hmm-source Bacteria_71 \
      -C default \
      --gene-names Ribosomal_L1,Ribosomal_L2,Ribosomal_L3,Ribosomal_L4,Ribosomal_L5,Ribosomal_L6
anvi-get-sequences-for-hmm-hits -c CONTIGS.db \
      -p PROFILE.db \
      -o seqs-for-phylogenomics.fa \
      --hmm-source Bacteria_71 \
      -C default \
      --gene-names Ribosomal_L1,Ribosomal_L2,Ribosomal_L3,Ribosomal_L4,Ribosomal_L5,Ribosomal_L6 \
      --concatenate-genes \
      --return-best-hit \
      --get-aa-sequences

anvi-gen-phylogenomic-tree -f seqs-for-phylogenomics.fa \
      -o phylogenomic-tree.txt

anvi-interactive --tree phylogenomic-tree.txt \
      -p temp-profile.db \
      --title "Pylogenomics of IGD Bins" \
      --manual


anvi-interactive -p PROFILE.db \
      -c CONTIGS.db \
      -C default \
      --tree phylogenomic-tree.txt

#Step Ten- Comparative Genomics

anvio

#Step Eleven- Pangenomics

anvio
