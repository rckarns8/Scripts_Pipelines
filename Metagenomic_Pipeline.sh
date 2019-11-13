#Metagenomics pipeline
#Initiated: October 30, 2019 by Rachael Storo
#Last Edit: November 13, 2019 by Rachael Storo
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
#Quality filtering will remove adapters and primers to give you clean fasta files.
java -jar trimmomatic-0.39.jar PE input_forward.fq.gz input_reverse.fq.gz output_forward_paired.fq.gz output_forward_unpaired.fq.gz output_reverse_paired.fq.gz output_reverse_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:keepBothReads LEADING:3 TRAILING:3 MINLEN:36

#Step Three- No Assembly
#You can choose to do a read-based analysis depending on your question. This would be a great
#way to take a quick look at the data before doing an assembly-based analysis, and sometimes is sufficient
#for answering your scientific question. The read-based analysis includes steps three through five.

metaphlan2.py SRS014476-Supragingival_plaque.fasta.gz  --input_type fasta > SRS014476-Supragingival_plaque_profile.txt

#Step Four- Gene Calling and Taxonomic Profiling
#In addition to being part of the read based analysis, you can also do this step after Recovering
#MAGs from your metagenomes (Step eight below).

prokka contigs.fa

#Step Five- Read Based Analysis

phyloseq??

#Step Six- (co)Assembly
#You can choose to assemble each metagenome individually or to co-assemble them all together,
#and then map individual reads to that co-assembly. Read up on the theory before deciding which is right for your data.

R1s=`ls 01_QC/*QUALITY_PASSED_R1* | python -c 'import sys; print ",".join([x.strip() for x in sys.stdin.readlines()])'`
R2s=`ls 01_QC/*QUALITY_PASSED_R2* | python -c 'import sys; print ",".join([x.strip() for x in sys.stdin.readlines()])'`
megahit -1 $R1s -2 $R2s --min-contig-len 1000 -m 0.85 -o 02_ASSEMBLY/ -t 40
mkdir 03_CONTIGS
anvi-script-reformat-fasta 02_ASSEMBLY/final.contigs.fa -o 03_CONTIGS/contigs.fa --min-len 2500 --simplify-names --report name_conversions.txt

#Step Seven- Mapping
#Map the individual sample reads to the co-assembly to generate covereage information.

mkdir 04_MAPPING
bowtie2-build 03_CONTIGS/contigs.fa 04_MAPPING/contigs
bowtie2 --threads $NUM_THREADS -x 04_MAPPING/contigs -1 01_QC/Sample_01-QUALITY_PASSED_R1.fastq -2 01_QC/Sample_01-QUALITY_PASSED_R2.fastq -S 04_MAPPING/Sample_01.sam
samtools view -F 4 -bS 04_MAPPING/Sample_01.sam > 04_MAPPING/Sample_01-RAW.bam
anvi-init-bam 04_MAPPING/Sample_01-RAW.bam -o 04_MAPPING/Sample_01.bam
rm 04_MAPPING/Sample_01.sam 04_MAPPING/Sample_01-RAW.bam

#Step Eight- Recovering Metagenome Assembled Genomes
#A note on MAGs- the accepted quality of MAGs has been put forward in the Woyke et al 2018 paper
#High-quality MAGs have greater than 90% completion, and <5% contamination. These are not to be
#considered the same as isolate genomes, but are a pretty reasonable representative genome of closely
#related genomic lineages.

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

#Step Ten- Comparative (Pan) Genomics

anvi-import-collection additional-files/collections/e-faecalis.txt \
    --bins-info additional-files/collections/e-faecalis-info.txt \
    -p PROFILE.db \
    -c CONTIGS.db \
    -C E_faecalis

anvi-gen-genomes-storage -i additional-files/pangenomics/internal-genomes.txt \
    -e additional-files/pangenomics/external-genomes.txt \
    -o Enterococcus-GENOMES.db

anvi-pan-genome -g Enterococcus-GENOMES.db \
    -n Enterococcus \
    -o PAN \
    --num-threads 10

anvi-display-pan -g Enterococcus-GENOMES.db \
    -p PAN/Enterococcus-PAN.db \
    --title "Enterococccus Pan"


anvi-import-state -p PAN/Enterococcus-PAN.db \
    --state additional-files/state-files/state-pan.json \
    --name default

anvi-display-pan -g Enterococcus-GENOMES.db \
    -p PAN/Enterococcus-PAN.db \
    --title "Enterococccus Pan"


anvi-compute-genome-similarity -e additional-files/pangenomics/external-genomes.txt \
    -i additional-files/pangenomics/internal-genomes.txt \
    --program pyANI \
    -o ANI \
    -T 6 \
    --pan-db PAN/Enterococcus-PAN.db

anvi-display-pan -g Enterococcus-GENOMES.db \
    -p PAN/Enterococcus-PAN.db \
    --title "Enterococccus Pan"

anvi-import-misc-data -p PAN/Enterococcus-PAN.db \
   --target-data-table layers \
   additional-files/pangenomics/additional-layers-data.txt


anvi-display-pan -g Enterococcus-GENOMES.db \
   -p PAN/Enterococcus-PAN.db \
   --title "Enterococccus Pan"

anvi-import-collection additional-files/pangenomics/pan-collection.txt \
    --bins-info additional-files/pangenomics/pan-collection-info.txt \
    -p PAN/Enterococcus-PAN.db \
    -C default

anvi-display-pan -g Enterococcus-GENOMES.db \
    -p PAN/Enterococcus-PAN.db \
    --title "Enterococccus Pan"

anvi-summarize -p PAN/Enterococcus-PAN.db \
    -g Enterococcus-GENOMES.db \
    -C default \
    -o PAN_SUMMARY

open PAN_SUMMARY/index.html
gzip -d PAN_SUMMARY/Enterococcus_protein_clusters_summary.txt.gz


#Step 11- SNV analysis
#Profiling SNVs allows for you to examine the microbial population genetics in your metagenomes

anvi-import-collection additional-files/collections/merens.txt \
      --bins-info additional-files/collections/merens-info.txt \
      -p PROFILE.db \
      -c CONTIGS.db \
      -C default

# importing taxonomy for gene calls
anvi-import-taxonomy-for-genes -c CONTIGS.db \
      -i additional-files/centrifuge-files/centrifuge_report.tsv \
      additional-files/centrifuge-files/centrifuge_hits.tsv \
      -p centrifuge

# importing the state file so things look pretty
anvi-import-state --state additional-files/state-files/state-merged.json \
      --name default \
      -p PROFILE.db

anvi-split -p PROFILE.db \
      -c CONTIGS.db \
      -C default \
      -b E_facealis \
      -o MAGs

anvi-estimate-genome-completeness -p MAGs/E_facealis/PROFILE.db \
      -c MAGs/E_facealis/CONTIGS.db \
      -C DEFAULT

anvi-interactive -p MAGs/E_facealis/PROFILE.db \
      -c MAGs/E_facealis/CONTIGS.db

anvi-gen-variability-profile -c MAGs/E_facealis/CONTIGS.db \
      -p MAGs/E_facealis/PROFILE.db \
      -C DEFAULT \
      -b ALL_SPLITS \
      --samples-of-interest additional-files/samples.txt \
      --min-coverage-in-each-sample 20 \
      --min-occurrence 3 \
      --include-split-names \
      --quince-mode \
      -o E-faecalis-SNVs.txt

cat additional-files/samples.txt
anvi-script-snvs-to-interactive E-faecalis-SNVs.txt -o e_faecalis_snvs
anvi-interactive --profile e_faecalis_snvs/profile.db \
                 --tree e_faecalis_snvs/tree.txt \
                 --view-data e_faecalis_snvs/view.txt \
                 --title "SNV Profile for the E. faecalis bin" \
                 --manual

anvi-import-state -p e_faecalis_snvs/profile.db \
                 --state additional-files/state-files/state-snvs.json \
                 --name default

anvi-interactive -d e_faecalis_snvs/view.txt \
                 -t e_faecalis_snvs/tree.txt \
                 -p e_faecalis_snvs/profile.db \
                 --title "SNV Profile for the E. faecalis bin" \
                 --manual

anvi-gen-fixation-index-matrix --variability-profile E-faecalis-SNVs.txt \
               --output-file FST_E_facealis.txt

anvi-matrix-to-newick FST_E_facealis.txt \
               --output-file FST_E_facealis.newick

anvi-interactive -t FST_E_facealis.newick \
                 -p FST_E_facealis.db \
                 --manual
