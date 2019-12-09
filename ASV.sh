#This code is for use in ecotyping 16S sequences (i.e. ASV analysis) as downloaded form a database (i.e. NCBI) in QIIME to make
#a phylogenetic tree which is usable in Interactive Tree of Life for figure making.
#Last edit: October 28, 2019, Rachael Storo


#Note: if running on the sapelo2 cluster, here is the job script template:
#!/bin/bash

#PBS -N Ecotyping_Overholt
#PBS -q highmem_q
#PBS -l nodes=1:ppn=24
#PBS -l walltime=35:00:00
#PBS -l mem=300gb

BASEDIR="/work/sbjlab/rck/Ecotyping_Overholt"
cd $BASEDIR

module load qiime2/2019.7_conda
echo "backend: Agg" > ~/.config/matplotlib/matplotlibrc
export LC_ALL=en_US.utf-8
export LANG=en_US.utf-8



### In command line/BASH:
# Import data file
# this is a file with all 16S sequences, each with a unique identifier
qiime tools import \
  --type 'SampleData[SequencesWithQuality]' \
  --input-path $BASEDIR/manifest-file.txt \
  --output-path single-end-demux_Overholt.qza \
  --input-format SingleEndFastqManifestPhred33V2

#Sequence Quality control and feature table construction
#DADA2- ie. ASVs
qiime dada2 denoise-single \
  --i-demultiplexed-seqs single-end-demux_Overholt.qza \
  --p-trim-left 0 \
  --p-trunc-len 120 \
  --o-representative-sequences $BASEDIR/rep-seqs-dada2.qza \
  --o-table $BASEDIR/table-dada2.qza \
  --o-denoising-stats $BASEDIR/stats-dada2.qza

qiime metadata tabulate \
  --m-input-file stats-dada2.qza \
  --o- $BASEDIR/visualization stats-dada2.qzv

mv rep-seqs-dada2.qza rep-seqs.qza
mv table-dada2.qza table.qza

#Summarize feature table and feature data
qiime feature-table summarize \
  --i-table table.qza \
  --o-visualization table.qzv \
  --m-sample-metadata-file Overholt-metadata.txt
qiime feature-table tabulate-seqs \
  --i-data rep-seqs.qza \
  --o-visualization rep-seqs.qzv

#Building a phylogenetic Tree
#Note: Qiime2 does have a workflow called align-to-tree-mafft-fasttree which
#does all of the next 3 steps, but it doesn't like large data inputs.

# Align sequences
qiime alignment mafft \
  --i-sequences rep-seqs.qza  \
  --p-n-threads 24 \
  --p-parttree \
  --o-alignment aligned-ep-seqs.qza \
  --verbose

#Build an unrooted tree
qiime phylogeny fasttree \
  --i-alignment aligned-ep-seqs.qza  \
  --p-n-threads 24 \
  --verbose  \
  --o-tree Overholt_unrooted_tree.qza

#Root the unrooted tree above
qiime phylogeny midpoint-root \
  --i-tree Overholt_unrooted_tree.qza \
  --o-rooted-tree Overholt_rooted_tree.qza \
  --verbose



#Alpha and Beta Diversity analyses
#This includes a list of statistics seen in this tutorial: https://docs.qiime2.org/2019.7/tutorials/moving-pictures/
#For more robust analyses, work with various diversity indecies ie inverse simpson
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny Overholt_rooted_tree.qza \
  --i-table table.qza \
  --p-sampling-depth 1103 \
  --m-metadata-file Overholt-metadata.txt \
  --output-dir core-metrics-results


qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/faith_pd_vector.qza \
  --m-metadata-file Overholt-metadata.txt \
  --o-visualization core-metrics-results/faith-pd-group-significance.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/evenness_vector.qza \
  --m-metadata-file Overholt-metadata.txt \
  --o-visualization core-metrics-results/evenness-group-significance.qzv

qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file Overholt-metadata.txt \
  --m-metadata-column body-site \
  --o-visualization core-metrics-results/unweighted-unifrac-body-site-significance.qzv \
  --p-pairwise

qiime diversity beta-group-significance \
  --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza \
  --m-metadata-file Overholt-metadata.txt \
  --m-metadata-column subject \
  --o-visualization core-metrics-results/unweighted-unifrac-subject-group-significance.qzv \
  --p-pairwise

#Generate emperor plots/pca visualizations
qiime emperor plot \
  --i-pcoa core-metrics-results/unweighted_unifrac_pcoa_results.qza \
  --m-metadata-file Overholt-metadata.txt \
  --o-visualization core-metrics-results/unweighted-unifrac-emperor-days-since-experiment-start.qzv

qiime emperor plot \
  --i-pcoa core-metrics-results/bray_curtis_pcoa_results.qza \
  --m-metadata-file Overholt-metadata.txt \
  --o-visualization core-metrics-results/bray-curtis-emperor-days-since-experiment-start.qzv


# Classify Taxonomy of Sequences
# You need to download the classifier of your preference at https://github.com/qiime2/docs/blob/master/source/data-resources.rst
# alternatively, see https://docs.qiime2.org/2019.4/tutorials/feature-classifier/ to train your own classifier to fit your needs (i.e. for a different gene of interest)

qiime feature-classifier classify-sklearn \
  --i-classifier silva-132-99-nb-classifier.qza \
  --i-reads rep-seqs.qza \
  --o-classification taxonomy.qza \
  --verbose

qiime taxa barplot \
  --i-table table.qza \
  --i-taxonomy taxonomy.qza \
  --m-metadata-file Overholt-metadata.txt \
  --o-visualization taxa-bar-plots.qzv



####qiime feature-table filter-samples \
  #--i-table table.qza \
  #--m-metadata-file sample-metadata.tsv \
  #--p-where "[body-site]='gut'" \
  #--o-filtered-table gut-table.qza

#qiime composition add-pseudocount \
  #--i-table gut-table.qza \
  #--o-composition-table comp-gut-table.qz

  #qiime composition ancom \
  #  --i-table comp-gut-table.qza \
  #  --m-metadata-file sample-metadata.tsv \
  #  --m-metadata-column subject \
  #  --o-visualization ancom-subject.qzv



  #qiime taxa collapse \
#  --i-table gut-table.qza \
#  --i-taxonomy taxonomy.qza \
#  --p-level 6 \
#  --o-collapsed-table gut-table-l6.qza

#qiime composition add-pseudocount \
#  --i-table gut-table-l6.qza \
#  --o-composition-table comp-gut-table-l6.qza

#qiime composition ancom \
#  --i-table comp-gut-table-l6.qza \
#  --m-metadata-file sample-metadata.tsv \
#  --m-metadata-column subject \
#  --o-visualization l6-ancom-subject.qzv

# Make a visual taxonomy table
qiime metadata tabulate \
  --m-input-file taxonomy.qza \
  --o-visualization taxonomy.qzv

# Export Tree and Taxa files from Qiime2

qiime tools export --input-path taxonomy.qza --output-path ~/Desktop/FastqReview
qiime tools export --input-path rooted-tree.qza --output-path ~/Desktop/FastqReview


# Format for ITOL
# be sure to take a look and make sure your identifiers match in both of these files. If they don't, you've done something wrong and won't
# see taxonomy on your tree.

echo $'LABELS\nSEPARATOR TAB\nDATA' > taxonomy.txt
sed "1d" taxonomy.tsv | cut -f1,2 >> taxonomy.txt
