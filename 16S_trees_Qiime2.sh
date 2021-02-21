#This code is for use in processing 16S sequences as downloaded from a database (i.e. NCBI) or extracted using REAGO in QIIME to make
#a phylogenetic tree which is usable in Interactive Tree of Life for figure making.
#You can do all of these steps individually without Qiime2.
#Last edit: February 21, 2021 Rachael Storo


### In command line/BASH:
# Import data file
# this is a file with all 16S sequences, each with a unique identifier

qiime tools import \
  --input-path Input_seqs.fasta \
  --output-path seqs.qza \
  --type 'FeatureData[Sequence]'

####If you have too many representative sequences and you are getting a Maaft error, align your sequences in maaft before importing to Qiime as follows:

#qiime tools import \
  #--input-path aligned-sequences.fna \
  #--output-path aligned-sequences.qza \
  #--type 'FeatureData[AlignedSequence]'

# Build Tree

 qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences seqs.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza


# Classify Taxonomy of Sequences
# You need to download the classifier of your preference at https://github.com/qiime2/docs/blob/master/source/data-resources.rst
# alternatively, see https://docs.qiime2.org/2019.4/tutorials/feature-classifier/ to train your own classifier to fit your needs (i.e. for a different gene of interest)

 qiime feature-classifier classify-sklearn \
  --i-classifier silva-132-99-nb-classifier.qza \
  --i-reads seqs.qza \
  --o-classification taxonomy.qza \
  --verbose

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
