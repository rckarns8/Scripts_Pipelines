#This code is for use in processing 16S sequences as downloaded from a database (i.e. NCBI) or extracted using REAGO in QIIME to make
#a phylogenetic tree which is usable in Interactive Tree of Life for figure making.
#You can do all of these steps individually without Qiime2.
#Last edit: February 21, 2021 Rachael Storo


### In command line/BASH:
# Import data file
# this is a file with all 16S sequences, each with a unique identifier

#!/bin/bash

#SBATCH --partition=joye_p
#SBATCH --job-name=wc98
#SBATCH --ntasks=1
#SBATCH --time=10:00:00
#SBATCH --mem=200G

module load QIIME2/2020.11

qiime tools import \
  --input-path /scratch/rck80079/BOEM/REAGO_out/wc98/unique_wc98_all.fasta \
  --output-path /scratch/rck80079/BOEM/REAGO_out/wc98/seqs.qza \
  --type 'FeatureData[Sequence]'



####If you have too many representative sequences and you are getting a Maaft error, align your sequences in maaft before importing to Qiime as follows:

#qiime tools import \
  #--input-path aligned-sequences.fna \
  #--output-path aligned-sequences.qza \
  #--type 'FeatureData[AlignedSequence]'



# Classify Taxonomy of Sequences
# You need to download the classifier of your preference at https://github.com/qiime2/docs/blob/master/source/data-resources.rst
# alternatively, see https://docs.qiime2.org/2019.4/tutorials/feature-classifier/ to train your own classifier to fit your needs (i.e. for a different gene of interest)


#!/bin/bash

#SBATCH --partition=joye_p
#SBATCH --job-name=wc98
#SBATCH --ntasks=1
#SBATCH --time=10:00:00
#SBATCH --mem=200G

module load QIIME2/2020.11

 qiime feature-classifier classify-sklearn \
  --i-classifier /scratch/rck80079/BOEM/REAGO_out/silva-138-99-nb-classifier.qza\
  --i-reads /scratch/rck80079/BOEM/REAGO_out/wc98/seqs.qza \
  --o-classification /scratch/rck80079/BOEM/REAGO_out/wc98/taxonomy.qza \
  --verbose

# Make a visual taxonomy table
#!/bin/bash

#SBATCH --partition=joye_p
#SBATCH --job-name=wc98
#SBATCH --ntasks=1
#SBATCH --time=10:00:00
#SBATCH --mem=200G

module load QIIME2/2020.11
 qiime metadata tabulate \
  --m-input-file /scratch/rck80079/BOEM/REAGO_out/wc98/taxonomy.qza \
  --o-visualization /scratch/rck80079/BOEM/REAGO_out/wc98/taxonomy.qzv


# Build Tree

#!/bin/bash

#SBATCH --partition=joye_p
#SBATCH --job-name=wc98
#SBATCH --ntasks=1
#SBATCH --time=10:00:00
#SBATCH --mem=200G
module load QIIME2/2020.11

qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences /scratch/rck80079/BOEM/REAGO_out/wc98/seqs.qza \
  --o-alignment /scratch/rck80079/BOEM/REAGO_out/wc98/aligned-rep-seqs.qza \
  --o-masked-alignment /scratch/rck80079/BOEM/REAGO_out/wc98/masked-aligned-rep-seqs.qza \
  --o-tree /scratch/rck80079/BOEM/REAGO_out/wc98/unrooted-tree.qza \
  --o-rooted-tree /scratch/rck80079/BOEM/REAGO_out/wc98/rooted-tree.qza


# Export Tree and Taxa files from Qiime2
#!/bin/bash

#SBATCH --partition=joye_p
#SBATCH --job-name=wc98
#SBATCH --ntasks=1
#SBATCH --time=10:00:00
#SBATCH --mem=200G

module load QIIME2/2020.11
qiime tools export --input-path /scratch/rck80079/BOEM/REAGO_out/wc98/taxonomy.qza --output-path /scratch/rck80079/BOEM/REAGO_out/wc98
qiime tools export --input-path /scratch/rck80079/BOEM/REAGO_out/wc98/rooted-tree.qza --output-path /scratch/rck80079/BOEM/REAGO_out/wc98



# Optional- if you want to keep only unique taxa, make a file using excel of only the unique IDs:
# i.e.
#IDs  Taxa
#gene_1 Odinarchaeota

#Remove duplicates based on the taxa column and then you have only the gene IDs you want to keep.
# From there, run these lines of code to keep only those unique IDs in the fasta file, then you can go back and rebuild the tree using the new fasta file.

grep -A1 -if <(tr -d '\r' < taxonomy_unique_ids.ids) wc98.fasta >> unique_wc98_all.fasta
sed -i '/--/d' unique_wc98_all.fasta




# Format for ITOL
# be sure to take a look and make sure your identifiers match in both of these files. If they don't, you've done something wrong and won't
# see taxonomy on your tree.

echo $'LABELS\nSEPARATOR TAB\nDATA' > /scratch/rck80079/BOEM/REAGO_out/wc98/taxonomy.txt
sed "1d" /scratch/rck80079/BOEM/REAGO_out/wc98/taxonomy.tsv | cut -f1,2 >> /scratch/rck80079/BOEM/REAGO_out/wc98/taxonomy.txt
