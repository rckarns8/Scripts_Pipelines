cd OC26_Bins

mkdir AT18-02_MUC19-8_16-21_Bin_1
cd AT18-02_MUC19-8_16-21_Bin_1
mkdir coverage
cut -f1,4 ~/Desktop/depth/AT18-02_MUC19-8_16-21_Bin_1-sorted.bam-depth.txt > coverage/AT18-02_MUC19_8_16-21_Bin_1_4_1_cov
rm(list=ls())

library(mmgenome2)
setwd("~/AT18-02_MUC19-8_16-21_Bin_1")

mm <- mmload(
  assembly = "/Users/Rachael Karns/Desktop/FNA/AT18_02_MUC19_8_16_21_Bin_1.fna.fna",
  coverage = "coverage/",
  verbose = TRUE,
  kmer_pca = FALSE,
  kmer_BH_tSNE = FALSE
)
mm


mmplot_pairs(mm,
             variables = c("cov_AT1802_MUC19_8_16_21_Bin_1_4_1",
                           "cov_test",
                           "cov_test2",
                           "cov_test3"),
             x_scale = "log10",
             y_scale = "log10",
             alpha = 0.1,
             size_scale = 0.4,
             textsize = 2)

