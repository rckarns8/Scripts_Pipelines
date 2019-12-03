#PBS -S /bin/bash
#PBS -q batch
#PBS -N Step_3
#PBS -l nodes=1:ppn=1
#PBS -l walltime=20:00:00
#PBS -l mem=10g
#PBS -j oe
BASEDIR=/work/sbjlab/rck/Sed_Assemblies
cd $BASEDIR


qlogin

module load R/3.4.4-foss-2016b-X11-20160819-GACRC
R
# at this point, be sure to install the packages d3heatmap and htmlwidgets. Instructions for install are here: https://blog.rstudio.com/2015/06/24/d3heatmap/
install.packages("d3heatmap")
install.packages("htmlwidgets")
R --no-save < /work/sbjlab/rck/Scripts_Pipelines/heatmap.r
scp *.html /work/sbjlab/rck/Scripts_Pipelines/
