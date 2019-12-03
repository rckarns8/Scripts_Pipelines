#PBS -S /bin/bash
#PBS -q highmem_q
#PBS -N Step_3
#PBS -l nodes=1:ppn=1
#PBS -l walltime=240:00:00
#PBS -l mem=100g
#PBS -j oe


BASEDIR=/work/sbjlab/rck/Sed_Assemblies/AT26-13-91/
cd $BASEDIR

module load MetaPhlAn2/2.7.8-foss-2016b-Python-2.7.14
metaphlan2.py AT26-13_91_Sed_S51_L003_R1_001.fa  --input_type fasta > AT26-13-91_R1_profile.txt
metaphlan2.py AT26-13_91_Sed_S51_L003_R2_001.fa  --input_type fasta > AT26-13-91_R2_profile.txt
cat AT26-13-91_R1_profile.txt AT26-13-91_R2_profile.txt > 91_Metaphlan_profile.txt



awk '!h[$4,$NF]++ { print $4, $NF }' FS="|" 91_Metaphlan_profile.txt > 91_temp.txt
awk '!h[$1,$NF]++ { print $1, $NF }' FS=" " 91_temp.txt  > 91_clean_profile.txt
ls -1 | grep "^o" 91_clean_profile.txt > 91_order_profile.txt
sed 's/ \+/,/g' 91_order_profile.txt > 91_order_profiles.csv
awk -F, '!seen[$1]++' 91_order_profiles.csv > 91_order_profile.csv
rm 91_temp.txt 91_clean_profile.txt  91_order_profiles.csv
