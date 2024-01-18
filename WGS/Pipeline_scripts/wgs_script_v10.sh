#!/bin/sh

#SBATCH --ntasks=1
#SBATCH --job-name=covWGS
#SBATCH --time=3:00:00
#SBATCH --cpus-per-task=24
#SBATCH --mem=24G
#SBATCH --output=covWGS_%A_%a.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hartman.l@wehi.edu.au

export TMPDIR=$(mktemp -d --tmpdir=/vast/scratch/users/$USER)

module load samtools
module load anaconda3

FILEPATH=`pwd -P`
wait
rename - _ $FILEPATH/raw_seqs/*.gz
rename - _ $FILEPATH/raw_seqs/*.gz
rename - _ $FILEPATH/raw_seqs/*.gz
rm $FILEPATH/sample_list.txt
ls $FILEPATH/raw_seqs/ | awk 'NR % 2 == 0' | sed 's/-/_/g' | sed 's/_S.*$//' > $FILEPATH/sample_list.txt
wait
readarray -t samples < $FILEPATH/sample_list.txt
sed 's/\r$//' -i /stornext/General/data/academic/lab_jex/QIASEQ_WGS_MINISEQ/mapping_table2.txt
rm $FILEPATH/wgs_run_summary.txt

# -----------------------------------------

for j in ${samples[@]}
do

# -----------------------------------------

rm -r $FILEPATH/${j}
mkdir $FILEPATH/${j}
mkdir $FILEPATH/${j}/inter
mkdir $FILEPATH/${j}/output
mkdir $FILEPATH/${j}/binning

# -----------------------------------------
# GET LINEAGE DATA
# -----------------------------------------

conda activate freyja

# trim & merge.
fastp --thread 16 -i $FILEPATH/raw_seqs/*${j}_*R1*.fastq.gz -I $FILEPATH/raw_seqs/*${j}_*R2*.fastq.gz -q 15 -u 40 -l 30 --trim_front1 15 --trim_front2 15 --cut_right --cut_window_size 20 --cut_mean_quality 25 --correction --trim_poly_g --merge --merged_out $FILEPATH/${j}/inter/${j}_paired.fastq

wait

# save read count data.
mv $FILEPATH/fastp.html $FILEPATH/${j}/output/${j}_fastp.html
mv $FILEPATH/fastp.json $FILEPATH/${j}/output/${j}_fastp.json
head -n 6 $FILEPATH/${j}/output/${j}_fastp.json | tail -n 1 | sed -e 's/.*:\(.*\),.*/\1/' > $FILEPATH/${j}/output/${j}_raw_reads.txt
head -n 17 $FILEPATH/${j}/output/${j}_fastp.json | tail -n 1 | sed -e 's/.*:\(.*\),.*/\1/' > $FILEPATH/${j}/output/${j}_filtered_reads.txt

wait

# align to ref & sort.
bwa-mem2 mem -t 24 /stornext/General/data/academic/lab_jex/QIASEQ_WGS_MINISEQ/cov2_genome/cov2.fa $FILEPATH/${j}/inter/${j}_paired.fastq | samtools sort - > $FILEPATH/${j}/inter/${j}.bam

wait

# call variants with ivar (requires piped input from samtools mpileup).
samtools mpileup -aa -d 0 -E -Q 20 --reference /stornext/General/data/academic/lab_jex/QIASEQ_WGS_MINISEQ/cov2_genome/cov2.fa $FILEPATH/${j}/inter/${j}.bam | ivar variants -t 0.02 -q 20 -p $FILEPATH/${j}/output/${j}_variants.tsv -r /stornext/General/data/academic/lab_jex/QIASEQ_WGS_MINISEQ/cov2_genome/cov2.fa -g /stornext/General/data/academic/lab_jex/QIASEQ_WGS_MINISEQ/cov2_genome/Sars_cov_2.ASM985889v3.101.gff3

wait

# create depth file for freyja.
samtools mpileup -aa -d 0 -E -Q 20 --reference /stornext/General/data/academic/lab_jex/QIASEQ_WGS_MINISEQ/cov2_genome/cov2.fa $FILEPATH/${j}/inter/${j}.bam | cut -f1-4 > $FILEPATH/${j}/inter/${j}_depth.txt

wait

# identify lineages with freyja.
freyja demix $FILEPATH/${j}/output/${j}_variants.tsv $FILEPATH/${j}/inter/${j}_depth.txt --output $FILEPATH/${j}/output/${j}_demix.tsv

wait

conda deactivate

# get lineage names & percentages.
cat $FILEPATH/${j}/output/${j}_demix.tsv | head -n -3 | tail -n -1 | cut -f2- > $FILEPATH/${j}/inter/${j}_lineages1.txt
cat $FILEPATH/${j}/output/${j}_demix.tsv | head -n -2 | tail -n -1 | cut -f2- > $FILEPATH/${j}/inter/${j}_percents1.txt

wait

# add line breaks.
sed -i 's/ /\n/g' $FILEPATH/${j}/inter/${j}_lineages1.txt
sed -i 's/ /\n/g' $FILEPATH/${j}/inter/${j}_percents1.txt

wait

# multiply percentages by 100 & round to one decimal place.
awk '{print $1*100}' $FILEPATH/${j}/inter/${j}_percents1.txt | xargs printf "%.1f\n" > $FILEPATH/${j}/inter/${j}_percents2.txt

wait

# count the number of percentages ≥1.0
gte1=`awk '($1>0.9){ ++count } END{ print count }' $FILEPATH/${j}/inter/${j}_percents2.txt`

wait 

# add parentheses around percentages.
sed -i 's/^/(/g' $FILEPATH/${j}/inter/${j}_percents2.txt
sed -i 's/$/)/g' $FILEPATH/${j}/inter/${j}_percents2.txt

wait

# join names and percentages.
pr -tm -s' ' $FILEPATH/${j}/inter/${j}_lineages1.txt $FILEPATH/${j}/inter/${j}_percents2.txt > $FILEPATH/${j}/inter/${j}_lineages2.txt

wait

# concatenate names & percentages.
cat $FILEPATH/${j}/inter/${j}_lineages2.txt | head -n $gte1 | sed ':a;N;$!ba;s/\n/, /g' > $FILEPATH/${j}/output/${j}_lineages3.txt

wait

# -----------------------------------------
# BIN THE IDENTIFIED LINEAGES
# -----------------------------------------

# replace the identified lineage names with the bin names.
while read x ; do
while read a b ; do
if [ $a = $x ] ; then
x=$b
break
fi
done < /stornext/General/data/academic/lab_jex/QIASEQ_WGS_MINISEQ/mapping_table2.txt
echo $x
done < $FILEPATH/${j}/inter/${j}_lineages1.txt > $FILEPATH/${j}/binning/${j}_bin1.txt

wait

# join bin names and percentages.
pr -tm -s' ' $FILEPATH/${j}/binning/${j}_bin1.txt $FILEPATH/${j}/inter/${j}_percents1.txt > $FILEPATH/${j}/binning/${j}_bin2.txt

wait

# get totals for each bin (note the $ anchor reqd on some search patterns due to naming conflicts).
# AY
awk '$1 ~ /AY/ {sum += $2} END {print sum}' $FILEPATH/${j}/binning/${j}_bin2.txt | awk '{print $1*100}' | xargs printf "%.1f\n" > $FILEPATH/${j}/binning/${j}_AY.txt
# BA.1
awk '$1 ~ /BA1/ {sum += $2} END {print sum}' $FILEPATH/${j}/binning/${j}_bin2.txt | awk '{print $1*100}' | xargs printf "%.1f\n" > $FILEPATH/${j}/binning/${j}_BA1.txt
# BA.2
awk '$1 ~ /BA2$/ {sum += $2} END {print sum}' $FILEPATH/${j}/binning/${j}_bin2.txt | awk '{print $1*100}' | xargs printf "%.1f\n" > $FILEPATH/${j}/binning/${j}_BA2.txt
# BA.3
awk '$1 ~ /BA3/ {sum += $2} END {print sum}' $FILEPATH/${j}/binning/${j}_bin2.txt | awk '{print $1*100}' | xargs printf "%.1f\n" > $FILEPATH/${j}/binning/${j}_BA3.txt
# BA.4
awk '$1 ~ /BA4/ {sum += $2} END {print sum}' $FILEPATH/${j}/binning/${j}_bin2.txt | awk '{print $1*100}' | xargs printf "%.1f\n" > $FILEPATH/${j}/binning/${j}_BA4.txt
# BA.5
awk '$1 ~ /BA5/ {sum += $2} END {print sum}' $FILEPATH/${j}/binning/${j}_bin2.txt | awk '{print $1*100}' | xargs printf "%.1f\n" > $FILEPATH/${j}/binning/${j}_BA5.txt
# BA.2.75.x
awk '$1 ~ /BA275x/ {sum += $2} END {print sum}' $FILEPATH/${j}/binning/${j}_bin2.txt | awk '{print $1*100}' | xargs printf "%.1f\n" > $FILEPATH/${j}/binning/${j}_BA275x.txt
# BN
awk '$1 ~ /BN/ {sum += $2} END {print sum}' $FILEPATH/${j}/binning/${j}_bin2.txt | awk '{print $1*100}' | xargs printf "%.1f\n" > $FILEPATH/${j}/binning/${j}_BN.txt
# BQ.1.x
awk '$1 ~ /BQ1x/ {sum += $2} END {print sum}' $FILEPATH/${j}/binning/${j}_bin2.txt | awk '{print $1*100}' | xargs printf "%.1f\n" > $FILEPATH/${j}/binning/${j}_BQ1x.txt
# BR.2
awk '$1 ~ /BR2/ {sum += $2} END {print sum}' $FILEPATH/${j}/binning/${j}_bin2.txt | awk '{print $1*100}' | xargs printf "%.1f\n" > $FILEPATH/${j}/binning/${j}_BR2.txt
# CH.1.1
awk '$1 ~ /CH11/ {sum += $2} END {print sum}' $FILEPATH/${j}/binning/${j}_bin2.txt | awk '{print $1*100}' | xargs printf "%.1f\n" > $FILEPATH/${j}/binning/${j}_CH11.txt
# XBB
awk '$1 ~ /XBB$/ {sum += $2} END {print sum}' $FILEPATH/${j}/binning/${j}_bin2.txt | awk '{print $1*100}' | xargs printf "%.1f\n" > $FILEPATH/${j}/binning/${j}_XBB.txt
# XBB.1.5
awk '$1 ~ /XBB15/ {sum += $2} END {print sum}' $FILEPATH/${j}/binning/${j}_bin2.txt | awk '{print $1*100}' | xargs printf "%.1f\n" > $FILEPATH/${j}/binning/${j}_XBB15.txt
# XBB.1.9.1
awk '$1 ~ /XBB191/ {sum += $2} END {print sum}' $FILEPATH/${j}/binning/${j}_bin2.txt | awk '{print $1*100}' | xargs printf "%.1f\n" > $FILEPATH/${j}/binning/${j}_XBB191.txt
# XBB.1.16
awk '$1 ~ /XBB116/ {sum += $2} END {print sum}' $FILEPATH/${j}/binning/${j}_bin2.txt | awk '{print $1*100}' | xargs printf "%.1f\n" > $FILEPATH/${j}/binning/${j}_XBB116.txt
# XBC
awk '$1 ~ /XBC/ {sum += $2} END {print sum}' $FILEPATH/${j}/binning/${j}_bin2.txt | awk '{print $1*100}' | xargs printf "%.1f\n" > $FILEPATH/${j}/binning/${j}_XBC.txt
# XBF
awk '$1 ~ /XBF/ {sum += $2} END {print sum}' $FILEPATH/${j}/binning/${j}_bin2.txt | awk '{print $1*100}' | xargs printf "%.1f\n" > $FILEPATH/${j}/binning/${j}_XBF.txt
# Recombinant
awk '$1 ~ /Recombinant/ {sum += $2} END {print sum}' $FILEPATH/${j}/binning/${j}_bin2.txt | awk '{print $1*100}' | xargs printf "%.1f\n" > $FILEPATH/${j}/binning/${j}_Recombinant.txt
# XBB.1.9.2
awk '$1 ~ /XBB192/ {sum += $2} END {print sum}' $FILEPATH/${j}/binning/${j}_bin2.txt | awk '{print $1*100}' | xargs printf "%.1f\n" > $FILEPATH/${j}/binning/${j}_XBB192.txt
# XCC
awk '$1 ~ /XCC/ {sum += $2} END {print sum}' $FILEPATH/${j}/binning/${j}_bin2.txt | awk '{print $1*100}' | xargs printf "%.1f\n" > $FILEPATH/${j}/binning/${j}_XCC.txt
# XBB.2.3
awk '$1 ~ /XBB23/ {sum += $2} END {print sum}' $FILEPATH/${j}/binning/${j}_bin2.txt | awk '{print $1*100}' | xargs printf "%.1f\n" > $FILEPATH/${j}/binning/${j}_XBB23.txt
# EG.5.1
awk '$1 ~ /EG51/ {sum += $2} END {print sum}' $FILEPATH/${j}/binning/${j}_bin2.txt | awk '{print $1*100}' | xargs printf "%.1f\n" > $FILEPATH/${j}/binning/${j}_EG51.txt
# HK.3
awk '$1 ~ /HK3/ {sum += $2} END {print sum}' $FILEPATH/${j}/binning/${j}_bin2.txt | awk '{print $1*100}' | xargs printf "%.1f\n" > $FILEPATH/${j}/binning/${j}_HK3.txt
# HV.1
awk '$1 ~ /HV1/ {sum += $2} END {print sum}' $FILEPATH/${j}/binning/${j}_bin2.txt | awk '{print $1*100}' | xargs printf "%.1f\n" > $FILEPATH/${j}/binning/${j}_HV1.txt
# DV.7.1
awk '$1 ~ /DV71/ {sum += $2} END {print sum}' $FILEPATH/${j}/binning/${j}_bin2.txt | awk '{print $1*100}' | xargs printf "%.1f\n" > $FILEPATH/${j}/binning/${j}_DV71.txt
# BA.2.86
awk '$1 ~ /BA286/ {sum += $2} END {print sum}' $FILEPATH/${j}/binning/${j}_bin2.txt | awk '{print $1*100}' | xargs printf "%.1f\n" > $FILEPATH/${j}/binning/${j}_BA286.txt
# JN.1
awk '$1 ~ /JN1/ {sum += $2} END {print sum}' $FILEPATH/${j}/binning/${j}_bin2.txt | awk '{print $1*100}' | xargs printf "%.1f\n" > $FILEPATH/${j}/binning/${j}_JN1.txt
# JD.1.1
awk '$1 ~ /JD11/ {sum += $2} END {print sum}' $FILEPATH/${j}/binning/${j}_bin2.txt | awk '{print $1*100}' | xargs printf "%.1f\n" > $FILEPATH/${j}/binning/${j}_JD11.txt
# GW.5.1.1
awk '$1 ~ /GW511/ {sum += $2} END {print sum}' $FILEPATH/${j}/binning/${j}_bin2.txt | awk '{print $1*100}' | xargs printf "%.1f\n" > $FILEPATH/${j}/binning/${j}_GW511.txt


wait

# calculate other ("Unassigned").
cat \
$FILEPATH/${j}/binning/${j}_AY.txt \
$FILEPATH/${j}/binning/${j}_BA1.txt \
$FILEPATH/${j}/binning/${j}_BA2.txt \
$FILEPATH/${j}/binning/${j}_BA3.txt \
$FILEPATH/${j}/binning/${j}_BA4.txt \
$FILEPATH/${j}/binning/${j}_BA5.txt \
$FILEPATH/${j}/binning/${j}_BA275x.txt \
$FILEPATH/${j}/binning/${j}_BN.txt \
$FILEPATH/${j}/binning/${j}_BQ1x.txt \
$FILEPATH/${j}/binning/${j}_BR2.txt \
$FILEPATH/${j}/binning/${j}_CH11.txt \
$FILEPATH/${j}/binning/${j}_XBB.txt \
$FILEPATH/${j}/binning/${j}_XBB15.txt \
$FILEPATH/${j}/binning/${j}_XBB191.txt \
$FILEPATH/${j}/binning/${j}_XBB116.txt \
$FILEPATH/${j}/binning/${j}_XBC.txt \
$FILEPATH/${j}/binning/${j}_XBF.txt \
$FILEPATH/${j}/binning/${j}_Recombinant.txt \
$FILEPATH/${j}/binning/${j}_XBB192.txt \
$FILEPATH/${j}/binning/${j}_XCC.txt \
$FILEPATH/${j}/binning/${j}_XBB23.txt \
$FILEPATH/${j}/binning/${j}_EG51.txt \
$FILEPATH/${j}/binning/${j}_HK3.txt \
$FILEPATH/${j}/binning/${j}_HV1.txt \
$FILEPATH/${j}/binning/${j}_DV71.txt \
$FILEPATH/${j}/binning/${j}_BA286.txt \
$FILEPATH/${j}/binning/${j}_JN1.txt \
$FILEPATH/${j}/binning/${j}_JD11.txt \
$FILEPATH/${j}/binning/${j}_GW511.txt \
| awk '{a+=$0}END{print 100-a}' > $FILEPATH/${j}/binning/${j}_other.txt

# -----------------------------------------
# GET DEPTH & COVERAGE DATA
# -----------------------------------------

# get mean depth & coverage & round to one decimal place.
samtools coverage $FILEPATH/${j}/inter/${j}.bam > $FILEPATH/${j}/inter/${j}_coverage.txt
awk '{print $6}' $FILEPATH/${j}/inter/${j}_coverage.txt | tail -n 1 | xargs printf "%.1f\n" > $FILEPATH/${j}/output/${j}_cov_percent.txt
awk '{print $7}' $FILEPATH/${j}/inter/${j}_coverage.txt | tail -n 1 | xargs printf "%.0f\n" > $FILEPATH/${j}/output/${j}_mean_depth.txt

# -----------------------------------------
# SAVE SAMPLE SUMMARY FILES
# -----------------------------------------

echo RNA-${j} > $FILEPATH/${j}/inter/${j}_sample_name.txt

paste \
$FILEPATH/${j}/inter/${j}_sample_name.txt \
$FILEPATH/${j}/output/${j}_filtered_reads.txt \
$FILEPATH/${j}/binning/${j}_AY.txt \
$FILEPATH/${j}/binning/${j}_BA1.txt \
$FILEPATH/${j}/binning/${j}_BA2.txt \
$FILEPATH/${j}/binning/${j}_BA3.txt \
$FILEPATH/${j}/binning/${j}_BA4.txt \
$FILEPATH/${j}/binning/${j}_BA5.txt \
$FILEPATH/${j}/binning/${j}_BA275x.txt \
$FILEPATH/${j}/binning/${j}_BN.txt \
$FILEPATH/${j}/binning/${j}_BQ1x.txt \
$FILEPATH/${j}/binning/${j}_BR2.txt \
$FILEPATH/${j}/binning/${j}_CH11.txt \
$FILEPATH/${j}/binning/${j}_XBB.txt \
$FILEPATH/${j}/binning/${j}_XBB15.txt \
$FILEPATH/${j}/binning/${j}_XBB191.txt \
$FILEPATH/${j}/binning/${j}_XBB116.txt \
$FILEPATH/${j}/binning/${j}_XBC.txt \
$FILEPATH/${j}/binning/${j}_XBF.txt \
$FILEPATH/${j}/binning/${j}_Recombinant.txt \
$FILEPATH/${j}/binning/${j}_XBB192.txt \
$FILEPATH/${j}/binning/${j}_XCC.txt \
$FILEPATH/${j}/binning/${j}_XBB23.txt \
$FILEPATH/${j}/binning/${j}_EG51.txt \
$FILEPATH/${j}/binning/${j}_HK3.txt \
$FILEPATH/${j}/binning/${j}_HV1.txt \
$FILEPATH/${j}/binning/${j}_DV71.txt \
$FILEPATH/${j}/binning/${j}_BA286.txt \
$FILEPATH/${j}/binning/${j}_JN1.txt \
$FILEPATH/${j}/binning/${j}_JD11.txt \
$FILEPATH/${j}/binning/${j}_GW511.txt \
$FILEPATH/${j}/binning/${j}_other.txt \
$FILEPATH/${j}/output/${j}_lineages3.txt \
$FILEPATH/${j}/output/${j}_cov_percent.txt \
$FILEPATH/${j}/output/${j}_mean_depth.txt \
> $FILEPATH/${j}/output/${j}_summary.txt

wait

cat $FILEPATH/${j}/output/${j}_summary.txt >> $FILEPATH/wgs_run_summary.txt

wait

# -----------------------------------------

done

# -----------------------------------------

# add header.
sed -i '1i Accession\tFiltered\tMapping_AY\tMapping_BA1\tMapping_BA2\tMapping_BA3\tMapping_BA4\tMapping_BA5\tMapping_BA275x\tMapping_BN\tMapping_BQ1x\tMapping_BR2\tMapping_CH11\tMapping_XBB\tMapping_XBB15\tMapping_XBB191\tMapping_XBB116\tMapping_XBC\tMapping_XBF\tMapping_Recombinant\tMapping_XXB192\tMapping_XCC\tMapping_XBB23\tMapping_EG51\tMapping_HK3\tMapping_HV1\tMapping_DV71\tMapping_BA286\tMapping_JN1\tMappingJD11\tMappingGW511\tMapping_Unassigned\tLineagesDetected\tCoveragePercent\tDepth' $FILEPATH/wgs_run_summary.txt

# -----------------------------------------
# MAKE LIST OF NON-MAPPED LINEAGES
# -----------------------------------------

rm $FILEPATH/all_lineages.txt

# make list of all lineages detected by freyja.
for j in ${samples[@]}
do
cat $FILEPATH/${j}/binning/${j}_bin1.txt >> $FILEPATH/all_lineages.txt
done

wait

# keep only one entry for each detected lineage.
sort -u $FILEPATH/all_lineages.txt > $FILEPATH/all_unique_lineages.txt

wait

# make list of CURRENT mapping table target names.
awk '{print $1}' /stornext/General/data/academic/lab_jex/QIASEQ_WGS_MINISEQ/mapping_table2.txt > $FILEPATH/mapable_lineages.txt

wait

# identify lineages that are missing from the mapping table.
sort $FILEPATH/all_unique_lineages.txt > $FILEPATH/all_unique_lineages_sorted.txt
sort $FILEPATH/mapable_lineages.txt > $FILEPATH/mapable_lineages_sorted.txt
comm -13 $FILEPATH/mapable_lineages_sorted.txt $FILEPATH/all_unique_lineages_sorted.txt > $FILEPATH/lineages_to_add.txt

# Note – Our system will not run the above task using command substitution :(
# comm -13 <(sort $FILEPATH/mapable_lineages.txt) <(sort $FILEPATH/all_unique_lineages.txt) > $FILEPATH/lineages_to_add.txt

wait

rm $FILEPATH/all_lineages.txt
rm $FILEPATH/all_unique_lineages.txt
rm $FILEPATH/all_unique_lineages_sorted.txt
rm $FILEPATH/mapable_lineages.txt
rm $FILEPATH/mapable_lineages_sorted.txt

# -----------------------------------------
# MAKE LIST OF SUB-LINEAGE BIN ASSIGNMENTS
# -----------------------------------------

sed 's/\r$//' -i /stornext/General/data/academic/lab_jex/QIASEQ_WGS_MINISEQ/mapping_table_dev3.txt

rm $FILEPATH/pre_binned_lineages_FINAL.txt

# -----------------------------------------

for j in ${samples[@]}
do

# -----------------------------------------

# count the sub-lineages that are ≥1.0%

gte1=`awk '($1>=0.01){ ++count } END{ print count }' $FILEPATH/${j}/inter/${j}_percents1.txt`


# write the sub-lineages that are ≥1.0% to a file - assumes they are in the descending order of abundance
cat $FILEPATH/${j}/inter/${j}_lineages1.txt | head -n $gte1 > $FILEPATH/${j}/inter/${j}_lineages_gte1.txt

# write the percents for the sub-lineages that are ≥1.0% to a file - assumes their order is the same as the sub-lineages
cat $FILEPATH/${j}/inter/${j}_percents1.txt | head -n $gte1 > $FILEPATH/${j}/inter/${j}_percents_gte1.txt

# map the sub-lineage names to their parental clades and lineages.
while read x ; do
while read a b ; do
if [ $a = $x ] ; then
x=$b
break
fi
done < /stornext/General/data/academic/lab_jex/QIASEQ_WGS_MINISEQ/mapping_table_dev3.txt
echo $x
done < $FILEPATH/${j}/inter/${j}_lineages_gte1.txt > $FILEPATH/${j}/binning/${j}_bin_lin1.txt

wait

# join clade/lineage and percent.
pr -tm -s' ' $FILEPATH/${j}/inter/${j}_percents_gte1.txt $FILEPATH/${j}/binning/${j}_bin_lin1.txt > $FILEPATH/${j}/binning/${j}_bin_lin2.txt

wait

# add sample name.
awk -v accession=${j} '{print accession,$0}' $FILEPATH/${j}/binning/${j}_bin_lin2.txt > $FILEPATH/${j}/binning/${j}_bin_lin3.txt

wait

# combine all results.
cat $FILEPATH/${j}/binning/${j}_bin_lin3.txt >> $FILEPATH/pre_binned_lineages_FINAL.txt

wait

# -----------------------------------------

done

# -----------------------------------------

# replace "z" separators with spaces.
sed -i 's/z/ /g' $FILEPATH/pre_binned_lineages_FINAL.txt

wait

# trim sub-lineage names after underscores in the gpLineage and pLineage columns
awk '{sub(/_.*/,"",$7)} 1' $FILEPATH/pre_binned_lineages_FINAL.txt | awk '{sub(/_.*/,"",$8)} 1' > $FILEPATH/pre_binned_lineages_FINAL2.txt

wait

mv $FILEPATH/pre_binned_lineages_FINAL2.txt $FILEPATH/pre_binned_lineages_FINAL.txt

# -----------------------------------------

rm -rf $TMPDIR

# -----------------------------------------
# THE END
# -----------------------------------------
