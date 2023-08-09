#!/bin/sh

#SBATCH --ntasks=1
#SBATCH --job-name=covWGS
#SBATCH --time=12:00:00
#SBATCH --cpus-per-task=24
#SBATCH --mem=24G
#SBATCH --output=covWGS_%A_%a.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hartman.l@wehi.edu.au

export TMPDIR=$(mktemp -d --tmpdir=/vast/scratch/users/$USER)

module load samtools
module load anaconda3

FILEPATH=`pwd -P`
sed 's/\r$//' -i /stornext/General/data/academic/lab_jex/QIASEQ_WGS_MINISEQ/wgs_codons.txt
sed 's/\r$//' -i /stornext/General/data/academic/lab_jex/QIASEQ_WGS_MINISEQ/wgs_lineages.txt
sed 's/\r$//' -i $FILEPATH/meta/sample_list.txt
readarray -t samples < $FILEPATH/meta/sample_list.txt

# -----------------------------------------

for j in ${samples[@]}
do

# -----------------------------------------

rm -r $FILEPATH/${j}
mkdir $FILEPATH/${j}
rm -r $FILEPATH/${j}/inter
mkdir $FILEPATH/${j}/inter
rm -r $FILEPATH/${j}/output
mkdir $FILEPATH/${j}/output
rm -r $FILEPATH/${j}/snps
mkdir $FILEPATH/${j}/snps
rm -r $FILEPATH/${j}/binning
mkdir $FILEPATH/${j}/binning

# -----------------------------------------
# GET LINEAGE DATA
# -----------------------------------------

conda activate lofreq

# trim & merge.
fastp --thread 24 -i $FILEPATH/raw_seqs/*${j}*R1*.fastq.gz -I $FILEPATH/raw_seqs/*${j}*R2*.fastq.gz -q 15 -u 40 -l 30 --trim_front1 15 --trim_front2 15 --cut_right --cut_window_size 20 --cut_mean_quality 25 --correction --trim_poly_g --merge --merged_out $FILEPATH/${j}/inter/${j}_paired.fastq

wait

# save read count data.
mv $FILEPATH/fastp.html $FILEPATH/${j}/output/${j}_fastp.html
mv $FILEPATH/fastp.json $FILEPATH/${j}/output/${j}_fastp.json
head -n 6 $FILEPATH/${j}/output/${j}_fastp.json | tail -n 1 | sed -e 's/.*:\(.*\),.*/\1/' > $FILEPATH/${j}/output/${j}_raw_reads.txt
head -n 17 $FILEPATH/${j}/output/${j}_fastp.json | tail -n 1 | sed -e 's/.*:\(.*\),.*/\1/' > $FILEPATH/${j}/output/${j}_filtered_reads.txt

conda deactivate
wait
conda activate freyja

# align to ref & sort.
bwa-mem2 mem -t 24 /stornext/General/data/academic/lab_jex/QIASEQ_WGS_MINISEQ/cov2_genome/cov2.fa $FILEPATH/${j}/inter/${j}_paired.fastq | samtools sort - > $FILEPATH/${j}/inter/${j}.bam

wait

# call variants with ivar (requires piped input from samtools mpileup).
samtools mpileup -aa -d 0 -E -Q 20 --reference /stornext/General/data/academic/lab_jex/QIASEQ_WGS_MINISEQ/cov2_genome/cov2.fa $FILEPATH/${j}/inter/${j}.bam | ivar variants -t 0.02 -q 20 -p $FILEPATH/${j}/output/${j}_variants.tsv -r /stornext/General/data/academic/lab_jex/QIASEQ_WGS_MINISEQ/cov2_genome/cov2.fa -g /stornext/General/data/academic/lab_jex/QIASEQ_WGS_MINISEQ/cov2_genome/Sars_cov_2.ASM985889v3.101.gff3

wait

# create depth file for freyja.
samtools mpileup -aa -d 0 -E -Q 20 --reference /stornext/General/data/academic/lab_jex/QIASEQ_WGS_MINISEQ/cov2_genome/cov2.fa $FILEPATH/${j}/inter/${j}.bam | cut -f1-4 > $FILEPATH/${j}/inter/${j}_depth.txt

wait

# identify lneages with freyja.
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
done < /stornext/General/data/academic/lab_jex/QIASEQ_WGS_MINISEQ/wgs_lineages.txt
echo $x
done < $FILEPATH/${j}/inter/${j}_lineages1.txt > $FILEPATH/${j}/binning/${j}_bin1.txt

wait

# join bin names and percentages.
pr -tm -s' ' $FILEPATH/${j}/binning/${j}_bin1.txt $FILEPATH/${j}/inter/${j}_percents1.txt > $FILEPATH/${j}/binning/${j}_bin2.txt

wait

# get totals for each bin (note the $ anchor reqd on some search patterns).
awk '$1 ~ /BA2$/ {sum += $2} END {print sum}' $FILEPATH/${j}/binning/${j}_bin2.txt | awk '{print $1*100}' | xargs printf "%.1f\n" > $FILEPATH/${j}/binning/${j}_BA2.txt
awk '$1 ~ /BA45/ {sum += $2} END {print sum}' $FILEPATH/${j}/binning/${j}_bin2.txt | awk '{print $1*100}' | xargs printf "%.1f\n" > $FILEPATH/${j}/binning/${j}_BA45.txt
awk '$1 ~ /BA2320/ {sum += $2} END {print sum}' $FILEPATH/${j}/binning/${j}_bin2.txt | awk '{print $1*100}' | xargs printf "%.1f\n" > $FILEPATH/${j}/binning/${j}_BA2320.txt
awk '$1 ~ /BA275x/ {sum += $2} END {print sum}' $FILEPATH/${j}/binning/${j}_bin2.txt | awk '{print $1*100}' | xargs printf "%.1f\n" > $FILEPATH/${j}/binning/${j}_BA275x.txt
awk '$1 ~ /BN/ {sum += $2} END {print sum}' $FILEPATH/${j}/binning/${j}_bin2.txt | awk '{print $1*100}' | xargs printf "%.1f\n" > $FILEPATH/${j}/binning/${j}_BN.txt
awk '$1 ~ /BQ1x/ {sum += $2} END {print sum}' $FILEPATH/${j}/binning/${j}_bin2.txt | awk '{print $1*100}' | xargs printf "%.1f\n" > $FILEPATH/${j}/binning/${j}_BQ1x.txt
awk '$1 ~ /BR2/ {sum += $2} END {print sum}' $FILEPATH/${j}/binning/${j}_bin2.txt | awk '{print $1*100}' | xargs printf "%.1f\n" > $FILEPATH/${j}/binning/${j}_BR2.txt
awk '$1 ~ /CH11/ {sum += $2} END {print sum}' $FILEPATH/${j}/binning/${j}_bin2.txt | awk '{print $1*100}' | xargs printf "%.1f\n" > $FILEPATH/${j}/binning/${j}_CH11.txt
awk '$1 ~ /XBB$/ {sum += $2} END {print sum}' $FILEPATH/${j}/binning/${j}_bin2.txt | awk '{print $1*100}' | xargs printf "%.1f\n" > $FILEPATH/${j}/binning/${j}_XBB.txt
awk '$1 ~ /XBB15/ {sum += $2} END {print sum}' $FILEPATH/${j}/binning/${j}_bin2.txt | awk '{print $1*100}' | xargs printf "%.1f\n" > $FILEPATH/${j}/binning/${j}_XBB15.txt
awk '$1 ~ /XBC/ {sum += $2} END {print sum}' $FILEPATH/${j}/binning/${j}_bin2.txt | awk '{print $1*100}' | xargs printf "%.1f\n" > $FILEPATH/${j}/binning/${j}_XBC.txt
awk '$1 ~ /XBF/ {sum += $2} END {print sum}' $FILEPATH/${j}/binning/${j}_bin2.txt | awk '{print $1*100}' | xargs printf "%.1f\n" > $FILEPATH/${j}/binning/${j}_XBF.txt

wait

# calculate other.
paste \
$FILEPATH/${j}/binning/${j}_BA2.txt \
$FILEPATH/${j}/binning/${j}_BA45.txt \
$FILEPATH/${j}/binning/${j}_BA2320.txt \
$FILEPATH/${j}/binning/${j}_BA275x.txt \
$FILEPATH/${j}/binning/${j}_BN.txt \
$FILEPATH/${j}/binning/${j}_BQ1x.txt \
$FILEPATH/${j}/binning/${j}_BR2.txt \
$FILEPATH/${j}/binning/${j}_CH11.txt \
$FILEPATH/${j}/binning/${j}_XBB.txt \
$FILEPATH/${j}/binning/${j}_XBB15.txt \
$FILEPATH/${j}/binning/${j}_XBC.txt \
$FILEPATH/${j}/binning/${j}_XBF.txt \
| awk '{print 100-($1+$2+$3+$4+$5+$6+$7+$8+$9+$10+$11+$12)}' > $FILEPATH/${j}/binning/${j}_other.txt

# -----------------------------------------
# GET DEPTH & COVERAGE DATA
# -----------------------------------------

# get mean depth & coverage & round to one decimal place.
samtools coverage $FILEPATH/${j}/inter/${j}.bam > $FILEPATH/${j}/inter/${j}_coverage.txt
awk '{print $6}' $FILEPATH/${j}/inter/${j}_coverage.txt | tail -n 1 | xargs printf "%.1f\n" > $FILEPATH/${j}/output/${j}_cov_percent.txt
awk '{print $7}' $FILEPATH/${j}/inter/${j}_coverage.txt | tail -n 1 | xargs printf "%.0f\n" > $FILEPATH/${j}/output/${j}_mean_depth.txt

# -----------------------------------------
# GET SNP DATA
# -----------------------------------------

# extract variant calls (setting a min depth threshold will exclude some deletions)
awk -F"\t" '{if ($12>0) print $2,$3,$4,$11}' $FILEPATH/${j}/output/${j}_variants.tsv | awk '{gsub (/-/, $2, $3); print $1,$2,$3,$4}' | awk  '{if (length($3) > 1) {print $1$3">"$2,$4} else {print $1$2">"$3,$4}}' > $FILEPATH/${j}/snps/${j}_snps1.tsv

wait

# remove header
sed '1d' -i $FILEPATH/${j}/snps/${j}_snps1.tsv

# retain one occcurence of any duplicates.
cat -n $FILEPATH/${j}/snps/${j}_snps1.tsv | sort -uk2 | sort -nk1 | cut -f2- > $FILEPATH/${j}/snps/${j}_snps2.tsv

wait

# rename mnps (only works if SNPs are on adjacent lines)
sed -zi 's/22577G>C\n22578G>A/s:G339H/ ; s/27382G>C\n27383A>T\n27384T>C/orf6:D61L/ ; s/22598A>G\n22599G>A/s:R346E/ ; s/22895G>C\n22896T>C/s:V445P/ ; s/23018T>C\n23019T>C/s:F486P/ ; s/10568C>T\n10569A>T/orf1b:H172F*/ ; s/10628C>T\n10629A>C/orf1b:Q192S*/ ; s/10628C>G\n10629A>T/orf1b:Q192V*/ ; s/10544C>A\n10545A>C/orf1b:H164N*/' $FILEPATH/${j}/snps/${j}_snps2.tsv

wait

# rename snps
unset codons
declare -A codons

wait

while IFS=';' read a b; do codons["$a"]="$b"; done < /stornext/General/data/academic/lab_jex/QIASEQ_WGS_MINISEQ/wgs_codons.txt

wait

while IFS=' ' read line; do
for i in ${!codons[@]}; do
line=${line//$i/${codons[$i]}}
done
echo $line
done < $FILEPATH/${j}/snps/${j}_snps2.tsv > $FILEPATH/${j}/snps/${j}_snps3.tsv

wait

# remove un-renamed snps
grep -v ">" $FILEPATH/${j}/snps/${j}_snps3.tsv > $FILEPATH/${j}/snps/${j}_snps4.tsv

wait

# multiply allele frequency values by 100
awk '{printf $1}; {printf " (%.1f)\n", $2*100}' $FILEPATH/${j}/snps/${j}_snps4.tsv > $FILEPATH/${j}/snps/${j}_snps5.tsv

wait

# separate snps by gene
unset genes
declare -a genes=(orf1a orf1b s orf3a e m orf6 orf7a orf7b orf8 orf10 n)

wait

for i in ${genes[@]}
do
grep '^'${i}':' $FILEPATH/${j}/snps/${j}_snps5.tsv > $FILEPATH/${j}/snps/${j}_${i}_snps.tsv
done

wait

# remove prefixes
for i in ${genes[@]}
do
sed -i 's/'${i}'://' $FILEPATH/${j}/snps/${j}_${i}_snps.tsv
done

wait

# replace line breaks with ", "
for i in ${genes[@]}
do
sed -i ':a;N;$!ba;s/\n/, /g' $FILEPATH/${j}/snps/${j}_${i}_snps.tsv
done

wait

# replace "+"
for i in ${genes[@]}
do
sed -i 's/+/-/g' $FILEPATH/${j}/snps/${j}_${i}_snps.tsv
done

wait

# # calculate SNP codon numbers.
# # awk -F"\t" '{printf "%.0f\n", ((($2-265)/3)+0.3)}' $FILEPATH/snps/${j}_snps2.tsv > $FILEPATH/snps/${j}_snps3.tsv

# # recalculate SNP codon numbers & add gene description.
# # awk '{ if ($1 > 0 && $1 < 4407) print $1,"orf1a"; \
# # else if ($1 >= 4407 && $1 < 7186) print $1-4406,"orf1b"; \
# # else if ($1 >= 7186 && $1 < 8462) print $1-7099,"s"}' \
# # $FILEPATH/snps/${j}_snps3.tsv

# -----------------------------------------
# SAVE SAMPLE SUMMARY FILES
# -----------------------------------------

echo RNA-${j} > $FILEPATH/${j}/inter/${j}_sample_name.txt

paste \
$FILEPATH/${j}/inter/${j}_sample_name.txt \
$FILEPATH/${j}/output/${j}_filtered_reads.txt \
$FILEPATH/${j}/binning/${j}_BA2.txt \
$FILEPATH/${j}/binning/${j}_BA2320.txt \
$FILEPATH/${j}/binning/${j}_BA275x.txt \
$FILEPATH/${j}/binning/${j}_BN.txt \
$FILEPATH/${j}/binning/${j}_CH11.txt \
$FILEPATH/${j}/binning/${j}_BQ1x.txt \
$FILEPATH/${j}/binning/${j}_BR2.txt \
$FILEPATH/${j}/binning/${j}_XBB.txt \
$FILEPATH/${j}/binning/${j}_XBC.txt \
$FILEPATH/${j}/binning/${j}_XBF.txt \
$FILEPATH/${j}/binning/${j}_XBB15.txt \
$FILEPATH/${j}/binning/${j}_BA45.txt \
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
sed -i '1i Accession\tFiltered\tMapping_BA2\tMapping_BA2320\tMapping_BA275x\tMapping_BNx\tMapping_CH11\tMapping_BQ1x\tMapping_BR2\tMapping_XBB\tMapping_XBC\tMapping_XBF\tMapping_XBB15\tMapping_BA45\tMapping_Other\tLineagesDetected\tCoveragePercent\tDepth' $FILEPATH/wgs_run_summary.txt

rm -r $FILEPATH/*/inter

rm -rf $TMPDIR

# -----------------------------------------
# THE END
# -----------------------------------------
