#!/bin/sh

#SBATCH --ntasks=1
#SBATCH --job-name=AmpSeq
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=24
#SBATCH --mem=24G
#SBATCH --output=AmpSeq_%A_%a.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hartman.l@wehi.edu.au

## ----------------------------- ##

module load cutadapt
module load trimmomatic
module load anaconda3

# redirect temp directory path for fast read/write processing.
export TMPDIR=$(mktemp -d --tmpdir=/vast/scratch/users/$USER)

# assign job directory path variable.
FILEPATH=`pwd -P`

# delete carriage returns in metadata files.
sed 's/\r$//' -i $FILEPATH/meta/oh2_meta.txt
sed 's/\r$//' -i /wehisan/general/academic/lab_jex/MINISEQ_AMPSEQ_DATA/sgene_grep_refs/ampseq_grep_MASTER_FILE_194bp.txt
sed 's/\r$//' -i /wehisan/general/academic/lab_jex/MINISEQ_AMPSEQ_DATA/sgene_grep_refs/ampseq_440_503_codons_REF.txt

# create files of reference variables.
sort -t ',' -k2 $FILEPATH/meta/oh2_meta.txt | awk -F, '{print $1}' > $FILEPATH/meta/oh2_index.txt
sort -t ',' -k2 $FILEPATH/meta/oh2_meta.txt | awk -F, '{print $2}' > $FILEPATH/meta/oh2_names.txt
awk -F '=' '{print $1}' /wehisan/general/academic/lab_jex/MINISEQ_AMPSEQ_DATA/sgene_grep_refs/ampseq_grep_MASTER_FILE_194bp.txt > $FILEPATH/meta/variants.txt

# read reference variables into memory.
readarray -t oh2 < $FILEPATH/meta/oh2_index.txt
readarray -t names < $FILEPATH/meta/oh2_names.txt
readarray -t variants < $FILEPATH/meta/variants.txt
source /wehisan/general/academic/lab_jex/MINISEQ_AMPSEQ_DATA/sgene_grep_refs/ampseq_grep_MASTER_FILE_194bp.txt
readarray -t codons < /wehisan/general/academic/lab_jex/MINISEQ_AMPSEQ_DATA/sgene_grep_refs/ampseq_440_503_codons_REF.txt

rm -r $FILEPATH/output
mkdir $FILEPATH/output

## ----------------------------- ##
##          DEMULTIPLEX          ##
## ----------------------------- ##

rm -r $FILEPATH/demulti
mkdir $FILEPATH/demulti

cutadapt \
-j 24 \
-e 0 \
--no-indels \
-g file:$FILEPATH/meta/fwd_indexes.fasta \
-G file:$FILEPATH/meta/rev_indexes.fasta \
-o $FILEPATH/demulti/{name1}_{name2}.fastq \
-p $FILEPATH/demulti/{name1}_{name2}_rev.fastq \
$FILEPATH/raw_seqs/*R1*.fastq.gz \
$FILEPATH/raw_seqs/*R2*.fastq.gz

# delete redundant files.
rm -r $FILEPATH/demulti/*_rev.fastq
rm -r $FILEPATH/demulti/*unknown*


wait

## ----------------------------- ##
##         RENAME FILES          ##
## ----------------------------- ##

awk -F, 'system("mv '$FILEPATH'/demulti/" $1 " '$FILEPATH'/demulti/" $2".fastq")' $FILEPATH/meta/oh2_meta.txt

wait

## ----------------------------- ##
##        REMOVE OVERHANGS       ##
## ----------------------------- ##

rm -r $FILEPATH/trim
mkdir $FILEPATH/trim

for i in ${names[@]}
do
cutadapt -g GTGACCTATGAACTCAGGAGTC $FILEPATH/demulti/${i}.fastq | cutadapt -a GCTGCGATGTGCAAGTCTCAG - -o $FILEPATH/trim/${i}.fastq
done
wait

## ----------------------------- ##
##   TRIM & FILTER OH2 3' READS  ##
## ----------------------------- ##

rm -r $FILEPATH/filter
mkdir $FILEPATH/filter

for i in ${names[@]}
do
trimmomatic SE -threads 24 -phred33 $FILEPATH/trim/${i}.fastq $FILEPATH/filter/${i}.fastq HEADCROP:23 CROP:198 SLIDINGWINDOW:5:25 MINLEN:198
done
wait

## ----------------------------- ##
##     RETAIN BEST 75k READS     ##
## ----------------------------- ##

rm -r $FILEPATH/best_75k
mkdir $FILEPATH/best_75k

conda activate lofreq

for i in ${names[@]}
do
filtlong --window_size 50 --target_bases 15000000 $FILEPATH/filter/${i}.fastq > $FILEPATH/best_75k/${i}.fastq
done
wait

## ------------------------------##
##         REV COMP READS        ##
## ------------------------------##

rm -r $FILEPATH/rev_comp
mkdir $FILEPATH/rev_comp

for i in ${names[@]}
do
seqkit seq --threads 24 $FILEPATH/best_75k/${i}.fastq -r -p --out-file $FILEPATH/rev_comp/${i}.fastq
done
wait

## ------------------------------##
##     CONVERT TO AA STRINGS     ##
## ------------------------------##

rm -r $FILEPATH/aa_strings
mkdir $FILEPATH/aa_strings

for i in ${names[@]}
do
seqkit translate --threads 24 --frame 1 --line-width 0 --clean $FILEPATH/rev_comp/${i}.fastq | awk '!/^>/' > $FILEPATH/aa_strings/${i}.txt
done
wait

conda deactivate

## ----------------------------- ##
##   REMOVE POORLY PRIMED READS  ##
## ----------------------------- ##

# remove strings that don't have Y505Y or Y505H because the 513_Rev primer has a 'Y' amino acid (matching C or T) at position 23075,
# therefore strings without the Y or H amino acid are unreliable.

rm -r $FILEPATH/primed_ok
mkdir $FILEPATH/primed_ok

for i in ${names[@]}
do
grep -E 'Y$|H$' $FILEPATH/aa_strings/${i}.txt > $FILEPATH/primed_ok/${i}.txt
done
wait

## ----------------------------- ##
##       TRIM TWO 3' CODONS      ##
## ----------------------------- ##

# remove the two 3' codons because they sit within the primer.

rm -r $FILEPATH/codon_trim
mkdir $FILEPATH/codon_trim

for i in ${names[@]}
do
cut -c-64 $FILEPATH/primed_ok/${i}.txt > $FILEPATH/codon_trim/${i}.txt
done
wait

## ----------------------------- ##
##   REMOVE READS WITH INDELS    ##
## ----------------------------- ##

rm -r $FILEPATH/no_indels
mkdir $FILEPATH/no_indels

for i in ${names[@]}
do
grep -E '^N|^K' $FILEPATH/codon_trim/${i}.txt > $FILEPATH/no_indels/${i}.txt
done
wait

## ----------------------------- ##
##    REMOVE LOW COUNT READS     ##
## ----------------------------- ##

rm -r $FILEPATH/no_low
mkdir $FILEPATH/no_low

# get and count all unique strings.
cat $FILEPATH/no_indels/*.txt | sort | uniq -c | sort -g -r > $FILEPATH/no_low/all_unique.txt

# save reference copy of all_unique.txt without the low count strings (i.e. without the strings detected ≤5 times across the dataset).
awk '!/      1/ && !/      2/ && !/      3/ && !/      4/ && !/      5/' $FILEPATH/no_low/all_unique.txt > $FILEPATH/output/all_unique.txt

# get the low count strings.
awk '/      1/ || /      2/ || /      3/ || /      4/ || /      5/' $FILEPATH/no_low/all_unique.txt > $FILEPATH/no_low/low_counts.txt

# remove first 8 characters (i.e. the count data).
sed -i 's/^.\{8\}//g' $FILEPATH/no_low/low_counts.txt

# remove low count strings from each sample.
for i in ${names[@]}
do
awk 'FNR==NR {hash[$0]; next} !($0 in hash)' $FILEPATH/no_low/low_counts.txt $FILEPATH/no_indels/${i}.txt > $FILEPATH/no_low/${i}.txt
done
wait

## ----------------------------- ##
## IDENTIFY TOP 10 'OTHER' READS ##
## ----------------------------- ##

rm -r $FILEPATH/other
mkdir $FILEPATH/other

# read variants2 into memory (excluding FLiP) so we can identify the most abundant non-target strings.
awk '{if ($1 != "FLiP") print}' $FILEPATH/meta/variants.txt > $FILEPATH/other/variants2.txt
readarray -t variants2 < $FILEPATH/other/variants2.txt

# get the line numbers ('indexes') of the pattern-matching unique reads.
for i in ${variants2[@]}
do
grep -n "${!i}" $FILEPATH/no_low/all_unique.txt | cut -f1 -d: >> $FILEPATH/other/variant_indexes.txt
done
wait

# extract the 10 most abundant non pattern-matching (i.e. "Other") reads.
awk 'NR==FNR { x[$1]; next } !(FNR in x)' $FILEPATH/other/variant_indexes.txt $FILEPATH/no_low/all_unique.txt | head -n 10 > $FILEPATH/other/top10_other.txt

# write each 'Other' read to a separate file.
split -d -l 1 $FILEPATH/other/top10_other.txt $FILEPATH/other/wide_

## ----------------------------- ##
##     GET 'OTHER' READ SNPS     ##
## ----------------------------- ##

# create AA ref file.
echo ${codons[@]} | sed 's/ /\n/g' > $FILEPATH/other/ref.txt

# create codon number file.
for p in $(seq 440 503); do echo ${p}; done | sed 's/ /\n/g' > $FILEPATH/other/num.txt

# save 'Other' abundance values.
for i in $(seq 0 9)
do
awk '{print $1}' $FILEPATH/other/wide_0${i} > $FILEPATH/other/abund_0${i}.txt
done

# convert from wide to long format.
for i in $(seq 0 9)
do
cat $FILEPATH/other/wide_0${i} | rev | cut -c-64 | rev | sed 's/./\0\n/g' > $FILEPATH/other/long_0${i}.txt
done

# paste snp fields.
for i in $(seq 0 9)
do
paste $FILEPATH/other/ref.txt $FILEPATH/other/num.txt $FILEPATH/other/long_0${i}.txt > $FILEPATH/other/join_0${i}.txt
done

# extract alt snp info.
for i in $(seq 0 9)
do
awk '{ if($1 != $3) {print "S:"$1$2$3} }' $FILEPATH/other/join_0${i}.txt >> $FILEPATH/other/snps_0${i}.txt
done

# convert from long to wide format.
for i in $(seq 0 9)
do
sed -zi 's/\n/, /g' $FILEPATH/other/snps_0${i}.txt
done

# replace last ", " with newline.
for i in $(seq 0 9)
do
sed -zi 's/, $/\n/g' $FILEPATH/other/snps_0${i}.txt
done

# insert abundance values.
for i in $(seq 0 9)
do
abund=`cat $FILEPATH/other/abund_0${i}.txt`
sed -i "1s/^/$abund\t/" $FILEPATH/other/snps_0${i}.txt
done

# combine all 'Other' snp info.
for i in $(seq 0 9)
do
cat $FILEPATH/other/snps_0${i}.txt >> $FILEPATH/output/top10_other_snps.txt
done

## ----------------------------- ##
## COUNT VARIANT-MATCHING READS  ##
## ----------------------------- ##

rm -r $FILEPATH/inter
mkdir $FILEPATH/inter

for i in ${variants[@]}
do
for j in ${names[@]}
do
grep -xc "${!i}" $FILEPATH/no_low/${j}.txt >> $FILEPATH/inter/${i}.txt
done
wait
done
wait

## ----------------------------- ##

# get total number of filtered reads per sample.
for i in ${names[@]}
do
wc -l $FILEPATH/no_low/${i}.txt | awk '{print $1}' >> $FILEPATH/inter/filtered_reads.txt
done
wait

# combine data.
paste \
$FILEPATH/meta/oh2_names.txt \
$FILEPATH/inter/filtered_reads.txt \
$FILEPATH/inter/BA2.txt \
$FILEPATH/inter/BA275.txt \
$FILEPATH/inter/CH11.txt \
$FILEPATH/inter/XBC.txt \
$FILEPATH/inter/XBB19.txt \
$FILEPATH/inter/XBB116.txt \
$FILEPATH/inter/EG51.txt \
$FILEPATH/inter/HK3.txt \
$FILEPATH/inter/BA286.txt \
$FILEPATH/inter/FLiP.txt \
$FILEPATH/inter/HV1.txt \
$FILEPATH/inter/DV71.txt \
$FILEPATH/inter/GW511.txt \
$FILEPATH/inter/JN1.txt \
> $FILEPATH/output/read_counts.txt

wait

# add header.
sed -i '1i sample\tfiltered\tBA2\tBA275\tCH11\tXBC\tXBB19\tXBB116\tEG51\tHK3\tBA286\tFLiP\tHV1\tDV71\tGW511\tJN1' $FILEPATH/output/read_counts.txt

wait

## ----------------------------- ##
##     COUNT REF & ALT CODONS    ##
## ----------------------------- ##

for i in ${names[@]}
do
for k in {0..63}
do

# get total filtered count.
j=`wc -l $FILEPATH/no_low/${i}.txt | awk '{print $1}'`

# adjust codon index for non-array use.
let "q=$k+1"

# make list of AAs detected at codon q.
cut -b ${q} $FILEPATH/no_low/${i}.txt | sort | uniq -c | awk '{print $2}' > $FILEPATH/inter/AAs.txt

# get number of AAs detected.
m=`wc -l $FILEPATH/inter/AAs.txt | awk '{print $1}'`

# get % of each AA - if counts are below 200, a small placeholder value is output
cut -b ${q} $FILEPATH/no_low/${i}.txt | sort | uniq -c | awk '{ if ($1 > 200) {print $1/'${j}'*100} else {print 0.1} }' > $FILEPATH/inter/per1.txt

# round % to 1 decimal place.
awk '{for (i=1; i<=NF; i++) printf "%.1f %s", $i, (i==NF?RS:FS)}' $FILEPATH/inter/per1.txt > $FILEPATH/inter/per2.txt

# get codon number.
let "l=$k+440"

# create codon number file with the required number of rows ($m).
for p in $(seq 1 $m); do echo ${l}; done > $FILEPATH/inter/num.txt

# create ref AA file with the required number of rows ($m).
for p in $(seq 1 $m); do echo ${codons[${k}]}; done > $FILEPATH/inter/ref.txt

# paste info 1.
paste -d "" \
$FILEPATH/inter/ref.txt \
$FILEPATH/inter/num.txt \
$FILEPATH/inter/AAs.txt \
> $FILEPATH/inter/paste1.txt

wait

# paste info 2.
paste \
$FILEPATH/inter/paste1.txt \
$FILEPATH/inter/ref.txt \
$FILEPATH/inter/AAs.txt \
$FILEPATH/inter/per2.txt \
> $FILEPATH/inter/paste2.txt

# concatenate info.
cat $FILEPATH/inter/paste2.txt >> $FILEPATH/inter/${i}_all_AAs.txt

done
wait
done
wait

## ----------------------------- ##
##        FILTER SNP LISTS       ##
## ----------------------------- ##

for i in ${names[@]}
do
# filter characterising SNPs and SNPs caused by V483 deletion in BA286/JN1 to leave only non-characterising SNPs
awk '{ if($2 != $3 && $4 >= 2 \
&& ($1!="N440K" \
&& $1!="L441K" \
&& $1!="D442L" \
&& $1!="S443D" \
&& $1!="K444S" \
&& $1!="V445K" \
&& $1!="G446H" \
&& $1!="G447S" \
&& $1!="N448G" \
&& $1!="Y449N" \
&& $1!="N450Y" \
&& $1!="Y451D" \
&& $1!="L452Y" \
&& $1!="Y453W" \
&& $1!="R454Y" \
&& $1!="L455R" \
&& $1!="F456S" \
&& $1!="R457F" \
&& $1!="K458R" \
&& $1!="S459K" \
&& $1!="N460S" \
&& $1!="L461K" \
&& $1!="K462L" \
&& $1!="P463K" \
&& $1!="F464P" \
&& $1!="E465F" \
&& $1!="R466E" \
&& $1!="D467R" \
&& $1!="I468D" \
&& $1!="S469I" \
&& $1!="T470S" \
&& $1!="E471T" \
&& $1!="I472E" \
&& $1!="Y473I" \
&& $1!="Q474Y" \
&& $1!="A475Q" \
&& $1!="G476A" \
&& $1!="S477G" \
&& $1!="T478N" \
&& $1!="P479K" \
&& $1!="C480P" \
&& $1!="N481C" \
&& $1!="G482K" \
&& $1!="V483G" \
&& $1!="K444T" \
&& $1!="V445P" \
&& $1!="V445H" \
&& $1!="G446S" \
&& $1!="N450D" \
&& $1!="L452R" \
&& $1!="L452W" \
&& $1!="L455F" \
&& $1!="F456L" \
&& $1!="N460K" \
&& $1!="A475V" \
&& $1!="S477N" \
&& $1!="T478K" \
&& $1!="T478I" \
&& $1!="T478R" \
&& $1!="N481K" \
&& $1!="E484A" \
&& $1!="E484K" \
&& $1!="F486S" \
&& $1!="F486P" \
&& $1!="F490S" \
&& $1!="Q493R" \
&& $1!="Q498R" \
&& $1!="N501Y")) \
print $1, $4}' $FILEPATH/inter/${i}_all_AAs.txt > $FILEPATH/inter/${i}_non_class_AAs.txt

# filter SNPs caused by V483 deletion in BA286/JN1 to leave only characterising SNPs and non-characterising SNPs 
awk '{ if($2 != $3 && $4 >= 2 \
&& ($1!="L441K" \
&& $1!="D442L" \
&& $1!="S443D" \
&& $1!="K444S" \
&& $1!="V445K" \
&& $1!="G446H" \
&& $1!="G447S" \
&& $1!="N448G" \
&& $1!="Y449N" \
&& $1!="N450Y" \
&& $1!="Y451D" \
&& $1!="L452Y" \
&& $1!="Y453W" \
&& $1!="R454Y" \
&& $1!="L455R" \
&& $1!="F456S" \
&& $1!="R457F" \
&& $1!="K458R" \
&& $1!="S459K" \
&& $1!="N460S" \
&& $1!="L461K" \
&& $1!="K462L" \
&& $1!="P463K" \
&& $1!="F464P" \
&& $1!="E465F" \
&& $1!="R466E" \
&& $1!="D467R" \
&& $1!="I468D" \
&& $1!="S469I" \
&& $1!="T470S" \
&& $1!="E471T" \
&& $1!="I472E" \
&& $1!="Y473I" \
&& $1!="Q474Y" \
&& $1!="A475Q" \
&& $1!="G476A" \
&& $1!="S477G" \
&& $1!="T478N" \
&& $1!="P479K" \
&& $1!="C480P" \
&& $1!="N481C" \
&& $1!="G482K" \
&& $1!="V483G")) \
print $1}' $FILEPATH/inter/${i}_all_AAs.txt > $FILEPATH/inter/${i}_all_reportable_AAs.txt

done
wait

## ----------------------------- ##
##   EXTRACT ALT CODON COUNTS    ##
## ----------------------------- ##

for i in ${names[@]}
do

# get non-classified SNP names.
awk -v RNA="$i" 'BEGIN {print RNA}; {print $1}' $FILEPATH/inter/${i}_non_class_AAs.txt | sed -z 's/\n/, /g' | sed -r 's/, $/\n/' | sed -r 's/, /\t/' >> $FILEPATH/inter/NonClassSNPs.txt

# get non-classified SNP percentages.
awk '{print $2}; END {if (NR == 0) {print "\t"}}' $FILEPATH/inter/${i}_non_class_AAs.txt | sed -z 's/\n/, /g' | sed -r 's/, $/\n/' >> $FILEPATH/inter/PercNonClassSNPs.txt

# get all detected SNPs.
awk '{print $1}; END {if (NR == 0) {print "\t"}}' $FILEPATH/inter/${i}_all_reportable_AAs.txt | sed -z 's/\n/, /g' | sed -r 's/, $/\n/' >> $FILEPATH/inter/DetectedSNPs.txt

done
wait

# add headers.
sed -i '1i sample\tNonClassSNPs' $FILEPATH/inter/NonClassSNPs.txt
sed -i '1i PercentNonClassSNPs' $FILEPATH/inter/PercNonClassSNPs.txt
sed -i '1i DetectedSNPs' $FILEPATH/inter/DetectedSNPs.txt

wait

# combine data.
paste \
$FILEPATH/inter/NonClassSNPs.txt \
$FILEPATH/inter/PercNonClassSNPs.txt \
$FILEPATH/inter/DetectedSNPs.txt \
> $FILEPATH/output/snp_counts.txt

wait

## ----------------------------- ##
## COUNT TOP 10 OTHER PER SAMPLE ##
## ----------------------------- ##

# get search strings.
awk '{print $2}' $FILEPATH/other/top10_other.txt > $FILEPATH/other/top10_other_strings.txt
readarray -t other < $FILEPATH/other/top10_other_strings.txt

# count search strings in each sample.
for i in $(seq 0 9)
do
for j in ${names[@]}
do
grep -xc "${other[$i]}" $FILEPATH/no_low/${j}.txt >> $FILEPATH/other/other_counts_0${i}.txt
done
wait
done
wait

# combine data.
paste \
$FILEPATH/meta/oh2_names.txt \
$FILEPATH/inter/filtered_reads.txt \
$FILEPATH/other/other_counts_00.txt \
$FILEPATH/other/other_counts_01.txt \
$FILEPATH/other/other_counts_02.txt \
$FILEPATH/other/other_counts_03.txt \
$FILEPATH/other/other_counts_04.txt \
$FILEPATH/other/other_counts_05.txt \
$FILEPATH/other/other_counts_06.txt \
$FILEPATH/other/other_counts_07.txt \
$FILEPATH/other/other_counts_08.txt \
$FILEPATH/other/other_counts_09.txt \
> $FILEPATH/output/other_counts.txt

# add header.
sed -i '1i sample\tfiltered\tA\tB\tC\tD\tE\tF\tG\tH\tI\tJ' $FILEPATH/output/other_counts.txt

## ----------------------------- ##
##           CLEAN UP            ##
## ----------------------------- ##

rm -r $FILEPATH/demulti
rm -r $FILEPATH/trim
rm -r $FILEPATH/filter
rm -r $FILEPATH/rev_comp
rm -r $FILEPATH/aa_strings
rm -r $FILEPATH/primed_ok
rm -r $FILEPATH/codon_trim
rm -r $FILEPATH/no_indels
rm -r $FILEPATH/no_low
rm -r $FILEPATH/other
rm -r $FILEPATH/inter
rm -r $FILEPATH/best_75k

rm -rf $TMPDIR

## ----------------------------- ##
## ---------- THE END ---------- ##
## ----------------------------- ##
