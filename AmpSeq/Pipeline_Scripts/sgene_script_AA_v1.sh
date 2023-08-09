#!/bin/sh

#SBATCH --ntasks=1
#SBATCH --job-name=AmpSeq
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=24
#SBATCH --mem=24G
#SBATCH --output=AmpSeq_%A_%a.out
#SBATCH --mail-type=ALL
#SBATCH --mail-user=hartman.l@wehi.edu.au

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
sed 's/\r$//' -i /wehisan/general/academic/lab_jex/MINISEQ_AMPSEQ_DATA/sgene_grep_refs/ampseq_variants_MASTER_LIST_194bp.txt
sed 's/\r$//' -i /wehisan/general/academic/lab_jex/MINISEQ_AMPSEQ_DATA/sgene_grep_refs/ampseq_440_503_codons_REF.txt

# create sample and file name reference files.
sort -t ',' -k2 $FILEPATH/meta/oh2_meta.txt | awk -F, '{print $1}' > $FILEPATH/meta/oh2_index.txt
sort -t ',' -k2 $FILEPATH/meta/oh2_meta.txt | awk -F, '{print $2}' > $FILEPATH/meta/oh2_names.txt
readarray -t oh2 < $FILEPATH/meta/oh2_index.txt
readarray -t names < $FILEPATH/meta/oh2_names.txt

# get data from master files.
source /wehisan/general/academic/lab_jex/MINISEQ_AMPSEQ_DATA/sgene_grep_refs/ampseq_grep_MASTER_FILE_194bp.txt
readarray -t variants < /wehisan/general/academic/lab_jex/MINISEQ_AMPSEQ_DATA/sgene_grep_refs/ampseq_variants_MASTER_LIST_194bp.txt
readarray -t codons < /wehisan/general/academic/lab_jex/MINISEQ_AMPSEQ_DATA/sgene_grep_refs/ampseq_440_503_codons_REF.txt

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
trimmomatic SE -threads 24 -phred33 $FILEPATH/trim/${i}.fastq $FILEPATH/filter/${i}.fastq HEADCROP:23 CROP:198 SLIDINGWINDOW:5:30 MINLEN:198
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

# remove strings that don't have Y505Y or Y505H because the 513_Rev primer has a 'Y' nucleotide (matching C or T) at position 23075,
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

# get low counts strings.
awk '/      1/ || /      2/ || /      3/ || /      4/ || /      5/' $FILEPATH/no_low/all_unique.txt > $FILEPATH/no_low/low_counts.txt

# remove first 8 characters (count data).
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

# get the line numbers ('indexes') of the pattern-matching unique reads.
for i in ${variants[@]}
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

# convert from wide to long.
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

# convert from long to wide.
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
cat $FILEPATH/other/snps_0${i}.txt >> $FILEPATH/top10_other_snps.txt
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
$FILEPATH/inter/AYX.txt \
$FILEPATH/inter/BA1.txt \
$FILEPATH/inter/BA2.txt \
$FILEPATH/inter/BA4X_BA5X.txt \
$FILEPATH/inter/BA2121.txt \
$FILEPATH/inter/BA213.txt \
$FILEPATH/inter/BA275.txt \
$FILEPATH/inter/BL13.txt \
$FILEPATH/inter/BL14.txt \
$FILEPATH/inter/BM1.txt \
$FILEPATH/inter/BM111.txt \
$FILEPATH/inter/BM21.txt \
$FILEPATH/inter/BM23.txt \
$FILEPATH/inter/BA2320.txt \
$FILEPATH/inter/CM21.txt \
$FILEPATH/inter/CM3.txt \
$FILEPATH/inter/CM8.txt \
$FILEPATH/inter/BH1.txt \
$FILEPATH/inter/BJ1.txt \
$FILEPATH/inter/BN1.txt \
$FILEPATH/inter/BN11.txt \
$FILEPATH/inter/BN3.txt \
$FILEPATH/inter/BN31.txt \
$FILEPATH/inter/BQ1.txt \
$FILEPATH/inter/BQ119.txt \
$FILEPATH/inter/BQ1125.txt \
$FILEPATH/inter/BQ2.txt \
$FILEPATH/inter/BR1.txt \
$FILEPATH/inter/BR12.txt \
$FILEPATH/inter/BR2.txt \
$FILEPATH/inter/BR3.txt \
$FILEPATH/inter/BR4.txt \
$FILEPATH/inter/BS1.txt \
$FILEPATH/inter/BY11.txt \
$FILEPATH/inter/BY121.txt \
$FILEPATH/inter/CH1.txt \
$FILEPATH/inter/CH11.txt \
$FILEPATH/inter/CH3.txt \
$FILEPATH/inter/CH31.txt \
$FILEPATH/inter/XAY.txt \
$FILEPATH/inter/XBB.txt \
$FILEPATH/inter/XBB13.txt \
$FILEPATH/inter/XBB41.txt \
$FILEPATH/inter/XBB4.txt \
$FILEPATH/inter/XBC.txt \
$FILEPATH/inter/XBC1.txt \
$FILEPATH/inter/XBC2.txt \
$FILEPATH/inter/BA2104.txt \
$FILEPATH/inter/XBF.txt \
$FILEPATH/inter/XBB15.txt \
$FILEPATH/inter/BU1.txt \
$FILEPATH/inter/BU2.txt \
$FILEPATH/inter/BU3.txt \
$FILEPATH/inter/BR21_F490S.txt \
$FILEPATH/inter/BN1_E484V.txt \
$FILEPATH/inter/BQ1_T478R.txt \
$FILEPATH/inter/BQ1_T478E.txt \
$FILEPATH/inter/BQ1_G476S.txt \
$FILEPATH/inter/XBC1_F490S.txt \
$FILEPATH/inter/XBF_F490P.txt \
$FILEPATH/inter/XBF_L452R.txt \
$FILEPATH/inter/XBF_K444T.txt \
$FILEPATH/inter/XBF_K444T_L452R.txt \
$FILEPATH/inter/XBF_K444T_G446G_L452R.txt \
$FILEPATH/inter/XBF_K444T_G446G.txt \
$FILEPATH/inter/XBF_G446N.txt \
$FILEPATH/inter/BQ1_F486I.txt \
$FILEPATH/inter/BQ1_G446S.txt \
$FILEPATH/inter/BQ1_G446S_F486I.txt \
$FILEPATH/inter/CB1_L452R.txt \
$FILEPATH/inter/XBB_L452R.txt \
> $FILEPATH/read_counts.txt

wait

# add header.
sed -i '1i sample\tfiltered\tAYX\tBA1\tBA2\tBA45\tBA2121\tBA213\tBA275\tBL13\tBL14\tBM1\tBM111\tBM21\tBM23\tBA2320\tCM21\tCM3\tCM8\tBH1\tBJ1\tBN1\tBN11\tBN3\tBN31\tBQ1\tBQ119\tBQ1125\tBQ2\tBR1BR3\tBR12\tBR2\tBR3\tBR4\tBS1\tCA1\tCA2\tCH1\tCH11\tCH3\tCH31\tXAY\tXBB\tXBB13\tXBB41\tXBB4\tXBC\tXBC1\tXBC2\tBA2104\tXBF\tXBB15\tBU1\tBU2\tBU3\tBR21_F490S\tBN1_E484V\tBQ1_T478R\tBQ1_T478E\tBQ1_G476S\tXBC1_F490S\tXBF_F490P\tXBF_L452R\tXBF_K444T\tXBF_K444T_L452R\tXBF_K444T_G446G_L452R\tXBF_K444T_G446G_L452R\tXBF_G446N\tBQ1_F486I\tBQ1_G446S\tBQ1_G446S_F486I\tCB1_L452R\tXBB_L452R' $FILEPATH/read_counts.txt

wait

## ----------------------------- ##
## COUNT VARIANT-MATCHING READS  ##
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

# get % of each AA - if counts are below 200, a small default value is output to prevent and reporting
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
##     COUNT REF & ALT CODONS    ##
## ----------------------------- ##

for i in ${names[@]}
do

awk '{ if($2 != $3 && $4 >= 2 \
&& ($1!="N440K" \
&& $1!="K444R" \
&& $1!="V445P" \
&& $1!="G446S" \
&& $1!="N450D" \
&& $1!="L452L" \
&& $1!="N460K" \
&& $1!="T470N" \
&& $1!="G476S" \
&& $1!="S477N" \
&& $1!="E484A" \
&& $1!="G485D" \
&& $1!="F486V" \
&& $1!="F490L" \
&& $1!="Q493R" \
&& $1!="S494P" \
&& $1!="G496S" \
&& $1!="Q498R" \
&& $1!="N501Y" \
&& $1!="K444T" \
&& $1!="G446D" \
&& $1!="L452Q" \
&& $1!="E484R" \
&& $1!="F486S" \
&& $1!="F490V" \
&& $1!="K444M" \
&& $1!="L452M" \
&& $1!="E484T" \
&& $1!="F486I" \
&& $1!="F490S" \
&& $1!="E484V" \
&& $1!="F490P" \
&& $1!="F486P" \
&& $1!="F486L" \
&& $1!="T478K" \
&& $1!="L452R")) \
print $1, $4}' $FILEPATH/inter/${i}_all_AAs.txt > $FILEPATH/inter/${i}_non_class_AAs.txt

awk '{ if($2 != $3 && $4 >= 2) print $1}' $FILEPATH/inter/${i}_all_AAs.txt > $FILEPATH/inter/${i}_all_reportable_AAs.txt

done
wait

## ----------------------------- ##
##   EXTRACT CODON COUNT DATA    ##
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
> $FILEPATH/snp_counts.txt

wait

## ----------------------------- ##
##           CLEAN UP            ##
## ----------------------------- ##

rm -r $FILEPATH/demulti
rm -r $FILEPATH/trim
rm -r $FILEPATH/filter
rm -r $FILEPATH/best_75k
rm -r $FILEPATH/rev_comp
rm -r $FILEPATH/aa_strings
rm -r $FILEPATH/primed_ok
rm -r $FILEPATH/codon_trim
rm -r $FILEPATH/no_indels
rm -r $FILEPATH/no_low
rm -r $FILEPATH/other
rm -r $FILEPATH/inter

rm -rf $TMPDIR

## ----------------------------- ##
## ---------- THE END ---------- ##
## ----------------------------- ##

