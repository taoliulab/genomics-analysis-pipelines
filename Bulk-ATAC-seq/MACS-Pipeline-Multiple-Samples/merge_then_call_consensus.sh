#!/bin/bash
#
# This script will find the consensus peak regions from peak files (in
# BED format) of multiple samples by:
#
# 1. Converting the peak file of each sample into non-overlapping 3
# cols BED file and concatenating them;
#
# 2. Sorting the concatenated file and Building a genome coverage
# track in BedGraph, of which the value (the 3rd col) indicates the
# number of samples with a peak at a particular location;
# 
# 3. Using MACS to call regions being covered by peaks from more than
# a certain number of samples.

# ------------------------------------------------------------------
# Modify the following parameters
#
# define the peaks from multiple samples
SAMPLE_PEAKS=`ls *.narrowPeak`

# define the CUTOFF
CUTOFF=3

# define the minlen and maxgap for peak calling (arbitrary)
MINLEN=200
MAXGAP=30

# define the genome file with the 1st column as chromosome name and
# the 2nd column as chromosome length
GENOME=~/db/genome/GRCh38/chrom.len

# define the OUTPUT_PREFIX
OUTPUT_PREFIX="consensus_peaks"

# ------------------------------------------------------------------

# 1
I=$SAMPLE_PEAKS
O=${OUTPUT_PREFIX}.all.bed

rm -f $O
touch $O

for f in $I; do
    bedtools sort -i $f | bedtools merge -i - | cut -f 1,2,3 >> $O;
done

# 2
I=${OUTPUT_PREFIX}.all.bed
O=${OUTPUT_PREFIX}.bdg

bedtools sort -i $I | bedtools genomecov -bga -i -  -g $GENOME > $O

#3 
I=${OUTPUT_PREFIX}.bdg
O=${OUTPUT_PREFIX}.consensus.bed

macs2 bdgpeakcall -i $I -o $O --no-trackline -c $CUTOFF -g $MAXGAP -l $MINLEN

# end
echo "All done. Check ${O}"






