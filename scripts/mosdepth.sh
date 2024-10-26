#!/bin/bash

DIR=$1
OUTFILE="$DIR"/summary_align_stat.txt
touch $OUTFILE

for file in "$DIR"/*
do
    if [[ $file == *.bam ]]; then
        echo $file
        /data1/biosoft/mosdepth/mosdepth -b 1000 "${file%.bam}" "$file"
        samtools stats $file > "${file%.bam}"_coverage_stat.txt
        a=$(grep -w SN "${file%.bam}"_coverage_stat.txt | grep "raw total sequences" | awk '{print $5}')
        b=$(grep -w SN "${file%.bam}"_coverage_stat.txt | grep "reads mapped:" | awk '{print $4}')
        echo $file >> $OUTFILE
        echo "Number of reads:" >> $OUTFILE
        echo $a >> $OUTFILE
        echo "Number of mapped reads:" >> $OUTFILE
        echo $b >> $OUTFILE
        echo $(python3 -c "print(($b / $a)*100)") >> $OUTFILE
        echo "Mean coverage:" >> $OUTFILE
        cat "${file%.bam}".mosdepth.summary.txt | grep -w NC_063383.1 | awk '{print $4}' >> $OUTFILE
    fi
done

python3 /data1/lyskovaa/mpxv/scripts/coverage_mosdepth.py $DIR
