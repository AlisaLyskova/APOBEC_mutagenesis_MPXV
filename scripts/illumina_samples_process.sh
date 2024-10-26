dir=data/
dir_genome=data/
bwa_dir=
freebayes_dir=

#download - sratoolkit, fastq-dump : 3.1.1
for SRR_id in <SRR list>
do
     echo $SRR_id
     fastq-dump $SRR_id --gzip -O $dir
done

#genome indexing
##bwa-0.7.17
$bwa_dir/./bwa index $dir_genome/NC_063383.1.fa

#alignment
##bwa-0.7.17
##samtools-1.20
for f in $dir/*.fastq.gz
do
    echo $f
    bwa_dir/./bwa mem $dir_genome/NC_063383.1.fa $f | samtools sort | samtools view -F 4 --write-index -b -o ${f%.fastq.gz}.bam
done

#SNP
##freebayes-1.3.6
for f in $dir/*sorted.bam
do
    echo $f
    freebayes_dir/./freebayes-1.3.6-linux-amd64-static -f $dir_genome/NC_063383.1.fasta $f > ${f%_sorted.bam}.vcf
done
