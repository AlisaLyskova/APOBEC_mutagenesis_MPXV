#!/bin/bash

# download fastq files (SRR)
for fastq in ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR105/081/ERR10543481/ERR10543481.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR105/074/ERR10513574/ERR10513574.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR105/019/ERR10550019/ERR10550019.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR105/021/ERR10550021/ERR10550021.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR105/022/ERR10550022/ERR10550022.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR105/024/ERR10550024/ERR10550024.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR105/023/ERR10550023/ERR10550023.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR105/031/ERR10550031/ERR10550031.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR105/025/ERR10550025/ERR10550025.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR105/033/ERR10550033/ERR10550033.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR109/007/ERR10963107/ERR10963107.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR105/028/ERR10550028/ERR10550028.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR109/010/ERR10963110/ERR10963110.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR109/012/ERR10963112/ERR10963112.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR109/014/ERR10963114/ERR10963114.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR105/034/ERR10550034/ERR10550034.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR109/016/ERR10963116/ERR10963116.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR109/022/ERR10963122/ERR10963122.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR105/035/ERR10550035/ERR10550035.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR109/026/ERR10963126/ERR10963126.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR109/027/ERR10963127/ERR10963127.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR109/009/ERR10963109/ERR10963109.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR109/028/ERR10963128/ERR10963128.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR109/013/ERR10963113/ERR10963113.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR109/021/ERR10963121/ERR10963121.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR109/025/ERR10963125/ERR10963125.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR105/079/ERR10543479/ERR10543479.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR105/080/ERR10543480/ERR10543480.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR105/020/ERR10550020/ERR10550020.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR105/026/ERR10550026/ERR10550026.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR105/027/ERR10550027/ERR10550027.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR105/029/ERR10550029/ERR10550029.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR105/030/ERR10550030/ERR10550030.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR105/032/ERR10550032/ERR10550032.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR105/036/ERR10550036/ERR10550036.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR109/008/ERR10963108/ERR10963108.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR109/011/ERR10963111/ERR10963111.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR109/015/ERR10963115/ERR10963115.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR109/017/ERR10963117/ERR10963117.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR109/018/ERR10963118/ERR10963118.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR109/019/ERR10963119/ERR10963119.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR109/020/ERR10963120/ERR10963120.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR109/023/ERR10963123/ERR10963123.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR109/024/ERR10963124/ERR10963124.fastq.gz
do
     wget $fastq -P /PRJEB56841/
done


# minimap2 index genome
./minimap2 -x ava-ont -d NC_063383.1.mmi NC_063383.1.fna


dir=/PRJEB56841
for f in $dir/*.fastq.gz
do
    # mapping reads to MPXV genome
    ./minimap2 -a -x ava-ont -t 10 NC_063383.1.mmi $f | samtools view --bam --with-header | samtools sort -o ${f%.fastq.gz}.bam
    samtools index ${f%.fastq.gz}.bam
done
