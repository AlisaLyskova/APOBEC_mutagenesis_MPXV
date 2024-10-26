INPUT_DIR="/data1/lyskovaa/mpxv"
PROJECT_DIR="PRJNA864638_nanopore"
GENOME_DIR=/data1/lyskovaa/mpxv/genome
SRR_acc_list=$(cat "$INPUT_DIR"/"$PROJECT_DIR"/SRR_Acc_List.txt)

#download - sratoolkit, fastq-dump : 3.1.1
for line in $SRR_acc_list
do
    SRR_id=$(echo -n "$line" | tr -d '\n')
    echo $SRR_id
    fastq-dump $SRR_id --gzip -O "$INPUT_DIR"/"$PROJECT_DIR"
done


#minimap2-2.28-r1209 index genome
/data1/biosoft/minimap2/./minimap2 -x ava-ont -d $GENOME_DIR/NC_063383.1.mmi $GENOME_DIR/NC_063383.1.fna

#mapping reads to MPXV genome
##samtools-1.20
for f in "$INPUT_DIR"/"$PROJECT_DIR"/*.fastq.gz
do
    /data1/biosoft/minimap2/./minimap2 -t 20 -a -x splice -Y -C5 --cs --MD -un -G 10000 $GENOME_DIR/NC_063383.1.fasta $f | samtools sort | samtools view --write-index -b -o ${f%.fastq.gz}_minimap.bam
done

#searching for SNP using clair3-v1.0.10
##https://github.com/HKU-BAL/Clair3
THREADS=10
MODEL_NAME="r941_prom_sup_g5014"

for file in ${INPUT_DIR}/${PROJECT_DIR}/*.bam
do
    mkdir "$INPUT_DIR"/"$PROJECT_DIR"/res_clair/${file##*/}
    OUTPUT_DIR="$INPUT_DIR"/"$PROJECT_DIR"/res_clair/"${file##*/}"
    sudo docker run -it \
  -v ${INPUT_DIR}:${INPUT_DIR} \
  -v ${OUTPUT_DIR}:${OUTPUT_DIR} \
  hkubal/clair3:latest \
  /opt/bin/run_clair3.sh \
  --bam_fn="${file}" \
  --ref_fn="${GENOME_DIR}"/NC_063383.1.fna \
  --threads=${THREADS} \
  --platform="ont" \
  --model_path="/opt/models/${MODEL_NAME}" \
  --output="${OUTPUT_DIR}" \
  --include_all_ctgs \
  --no_phasing_for_fa \
  --chunk_size=500000
done
