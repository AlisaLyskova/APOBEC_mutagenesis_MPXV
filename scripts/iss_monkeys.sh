export picard=
export bwa=

ref_dir=
ref_MPXV=NC_063383.1.fa
ref_monkey=GCF_015252025.1_Vero_WHO_p1.0_genomic.fna
ref_monkey_MPXV=monkey_mpox.fa
ref_MacFas=MacacaFas.fna
ref_rheMac=rheMac.fna
ref_MacFas_MPXV=MacacaFas_MPXV.fna
ref_rheMac_MPXV=rheMac_MPXV.fna

WORK_DIR=
OUT_DIR="${WORK_DIR}/bam_files"
mkdir $OUT_DIR
export SUMMARY_FILE="${WORK_DIR}/Nreads_iss_monkeys.txt"

function TO_FASTQ()
{
  PARAMS=("$@")
  FILE_PATH="${PARAMS[0]}"
  R1_out="${FILE_PATH%.bam}_R1.fastq"
  R2_out="${FILE_PATH%.bam}_R2.fastq"
  java -jar $picard SamToFastq -I $FILE_PATH -F $R1_out -F2 $R2_out --VALIDATION_STRINGENCY SILENT
}
export -f TO_FASTQ

function COUNT_READS()
{
  PARAMS=("$@")
  FILE_PATH="${PARAMS[0]}"
  file_name=$(basename $FILE_PATH)
  nreads=$(( $(cat $FILE_PATH | wc -l) / 4 ))
  nread_mpxv=$(cat $FILE_PATH | grep '@NC_063383.1' | wc -l)
  printf '%s\t%s\t%s\n' $file_name $nreads $nread_mpxv >> $SUMMARY_FILE
}
export -f COUNT_READS

echo Generating MPXV reads, model NovaSeq
iss generate --genomes "${ref_dir}/${ref_MPXV}" --model NovaSeq -n 0.39m --compress -p 10 --output "${WORK_DIR}"/mpxv_generated_iss_reads

echo Generating rheMacaca reads, model NovaSeq
iss generate --genomes "${ref_dir}/${ref_rheMac}" --model NovaSeq -n 580m --compress -p 10 --output "${WORK_DIR}"/rheMac_generated_iss_reads

echo Merging reads
cat "${WORK_DIR}/rheMac_generated_iss_reads_R1.fastq.gz" "${WORK_DIR}/mpxv_generated_iss_reads_R1.fastq.gz" > "${WORK_DIR}/rheMac_mpxv_combined_generated_iss_reads_R1.fastq.gz"
cat "${WORK_DIR}/rheMac_generated_iss_reads_R2.fastq.gz" "${WORK_DIR}/mpxv_generated_iss_reads_R2.fastq.gz" > "${WORK_DIR}/rheMac_mpxv_combined_generated_iss_reads_R2.fastq.gz"
rm "${WORK_DIR}/rheMac_generated_iss_reads_R1.fastq.gz"
rm "${WORK_DIR}/rheMac_generated_iss_reads_R2.fastq.gz"


echo Alignment for rheMac_MPXV reads
echo rheMac_MPXV reads, model NovaSeq >> $SUMMARY_FILE
echo 390000 reads from MPXV >> $SUMMARY_FILE
echo Align combined generated reads to MPXV
$bwa mem -t 10 "${ref_dir}/${ref_MPXV}" "${WORK_DIR}/rheMac_mpxv_combined_generated_iss_reads_R1.fastq.gz" "${WORK_DIR}/rheMac_mpxv_combined_generated_iss_reads_R2.fastq.gz" \
	| samtools sort --write-index -o "${OUT_DIR}/rheMacMPXV_generated_reads_alignMPXV.bam"
echo -e "Select only mapped reads: exclude read unmapped, not primary alignment, supplementary alignment and with quality more or equal 1"
samtools view -F 4 -F 2048 -F 256 "${OUT_DIR}/rheMacMPXV_generated_reads_alignMPXV.bam" --min-MQ 1 --write-index -b -o "${OUT_DIR}/rheMacMPXV_generated_reads_alignMPXV_OnlyMapped_MQ1.bam"
echo -e "Select only mapped reads: exclude read unmapped, not primary alignment, supplementary alignment and with quality 0"
samtools view -F 4 -F 2048 -F 256 "${OUT_DIR}/rheMacMPXV_generated_reads_alignMPXV.bam" --min-MQ 0 --write-index -b -o "${OUT_DIR}/rheMacMPXV_generated_reads_alignMPXV_OnlyMapped_MQ0.bam"
TO_FASTQ "${OUT_DIR}/rheMacMPXV_generated_reads_alignMPXV_OnlyMapped_MQ1.bam"
TO_FASTQ "${OUT_DIR}/rheMacMPXV_generated_reads_alignMPXV_OnlyMapped_MQ0.bam"
for file in "${OUT_DIR}"/*
do
  if [[ $file == *.fastq ]]; then
    COUNT_READS $file
  fi
done
rm $OUT_DIR/*.bam*
rm $OUT_DIR/*.fastq


echo Align combined generated reads to host
$bwa mem -t 10 "${ref_dir}/${ref_rheMac}" "${WORK_DIR}/rheMac_mpxv_combined_generated_iss_reads_R1.fastq.gz" "${WORK_DIR}/rheMac_mpxv_combined_generated_iss_reads_R2.fastq.gz" \
	| samtools sort --write-index -o "${OUT_DIR}/rheMacMPXV_generated_reads_alignHost.bam" 
echo Select all unmapped reads
samtools view -f 4 "${OUT_DIR}/rheMacMPXV_generated_reads_alignHost.bam" --write-index -b -o "${OUT_DIR}/rheMacMPXV_generated_reads_alignHost_unmapped.bam"
TO_FASTQ "${OUT_DIR}/rheMacMPXV_generated_reads_alignHost_unmapped.bam"
COUNT_READS "${OUT_DIR}/rheMacMPXV_generated_reads_alignHost_unmapped_R1.fastq"
COUNT_READS "${OUT_DIR}/rheMacMPXV_generated_reads_alignHost_unmapped_R2.fastq"
rm $OUT_DIR/*.bam*
rm $OUT_DIR/*.fastq


echo Align combined generated reads to host+MPXV
$bwa mem -t 10 "${ref_dir}/${ref_rheMac_MPXV}" "${WORK_DIR}/rheMac_mpxv_combined_generated_iss_reads_R1.fastq.gz" "${WORK_DIR}/rheMac_mpxv_combined_generated_iss_reads_R2.fastq.gz" \
	| samtools sort --write-index -o "${OUT_DIR}/rheMacMPXV_generated_reads_alignHostMPXV.bam"
echo -e "Selet only mapped reads: exclude read unmapped, not primary alignment, supplementary alignment; with quality more or equal 1 and only NC_063383.1 chromosome"
samtools view -F 4 -F 2048 -F 256 "${OUT_DIR}/rheMacMPXV_generated_reads_alignHostMPXV.bam" NC_063383.1 --min-MQ 1 --write-index -b -o "${OUT_DIR}/rheMacMPXV_generated_reads_alignHostMPXV_MPXVchr_MQ1.bam"
echo -e "Selet only mapped reads: exclude read unmapped, not primary alignment, supplementary alignment; with quality 0 and only NC_063383.1 chromosome"
samtools view -F 4 -F 2048 -F 256 "${OUT_DIR}/rheMacMPXV_generated_reads_alignHostMPXV.bam" NC_063383.1 --min-MQ 0 --write-index -b -o "${OUT_DIR}/rheMacMPXV_generated_reads_alignHostMPXV_MPXVchr_MQ0.bam"
TO_FASTQ "${OUT_DIR}/rheMacMPXV_generated_reads_alignHostMPXV_MPXVchr_MQ1.bam"
TO_FASTQ "${OUT_DIR}/rheMacMPXV_generated_reads_alignHostMPXV_MPXVchr_MQ0.bam"
for file in "${OUT_DIR}"/*
do
  if [[ $file == *.fastq ]]; then
    COUNT_READS $file
  fi
done
rm $OUT_DIR/*.bam*
rm $OUT_DIR/*.fastq

rm "${WORK_DIR}/rheMac_mpxv_combined_generated_iss_reads_R1.fastq.gz"
rm "${WORK_DIR}/rheMac_mpxv_combined_generated_iss_reads_R2.fastq.gz"


echo Generating monkey reads, model NovaSeq
iss generate --genomes "${ref_dir}/${ref_monkey}" --model NovaSeq -n 580m --compress -p 10 --output "${WORK_DIR}"/monkey_generated_iss_reads

echo Merging reads
cat "${WORK_DIR}/monkey_generated_iss_reads_R1.fastq.gz" "${WORK_DIR}/mpxv_generated_iss_reads_R1.fastq.gz" > "${WORK_DIR}/monkey_mpxv_combined_generated_iss_reads_R1.fastq.gz"
cat "${WORK_DIR}/monkey_generated_iss_reads_R2.fastq.gz" "${WORK_DIR}/mpxv_generated_iss_reads_R2.fastq.gz" > "${WORK_DIR}/monkey_mpxv_combined_generated_iss_reads_R2.fastq.gz"
rm "${WORK_DIR}/monkey_generated_iss_reads_R1.fastq.gz"
rm "${WORK_DIR}/monkey_generated_iss_reads_R2.fastq.gz"

echo Alignment for monkey_MPXV reads
echo monkey_MPXV reads, model NovaSeq >> $SUMMARY_FILE
echo 390000 reads from MPXV >> $SUMMARY_FILE
echo Align combined generated reads to MPXV
$bwa mem -t 10 "${ref_dir}/${ref_MPXV}" "${WORK_DIR}/monkey_mpxv_combined_generated_iss_reads_R1.fastq.gz" "${WORK_DIR}/monkey_mpxv_combined_generated_iss_reads_R2.fastq.gz" \
	| samtools sort --write-index -o "${OUT_DIR}/monkeyMPXV_generated_reads_alignMPXV.bam"
echo -e "Select only mapped reads: exclude read unmapped, not primary alignment, supplementary alignment and with quality more or equal 1"
samtools view -F 4 -F 2048 -F 256 "${OUT_DIR}/monkeyMPXV_generated_reads_alignMPXV.bam" --min-MQ 1 --write-index -b -o "${OUT_DIR}/monkeyMPXV_generated_reads_alignMPXV_OnlyMapped_MQ1.bam"
echo -e "Select only mapped reads: exclude read unmapped, not primary alignment, supplementary alignment and with quality 0"
samtools view -F 4 -F 2048 -F 256 "${OUT_DIR}/monkeyMPXV_generated_reads_alignMPXV.bam" --min-MQ 0 --write-index -b -o "${OUT_DIR}/monkeyMPXV_generated_reads_alignMPXV_OnlyMapped_MQ0.bam"
TO_FASTQ "${OUT_DIR}/monkeyMPXV_generated_reads_alignMPXV_OnlyMapped_MQ1.bam"
TO_FASTQ "${OUT_DIR}/monkeyMPXV_generated_reads_alignMPXV_OnlyMapped_MQ0.bam"
for file in "${OUT_DIR}"/*
do
  if [[ $file == *.fastq ]]; then
    COUNT_READS $file
  fi
done
rm $OUT_DIR/*.bam*
rm $OUT_DIR/*.fastq


echo Align combined generated reads to host
$bwa mem -t 10 "${ref_dir}/${ref_rheMac}" "${WORK_DIR}/monkey_mpxv_combined_generated_iss_reads_R1.fastq.gz" "${WORK_DIR}/monkey_mpxv_combined_generated_iss_reads_R2.fastq.gz" \
	| samtools sort --write-index -o "${OUT_DIR}/monkeyMPXV_generated_reads_alignHost.bam" 
echo Select all unmapped reads
samtools view -f 4 "${OUT_DIR}/monkeyMPXV_generated_reads_alignHost.bam" --write-index -b -o "${OUT_DIR}/monkeyMPXV_generated_reads_alignHost_unmapped.bam"
TO_FASTQ "${OUT_DIR}/monkeyMPXV_generated_reads_alignHost_unmapped.bam"
COUNT_READS "${OUT_DIR}/monkeyMPXV_generated_reads_alignHost_unmapped_R1.fastq"
COUNT_READS "${OUT_DIR}/monkeyMPXV_generated_reads_alignHost_unmapped_R2.fastq"
rm $OUT_DIR/*.bam*
rm $OUT_DIR/*.fastq


echo Align combined generated reads to host+MPXV
$bwa mem -t 10 "${ref_dir}/${ref_rheMac_MPXV}" "${WORK_DIR}/monkey_mpxv_combined_generated_iss_reads_R1.fastq.gz" "${WORK_DIR}/monkey_mpxv_combined_generated_iss_reads_R2.fastq.gz" \
	| samtools sort --write-index -o "${OUT_DIR}/monkeyMPXV_generated_reads_alignHostMPXV.bam"
echo -e "Selet only mapped reads: exclude read unmapped, not primary alignment, supplementary alignment; with quality more or equal 1 and only NC_063383.1 chromosome"
samtools view -F 4 -F 2048 -F 256 "${OUT_DIR}/monkeyMPXV_generated_reads_alignHostMPXV.bam" NC_063383.1 --min-MQ 1 --write-index -b -o "${OUT_DIR}/monkeyMPXV_generated_reads_alignHostMPXV_MPXVchr_MQ1.bam"
echo -e "Selet only mapped reads: exclude read unmapped, not primary alignment, supplementary alignment; with quality 0 and only NC_063383.1 chromosome"
samtools view -F 4 -F 2048 -F 256 "${OUT_DIR}/monkeyMPXV_generated_reads_alignHostMPXV.bam" NC_063383.1 --min-MQ 0 --write-index -b -o "${OUT_DIR}/monkeyMPXV_generated_reads_alignHostMPXV_MPXVchr_MQ0.bam"
TO_FASTQ "${OUT_DIR}/monkeyMPXV_generated_reads_alignHostMPXV_MPXVchr_MQ1.bam"
TO_FASTQ "${OUT_DIR}/monkeyMPXV_generated_reads_alignHostMPXV_MPXVchr_MQ0.bam"
for file in "${OUT_DIR}"/*
do
  if [[ $file == *.fastq ]]; then
    COUNT_READS $file
  fi
done
rm $OUT_DIR/*.bam*
rm $OUT_DIR/*.fastq
rm $WORK_DIR/*.fastq.gz
rm $WORK_DIR/*.vcf
rm $WORK_DIR/*.txt


echo Generating MPXV reads, model HiSeq
iss generate --genomes "${ref_dir}/${ref_MPXV}" --model hiseq -n 0.47m --compress -p 10 --output "${WORK_DIR}"/mpxv_generated_iss_reads

echo Generating rheMacaca reads, model HiSeq
iss generate --genomes "${ref_dir}/${ref_MacFas}" --model hiseq -n 720m --compress -p 10 --output "${WORK_DIR}"/MacFas_generated_iss_reads

echo Merging reads
cat "${WORK_DIR}/MacFas_generated_iss_reads_R1.fastq.gz" "${WORK_DIR}/mpxv_generated_iss_reads_R1.fastq.gz" > "${WORK_DIR}/MacFas_mpxv_combined_generated_iss_reads_R1.fastq.gz"
cat "${WORK_DIR}/MacFas_generated_iss_reads_R2.fastq.gz" "${WORK_DIR}/mpxv_generated_iss_reads_R2.fastq.gz" > "${WORK_DIR}/MacFas_mpxv_combined_generated_iss_reads_R2.fastq.gz"
rm "${WORK_DIR}/MacFas_generated_iss_reads_R1.fastq.gz"
rm "${WORK_DIR}/MacFas_generated_iss_reads_R2.fastq.gz"

echo Alignment for MacFas_MPXV reads
echo MacFas_MPXV reads, model HiSeq >> $SUMMARY_FILE
echo 470000 reads from MPXV >> $SUMMARY_FILE
echo Align combined generated reads to MPXV
$bwa mem -t 10 "${ref_dir}/${ref_MPXV}" "${WORK_DIR}/MacFas_mpxv_combined_generated_iss_reads_R1.fastq.gz" "${WORK_DIR}/MacFas_mpxv_combined_generated_iss_reads_R2.fastq.gz" \
	| samtools sort --write-index -o "${OUT_DIR}/MacFasMPXV_generated_reads_alignMPXV.bam"
echo -e "Select only mapped reads: exclude read unmapped, not primary alignment, supplementary alignment and with quality more or equal 1"
samtools view -F 4 -F 2048 -F 256 "${OUT_DIR}/MacFasMPXV_generated_reads_alignMPXV.bam" --min-MQ 1 --write-index -b -o "${OUT_DIR}/MacFasMPXV_generated_reads_alignMPXV_OnlyMapped_MQ1.bam"
echo -e "Select only mapped reads: exclude read unmapped, not primary alignment, supplementary alignment and with quality 0"
samtools view -F 4 -F 2048 -F 256 "${OUT_DIR}/MacFasMPXV_generated_reads_alignMPXV.bam" --min-MQ 0 --write-index -b -o "${OUT_DIR}/MacFasMPXV_generated_reads_alignMPXV_OnlyMapped_MQ0.bam"
TO_FASTQ "${OUT_DIR}/MacFasMPXV_generated_reads_alignMPXV_OnlyMapped_MQ1.bam"
TO_FASTQ "${OUT_DIR}/MacFasMPXV_generated_reads_alignMPXV_OnlyMapped_MQ0.bam"
for file in "${OUT_DIR}"/*
do
  if [[ $file == *.fastq ]]; then
    COUNT_READS $file
  fi
done
rm $OUT_DIR/*.bam*
rm $OUT_DIR/*.fastq


echo Align combined generated reads to host
$bwa mem -t 10 "${ref_dir}/${ref_rheMac}" "${WORK_DIR}/MacFas_mpxv_combined_generated_iss_reads_R1.fastq.gz" "${WORK_DIR}/MacFas_mpxv_combined_generated_iss_reads_R2.fastq.gz" \
	| samtools sort --write-index -o "${OUT_DIR}/MacFasMPXV_generated_reads_alignHost.bam" 
echo Select all unmapped reads
samtools view -f 4 "${OUT_DIR}/MacFasMPXV_generated_reads_alignHost.bam" --write-index -b -o "${OUT_DIR}/MacFasMPXV_generated_reads_alignHost_unmapped.bam"
TO_FASTQ "${OUT_DIR}/MacFasMPXV_generated_reads_alignHost_unmapped.bam"
COUNT_READS "${OUT_DIR}/MacFasMPXV_generated_reads_alignHost_unmapped_R1.fastq"
COUNT_READS "${OUT_DIR}/MacFasMPXV_generated_reads_alignHost_unmapped_R2.fastq"
rm $OUT_DIR/*.bam*
rm $OUT_DIR/*.fastq


echo Align combined generated reads to host+MPXV
$bwa mem -t 10 "${ref_dir}/${ref_rheMac_MPXV}" "${WORK_DIR}/MacFas_mpxv_combined_generated_iss_reads_R1.fastq.gz" "${WORK_DIR}/MacFas_mpxv_combined_generated_iss_reads_R2.fastq.gz" \
	| samtools sort --write-index -o "${OUT_DIR}/MacFasMPXV_generated_reads_alignHostMPXV.bam"
echo -e "Selet only mapped reads: exclude read unmapped, not primary alignment, supplementary alignment; with quality more or equal 1 and only NC_063383.1 chromosome"
samtools view -F 4 -F 2048 -F 256 "${OUT_DIR}/MacFasMPXV_generated_reads_alignHostMPXV.bam" NC_063383.1 --min-MQ 1 --write-index -b -o "${OUT_DIR}/MacFasMPXV_generated_reads_alignHostMPXV_MPXVchr_MQ1.bam"
echo -e "Selet only mapped reads: exclude read unmapped, not primary alignment, supplementary alignment; with quality 0 and only NC_063383.1 chromosome"
samtools view -F 4 -F 2048 -F 256 "${OUT_DIR}/MacFasMPXV_generated_reads_alignHostMPXV.bam" NC_063383.1 --min-MQ 0 --write-index -b -o "${OUT_DIR}/MacFasMPXV_generated_reads_alignHostMPXV_MPXVchr_MQ0.bam"
TO_FASTQ "${OUT_DIR}/MacFasMPXV_generated_reads_alignHostMPXV_MPXVchr_MQ1.bam"
TO_FASTQ "${OUT_DIR}/MacFasMPXV_generated_reads_alignHostMPXV_MPXVchr_MQ0.bam"
for file in "${OUT_DIR}"/*
do
  if [[ $file == *.fastq ]]; then
    COUNT_READS $file
  fi
done
rm $OUT_DIR/*.bam*
rm $OUT_DIR/*.fastq
rm $WORK_DIR/*.fastq.gz
rm $WORK_DIR/*.vcf
rm $WORK_DIR/*.txt
