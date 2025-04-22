#PATHS for work directory, script and tool
bwa=
picard=

ref_dir=
ref_MPXV=NC_063383.1.fa
ref_human=GRCh38.d1.vd1.fa
ref_human_MPXV=GRCh38.d1.vd1_combined_NC_063383.1.fa

WORK_DIR=
OUT_DIR="${WORK_DIR}/bam_files"
SUMMARY_FILE="${OUT_DIR}/Nreads_summary.txt"


#10 reads generation
#for n in {8..10};
for n in 10
do
  echo -e "Generation #$n"
  echo -e "Generation #$n" >> $SUMMARY_FILE
  iss generate --genomes "${ref_dir}/${ref_MPXV}" --model miseq -n 0.19m --compress -p 10 --output "${WORK_DIR}"/mpxv_generated_iss_reads
  iss generate --genomes "${ref_dir}/${ref_human}" --model miseq -n 300m --compress -p 10 --output "${WORK_DIR}"/human_generated_iss_reads
  echo "Merge generated reads for human and MPXV"
  cat "${WORK_DIR}/human_generated_iss_reads_R1.fastq.gz" "${WORK_DIR}/mpxv_generated_iss_reads_R1.fastq.gz" > "${WORK_DIR}/human_mpxv_combined_generated_iss_reads_R1.fastq.gz"
  cat "${WORK_DIR}/human_generated_iss_reads_R2.fastq.gz" "${WORK_DIR}/mpxv_generated_iss_reads_R2.fastq.gz" > "${WORK_DIR}/human_mpxv_combined_generated_iss_reads_R2.fastq.gz"
  rm "${WORK_DIR}/human_generated_iss_reads_R1.fastq.gz"
  rm "${WORK_DIR}/human_generated_iss_reads_R2.fastq.gz"

  #align on MPXV and save only mapped reads
  echo Align combined generated reads to MPXV
  $bwa mem -t 10 "${ref_dir}/${ref_MPXV}" "${WORK_DIR}/human_mpxv_combined_generated_iss_reads_R1.fastq.gz" "${WORK_DIR}/human_mpxv_combined_generated_iss_reads_R2.fastq.gz" \
	| samtools sort --write-index -o "${OUT_DIR}/HumanMPXV_generated_reads_alignMPXV.bam"
  echo -e "Select only mapped reads: exclude read unmapped, not primary alignment, supplementary alignment and with quality more or equal 1"
  samtools view -F 4 -F 2048 -F 256 "${OUT_DIR}/HumanMPXV_generated_reads_alignMPXV.bam" --min-MQ 1 --write-index -b -o "${OUT_DIR}/HumanMPXV_generated_reads_alignMPXV_OnlyMapped_MQ1.bam"
  echo -e "Select only mapped reads: exclude read unmapped, not primary alignment, supplementary alignment and with quality 0"
  samtools view -F 4 -F 2048 -F 256 "${OUT_DIR}/HumanMPXV_generated_reads_alignMPXV.bam" --min-MQ 0 --write-index -b -o "${OUT_DIR}/HumanMPXV_generated_reads_alignMPXV_OnlyMapped_MQ0.bam"
  rm "${OUT_DIR}/HumanMPXV_generated_reads_alignMPXV.bam"

  #align on human and save not mapped reads
  echo Align combined generated reads to Human
  $bwa mem -t 10 "${ref_dir}/${ref_human}" "${WORK_DIR}/human_mpxv_combined_generated_iss_reads_R1.fastq.gz" "${WORK_DIR}/human_mpxv_combined_generated_iss_reads_R2.fastq.gz" \
	| samtools sort --write-index -o "${OUT_DIR}/HumanMPXV_generated_reads_alignHuman.bam" 
  echo Select all unmapped reads
  samtools view -f 4 "${OUT_DIR}/HumanMPXV_generated_reads_alignHuman.bam" --write-index -b -o "${OUT_DIR}/HumanMPXV_generated_reads_alignHuman_unmapped.bam"
  rm "${OUT_DIR}/HumanMPXV_generated_reads_alignHuman.bam"

  #align on human+mpxv combined ref and save only virus chromosome
  echo Align combined generated reads to Human+MPXV
  $bwa mem -t 10 "${ref_dir}/${ref_human_MPXV}" "${WORK_DIR}/human_mpxv_combined_generated_iss_reads_R1.fastq.gz" "${WORK_DIR}/human_mpxv_combined_generated_iss_reads_R2.fastq.gz" \
	| samtools sort --write-index -o "${OUT_DIR}/HumanMPXV_generated_reads_alignHumanMPXV.bam"
  echo -e "Selet only mapped reads: exclude read unmapped, not primary alignment, supplementary alignment; with quality more or equal 1 and only NC_063383.1 chromosome"
  samtools view -F 4 -F 2048 -F 256 "${OUT_DIR}/HumanMPXV_generated_reads_alignHumanMPXV.bam" NC_063383.1 --min-MQ 1 --write-index -b -o "${OUT_DIR}/HumanMPXV_generated_reads_alignHumanMPXV_MPXVchr_MQ1.bam"
  echo -e "Selet only mapped reads: exclude read unmapped, not primary alignment, supplementary alignment; with quality 0 and only NC_063383.1 chromosome"
  samtools view -F 4 -F 2048 -F 256 "${OUT_DIR}/HumanMPXV_generated_reads_alignHumanMPXV.bam" NC_063383.1 --min-MQ 0 --write-index -b -o "${OUT_DIR}/HumanMPXV_generated_reads_alignHumanMPXV_MPXVchr_MQ0.bam"
  rm "${OUT_DIR}/HumanMPXV_generated_reads_alignHumanMPXV.bam"

  echo Creating fastq files from bam files
  java -jar $picard SamToFastq -I "${OUT_DIR}/HumanMPXV_generated_reads_alignMPXV_OnlyMapped_MQ1.bam" -F "${OUT_DIR}/HumanMPXV_generated_reads_alignMPXV_OnlyMapped_MQ1_R1.fastq" -F2 "${OUT_DIR}/HumanMPXV_generated_reads_alignMPXV_OnlyMapped_MQ1_R2.fastq" --VALIDATION_STRINGENCY SILENT
  java -jar $picard SamToFastq -I "${OUT_DIR}/HumanMPXV_generated_reads_alignMPXV_OnlyMapped_MQ0.bam" -F "${OUT_DIR}/HumanMPXV_generated_reads_alignMPXV_OnlyMapped_MQ0_R1.fastq" -F2 "${OUT_DIR}/HumanMPXV_generated_reads_alignMPXV_OnlyMapped_MQ0_R2.fastq" --VALIDATION_STRINGENCY SILENT
  java -jar $picard SamToFastq -I "${OUT_DIR}/HumanMPXV_generated_reads_alignHuman_unmapped.bam" -F "${OUT_DIR}/HumanMPXV_generated_reads_alignHuman_unmapped_R1.fastq" -F2 "${OUT_DIR}/HumanMPXV_generated_reads_alignHuman_unmapped_R2.fastq" --VALIDATION_STRINGENCY SILENT
  java -jar $picard SamToFastq -I "${OUT_DIR}/HumanMPXV_generated_reads_alignHumanMPXV_MPXVchr_MQ1.bam" -F "${OUT_DIR}/HumanMPXV_generated_reads_alignHumanMPXV_MPXVchr_MQ1_R1.fastq" -F2 "${OUT_DIR}/HumanMPXV_generated_reads_alignHumanMPXV_MPXVchr_MQ1_R2.fastq" --VALIDATION_STRINGENCY SILENT
  java -jar $picard SamToFastq -I "${OUT_DIR}/HumanMPXV_generated_reads_alignHumanMPXV_MPXVchr_MQ0.bam" -F "${OUT_DIR}/HumanMPXV_generated_reads_alignHumanMPXV_MPXVchr_MQ0_R1.fastq" -F2 "${OUT_DIR}/HumanMPXV_generated_reads_alignHumanMPXV_MPXVchr_MQ0_R2.fastq" --VALIDATION_STRINGENCY SILENT

  echo Count number of reads
  for file in "${OUT_DIR}"/*
  do
    if [[ $file == *.fastq ]]; then
      file_name=$(basename $file)
      nreads=$(( $(cat $file | wc -l) / 4 ))
      nread_mpxv=$(cat $file | grep '@NC_063383.1' | wc -l)
      printf '%s\t%s\t%s\n' $file_name $nreads $nread_mpxv >> $SUMMARY_FILE
    fi
  done

  rm $WORK_DIR/*.fastq.gz
  rm $WORK_DIR/*.vcf
  rm $WORK_DIR/*.txt
  rm $OUT_DIR/*.fastq
  rm $OUT_DIR/*.bam*

done