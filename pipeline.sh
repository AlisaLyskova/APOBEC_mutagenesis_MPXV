#PATHS for work directory, script and tools
export mosdepth=
export bwa=
export fastQC=
export minimap=
export freebayes=
export star=

export ref_MPXV="NC_063383.1.fa"
export ref_human="GRCh38.d1.vd1.fa"
export ref_human_MPXV="GRCh38.d1.vd1_combined_NC_063383.1.fa"
export ref_monkey_MPXV="monkey_mpox.fa"
export WORK_DIR=
export SCRIPTS_DIR="${WORK_DIR}/scripts"
export ref_dir="${WORK_DIR}/ref"
#table with samples description
export samples_list="${WORK_DIR}/data/samples_description.csv"

#paths to output directories
export FASTQ_DIR="${WORK_DIR}/fastq"
export FASTQC_DIR="${WORK_DIR}/fastqc_results"
mkdir "$FASTQC_DIR"
export fastqc_file="${WORK_DIR}/summary_fastqc_res.txt"
export BAM_DIR="${WORK_DIR}/bam"
mkdir "$BAM_DIR"
#directory for star temporary files
export tmp=
export COVERAGE_DIR="${BAM_DIR}/coverage"
mkdir "$COVERAGE_DIR"
export statistics_summary_file="${WORK_DIR}/statistics_summary_alignment.txt"
printf '%s\t%s\t%s\t%s\t%s\t%s\n' Sample Project Nreads N_mapped_reads Precent_mapped Mean_coverage > $statistics_summary_file
export VCF_DIR="${WORK_DIR}/vcf"
mkdir "$VCF_DIR"
export SIMULATION_DIR="${WORK_DIR}/simulation"
mkdir "$SIMULATION_DIR"
export PICTURES_DIR="${WORK_DIR}/pictures"
mkdir "$PICTURES_DIR"
export circos_dir="${WORK_DIR}/circos"
mkdir "$circos_dir"
export genes_type_file="${WORK_DIR}/data/ucsc_early_late.txt"
export gtf_file="${WORK_DIR}/data/GCF_014621545.1_ASM1462154v1_genomic.gff"
export repeats_file="${WORK_DIR}/data/repeats_MPXV.tsv"

##FASTQC
function FASTQC()
{
  PARAMS=("$@")
  SAMPLE_PATH="${PARAMS[0]}"
  SAMPLE=$(basename ${SAMPLE_PATH%.fastq.gz})
  echo $SAMPLE
  $fastQC $SAMPLE_PATH -o $FASTQC_DIR --memory 1200
  unzip "${FASTQC_DIR}/${SAMPLE}_fastqc.zip" -d $FASTQC_DIR
  text_file="${FASTQC_DIR}/${SAMPLE}_fastqc/summary.txt"
  res=$(cat "$text_file" | grep 'Per base sequence quality' | awk -F '\t' '{print $1}')
  echo $SAMPLE:$res >> $fastqc_file

  mv "${FASTQC_DIR}/${SAMPLE}_fastqc/Images/per_base_quality.png" "${FASTQC_DIR}/${SAMPLE}_per_base_quality.png"
  rm -r "${FASTQC_DIR}/${SAMPLE}_fastqc"
  rm $FASTQC_DIR/*.html
  rm $FASTQC_DIR/*.zip
}
export -f FASTQC

##ALIGNMENT
function ALIGNMENT_STAR()
{
  PARAMS=("$@")
  SAMPLE="${PARAMS[0]}"
  REF_name="${PARAMS[1]}"
  ref="${ref_dir}/star/${REF_name}"
  fastq1="${FASTQ_DIR}/${SAMPLE}_1.fastq.gz"
  fastq2="${FASTQ_DIR}/${SAMPLE}_2.fastq.gz"
  OUTFILE="${BAM_DIR}/${SAMPLE}.bam"
  OUTFILE2="${BAM_DIR}/${SAMPLE}_MPXVchr_MQ0.bam"
  "${star}" --genomeDir $ref --runThreadN 20 -readFilesIn $fastq1 $fastq2 --readFilesCommand "gunzip -c" --outSAMtype BAM SortedByCoordinate --twopassMode Basic --outFileNamePrefix "${BAM_DIR}" --outTmpKeep All --outStd BAM_SortedByCoordinate --limitBAMsortRAM 1034634632 --outTmpDir "${tmp}/STAR/${SAMPLE}" | samtools view --write-index -b -o "$OUTFILE"
  samtools view -F 4 -F 2048 -F 256 "$OUTFILE" NC_063383.1 --min-MQ 0 --write-index -b -o "$OUTFILE2"
}
export -f ALIGNMENT_STAR

function ALIGNMENT_MINIMAP()
{
  PARAMS=("$@")
  SAMPLE="${PARAMS[0]}"
  fastq="${FASTQ_DIR}/${SAMPLE}.fastq.gz"
  REF_name="${PARAMS[1]}"
  ref="${ref_dir}/${REF_name}"
  echo ALIGNMENT_MINIMAP for $SAMPLE
  OUTFILE="${BAM_DIR}/${SAMPLE}.bam"
  OUTFILE2="${BAM_DIR}/${SAMPLE}_MPXVchr_MQ0.bam"
  $minimap -t 10 -a -x splice -Y -C5 --cs --MD -un -G 10000 $ref $fastq | samtools sort | samtools view --write-index -b -o $OUTFILE
  samtools view -F 4 -F 2048 -F 256 $OUTFILE NC_063383.1 --min-MQ 0 --write-index -b -o $OUTFILE2
}
export -f ALIGNMENT_MINIMAP


##ALIGNMENT QUALITY
###INPUT: a bam file aligned to HybridRef and a filtered bam file with Virus
function ALIGNMENT_QUALITY()
{
  PARAMS=("$@")
  SAMPLE="${PARAMS[0]}"
  bam="${BAM_DIR}/${SAMPLE}.bam"
  bam_virus="${BAM_DIR}/${SAMPLE}_virus.bam"
  PROJECT=$(cat $samples_list | grep $SAMPLE | awk -F, '{print $2}') 
  echo -e "Calculate mosdepth for $SAMPLE"
  $mosdepth -t 2 "${COVERAGE_DIR}/${SAMPLE}" $bam
  samtools_outfile="${COVERAGE_DIR}/${SAMPLE}_CoverageSamtoolsStat.txt"
  samtools_outfile_virus="${COVERAGE_DIR}/${SAMPLE}_virus_CoverageSamtoolsStat.txt"
  echo -e "Calculate samtools stat for $SAMPLE"
  samtools stats $bam > $samtools_outfile
  samtools stats $bam_virus > $samtools_outfile_virus
  a=$(grep -w SN $samtools_outfile | grep "raw total sequences" | awk '{print $5}')
  a2=$(grep -w SN $samtools_outfile_virus | grep "raw total sequences" | awk '{print $5}')
  b=$(grep -w SN $samtools_outfile | grep "reads mapped:" | awk '{print $4}')
  #percent mapped reads
  c=$(python3 -c "print(($b / $a)*100)")
  c2=$(python3 -c "print(($a2 / $b)*100)")
  #mean coverage
  d=$(cat "${COVERAGE_DIR}/${SAMPLE}.mosdepth.summary.txt" | grep -w NC_063383.1 | awk '{print $4}')
  printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' $SAMPLE $PROJECT $a $b $c $d $a2 $c2 >> $statistics_summary_file
  
  files=( $(find "${COVERAGE_DIR}" ! -name "*.mosdepth.global.dist.txt") )
  rm "${files[@]}"
}
export -f ALIGNMENT_QUALITY


###change dictionary samples_features in python script mosdepth_distribution.py to change samples names, change hight if there are many samples (1200*500 for <20 samples)
function COVERAGE_PLOT()
{
  mosdepth_files=( $(find "${COVERAGE_DIR}" -name "*.mosdepth.global.dist.txt") )
  echo Create coverage distribution plot
  python3 "${SCRIPTS_DIR}/mosdepth_distribution.py" "${mosdepth_files[@]}" -o "${COVERAGE_DIR}/coverage_dist.html"
}
export -f COVERAGE_PLOT


##Variant calling
###clair3
function VARIANT_CLAIR3()
{
  PARAMS=("$@")
  SAMPLE="${PARAMS[0]}"
  MODEL_NAME="${PARAMS[1]}"
  PLATFORM="${PARAMS[2]}"
  THREADS=5

  BAM_FILE="${BAM}/${SAMPLE}_MPXVchr_MQ0.bam"

  OUT_DIR="${VCF_DIR}/clair3"
  mkdir -p "$OUT_DIR"

  sudo docker run -it \
          -v ${VCF_DIR}:${VCF_DIR} \
          -v ${OUT_DIR}:${OUT_DIR} \
          hkubal/clair3:latest \
          /opt/bin/run_clair3.sh \
          --bam_fn="${BAM_FILE}" \
          --ref_fn="${ref_dir}/${ref_MPXV}" \
          --threads=${THREADS} \
          --platform=$PLATFORM \
          --model_path="/opt/models/${MODEL_NAME}" \
          --output="${OUT_DIR}" \
          --include_all_ctgs \
          --no_phasing_for_fa \
          --chunk_size=500000 \
          --min_mq=0 \
          --snp_min_af=0.005

  mv "${OUT_DIR}/merge_output.vcf.gz" "${VCF_DIR}/${SAMPLE}_clair3.vcf.gz"

  rm -r "$OUT_DIR"
}
export -f VARIANT_CLAIR3


##Distribution of VAF in vcf files
###input - path to directory with vcf files and variant caller tool
###output - VAF plot
function VAF_DIST()
{
  PARAMS=("$@")
  DIR_VCF_FILES="${PARAMS[0]}"
  OUTFILE="${PARAMS[1]}"
  echo Create VAF distribution plot
  python3 "${SCRIPTS_DIR}/VAF_plot.py" "$DIR_VCF_FILES" "$OUTFILE"
}
export -f VAF_DIST


##Plots with all SPNs from VCF file
###input - path to directory with vcf files
###output - Heatmap with SNPs shares for all types of substitutions
function VCF_plot()
{
  PARAMS=("$@")
  VCF_DIR="${PARAMS[0]}"
  OUTDIR="${PARAMS[1]}"
  python3 "${SCRIPTS_DIR}/Heatmap_SNPs.py" "$VCF_DIR" "$OUTDIR" "$statistics_summary_file"
}
export -f VCF_plot

##Plots with 3nucl motives for all types of substitutions
###input - path to VCF file
function SUBST_MOTIVES()
{
  PARAMS=("$@")
  VCF_PATH="${PARAMS[0]}"
  OUTDIR="${PARAMS[1]}"
  python3 "${SCRIPTS_DIR}/SUBST_MOTIVES.py" "$VCF_PATH" "$OUTDIR" "${ref_dir}/${ref_MPXV}"
}
export -f SUBST_MOTIVES

function APOBEC_filter()
{
  PARAMS=("$@")
  VCF="${PARAMS[0]}"
  OUT="${PARAMS[1]}"
  python3 "${SCRIPTS_DIR}/APOBEC_vcf_filter.py" "$VCF" "$OUT" "${ref_dir}/${ref_MPXV}" 
}
export -f APOBEC_filter

function APOBEC_3nucl_filter()
{
  PARAMS=("$@")
  VCF="${PARAMS[0]}"
  OUT="${PARAMS[1]}"
  POS3_NUCL="${PARAMS[2]}"
  python3 "${SCRIPTS_DIR}/APOBEC_3nucl_motif_filter.py" "$VCF" "$OUT" "${ref_dir}/${ref_MPXV}" "$POS3_NUCL"
}
export -f APOBEC_3nucl_filter

#LOGO 3 nucl context for samples from one project
function APOBEC_LOGO()
{
  PARAMS=("$@")
  VCF_DIR="${PARAMS[0]}"
  OUTDIR="${PARAMS[1]}"
  python3 "${SCRIPTS_DIR}/LOGO_APOBEC.py" "$VCF_DIR" "$OUTDIR" "${ref_dir}/${ref_MPXV}" 
}
export -f APOBEC_LOGO

##Calculate number of reads for APOBEC 3 nucleotides motives and compare to all SNPs for all samples in directory
###input - path to directory with vcf files/ path to ref genome, path to outfile with statistics
###output - Number of reads and freq for 3 nucleotides context for all files with percents
function VCF_category_stat()
{
  PARAMS=("$@")
  PROJECT_DIR="${PARAMS[0]}"
  OUTFILE="${PROJECT_DIR}/SNPs_pos_by_category_stat_pc.txt"
  python3 "${SCRIPTS_DIR}/SNPs_by_category.py" "$PROJECT_DIR" "${ref_dir}/${ref_MPXV}" "$OUTFILE"
}
export -f VCF_category_stat

function VCF_category_plot()
{
  PARAMS=("$@")
  VCF_PATH="${PARAMS[0]}"
  OUTDIR="${PARAMS[1]}"
  python3 "${SCRIPTS_DIR}/SNPs_by_category_plot.py" "$VCF_PATH" "$OUTDIR" "${ref_dir}/${ref_MPXV}"
}
export -f VCF_category_plot

##Function to find given substitution
###input - path to directory with vcf files
###output - vcf files with only given SNP and substitution motif
###can change substitution and coordinates of motif in python function
function find_SNP_vcf()
{
  PARAMS=("$@")
  PROJECT_DIR="${PARAMS[0]}"
  python3 "${SCRIPTS_DIR}/find_SNP_vcf.py" $PROJECT_DIR "${ref_dir}/${ref_MPXV}"
}

##Random modulation procedure
###input: sample id, path to bam file, path to vcf file with APOBEC substitutions 
###output: plot with normalized number of reads and potential motif positions for SNPs and for simulated SNPs
###change: motives array
function random_modulation_procedure()
{
  #motives=(TCT TCA TCG TCC)
  #motif=TCA
  CHR="NC_063383.1"
  PARAMS=("$@")
  SAMPLE="${PARAMS[0]}"
  BAM_PATH="${PARAMS[1]}"
  VCF_PATH="${PARAMS[2]}"
  MOTIF="${PARAMS[3]}"
  OUTDIR="${SIMULATION_DIR}/${SAMPLE}"
  mkdir "$OUTDIR"
  
  #find potential targets
  targets_file="${SIMULATION_DIR}/${SAMPLE}/targets_potential_${MOTIF}.csv"
  python3 -c 'from Bio import SeqIO; from Bio.Seq import Seq; import re; import pandas as pd; import sys; CHR = "NC_063383.1"; genome_file = open(sys.argv[1]); motif = sys.argv[2]; genome_dict = SeqIO.to_dict(SeqIO.parse(genome_file, "fasta")); GENOME = str(genome_dict[CHR].seq); targets_coordinates = [i.start() for i in re.finditer(motif, GENOME, flags=0)] + [i.start() for i in re.finditer(str(Seq(motif).reverse_complement()), GENOME, flags=0)]; targets_coordinates.sort(); print(len(targets_coordinates)); out_df = pd.DataFrame({"target_coordinate":targets_coordinates}); out_df.to_csv(sys.argv[3], index=False, header=False)' "${ref_dir}/${ref_MPXV}" "$MOTIF" "$targets_file"
  ##filter targets with 1 or more reads from bam file
  targets=( $(cat "$targets_file") )
  nreads_max_sample=1
  targets_file_filtered="${SIMULATION_DIR}/${SAMPLE}/targets_potential_${MOTIF}_filtered.csv"
  for target in "${targets[@]}"
  do
    new_region=$(python3 -c 'import sys; region="{0}:{1}-{2}"; new_region=region.format(sys.argv[1], int(sys.argv[2])-2, int(sys.argv[2])+2); print(new_region)' "$CHR" "$target")
    nreads=$(samtools view -c "$BAM_PATH" "$new_region")
    if [ "$nreads" -ge 1 ]; then
      echo "$target" >> "$targets_file_filtered"
      if [ "$nreads" -gt "$nreads_max_sample" ]; then
        nreads_max_sample="$nreads"
      fi
    fi
  done
  echo After filtering number of reads is $(wc -l "$targets_file_filtered")
  echo Max N reads "$nreads_max_sample"
  
  python3 "${SCRIPTS_DIR}/simulation.py" "$SAMPLE" "$BAM_PATH" "$VCF_PATH" "$targets_file_filtered" "${ref_dir}/${ref_MPXV}" "$MOTIF" "$OUTDIR" "$nreads_max_sample"

}
export -f random_modulation_procedure

#Searching for codon and amino acids for variants and predicting amino acids for edited allele
##output - bar plot with comparing potential targets and APOBEC substitutions by aa mutation category and plot with aa changes by grantham score
function AA_changes()
{
  PARAMS=("$@")
  SAMPLE_VCF="${PARAMS[0]}"
  RES_DIR="${WORK_DIR}/aa_changes"
  mkdir "${RES_DIR}"
  OUT_SAMPLE="${RES_DIR}/APOBEC_aa.csv"
  OUT_TARGETS="${RES_DIR}/targets_aa.csv"
  
  #python3 "${SCRIPTS_DIR}/aa_subst.py" "$SAMPLE_VCF" "$OUT_SAMPLE" "$OUT_TARGETS" "${ref_dir}/${ref_MPXV}" "$gtf_file"

  #Rscript "${RES_DIR}/grantham_score_APOBEC.csv"

  python3 "${SCRIPTS_DIR}/sankey_aa_changes.py" "${RES_DIR}/grantham_score_APOBEC.csv" "${RES_DIR}/aa_changes.html"
}
export -f AA_changes


#I. RNA-seq data processing
<<"PRJEB56841"
echo PRJEB56841
PROJECT="PRJEB56841"
#samples=( $(cat $samples_list | grep "$PROJECT" | grep ERR109 | awk -F, '{print $1}') )
samples=(ERR10963117)
for sample in "${samples[@]}"
do
  echo $sample
  FASTQ="${FASTQ_DIR}/${sample}.fastq.gz"
  FASTQC $FASTQ
done

for sample in "${samples[@]}"
do
  echo $sample
  ALIGNMENT_MINIMAP "$sample" "$ref_monkey_MPXV"
  ALIGNMENT_QUALITY "$sample"
done

COVERAGE_PLOT
mv "${COVERAGE_DIR}/coverage_dist.html" "${PICTURES_DIR}/${PROJECT}_coverage_dist.html"
rm $COVERAGE_DIR/*

for sample in "${samples[@]}"
mkdir "${VCF_DIR}/${PROJECT}"
do
  echo $sample
  VARIANT_CLAIR3 "$sample" "r941_prom_sup_g5014" "ont"

  VCF="${VCF_DIR}/${sample}_clair3.vcf.gz"
  new_vcf="${VCF_DIR}/${PROJECT}/${sample}.vcf.gz"
  bcftools filter -e "FORMAT/VAF < 0.01" "$VCF" -O z -o "$new_vcf"
  rm "$VCF"
done
outfile="${PICTURES_DIR}/${PROJECT}_VAF.html"
VAF_DIST "${VCF_DIR}/${PROJECT}" "$outfile"
PRJEB56841


<<"PRJEB60728"
echo PRJEB60728
PROJECT="PRJEB60728"
samples=( $(cat $samples_list | grep "$PROJECT" | awk -F, '{print $1}') )
for sample in "${samples[@]}"
do
  FASTQ1="${FASTQ_DIR}/${sample}_1.fastq.gz"
  FASTQ2="${FASTQ_DIR}/${sample}_2.fastq.gz"
  FASTQC $FASTQ1
  FASTQC $FASTQ2
done

for sample in "${samples[@]}"
do
  echo $sample
  ALIGNMENT_STAR "$sample" "$ref_human_MPXV"
  ALIGNMENT_QUALITY "$sample"
done

COVERAGE_PLOT
mv "${COVERAGE_DIR}/coverage_dist.html" "${PICTURES_DIR}/${PROJECT}_coverage_dist.html"
rm $COVERAGE_DIR/*

for sample in "${samples[@]}"
mkdir "${VCF_DIR}/${PROJECT}"
do
  echo $sample
  VARIANT_CLAIR3 "$sample" "ilmn" "ilmn"

  VCF="${VCF_DIR}/${sample}_clair3.vcf.gz"
  new_vcf="${VCF_DIR}/${PROJECT}/${sample}.vcf.gz"
  bcftools filter -e "FORMAT/VAF < 0.01" "$VCF" -O z -o "$new_vcf"
  rm "$VCF"
done
outfile="${PICTURES_DIR}/${PROJECT}_VAF.html"
VAF_DIST "${VCF_DIR}/${PROJECT}" "$outfile"
PRJEB60728


<<"PRJNA1183318"
echo PRJNA1183318
samples=( $(cat $samples_list | grep PRJNA1183318 | awk -F, '{print $1}') )
PRJNA1183318


<<"PRJNA906618"
echo PRJNA906618
PRJNA906618


<<"PRJNA980137"
echo PRJNA980137
PRJNA980137


#II. Searching for APOBEC substitutions

echo Heatmap with all substitutions
VCF_plot "$VCF_DIR" "$PICTURES_DIR"

echo LOGO plot
APOBEC_LOGO "$VCF_DIR" "$PICTURES_DIR"

echo Create plot with 3 nucleotides motives for sample ERR10513574
SUBST_MOTIVES "${VCF_DIR}/${PROJECT}/ERR10513574.vcf.gz" "$PICTURES_DIR"

echo Create plot for categories
VCF_category_plot "${VCF_DIR}/${PROJECT}/ERR10513574.vcf.gz" "$PICTURES_DIR"

echo Calculate number of reads for categories
VCF_category_stat "$VCF_DIR"

echo Create VAF distribution for APOBEC substitutions
mkdir "${VCF_DIR}/APOBEC"
for file in "$VCF_DIR"/*.vcf.gz
do
  filename=$(basename "${file%.vc.gz}")
  APOBEC_vcf="${VCF_DIR}/APOBEC/${filename}_APOBEC.vcf.gz"
  APOBEC_filter "$new_vcf" "$APOBEC_vcf"
  bcftools index "$APOBEC_vcf"
done
outfile="${PICTURES_DIR}/VAF_APOBEC.html"
VAF_DIST "${VCF_DIR}/APOBEC" "$outfile"


<<"Simulation"
VCF_FILE="${VCF_DIR}/ERR10513574_APOBEC.vcf.gz"
BAM_FILE="${BAM_DIR}/ERR10513574_MPXVchr_MQ0.bam"

random_modulation_procedure "ERR10513574" "$BAM_FILE" "$VCF_FILE" "TCA"
random_modulation_procedure "ERR10513574" "$BAM_FILE" "$VCF_FILE" "TCC"
random_modulation_procedure "ERR10513574" "$BAM_FILE" "$VCF_FILE" "TCG"
random_modulation_procedure "ERR10513574" "$BAM_FILE" "$VCF_FILE" "TCT"
Simulation


#echo Processing file for circos map
#circos_dir="/data1/lyskovaa/mpxv_pipeline/circos"
#python3 "${SCRIPTS_DIR}/circos_preprocess_vcf.py" "${VCF_DIR}/${PROJECT}/ERR10513574.vcf.gz" "${circos_dir}/ERR10513574.txt" "${ref_dir}/${ref_MPXV}" "$gtf_file"
#change filenames inside R script
#Rscript "${SCRIPTS_DIR}/circ_genome.R"


echo "Calculate APOBEC density for sample ERR11026637 (24h infected)"
VCF1="${VCF_DIR}/APOBEC/ERR11026637_APOBEC.vcf.gz"
VCF_common_pos=
VCF_all_pos=

python3 "${SCRIPTS_DIR}/APOBEC_density.py" "$VCF1" "$VCF_common_pos" "$VCF_all_pos" "$genes_type_file" "${ref_dir}/${ref_MPXV}" "$gtf_file" "$repeats_file"


#echo Translation APOBEC substitution to coding strand
mkdir "${VCF_DIR}/APOBEC/to_coding_strand"
#python3 "${SCRIPTS_DIR}/APOBEC_subst_strand_translation.py" "${VCF_DIR}/APOBEC/to_coding_strand" "$PICTURES_DIR" "${ref_dir}/${ref_MPXV}" "$gtf_file"


#echo Find amino acids changes
VCF_all_pos=
#AA_changes "$VCF_all_pos"


<<"RNASSELEM"
echo "Secondary structure for APOBEC subatitutions and targets"
OUTDIR="${VCF_DIR}/rnasselem"
mkdir "$OUTDIR"
python3 "${SCRIPTS_DIR}/secondary_structure_predict.py" "$VCF_all_pos" "${ref_dir}/${ref_MPXV}" "$OUTDIR" "clair3"
RNASSELEM


#echo Antisense transcription
#TSS and TES data from https://journals.asm.org/doi/10.1128/msphere.00356-24
#TSS="${WORK_DIR}/data/TSS.xlsx"
#TES="${WORK_DIR}/data/TES.xlsx"
#VCF="${VCF_DIR}/all_APOBEC_subst_clair3.vcf.gz"
#python3 "${SCRIPTS_DIR}/antisense_transcription.py" "$VCF" "${ref_dir}/${ref_MPXV}" "$gtf_file" "$TSS" "$TES"


#echo Specificity matrix
mkdir "${WORK_DIR}/specificity_matrix"
#create file with APOBEC positions with header 'pos'
#echo 'pos' > "${WORK_DIR}/specificity_matrix/all_APOBEC_subst_clair3_pos.txt"
#bcftools view --no-header "$VCF" | awk '{print $2}' >> "${WORK_DIR}/specificity_matrix/all_APOBEC_subst_clair3_pos.txt"
#A3A_Gordenin="${WORK_DIR}/data/S1A.M11_A3A_Gordenin.txt"
#A3B_Gordenin="${WORK_DIR}/data/S1B.M11_A3B_Gordenin.txt"
#python3 "${SCRIPTS_DIR}/script_for_motives_specificity.py" -i "${WORK_DIR}/specificity_matrix/all_APOBEC_subst_clair3_pos.txt" -g "${ref_dir}/${ref_MPXV}" -m "$A3A_Gordenin" -o "${WORK_DIR}/specificity_matrix/matrix_A3A_gordenin.txt"
#python3 "${SCRIPTS_DIR}/script_for_motives_specificity.py" -i "${WORK_DIR}/specificity_matrix/all_APOBEC_subst_clair3_pos.txt" -g "${ref_dir}/${ref_MPXV}" -m "$A3B_Gordenin" -o "${WORK_DIR}/specificity_matrix/matrix_A3B_gordenin.txt"


#echo Genomic samples
#projects=(PRJNA845087 PRJNA981509)

echo Create Venn diagram
VCF_genomic=
VCF_transcriptomic=
n1=$(bcftools index -n "$VCF_transcriptomic")
n2=$(bcftools index -n "$VCF_genomic")
positions_expr=$(bcftools view --no-header "$VCF_genomic" | awk -F '\t' '{print $2}' | tr '\n' '|') 
n3=$(bcftools view --no-header "$VCF_transcriptomic" | awk -F '\t' '{print $2}' | grep -w -E "$positions_expr" | wc -l)
echo $n1
echo $n2
echo $n3
python3 -c 'import sys; from matplotlib_venn import venn2; from matplotlib import pyplot as plt; venn2(subsets = (int(sys.argv[1]), int(sys.argv[2]), int(sys.argv[3])), set_labels = (sys.argv[4],sys.argv[5]), set_colors=("#8AB6F9","#00246B"),alpha=0.7); plt.savefig(sys.argv[6], dpi=800)' "$n1" "$n2" "$n3" "APOBEC позиции из транскриптомных образцов" "APOBEC позиции из геномных образцов" "${PICTURES_DIR}/venn_clair3_MQ5_transcriptomic_genomic_VAF_0_3.png"
