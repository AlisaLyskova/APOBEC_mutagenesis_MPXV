INPUT_DIR=""
BAM_PATH="/PRJEB56841"
THREADS=10
MODEL_NAME="r941_prom_sup_g5014"

for file in ${INPUT_DIR}${BAM_PATH}/*.bam
do
    mkdir /res_clair/${file##*/}
    OUTPUT_DIR="/res_clair/${file##*/}"
    sudo docker run -it \
  -v ${INPUT_DIR}:${INPUT_DIR} \
  -v ${OUTPUT_DIR}:${OUTPUT_DIR} \
  hkubal/clair3:latest \
  /opt/bin/run_clair3.sh \
  --bam_fn="${file}" \
  --ref_fn="${INPUT_DIR}/NC_063383.1.fna" \
  --threads=${THREADS} \
  --platform="ont" \
  --model_path="/opt/models/${MODEL_NAME}" \
  --output="${OUTPUT_DIR}" \
  --include_all_ctgs \
  --no_phasing_for_fa \
  --chunk_size=500000
done


