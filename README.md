# transcriptformer
transcriptformer

conda activate transcritpformer
srun --partition=scu-gpu --gres=gpu:1 --cpus-per-task=8 --mem=32G --time=2:00:00 --pty bash
pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu118

----- installation done sucessfully 
----- try running on your own data


### dataset downloaded from (test data)
wget https://datasets.cellxgene.cziscience.com/ed440808-09de-4fc6-9a1a-5c15056609dd.h5ad


#!/bin/bash
#SBATCH --job-name=tf_inference
#SBATCH --partition=scu-gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=2:00:00
#SBATCH --exclude=g34
#SBATCH --output=tf_inference_%j.out
#SBATCH --error=tf_inference_%j.err

# --- Setup ---
set -euo pipefail

echo "Job started on $(hostname) at $(date)"
nvidia-smi

# Activate conda environment
source /home/sor4003/store_sor4003/anaconda3/etc/profile.d/conda.sh
conda activate transcriptformer

# --- Run inference ---
transcriptformer inference \
  --checkpoint-path /home/sor4003/store_sor4003/transcriptformer/tf_sapiens \
  --data-file /home/sor4003/store_sor4003/transcriptformer/data/human_custom/ed440808-09de-4fc6-9a1a-5c15056609dd.h5ad \
  --output-path /home/sor4003/store_sor4003/transcriptformer/inference_results \
  --batch-size 8 \
  --num-gpus 1

echo "Job finished at $(date)"


###  rs_metedata_extraction_script
It is for subsetting downloaded datasets and its metdata ---only downloaded folder 

### data reformatting 
