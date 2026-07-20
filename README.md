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


###### data restructuring 


################################################################################
# TranscriptFormer H5AD Preparation and Inference Documentation
#
# Problem:
# TranscriptFormer requires the gene identifier column specified by
# --gene-col-name (default: ensembl_id) to exist in adata.var columns.
#
# Cellxgene H5AD files often store Ensembl IDs as the AnnData var index:
#
# adata.var.index:
# ENSG00000243485
# ENSG00000237613
# ...
#
# but do NOT contain:
#
# adata.var["ensembl_id"]
#
# Therefore TranscriptFormer fails with:
#
# ValueError:
# Gene column 'ensembl_id' not found in var DataFrame columns
#
################################################################################


################################################################################
# 1. Check original H5AD structure
################################################################################

conda activate transcriptformer

python

import anndata

adata = anndata.read_h5ad(
    "ed440808-09de-4fc6-9a1a-5c15056609dd.h5ad"
)

print(adata)
print("\nVAR columns:")
print(adata.var.columns)

print("\nVAR index:")
print(adata.var.index[:10])


################################################################################
# Expected output before fixing:
#
# VAR columns:
# Index([
# 'feature_types',
# 'Symbol',
# 'feature_name',
# 'feature_reference',
# 'feature_biotype',
# 'feature_length',
# 'feature_type'
# ])
#
# VAR index:
# ENSG00000243485
# ENSG00000237613
#
# Problem:
# Ensembl IDs exist only in the index.
################################################################################



################################################################################
# 2. Create TranscriptFormer-compatible H5AD
#
# Add Ensembl IDs as a real column in adata.var
################################################################################

import anndata

input_file = (
"ed440808-09de-4fc6-9a1a-5c15056609dd.h5ad"
)

output_file = (
"transcriptformer_ready.h5ad"
)


adata = anndata.read_h5ad(input_file)


# Copy gene IDs from var index into a new column
adata.var["ensembl_id"] = adata.var.index.astype(str)


# Remove Ensembl version numbers if present
# Example:
# ENSG00000123456.5 -> ENSG00000123456

adata.var["ensembl_id"] = (
    adata.var["ensembl_id"]
    .str.split(".")
    .str[0]
)


# Verify
print("\nAfter modification:")
print(adata.var.columns)

print("\nExample Ensembl IDs:")
print(adata.var["ensembl_id"].head())


# Save new file
adata.write_h5ad(output_file)

print("\nSaved:")
print(output_file)



################################################################################
# 3. Validate final H5AD
################################################################################

adata = anndata.read_h5ad(
    "transcriptformer_ready.h5ad"
)

assert "ensembl_id" in adata.var.columns

assert len(adata.var["ensembl_id"]) == adata.shape[1]


print("H5AD validation successful")
print(adata)


################################################################################
# 4. Run TranscriptFormer inference
################################################################################

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


source /home/sor4003/store_sor4003/anaconda3/etc/profile.d/conda.sh

conda activate transcriptformer


echo "GPU information"
nvidia-smi


DATA_FILE="/home/sor4003/store_sor4003/transcriptformer/data/human_cancer_rs/transcriptformer_ready.h5ad"


ls -lh "$DATA_FILE"


transcriptformer inference \
 --checkpoint-path /home/sor4003/store_sor4003/transcriptformer/tf_sapiens \
 --data-file "$DATA_FILE" \
 --output-path /home/sor4003/store_sor4003/transcriptformer/inference_results \
 --batch-size 8 \
 --num-gpus 1 \
 --oom-dataloader \
 --n-data-workers 0 \
 --gene-col-name ensembl_id \
 --use-raw auto



################################################################################
# Expected successful log:
#
# Loading vocabulary file
# Building gene vocabulary
# Instantiating TranscriptFormer model
# Model instantiated successfully
# Loading model checkpoint
# Model weights loaded successfully
# Processing AnnData
# Running inference
#
################################################################################



################################################################################
# Important notes:
#
# 1. Do NOT modify X matrix.
#    TranscriptFormer expects raw count data.
#
# 2. Do NOT convert H5AD to Seurat/H5Seurat.
#    Direct AnnData workflow is preferred.
#
# 3. The only required modification was:
#
#    adata.var["ensembl_id"] = adata.var.index
#
# 4. The dataset description was correct:
#
#    "AnnData object in H5AD format containing raw count data"
#
#    The issue was not the counts.
#    The issue was only missing gene metadata column.
#
################################################################################
