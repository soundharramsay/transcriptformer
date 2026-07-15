# transcriptformer
transcriptformer

conda activate transcritpformer
srun --partition=scu-gpu --gres=gpu:1 --cpus-per-task=8 --mem=32G --time=2:00:00 --pty bash
pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu118

----- installation done sucessfully 
----- try running on your own data


### dataset downloaded from (test data)
wget https://datasets.cellxgene.cziscience.com/ed440808-09de-4fc6-9a1a-5c15056609dd.h5ad
