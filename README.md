# transcriptformer
transcriptformer

conda activate transcritpformer
srun --partition=scu-gpu --gres=gpu:1 --cpus-per-task=8 --mem=32G --time=2:00:00 --pty bash
pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cu118

