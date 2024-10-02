#!/bin/bash

#SBATCH --job-name=pmpnn
#SBATCH --partition=qgpu_3090
#SBATCH -N 1
#SBATCH -n 20
#SBATCH --mem=50G
#SBATCH --gres=gpu:1
#SBATCH --mail-type=end
#SBATCH --mail-user=YOU@MAIL.COM
#SBATCH --output=%j.out
#SBATCH --error=%j.err


singularity exec --nv \
  --bind /hpcfs/fhome/yangchh/tools_deployed/esm-if/input/5YH2.pdb:/workspace/esm-if/data/5YH2.pdb \
  --bind /hpcfs/fhome/yangchh/tools_deployed/esm-if/input/input.csv:/workspace/esm-if/data/input.csv \
  --bind /hpcfs/fhome/yangchh/tools_deployed/esm-if/tmp/results/:/tmp/results/ \
  --bind /hpcfs/fpublic/container/singularity/app/esm-if/checkpoint:/root/.cache/torch/hub/checkpoints \
  --bind /hpcfs/fpublic/container/singularity/app/esm-if/esm-if1/get_mute_fasta.py:/workspace/esm-if/get_mute_fasta.py \
  --bind /hpcfs/fpublic/container/singularity/app/esm-if/esm-if1/score_log_likelihoods.py:/workspace/esm-if/score_log_likelihoods.py \
  run_esm-if.sif \
  /opt/inverse/bin/python /workspace/esm-if/score_log_likelihoods.py \
  /workspace/esm-if/data/5YH2.pdb \
  --mute_csv /workspace/esm-if/data/input.csv \
  --chain C \
  --outpath /tmp/results/scores.csv