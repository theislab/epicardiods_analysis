#!/bin/bash
#SBATCH -o scjoint_%A.out
#SBATCH -e scjoint_%A.out
#SBATCH -J scjoint
#SBATCH -p gpu_p
#SBATCH --qos gpu
#SBATCH --gres=gpu:1
#SBATCH -c 4
#SBATCH --mem 80Gb
#SBATCH --time 12:00:00
#SBATCH --nice=10000

source $HOME/.bashrc
env_name=poisson_atac_new
path_to_conda_sh=$(conda info --base)/etc/profile.d/conda.sh
source $path_to_conda_sh

conda activate $env_name

python /storage/groups/ml01/code/laura.martens/github_repos/scJoint/main.py