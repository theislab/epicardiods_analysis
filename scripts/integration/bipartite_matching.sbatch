#!/bin/bash
#SBATCH -o matching_%A.out
#SBATCH -e matching_%A.out
#SBATCH -J matching
#SBATCH -p cpu_p
#SBATCH -c 4
#SBATCH --mem 40Gb
#SBATCH --time 01:00:00
#SBATCH --nice=10000

source $HOME/.bashrc
env_name=cellrank
path_to_conda_sh=$(conda info --base)/etc/profile.d/conda.sh
source $path_to_conda_sh

python /home/icb/laura.martens/organoid_regulomes/integration/get_bipartite_matches.py /lustre/groups/ml01/workspace/laura.martens/moretti_colab/signac/adata_xglue_embedding.h5ad -o  /lustre/groups/ml01/workspace/laura.martens/moretti_colab/scglue/matching_k15_99.tsv -s "tech" -k 15 -p 99
