#!/bin/bash
#SBATCH -o amulet_%A.out
#SBATCH -e amulet_%A.out
#SBATCH -J amulet
#SBATCH -p cpu_p
#SBATCH -c 32
#SBATCH --mem 80Gb
#SBATCH --time 72:00:00
#SBATCH --nice=10000

source $HOME/.bashrc
env_name=amulet
path_to_conda_sh=$(conda info --base)/etc/profile.d/conda.sh
source $path_to_conda_sh

/lustre/groups/ml01/code/laura.martens/github_repos/amulet/AMULET.sh /lustre/groups/ml01/workspace/laura.martens/moretti_colab/cellranger/merged/outs/fragments.tsv.gz /lustre/groups/ml01/workspace/laura.martens/moretti_colab/cellranger/merged/outs/singlecell.csv /lustre/groups/ml01/code/laura.martens/github_repos/amulet/human_autosomes.txt /lustre/groups/ml01/workspace/laura.martens/moretti_colab/amulet/repeats_hg19.bed /lustre/groups/ml01/workspace/laura.martens/moretti_colab/amulet/ /lustre/groups/ml01/code/laura.martens/github_repos/amulet/