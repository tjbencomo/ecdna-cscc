#!/bin/bash
#SBATCH --job-name=ecdna
#SBATCH --output=/home/users/tbencomo/cscc-ecdna/workflow/log
#SBATCH --nodes=1
#SBATCH --time=05-00:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=200
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=tbencomo@stanford.edu
#SBATCH --qos=long

set -e
cd /home/users/tbencomo/cscc-ecdna/workflow
snakemake --cluster-config cluster.json -j 499 \
    --rerun-incomplete \
    --rerun-triggers mtime \
    --use-singularity \
    --cluster 'sbatch -p {cluster.partition} -q {cluster.qos} -t {cluster.time} --mem {cluster.mem} -c {cluster.ncpus} -o {cluster.out}'
