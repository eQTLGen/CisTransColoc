#!/bin/bash

#SBATCH --time=01:00:00
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=6G
#SBATCH --mail-type=BEGIN
#SBATCH --mail-type=END
#SBATCH --mail-type=FAIL
#SBATCH --job-name="HyprColoc"

# Here load needed system tools (Java 1.8 is required, one of singularity or anaconda - python 2.7 are needed,
# depending on the method for dependancy management)

#module load jdk/16.0.1
#module load openjdk/11.0.2
#module load squashfs
#module load singularity

set -f

nextflow_path=[path to nextflow executable]

eqtls=[path to eQTL folder]
allele_info=[path to allele info .parquet]
pqtls=[path to pQTL folder]
pqtls_meta=[pQTL meta-table]

genes=[table with genes]
gtf=[path to annotation gtf]

output_folder=[output folder]

NXF_VER=23.04.1 ${nextflow_path}/nextflow run main.nf \
--outdir ${output_folder} \
--eqtl_files ${empirical}/eqtls \
--pqtl_files ${pqtls} \
--allele_info ${allele_info} \
--pqtl_meta_table ${pqtls_meta} \
--genes ${genes} \
--gtf ${gtf} \
--ManhattanPlots true \
-resume \
-profile singularity,slurm
