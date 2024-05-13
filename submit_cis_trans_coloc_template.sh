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

nextflow_path=[path to Nextflow executable]

empirical=[path to full eQTL files in .parquet format]
sig=[path to significant eQTL file in txt/txt.gz format]
ref=[path to SNP reference file in .parquet format]
gtf=[path to ENSEMBL .gtf file for cis/trans classification]
output_folder=[output folder]

NXF_VER=23.04.1 ${nextflow_path}/nextflow run main.nf \
--sig_eqtls ${sig} \
--eqtl_files ${empirical} \
--allele_info ${ref} \
--gtf ${gtf} \
--p_thresh 5e-8 \
--OutputDir ${output_folder} \
--minN_thresh 6000 \
--maxN_thresh 0.5 \
--i2_thresh 100 \
--leadvar_window  500000 \
--cis_window 500000 \
-profile docker \
-resume

