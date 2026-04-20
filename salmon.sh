#!/bin/bash

#SBATCH --job-name=salmon_quant
#SBATCH --mem=200G
#SBATCH --cpus-per-task=30
#SBATCH --account=5-39268
#SBATCH --time=7-00:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=mgaupp42@tntech.edu
#SBATCH -o /work/projects/5-39268/mgaupp42/out/salmon_quant_%a.out
#SBATCH -e /work/projects/5-39268/mgaupp42/error/salmon_quant_%a.err
#SBATCH --array=0-78

spack load /vacau76

#cd /work/projects/5-39268/mgaupp42/transcriptome_assembly/evigene_final/cdhit/
#salmon index -t ./Abar_trinity_spades.95.unitigs.fasta -i index

export MY_FOLDER=/work/projects/5-39268/mgaupp42/diff_express/salmon/
cd ${MY_FOLDER}

# array is zero-indexed, so ${file_array[0] is BOON_WALT_2_2
file_array=(WILS_FG_1_3 WILS_FG_1_8 WILS_FG_1_9 WILS_FG_2_9 WILS_FG_2_10 WILS_FG_2_15 WILS_FG_3_10 WILS_FG_3_11 WILS_FG_3_14 \
  WILS_RW_1_3 WILS_RW_1_4 WILS_RW_1_6 WILS_RW_2_3 WILS_RW_2_4 WILS_RW_2_5 WILS_FG_1_13 WILS_FG_1_15 WILS_FG_1_17 WILS_FG_2_1 \
  WILS_FG_2_2 WILS_FG_2_3 WILS_FG_3_2 WILS_FG_3_3 WILS_FG_3_4 WILS_RW_1_7 WILS_RW_1_8 WILS_RW_1_9 WILS_RW_2_12 WILS_RW_2_13 \
  WILS_RW_2_14 RUTH_LCC_2_2 RUTH_LCC_2_3 RUTH_LCC_2_4 RUTH_LCC_2_8 RUTH_LCC_1_3 RUTH_LCC_1_5 RUTH_LCC_1_8 RUTH_LCC_1_10 \
  RUTH_LCC_1_12 RUTH_LCC_1_14 RUTH_LCC_2_10 RUTH_LCC_2_11 RUTH_LCC_2_12 RUTH_LCC_2_13 BOON_WALT_2_2 BOON_WALT_2_4 \
  BOON_WALT_2_5 BOON_WALT_2_8 BOON_WALT_2_10 BOON_WALT_2_11 BOON_WALT_3_3 BOON_WALT_3_5 BOON_WALT_3_6 \
  BOON_WALT_3_7 BOON_WALT_3_8 BOON_WALT_3_10 SCOT_JAM_1_1 SCOT_JAM_1_2 SCOT_JAM_1_5 SCOT_JAM_1_7 SCOT_JAM_1_8 \
  SCOT_JAM_1_11 SCOT_JAM_2_2 SCOT_JAM_2_5 SCOT_JAM_2_6 SCOT_JAM_2_8 SCOT_JAM_2_10 SCOT_JAM_2_11 \
  PREB_HW_1_3 PREB_HW_1_5 PREB_HW_1_7 PREB_HW_1_8 PREB_HW_1_9 PREB_HW_2_3 PREB_HW_2_5 PREB_HW_2_2 \
  PREB_HW_2_8 PREB_HW_2_9 PREB_HW_2_10)

export sample=${file_array[${SLURM_ARRAY_TASK_ID}]}
export data=/work/projects/5-39268/mgaupp42/transcriptome_assembly/rcorrector/

salmon quant -i ${MY_FOLDER}/index -l A --gcBias -p 20 -o ${MY_FOLDER}/${sample} \
    -1 ${data}${sample}.R1.cor.fq -2 ${data}${sample}.R2.cor.fq 
