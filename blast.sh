#!/bin/bash

#SBATCH --job-name=blast
#SBATCH --mem=200G
#SBATCH --cpus-per-task=20
#SBATCH --account=5-39268
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=mgaupp42@tntech.edu
#SBATCH -o /work/projects/5-39268/mgaupp42/transcriptome_assembly/transdecoder/wb.out
#SBATCH -e /work/projects/5-39268/mgaupp42/transcriptome_assembly/transdecoder/wb.err

spack load blast-plus

cd /work/projects/5-39268/mgaupp42/transcriptome_assembly/transdecoder/

#makeblastdb -in Abar_trinity_spades.95.unitigs.fasta.transdecoder.pep -dbtype prot -parse_seqids -title "Abar_pep" -out Abar_transcriptome_pep

blastp -query wollenberg_proteins_final.txt \
       -db Abar_transcriptome_pep \
       -evalue 1e-5 -seg yes -comp_based_stats 2 -num_alignments 1 \
       -num_threads ${SLURM_CPUS_PER_TASK:-16} \
       -outfmt 6 \
       -out blastp_wollenberg_final.tsv
