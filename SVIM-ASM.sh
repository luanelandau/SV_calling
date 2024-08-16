#!/bin/bash
#SBATCH --qos=general-compute
#SBATCH --partition=general-compute
#SBATCH --account=omergokc
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=54
#SBATCH --mem=500G
#SBATCH --job-name=SVIM_Opossum
#SBATCH --output=SVIM_Opossum.out
#SBATCH --error=SVIM_Opossum.err
#SBATCH --constraint=AVX512
#SBATCH --export=NONE

eval "$(/projects/academic/omergokc/own_modules/2023.01/2023.01/software/Core/anaconda3/2023.03-1/bin/conda shell.bash hook)"
conda activate Assembly

#input
#hap1="/projects/academic/omergokc/Luane/Opossum1/ragtag/ragtag_output_phased1/ragtag.scaffold.fasta" #path to draft genome
#hap2="/projects/academic/omergokc/Luane/Opossum1/ragtag/ragtag_output/ragtag.scaffold.fasta" #path to draft genome

#HapDup
hap1="/projects/academic/omergokc/Luane/Opossum1/HapDup/hapdup/hapdup_dual_1.fasta" #path to draft genome
hap2="/projects/academic/omergokc/Luane/Opossum1/HapDup/hapdup/hapdup_dual_2.fasta"


num_threads=40 #number threads
sample="Opossum1" #sample name
Ref="/projects/academic/omergokc/Luane/Opossum1/ncbi_dataset/data/GCF_027887165.1/GCF_027887165.1_mMonDom1.pri_genomic.fna"


minimap2 -a -x asm5 --cs -r2k -t ${num_threads} $Ref $hap1 > alignments_hap1.sam
minimap2 -a -x asm5 --cs -r2k -t ${num_threads} $Ref $hap2 > alignments_hap2.sam
samtools sort -m4G -@4 -o alignments_hap1.sorted.bam alignments_hap1.sam
samtools sort -m4G -@4 -o alignments_hap2.sorted.bam alignments_hap2.sam
samtools index alignments_hap1.sorted.bam
samtools index alignments_hap2.sorted.bam
svim-asm diploid /projects/academic/omergokc/Luane/Opossum1/SVIM-ASM alignments_hap1.sorted.bam alignments_hap2.sorted.bam $Ref

conda deactivate
