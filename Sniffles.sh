#!/bin/bash
#SBATCH --qos=general-compute
#SBATCH --partition=general-compute
#SBATCH --account=omergokc
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=40
#SBATCH --mem=700G
###SBATCH #--exclusive
#SBATCH --job-name=_T2TSniff_HG1946
#SBATCH --output=_T2TSniff_HG1946.out
#SBATCH --error=_T2TSniff_HG1946.err
#SBATCH --reservation=ubhpc-future
#SBATCH --constraint=AVX512
#SBATCH --export=NONE

#AISO Script (ABOUT, INPUT, STEPS, OUTPUT) for Sniffles2 to identify genomic structural variations using Oxford Nanopre Long-reads and a reference genome.

#ABOUT
#This script will take Nanopore reads in a gzipped fastq file (file.fastq.gz) and a reference geneome in a fasta file (ref.fasta), and output filtered vcf and bed files containing filtered structrual variant calls. The filtering implemented in this script represent the current standard for retaining only high confidence SV calls.
#Christopher Osborne, 05/2023, University at Buffalo
# This version runs on the CCR
#Required software: Minimap2, Sniffles2, samtools, BCFTools, BedTools, and SURVIVOR


#Input
Ref="/projects/academic/omergokc/Kendra/Ancients/CONGA/T2T.fa" #path to reference genome
reads="/projects/academic/omergokc/Luane/HG01946/results_nanoq_HG01946_Aug30/HG01946_20kbreads.fastq.gz" #path to raw Nanopore reads
sample="HG01946_20kbreads" #sample ID
supp="5" #minimum number of reads supporting any given SV call, SVs with fewer reads support will be filtered out. The authors of Sniffles2 suggest 5 reads
maxsupp="60" #maximum number of reads supporting any given SV call. SVs with more reads support will be filtered out. The authors of Sniffles2 suggest that, depending on genome coverage (i.e., starting data amount), any SVs with over 60 reads supporting them are likely erronious. 
minSV=100 #minimum length of SV to be kept after filtering. Most studies use 50 or 100bp
maxSV=500000 #max length of SV to be kept after filtering. The authors of Sniffles2 suggest that SVs greater thatn 100kbp in length are likely false.
threads=40 #number of threads

#Load Sniffles2 Conda Environment (This is a conda environment has all the required software installed in it)
eval "$(/projects/academic/omergokc/own_modules/2023.01/2023.01/software/Core/anaconda3/2023.03-1/bin/conda shell.bash hook)"

#conda activate Sniffles_LL

#Map your long reads to the reference genome using Minimap2

#minimap2 -t ${threads} -ax map-ont --secondary=no ${Ref} ${reads} | samtools view -@ ${threads} -S -b | samtools sort -@ ${threads} -o ${sample}_sorted.bam
#
##No we will remove unmapped (4), non-primary (256) and supplemental (2048) alignments from the resultant bam file as these are likely to cause false detection of SVs.
#samtools view -@ ${threads} -F 2308 -o ${sample}_filtered.bam ${sample}_sorted.bam
#samtools index ${sample}_filtered.bam
##wont be needing this anymore and it's likely a large file. 
#rm ./${sample}_sorted.bam
#
##Now we will actually run Sniffles2 to detect SVs using the bam file we created
#
#sniffles --threads ${threads} \
#--cluster-merge-pos 150 \
#--minsvlen ${minSV} \
#--no-consensus \
#--input ./${sample}_filtered.bam \
#--output-rnames \
#--mosaic \
#--minsupport ${supp} \
#--reference ${Ref} \
#--vcf ${sample}_Sniffles_SV_sup${supp}.vcf \
#--allow-overwrite \
#
#Filter Sniffles Results
#Then next few steps will be used to remove low confidence SVs from the Sniffles output

#Here we use bcftools to filter the Sniffles SV
#It has been recomened that SVs over 100 Kbp cannont be reliably called with Sniffles2, so we get rid of those
#Similarly, it has been discoverd that SVs with over 60 reads supporting them are likely false SVs associated with high mapping depth in repeat regins. So we will get rid of those too. 

conda activate Assembly
module load gcc
module load samtools

bcftools view -i ' SVLEN<100000 && SVLEN>-100000 && SUPPORT<60' ${sample}_Sniffles_SV_sup${supp}.vcf >${sample}_Sniffles_SV_sup${supp}_filtered.vcf

#Next, we will identify regions of the genome where read mapping quality is low (below MP=5) as these have been shown to be associated with false SV calls. 
#Start by making a seperate bam file that is comprised of low quality reads 


samtools view -H ${sample}_filtered.bam > ${sample}_lowMQ.sam

samtools view ${sample}_filtered.bam | awk '$5<5 {print $0}' >>  ${sample}_lowMQ.sam

samtools view -S -b -h ${sample}_lowMQ.sam > ${sample}_lowMQ.bam

#Now we should have a sorted bam file of all the reads that have an MQ <5 (second line awk statement). Next we compute the bp coverage of low map quality reads accross the genome. 
samtools depth ${sample}_lowMQ.bam >  ${sample}_lowMQ.cov

#Then we use SURVIVOR to cluster the coverage track into a bed file for filtering that we can use to remove any SVs that ovelap problematic regions. 

SURVIVOR bincov ${sample}_lowMQ.cov 10 2 > ${sample}_lowMQ.bed

#wont be needign these anymore!
rm ./${sample}_lowMQ.sam
rm ./${sample}_lowMQ.bam

#Finally, we can use use BEDtools to remove any SVs that overlap the low mapping quality regions

bedtools subtract -a ${sample}_Sniffles_SV_sup${supp}_filtered.vcf -b ${sample}_lowMQ.bed -A > temp.vcf

#This step is required to reappend the vcf file header to the vcf file because it somehow glitches when you use the subtract function from bedtools.If you can finde a way to fix this please let me know! The number 4168 below is the number of lines from the first line in filtered sniffles vcf to the vcf header which looks like this "#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  SAMPLE"

###MODIFY THE NUMBER OF LINES

##this creates a dummy file with the header information.
#head -n 4168 ${sample}_Sniffles_SV_sup${supp}_filtered.vcf >temp.file
#
##now to crudly paste it onto the temporary vcf file to give us a fully functioning vcf file!
#cat temp.file temp.vcf >${sample}_Sniffles_SV_sup${supp}_filtered_FINAL.vcf
#
##Take out the rubbish!
#rm ./temp.file
#rm ./temp.vcf
#rm ${sample}_Sniffles_SV_sup${supp}_filtered.vcf
#rm ${sample}_filtered.bam
#
##Bed files are a wildly useful file format that makes comparisons between SV files and between a SV file and a genome annotation easy. 
##Here we convert our final, filtered SVs vcf into a bed file using survivor.
#
#SURVIVOR vcftobed ${sample}_Sniffles_SV_sup${supp}_filtered_FINAL.vcf ${minSV} ${maxSV} ${sample}_Sniffles_SV_sup${supp}_filtered_FINAL.bed
#
#

conda deactivate
