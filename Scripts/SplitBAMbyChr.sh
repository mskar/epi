#!/bin/sh

#  SplitBAMbyChr.sh
#  
#
#  Created by hari easwaran on 11/30/15.
#
################################################
#Set options
################################################
#Set resource requirements
#$ -l mem_free=20G
#$ -l h_vmem=22G

#Name of the job
#$ -N SplitBAMbyChr.sh

#Send email at end of job
#$ -m e
#$ -M heaswar2@jhmi.edu

#$ -e /amber2/scratch/baylin/Hari/TRAC/Scripts/error.txt
#$ -o /amber2/scratch/baylin/Hari/TRAC/Scripts/output.txt



######################################################################
# Download data from ENCODE
######################################################################
#(Do "command + / " for multiline commenting)
## Example of downloading data from ENCODE using command line (in linux)
#wget https://www.encodeproject.org/files/ENCFF000CRM/@@download/ENCFF000CRM.bam -P /amber2/scratch/baylin/Hari/TRAC/ENCODE_ChIP-seq_NHLF
#wget https://www.encodeproject.org/files/ENCFF000CRQ/@@download/ENCFF000CRQ.bam -P /amber2/scratch/baylin/Hari/TRAC/ENCODE_ChIP-seq_NHLF


# Load Samtools
module avail
module load samtools/1.1

######################################################################
# Extract chr18 data from the EZH2 and K27me3 ChIP-seq
######################################################################
directory="/amber2/scratch/baylin/Hari/TRAC/ENCODE_ChIP-seq_NHLF"
outputDirectory="$directory/chr18"
mkdir -vp $outputDirectory
cd $directory
for BAMfile in *.bam
do
echo "sorting and indexing downloaded $BAMfile"
# sort: .bam extension is added to the output file
samtools sort "$directory/$BAMfile" "$outputDirectory/${BAMfile}_sorted"
# index: .bai extension is added to the output file
samtools index "$outputDirectory/${BAMfile}_sorted.bam"

# Next extract/split the chr18 data from the sorted/indexed bam file and sort and index the resulting chr18 bam file
echo "splitting sorted bam file by chr18: $outputDirectory/${BAMfile}_sorted.bam"
samtools view -b "$outputDirectory/${BAMfile}_sorted.bam" chr18 | \
## sort: .bam extension is added to the output file
samtools sort - "$outputDirectory/${BAMfile}_chr18_sorted"
# index: .bai extension is added to the output file
samtools index "$outputDirectory/${BAMfile}_chr18_sorted.bam"

# Get stats
echo "Stats of chr18 data for $outputDirectory/${BAMfile}_chr18_sorted.bam is:"
samtools flagstat "$outputDirectory/${BAMfile}_chr18_sorted.bam"

echo "Removing files $outputDirectory/${BAMfile}_sorted.bam"
rm "$outputDirectory/${BAMfile}_sorted.bam"
rm "$outputDirectory/${BAMfile}_sorted.bam.bai"
done

samtools flagstat /amber2/scratch/baylin/Hari/TRAC/ENCODE_ChIP-seq_NHLF/ENCODE_ChIP-seq_NHLF_K27me3_ENCFF000CRQ.bam


######################################################################
# MACS on ENCODE EZH2 data
######################################################################
module load macs/2.1.0
cd /amber2/scratch/baylin/Hari/TRAC/ENCODE_ChIP-seq_NHLF/chr18
macsOutputDir="/amber2/scratch/baylin/Hari/TRAC/ENCODE_ChIP-seq_NHLF/MACS_EZH2_ENCODE"
mkdir -vp $macsOutputDir

chip1="ENCODE_ChIP-seq_NHLF_EZH2_ENCFF000CQR.bam_chr18_sorted.bam"
chip2="ENCODE_ChIP-seq_NHLF_EZH2_ENCFF000CQS.bam_chr18_sorted.bam"
control1="ENCODE_ChIP-seq_NHLF_Control_ENCFF000CQE.bam_chr18_sorted.bam"
control2="ENCODE_ChIP-seq_NHLF_Control_ENCFF000CQF.bam_chr18_sorted.bam"
macs2 callpeak -t $chip1 $chip2 -c $control1 $control2 -f BAM -g hs -n EZH2_ENCODE --outdir $macsOutputDir  -B -q 0.01

######################################################################
# MACS on ENCODE K27me3 data
######################################################################
macsOutputDir="/amber2/scratch/baylin/Hari/TRAC/ENCODE_ChIP-seq_NHLF/MACS_K27me3_ENCODE"
mkdir -vp $macsOutputDir

chip1="ENCODE_ChIP-seq_NHLF_K27me3_ENCFF000CRQ.bam_chr18_sorted.bam"
chip2="ENCODE_ChIP-seq_NHLF_K27me3_ENCFF000CRM.bam_chr18_sorted.bam"
control1="ENCODE_ChIP-seq_NHLF_Control_ENCFF000CQE.bam_chr18_sorted.bam"
control2="ENCODE_ChIP-seq_NHLF_Control_ENCFF000CQF.bam_chr18_sorted.bam"
macs2 callpeak -t $chip1 $chip2 -c $control1 $control2 -f BAM -g hs -n K27me3_ENCODE --outdir $macsOutputDir  -B -q 0.01
macs2 callpeak -t $chip1 $chip2 -c $control1 $control2 -f BAM --broad -g hs -n Broad_K27me3_ENCODE --outdir $macsOutputDir -B --broad-cutoff 0.1


######################################################################
# MACS on CSC data
######################################################################
#module load macs/2.1.0
#cd /amber2/scratch/baylin/Hari/TRAC/HE_ChIP_seq/chr18
#macsOutputDir="/amber2/scratch/baylin/Hari/TRAC/HE_ChIP_seq/MACS_HE_ChIP_seq"
#mkdir -vp $macsOutputDir
#chip1="C10DEZH2_Chr18_sorted.bam"
#control1="C10DInput_Chr18_sorted.bam"
#chip2="CSC10DEZH2_Chr18_sorted.bam"
#control2="CSC10DInput_Chr18_sorted.bam"
#macs2 callpeak -t $chip1 -c $control1 -f BAM -g hs -n C10DEZH2 --outdir $macsOutputDir  -B -q 0.01
#macs2 callpeak -t $chip2 -c $control2 -f BAM -g hs -n CSC10DEZH2 --outdir $macsOutputDir  -B -q 0.01
#


