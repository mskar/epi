############################################
# ChIP-seq analysis
## http://bioconductor.org/help/course-materials/2015/BioC2015/csaw_lab.html
############################################

bamDirectory <- "Z:/R_Working_Directory_onc-cbio/TRAC_lecture/ENCODE/chr18"

bamFiles <- file.path(bamDirectory, dir(bamDirectory)[c(1,3,5,7,9,11)])
macsPeakFiles <- file.path("Z:/R_Working_Directory_onc-cbio/TRAC_lecture/MACS_PEAKS", dir("Z:/R_Working_Directory_onc-cbio/TRAC_lecture/MACS_PEAKS"))

# Install and load csaw package
source("https://bioconductor.org/biocLite.R")
biocLite("csaw")
biocLite("GenomicRanges")
biocLite("PICS")
biocLite("Gviz")

#####################################################
# Peak calling by PICS
#####################################################
library("Rsamtools")
library("GenomicRanges")
library("PICS")
{
input1 <- readGAlignments(bamFiles[1])
seqlevels(input1) <- "chr18"
input1 <- as(input1, "GRanges")

input2 <- readGAlignments(bamFiles[2])
seqlevels(input2) <- "chr18"
input2 <- as(input2, "GRanges")

EZH2.1 <- readGAlignments(bamFiles[3])
seqlevels(EZH2.1) <- "chr18"
EZH2.1 <- as(EZH2.1, "GRanges")

EZH2.2 <- readGAlignments(bamFiles[4])
seqlevels(EZH2.2) <- "chr18"
EZH2.2 <- as(EZH2.2, "GRanges")

K27Me3.1 <- readGAlignments(bamFiles[5])
seqlevels(K27Me3.1) <- "chr18"
K27Me3.1 <- as(K27Me3.1, "GRanges")

K27Me3.2 <- readGAlignments(bamFiles[6])
seqlevels(K27Me3.2) <- "chr18"
K27Me3.2 <- as(K27Me3.2, "GRanges")

# Find peaks in EZH2.1
EZH2.1_seg<-segmentPICS(data=EZH2.1, dataC=input1, minReads=1)
library(parallel)
EZH2.1_pics<-PICS(EZH2.1_seg, nCores=2)
EZH2.1_segC<-segmentPICS(data=input1,dataC=EZH2.1,minReads=1)
EZH2.1_picsC<-PICS(EZH2.1_segC, nCores=2)
EZH2.1_fdr<-picsFDR(EZH2.1_pics,EZH2.1_picsC,filter=list(delta=c(50,Inf),se=c(0,50), sigmaSqF=c(0,22500),sigmaSqR=c(0,22500)))
plot(EZH2.1_pics, EZH2.1_picsC, xlim=c(2,17), ylim=c(0,.2), filter=list(delta=c(50,300), se=c(0,50), sigmaSqF=c(0,22500), sigmaSqR=c(0,22500)), type="l")
plot(EZH2.1_fdr[,c(3,1)], type="b", pch=19)

# Find peaks in K27Me3.1
K27Me3.1_seg<-segmentPICS(data=K27Me3.1, dataC=input1, minReads=1)
library(parallel)
K27Me3.1_pics<-PICS(K27Me3.1_seg, nCores=2)
K27Me3.1_segC<-segmentPICS(data=input1,dataC=K27Me3.1,minReads=1)
K27Me3.1_picsC<-PICS(K27Me3.1_segC, nCores=2)
K27Me3.1_fdr<-picsFDR(K27Me3.1_pics,K27Me3.1_picsC,filter=list(delta=c(50,Inf),se=c(0,50), sigmaSqF=c(0,22500),sigmaSqR=c(0,22500)))
plot(K27Me3.1_pics, K27Me3.1_picsC, xlim=c(2,8), ylim=c(0,.2), filter=list(delta=c(50,300), se=c(0,50), sigmaSqF=c(0,22500), sigmaSqR=c(0,22500)), type="l")
plot(K27Me3.1_fdr[,c(3,1)], type="b", pch=19)

# Export Ezh2.1 peaks as bed files and wig files
library(rtracklayer)
myFilter=list(delta=c(50,300),se=c(0,50),sigmaSqF=c(0,22500),sigmaSqR=c(0,22500))
EZH2.1_Bed<-makeRangedDataOutput(EZH2.1_pics,type="bed",filter=c(myFilter,list(score=c(1,Inf))))
export(EZH2.1_Bed,"Z:/R_Working_Directory_onc-cbio/TRAC_lecture/EZH2.1_peaks.bed")

EZH2.1_Wig<-makeRangedDataOutput(EZH2.1_pics,type="wig",filter=c(myFilter,list(score=c(1,Inf))))
export(EZH2.1_Wig,"Z:/R_Working_Directory_onc-cbio/TRAC_lecture/EZH2.1_peaks.wig")

# Export K27Me3.1 peaks as bed files and wig files
myFilter=list(delta=c(50,300),se=c(0,50),sigmaSqF=c(0,22500),sigmaSqR=c(0,22500))
K27Me3.1_Bed<-makeRangedDataOutput(K27Me3.1_pics,type="bed",filter=c(myFilter,list(score=c(1,Inf))))
export(K27Me3.1_Bed,"Z:/R_Working_Directory_onc-cbio/TRAC_lecture/K27Me3.1_peaks.bed")

K27Me3.1_Wig<-makeRangedDataOutput(K27Me3.1_pics,type="wig",filter=c(myFilter,list(score=c(1,Inf))))
export(K27Me3.1_Wig,"Z:/R_Working_Directory_onc-cbio/TRAC_lecture/K27Me3.1_peaks.wig")

}
#####################################################
# Some quality metrics (csaw package)
#####################################################
library("csaw")
# Count reads into windows
frag.len <- 300
window.width <- 10
spacing <- 50


param <- readParam(minq=50)
#which.param <- RangesList(IRanges(1, 78077248)) #defines which chromosome and what region is to be obtained
#names(which.param) <- "chr18"
#param.what <- c("rname", "strand", "pos", "qwidth")
#param.what <- c("qname", "flag")
#defines what fields from the bam file should be obtained
#param <- ScanBamParam(which=which.param, what=param.what) #creates instance of class param to be interpreted by scanBam

EZH2.1_data <- windowCounts(bamFiles[3], ext=frag.len, width=window.width, spacing=spacing, param=param)
EZH2.1_data
assay(EZH2.1_data)

K27me3_data <- windowCounts(bamFiles[3], ext=frag.len, width=window.width, spacing=spacing, param=param)
K27me3_data
assay(K27me3_data)

# Cross Correlation Plots
max.delay <- 500
dedup.on <- readParam(dedup=TRUE, minq=50)
x <- correlateReads(bamFiles[3], max.delay, param=dedup.on)
plot(0:max.delay, x, type="l", ylab="CCF", xlab="Delay (bp)")

#####################################################
# Annotation using ChIPpeakAnno
#####################################################
library("ChIPpeakAnno")
library("biomaRt")
library("RCurl")
human = useMart("ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl", host="www.ensembl.org")
#hg19TSSAnnotation <- getAnnotation(human, featureType="TSS")
hg19TSSAnnotation <- getAnnotation(human, featureType="Exon")
macsEZH2 <- read.table("Z:/R_Working_Directory_onc-cbio/TRAC_lecture/MACS_PEAKS/EZH2_ENCODE_peaks.narrowPeak.bed", header=F, sep="\t", stringsAsFactors=F)
colnames(macsEZH2) <- c('chr','start','end','id','score')
head(macsEZH2)
macsEZH2.Granges=GRanges( seqnames=Rle("chr18", nrow(macsEZH2)), ranges=IRanges( start=macsEZH2$start, end=macsEZH2$end, names=macsEZH2$id) )
annotatedPeak = annotatePeakInBatch(macsEZH2.Granges, AnnotationData=hg19TSSAnnotation)
as.data.frame(annotatedPeak)
write.table(as.data.frame(annotatedPeak), "Z:/R_Working_Directory_onc-cbio/TRAC_lecture/EZh2_MACSPeak_Annotation.txt", sep="\t", quote=F, row.names=F)



#####################################################
# Visualization using gviz
#####################################################
library("Gviz")

cur.region <- GRanges("chr18", IRanges(77806807, 77807165))
collected <- list()
for (i in 1:length(bamFiles)) {
  reads <- extractReads(bamFiles[i], cur.region)
  pcov <- as(coverage(reads[strand(reads)=="+"])/data$totals[i]*1e6, "GRanges")
  ncov <- as(coverage(reads[strand(reads)=="-"])/data$totals[i]*1e6, "GRanges")
  ptrack <- DataTrack(pcov, type="histogram", lwd=0, fill=rgb(0,0,1,.4), ylim=c(0,1),
                      name=bamFiles[i], col.axis="black", col.title="black")
  ntrack <- DataTrack(ncov, type="histogram", lwd=0, fill=rgb(1,0,0,.4), ylim=c(0,1))
  collected[[i]] <- OverlayTrack(trackList=list(ptrack,ntrack))
}
gax <- GenomeAxisTrack(col="black")
plotTracks(c(gax, collected), from=start(cur.region), to=end(cur.region))




#####################################################
#Miscallaneous
#####################################################
x="Z:/R_Working_Directory_onc-cbio/TRAC_lecture/HE_ChIP_seq/C10DEZH2_Chr18_sorted.bam"
frag.len <- 300
window.width <- 10
spacing <- 50
param <- readParam(minq=50)
x_data <- windowCounts(x, ext=frag.len, width=window.width, spacing=spacing, param=param)
max.delay <- 500
dedup.on <- readParam(dedup=TRUE, minq=50)
x <- correlateReads(x, max.delay, param=dedup.on)
plot(0:max.delay, x, type="l", ylab="CCF", xlab="Delay (bp)")

