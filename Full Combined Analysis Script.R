#combined script to QC, trim, map, count features and perform differential expression analysis of RNA-seq data

#STEP1 - Use rfastp to perfrom QC and trim reads
#install rfastp
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("Rfastp")
library(Rfastp)
library(ggplot2)

#load data and set prefix
outputPrefix <- tempfile(tmpdir = tempdir())
col_ctl_1_1 <- file.path("/Users/roryburke/Documents/RNAseq_raw_data_july22/col_ctl_1_1.fq.gz")
col_ctl_1_2 <- file.path("/Users/roryburke/Documents/RNAseq_raw_data_july22/col_ctl_1_2.fq.gz")

#run rfastp QC on paired reads files
pe_json_report <- rfastp(read1 = col_ctl_1_1, read2 = col_ctl_1_2,
                         outputFastq = paste0(outputPrefix, "_pe"))
#plots of col_ctl_2 pe_json_report
p1 <- curvePlot(pe_json_report)
p1
p2 <- curvePlot(pe_json_report, curve="content_curves")
p2
#dataframe of col_ctl_2 report
dfsummary <- qcSummary(pe_json_report)   

#custom QC with sliding window trimming with average phred of 20+ across 4 bases
outputPrefix <- file.path("/Users/roryburke/Documents/RNAseq_raw_data_july22/trimmed/paired_trimmed_reads/col_ctl_1")
clipr_json_report <- rfastp(read1 = col_ctl_1_1, read2 = col_ctl_1_2,
                            outputFastq = paste0(outputPrefix, '_clipr'),
                            disableTrimPolyG = TRUE,
                            cutLowQualFront = TRUE,
                            cutFrontWindowSize = 4,
                            cutFrontMeanQual = 20
)                        


#STEP2 - Use Rsubread to map and count reads

#install
BiocManager::install("Rsubread")
library("Rsubread")

#build index for align function
ref <- file.path("/Users/roryburke/Documents/RNAseq_raw_data_july22/Arabidopsis_thaliana.Ensembl_TAIR10_release23.genome.fa.gz].fasta.gz")
buildindex(
  
  # basic input/output options
  basename = "subread_index",
  reference = ref,
  
  # options for the details of the index
  gappedIndex = FALSE,
  indexSplit = FALSE,
  memory = 8000,
  TH_subread = 100,
  colorspace = FALSE)

#use align function to map reads
reads1 <- file.path("/Users/roryburke/Documents/RNAseq_raw_data_july22/trimmed/paired_trimmed_reads/col_ctl_1_clipr_R1.fastq.gz")
reads2 <- file.path("/Users/roryburke/Documents/RNAseq_raw_data_july22/trimmed/paired_trimmed_reads/col_ctl_1_clipr_R2.fastq.gz")

align(index = "subread_index", readfile1 = reads1, readfile2 = reads2, output_file="./col_ctl_1_Rsubread_alignment.BAM")

#count mapped reads using featurecounts
col_ctl_1_fc <- featureCounts(files = file.path("/Users/roryburke/Documents/DE_analysis_july22/col_ctl_1_Rsubread_alignment.BAM"),
                              annot.ext = file.path("/Users/roryburke/Documents/RNAseq_raw_data_july22/Arabidopsis_thaliana.Ensembl_TAIR10_release23.exon.gtf.gz.gff"),
                              isGTFAnnotationFile = TRUE,
                              strandSpecific = 0, isPairedEnd = TRUE)

#merge count matrices to one matrix
count_matrix <- cbind(col_ctl_1_fc$counts, col_ctl_1_fc$counts,col_ctl_3_fc$counts, col_HS_1_fc$counts, col_HS_2_fc$counts, col_HS_3_fc$counts)

#Tidy column names
colnames(count_matrix)[1] = "col_ctl_1"
colnames(count_matrix)[2] = "col_ctl_2"
colnames(count_matrix)[3] = "col_ctl_3"
colnames(count_matrix)[4] = "col_HS_1"
colnames(count_matrix)[5] = "col_HS_2"
colnames(count_matrix)[6] = "col_HS_3"

#create dataframe with column data for DESeq2
condition <- c("control", "control", "control", "HS", "HS", "HS") 
type <- c("PE","PE","PE","PE","PE","PE")
column_data <- data.frame(condition, type)
rownames(column_data) <- c("col_ctl_1","col_ctl_2","col_ctl_3","col_HS_1", "col_HS_2", "col_HS_3")

#install DESeq2 and apeglm
BiocManager::install("DESeq2")
library("DESeq2")

BiocManager::install("apeglm")
library("apeglm")

#create deseq2 dataset
dds <- DESeqDataSetFromMatrix(countData = count_matrix,
                              colData = column_data,
                              design = ~ condition)
dds

#trim rows with less tha  10 counts
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

#set factor levels
dds$condition <- factor(dds$condition, levels = c("control","HS"))

#run deseq2 on dds
dds <- DESeq(dds)
res <- results(dds)
res
summary(res)


#change alpha / pval adjusted cutoff to 0.05
res05 <- results(dds, alpha=0.05)
summary(res05)

# shrink by effect size (LFC) and construct MA plot 
resLFC <- lfcShrink(dds, coef="condition_HS_vs_control", type="apeglm")
plotMA(resLFC, ylim=c(-10,10))
