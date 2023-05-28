#First it is necessary to split the genome reference
fold -b1000 '/home/gloria/Documentos/2.Master_Bioinfomatica_i_Bioestadistica/TFM/RNA_seq/BF_referenceprova.fa' > BF_ref_split.fa

#FastQC =============================================
#This step is done in R studio

```{r}
library(Rsubread)
library(limma)
library(edgeR)
library(fastqcr)

# set base directory variable
base_dir <- "/home/gloria/Documentos/2.Master_Bioinfomatica_i_Bioestadistica/TFM/RNA_seq/data/"

# Using FastQC to analyze the raw sequencing data

#list forward files from the base directory
fastq1 <- list.files(path = file.path(base_dir), pattern = "*1.fq.gz$", full.names = TRUE)
#list reverse files from the base directory
fastq2 <- list.files(path = file.path(base_dir), pattern = "*2.fq.gz$", full.names = TRUE)

#verify if both lists (forward and reverse) lenght is equal
all.equal(length(fastq1),length(fastq2))

#set fastqc directotory variable. This will contain the results of fastqc() function
fastqc_dir <-"/home/gloria/Documentos/2.Master_Bioinfomatica_i_Bioestadistica/TFM/RNA_seq/BFrefgenome/"

#fastqc() function
fastqc(fq.dir = base_dir, qc.dir = fastqc_dir, threads = 4)

#Plotting results of fastaqc() function
qc <- qc_aggregate(fastqc_dir)
qc <- qc_read(file.path(fastqc_dir, "RNA_1_1_fastqc.zip"))
par(mfrow=c(1,2))
qc_plot(qc, "Per base sequence quality")
qc_plot(qc, "Per sequence quality scores")
qc_plot(qc, "Per base sequence content")
qc_plot(qc, "sequence length distribution")

```

#TrimmGalore ==============================================
#to trimm all the files 
ls *_1.fq.gz | xargs -P24 -I@ bash -c 'trim_galore -q 20 --paired -o trimmed "$1" ${1%_1.*.*}_2.fq.gz' _ @


#Genome Aligment ==========================================
#This step is done in R studio

```{r}
# set trimmed directory variable
trimmed_dir <- "/home/gloria/Documentos/2.Master_Bioinfomatica_i_Bioestadistica/TFM/RNA_seq/data/trimmed"

#list forward files from the trimmed directory
reads1 <- list.files(path = file.path(trimmed_dir), pattern = "*_1_val_1.fq.gz$", 
                     full.names = TRUE)

#list reverse files from the trimmed directory
reads2 <- list.files(path = file.path(trimmed_dir), pattern = "*_2_val_2.fq.gz$", 
                     full.names = TRUE)

#verify if both lists (forward and reverse) lenght is equal

all.equal(length(reads1),length(reads2))

#Build index for reference genome
Rsubread::buildindex(basename="BF_ref", 
                     reference= "/home/gloria/Documentos/2.Master_Bioinfomatica_i_Bioestadistica/TFM/RNA_seq/data/BF_ref_split.fa", memory = 3600)

#Align the files to reference genome
Rsubread::align(index = "BF_ref",
                readfile1 = reads1,
                readfile2 = reads2,
                type = "rna",
                input_format = "gzFASTQ",
                output_format = "BAM",
                PE_orientation = "rf",
                nthreads = 6 )
                
```

#Samtools ========================================================
#Use of Samtools to sort (an example)
samtools sort N_SR_C_1_val_1.fq.gz.subread.BAM > N_SR_C_1_val_1_sorted.fq.gz.subread.BAM

#Use of Samtools to index (an example)
samtools index N_SR_C_1_val_1_sorted.fq.gz.subread.BAM


#HT-sq ===========================================================
#Use htseq-count tool to count (an example)
htseq-count -f bam -r pos -i ID -a 10 --type=mRNA -s yes N_SR_C_1_val_1_sorted.fq.gz.subread.BAM /home/gloria/Documentos/2.Master_Bioinfomatica_i_Bioestadistica/TFM/RNA_seq/data/BF_ref_anotation.gff > N_SR_C_1
