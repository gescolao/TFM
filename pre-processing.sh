#First it is necessary to split the genome reference
fold -b1000 '/home/gloria/Documentos/2.Master_Bioinfomatica_i_Bioestadistica/TFM/RNA_seq/BF_referenceprova.fa' > BF_ref_split.fa

#FastQC =============================================
#This step is done in R studio

```{r}
base_dir <- "/home/gloria/Documentos/2.Master_Bioinfomatica_i_Bioestadistica/TFM/RNA_seq/data/"

# Using FastQC to analyze the raw sequencing data
fastq1 <- list.files(path = file.path(base_dir), pattern = "*1.fq.gz$", full.names = TRUE)
fastq2 <- list.files(path = file.path(base_dir), pattern = "*2.fq.gz$", full.names = TRUE)
all.equal(length(fastq1),length(fastq2))

fastqc_dir <-"/home/gloria/Documentos/2.Master_Bioinfomatica_i_Bioestadistica/TFM/RNA_seq/BFrefgenome/"
fastqc(fq.dir = base_dir, qc.dir = fastqc_dir, threads = 4)

qc <- qc_aggregate(fastqc_dir)
qc

qc <- qc_read(file.path(fastqc_dir, "RNA_1_1_fastqc.zip"))
names(qc)

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
trimmed_dir <- "/home/gloria/Documentos/2.Master_Bioinfomatica_i_Bioestadistica/TFM/RNA_seq/data/trimmed"

reads1 <- list.files(path = file.path(trimmed_dir), pattern = "*_1_val_1.fq.gz$", 
                     full.names = TRUE)

reads2 <- list.files(path = file.path(trimmed_dir), pattern = "*_2_val_2.fq.gz$", 
                     full.names = TRUE)

all.equal(length(reads1),length(reads2))


Rsubread::buildindex(basename="BF_ref", 
                     reference= "/home/gloria/Documentos/2.Master_Bioinfomatica_i_Bioestadistica/TFM/RNA_seq/data/BF_ref_split.fa", memory = 3600)


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
#SORT 
samtools sort RNA_33_1_val_1.fq.gz.subread.BAM > RNA_33_1_val_1_sorted.fq.gz.subread.BAM

#INDEX
samtools index RNA_33_1_val_1_sorted.fq.gz.subread.BAM


#HT-sq ===========================================================
htseq-count -f bam -r pos -i ID -a 10 --type=mRNA -s yes RNA_33_1_val_1_sorted.fq.gz.subread.BAM /home/gloria/Documentos/2.Master_Bioinfomatica_i_Bioestadistica/TFM/RNA_seq/data/BF_ref_anotation.gff > RNA_33_1_ok
