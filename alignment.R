# Once we have the fastq files downloaded it's time to do the alignment to the reference genome


## ----GENOME ALIGNMNENT--------------------------------------------------------

# We're going to use the Rsubread package
library(Rsubread)

# 1. First we need to tell Rsubread which files to use
# We set the fastq files directory
fastq_data <- "D:/RNAseq/fastq"


# List the forward files from the trimmed directory
reads1 <- list.files(path = file.path(fastq_data), pattern = "*_1.fastq.gz$", 
                     full.names = T)

# List the reverse files from the trimmed directory
reads2 <- list.files(path = file.path(fastq_data), pattern = "*_2.fastq.gz$", 
                     full.names = T)

# Verify that both lists length is equal
all.equal(length(reads1), length(reads2))



# 2. Now we're going to create an index of the reference genome, in this case the
#Solanum lycopersicum SL4.0 genome

buildindex(basename = "genome_index", reference = "/users/elise/OneDrive/Documentos/MÀSTER BIOINFORMÀTICA I BIOESTADÍSTICA/TFM/RNAseq_analysis/SL4.0.genome.fasta")


# 3. We use the align function to do the alignment
align(index = "genome_index",
      readfile1 = reads1,
      readfile2 = reads2,
      type = "rna",
      input_format = "gzFASTQ",
      output_format = "BAM",
      nthreads = 6)
 
     
 

# Sort the BAM file we got from the alignment using "Rsamtools"  
library(Rsamtools)

bam_in <- "D:/RNAseq/BAM files/SRR12026426_1.fastq.gz.subread.BAM"
bam_out <- "/users/elise/OneDrive/Documentos/MÀSTER BIOINFORMÀTICA I BIOESTADÍSTICA/TFM/RNAseq_analysis/sorted BAM/SRR12026426_1.sorted.fastq.gz.subread.BAM"

sortBam(bam_in, bam_out)


# NOw we create the index of the sorted BAM files

# Define the path to the sorted BAM file
sorted_bam_file <- "./sorted BAM/SRR12026426_1.sorted.fastq.gz.subread.BAM.bam"

# Indexing the BAM file
indexBam(sorted_bam_file)




## ----COUNTING OF THE READS----------------------------------------------------

# For obtaining the counts we will use the function "featureCounts"

# 1. Define the path of the sorted BAM files
bam_files <- list.files(path = "./sorted BAM", pattern = "*sorted.fastq.gz.subread.BAM.bam$", full.names = TRUE)

# 2. Define the path of the .gff file
annotation_file <- "./ITAG4.0_gene_models.gff"

# 3. Obtaining the counts for each gene
counts <- featureCounts(files = bam_files, 
                        annot.ext= annotation_file, 
                        isGTFAnnotationFile=TRUE,
                        isPairedEnd = TRUE,
                        GTF.featureType = "exon",
                        GTF.attrType = "ID",
                        useMetaFeatures = TRUE,
                        nthreads = 6)

head(counts$counts)

# 4. Extracting the counts matrix
write.csv(counts$counts, file = "./gene_counts.csv")      


names(counts)

# Data frame with info on the alignment: assigned reads, unassigned...
counts$stat

# Data frame with features of the annotation: geneID, chr, start, end, strand and length
head(counts$annotation)



# 5. Accessing and plotting the summary statistics
mapping_stats <- counts$stat

library(tidyr)
library(dplyr)
library(ggplot2)

# Reshape data
mapping_stats_long <- mapping_stats %>%
  pivot_longer(cols = -Status, names_to = "Sample", values_to = "Reads") %>%
  filter(Reads > 0)

#Shorten sample names
mapping_stats_long$Sample <- gsub("_.*", "", mapping_stats_long$Sample)



ggplot(mapping_stats_long, aes(x = Sample, y = Reads, fill = Status)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  coord_flip() +
  theme_minimal() +
  labs(title = "Mapped and Unmapped Reads",
       x = "Samples",
       y = "%",
       fill = "Read Status") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))



