##########################
## METAGENOMIC ANALYSIS ##
##########################

# Download function gdrive_download:
function gdrive_download () {
 CONFIRM=$(wget --quiet --save-cookies /tmp/cookies.txt 
--keep-session-cookies --no-check-certificate 
"https://docs.google.com/uc?export=download&id=$1" -O- | sed -rn 
's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')
 wget --load-cookies /tmp/cookies.txt 
"https://docs.google.com/uc?export=download&confirm=$CONFIRM&id=$1" -O $2
 rm -rf /tmp/cookies.txt
}

# Installation of khmer:
sudo yum install -y python3-devel gcc-c++ make
conda create --name khmerEnv python=3.6

# Open the terminal and activate the conda environment:
conda activate base

# Create a new conda environment for Khmer:
conda create --name khmerEnv

# Activate the conda environment
conda install -c bioconda khmer

#######
## 1 ##
#######

######################
## QUALITY ANALYSIS ##
######################

# A FastQC analysis is used to assess the quality of genomic sequencing data, such as those generated by platforms like Illumina. It evaluates base quality, 
# base composition, the presence of adapters, sequence duplication levels, read length distribution, and the overrepresentation of sequences.

# Installation of FASTQC 
conda install -c bioconda fastqc

# Quality analysis
conda activate khmerEnv
fastqc *.gz -o ~/metagenomic_analysis/1_FASTQC_RESULTS
conda deactivate

#######
## 2 ##
#######

################
## INTERLEAVE ##
################

# Interleaving in metagenomics is the process of combining two paired-end read files into a single file. In this interleaved file, 
# the forward and reverse reads of each pair are arranged alternately (i.e., the forward read of the pair is followed by its corresponding reverse read). 
# Interleaving is primarily used to facilitate data processing by bioinformatics tools that require paired-end reads to be stored in a single file.

# Unzip the fastq.gz in fastq
for i in {1..32}
do
    gunzip -c ${i}_R1_001.fastq.gz > ${i}_R1_001.fastq
    gunzip -c ${i}_R2_001.fastq.gz > ${i}_R2_001.fastq
done

# Interleaved
for file in *_R1_001.fastq
do
   sample=${file%%_R1_001.fastq}
   echo “interleave-reads.py ${sample}_R1_001.fastq ${sample}_R2_001.fastq -or ${sample}.pe.fq”
done > interleave.sh

cat interleave.sh | parallel

# Remove unnecessary files and organize them
rm -rf *.fastq
cd ...
mkdir 2_INTERLEAVED
cd 0_SAMPLES
mv *.pe.fq ../2_INTERLEAVED
cd ../2_INTERLEAVED

#######
## 3 ##
#######

#######################
## QUALITY FILTERING ##
#######################

# The purpose of this step is to remove low-quality reads from sequencing data, 
# thereby improving the reliability of downstream analyses such as assemblies or annotations. 
# This ensures that the reads used meet a minimum quality standard, reducing errors and artifacts in the final results.

# The filtering process employs the following parameters:
# -Q33: Specifies that the quality scores are encoded in the Phred+33 format, commonly used in Illumina sequencing platforms.
# -q 30: Filters out reads where the average base quality is below 30, corresponding to high-quality bases.
# -p 50: Retains only reads in which at least 50% of the bases meet or exceed the specified quality threshold.

for file in *.pe.fq 
do
  newfile=${file%%.pe.fq}   
  echo "fastq_quality_filter -i ${file} -Q33 -q 30 -p 50 -o ${newfile}.pe.qc.fq"
done > qual_filter.sh

cat qual_filter.sh | parallel

#######
## 4 ##
#######

############################
## REMOVE SHORT SEQUENCES ## 
############################

# Download function gdrive_download:
function gdrive_download () {
 CONFIRM=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 
"https://docs.google.com/uc?export=download&id=$1" -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')
 wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$CONFIRM&id=$1" -O $2
 rm -rf /tmp/cookies.txt
}

# Downlowad filter_fastq_by_length.py script
gdrive_download 1w-OyfdEuMi38utz4cN9g_ng-S9kNeOj9 filter_fastq_by_length.py

# Remove short sequences
for file in *pe.qc.fq
do
  echo "python2.7 filter_fastq_by_length.py ${file} ${file}.cut 50"
done > remove_short.sh

cat remove_short.sh | parallel

#######
## 5 ##
#######

#######################################################
## EXTRACT PAIRED ENDS, RENAME FILES AND MERGE FILES ## 
#######################################################

# In this step, paired-end sequence files are processed after quality cleaning. It includes three main steps: extracting paired reads, 
# removing unnecessary files, and renaming and organizing the output files to facilitate subsequent analysis.

# Extracting paired-ends
for file in *.pe.qc.fq.cut
do
   echo "extract-paired-reads.py ${file}"
done > extract_command.sh

cat extract_command.sh | parallel

# Remove unnecessary files
rm -rf *.tr.qc.fq.cut

# Rename files and merging pe and se files
for file in *.pe
do
   sample=${file%%.pe.qc.fq.cut.pe}
   mv ${file} ${sample}.pe.qc.fq
done

for file in *.se
do
   sample=${file%%.pe.qc.fq.cut.se}
   mv ${file} ${sample}.se.qc.fq
done

#######
## 6 ##
#######

##############################
## PREPARATION FOR ASSEMBLY ##
##############################

# In this step, paired-end sequencing data are prepared for assembly. The script split-paired-reads.py is used to split each paired-end file (*.pe.qc.fq) 
# into two separate files: one containing the forward reads (R1) and the other containing the reverse reads (R2). This step is necessary for assembly tools 
# such as MEGAHIT, which require paired-end reads to be provided in individual files.

for file in *.pe.qc.fq
do
   echo "split-paired-reads.py ${file}"
done > split_command.sh

cat split_command.sh | parallel

#######
## 7 ##
#######

##############
## ASSEMBLY ##
##############

# The assembly step aims to combine sequencing reads (forward, reverse, and unpaired) to reconstruct complete or contiguous genomic sequences 
# from smaller fragments (reads).

mkdir 3_FOR_ASSEMBLY

# Create a file with all forward sequences
cat *.1 > 3_FOR_ASSEMBLY/all.pe.qc.fq.1

# Create a file with all reverse sequences
cat *.2 > 3_FOR_ASSEMBLY/all.pe.qc.fq.2

# Create a file with all unpaired sequences
cat *.se.qc.fq > 3_FOR_ASSEMBLY/all.se.qc.fq

# Assembly using MEGAHIT
megahit -m 0.75 -t 120 -1 all.pe.qc.fq.1 -2 all.pe.qc.fq.2 -r all.se.qc.fq 
-o all.Megahit.assembly

#######
## 8 ##
#######

############################
## ASSEMBLY QUALITY CHECK ##
############################

# In this step, MetaQUAST is used, a tool designed to assess the quality of genomic assemblies. 
# The command takes the final contigs generated by the MEGAHIT assembly (final.contigs.fa) as input. Several evaluation options are specified: 
# --rna-finding identifies potential RNA regions. 
# --conserved-genes-finding searches for conserved genes 
# --max-ref-number 20 limits the maximum number of references for comparison. 

# This analysis allows for the verification of the assembly's quality and integrity.

conda activate quast

# Assembly quality check
metaquast /mnt/DATA/belen/3_FOR_ASSEMBLY/all.Megahit.assembly/final.contigs.fa -t 120 --rna-finding --conserved-genes-finding --max-ref-number 20

# This is the quality report of the samples:

# contigs 	5259264
# contigs (>= 0 bp)	12175284
# contigs (>= 1000 bp)	1377527
# contigs (>= 5000 bp)	46151
# contigs (>= 10000 bp)	10406
# contigs (>= 25000 bp)	1303
# contigs (>= 50000 bp)	229
# Largest contig	213010
# Total length	5241019665
# Total length (>= 0 bp)	7719201062
# Total length (>= 1000 bp)	2613105109
# Total length (>= 5000 bp)	417454729
# Total length (>= 10000 bp)	180508676
# Total length (>= 25000 bp)	52451342
# Total length (>= 50000 bp)	16828383
# N50	997
# N90	566
# auN	2360.7
# L50	1384946
# L90	4272456
# GC (%)	...

#######
## 9 ##
#######

#######################
## GENECALLING - FGS ##
#######################

# FragGeneScan is a program used to predict genes in DNA sequences. 
# The main objective of FragGeneScan is to identify coding sequences (CDS) in DNA sequences that may be fragmented or incomplete. 

mkdir 4_FGS
cd 4_FGS

# Link creation
ln -s /home/kdanielmorais/bioinformatics/tools/fraggenescan/FragGeneScan1.31/train/ ./

# FragGeneScan
FragGeneScan -s ~/metagenomic_analysis/3_FOR_ASSEMBLY/all.Megahit.assembly/final.contigs.fa -w 1 -o belen_MG_Megahit_genecalling_fgs -t complete -p 120

########
## 10 ##
########

#############
## MAPPING ##
#############

# Mapping is a key step in metagenomic analysis. It consists of mapping the DNA or RNA sequences obtained from the sample to a reference database, 
# usually a database of known sequences.

# In this analysis, sequencing reads were mapped against a reference assembly. The process involved several steps, 
# starting with the preparation of the reference and culminating in the generation of alignment statistics. 
# First, an index was created for the reference contigs file (final.contigs.fa) using Bowtie2, optimizing the alignment process. 
# Next, paired-end (.pe.qc.fq) and unpaired (.se.qc.fq) reads were combined into a single file for each sample to facilitate joint mapping. 
# The combined reads were then aligned against the reference using Bowtie2, producing alignment files in SAM format, 
# which were subsequently converted to compressed BAM format using Samtools. 
# Mapped and unmapped reads were counted to assess the quality and efficiency of the alignment. BAM files were sorted by reference position, 
# indexed, and finally, statistics on read distribution across contigs were generated.

REF=final.contigs.fa
reference=${REF%%.fa}
echo "reference is" ${reference}
mkdir ${reference}_build
bowtie2-build ~/metagenomic_analysis/3_FOR_ASSEMBLY/all.Megahit.assembly/${REF} ${reference}_build/${reference}.build

conda activate khmerEnv

for file in *.pe.qc.fq
do
 sample=${file%%.pe.qc.fq}
 cat ${sample}.pe.qc.fq ${sample}.se.qc.fq > ~/metagenomic_analysis/5_SAMPLE_MAPPING/${sample}.all.qc.fq
 echo "processing ${sample}...}"

 bowtie2 -p 70 -x ~/metagenomic_analysis/5_SAMPLE_MAPPING/final.contigs_build/final.contigs.build -q ~/metagenomic_analysis/5_SAMPLE_MAPPING/${sample}.all.qc.fq -S ~/metagenomic_analysis/5_SAMPLE_MAPPING/${sample}.sam 
 echo "sam file is done..."
 
 rm -rf ~/metagenomic_analysis/5_SAMPLE_MAPPING/${sample}.all.qc.fq

 samtools view -Sb ~/metagenomic_analysis/5_SAMPLE_MAPPING/${sample}.sam > ~/metagenomic_analysis/5_SAMPLE_MAPPING/${sample}.bam
 echo "bam file is done..."

 rm -rf ~/metagenomic_analysis/5_SAMPLE_MAPPING/${sample}.sam

 samtools view -c -f 4 ~/metagenomic_analysis/5_SAMPLE_MAPPING/${sample}.bam > ~/metagenomic_analysis/5_SAMPLE_MAPPING/${sample}.reads-unmapped.count.txt
 echo "unmapped reads info done..."

 samtools view -c -F 4 ~/metagenomic_analysis/5_SAMPLE_MAPPING/${sample}.bam > ~/metagenomic_analysis/5_SAMPLE_MAPPING/${sample}.reads-mapped.count.txt
 echo "mapped reads info done..."

 samtools sort -o ~/metagenomic_analysis/5_SAMPLE_MAPPING/${sample}.sorted.bam ~/metagenomic_analysis/5_SAMPLE_MAPPING/${sample}.bam
 echo "bam file was sorted..."

 rm -rf ~/metagenomic_analysis/5_SAMPLE_MAPPING/${sample}.bam

 samtools index ~/metagenomic_analysis/5_SAMPLE_MAPPING/${sample}.sorted.bam
 echo "soerted bam file was indexed..."

 samtools idxstats ~/metagenomic_analysis/5_SAMPLE_MAPPING/${sample}.sorted.bam > ~/metagenomic_analysis/5_SAMPLE_MAPPING/${sample}.reads.by.contigs.txt

 echo "${sample} is done..."
done

# Download count-up-mapped-from-results-txt-with-ctg-length.py script
gdrive_download 1HDB2EF-pq-EJxQxsI1uv1-iVjTl6tAlo count-up-mapped-from-results-txt-with-ctg-length.py

python2.7 count-up-mapped-from-results-txt-with-ctg-length.py *.reads.by.contigs.txt

# Validate the consistency between the assembled contigs and the data generated from the mapping
wc -l summary-count-mapped.tsv
12175286 summary-count-mapped.tsv

grep '>' /mnt/DATA/belen/4_FOR_ASSEMBLY/all.Megahit.assembly/final.contigs.fa | wc -l
12175284

# Coverage is a key metric in genomics, as it indicates how many times a genomic region has been sequenced, 
# providing insights into the reliability of the assembly and the relative abundance of the contigs. 
# The purpose of this step is to generate a coverage file that associates each assembled contig with its average coverage.

# Download function gdrive_download:
function gdrive_download () {
 CONFIRM=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate "https://docs.google.com/uc?export=download&id=$1" -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')
 wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$CONFIRM&id=$1" -O $2
 rm -rf /tmp/cookies.txt
}

# Download get_assembly_coverage.py script
gdrive_download 1S2AQHd2YIjnxZz2kIa2avo1RSAuj-pWT get_assembly_coverage.py

# Obtain coverage of our data
python get_assembly_coverage.py summary-count-mapped.tsv 151 belen_MG_Megahit_assembly_DN_coverage.txt

########
## 11 ##
######## 

######################################
## NORMALISE MAPPING TABLE PER BASE ##
######################################

# In this step, the aim is to normalize the mapping table to adjust coverage based on the length of sequencing reads and contigs. 
# The script normalize-mapping-table-by-read-length-and-ctg-length.py takes as input the mapping count file (summary-count-mapped.tsv) and the average read length 
# (in this case, 151 bases). It generates an output file (TABLE_normalised.txt) where the mapping values are adjusted to provide a more accurate comparative 
# measure of relative coverage, regardless of differences in contig or read lengths. 
# This procedure is essential to correct potential biases arising from variations in lengths and allows for fair comparisons between different contigs or 
# genomic regions.

# Download normalize-mapping-table-by-read-length-and-ctg-length.py script
gdrive_download 1w0bfttjXFZ64NHD8bP7UDaQcS1yd20qR normalize-mapping-table-by-read-length-and-ctg-length.py

python2.7 normalize-mapping-table-by-read-length-and-ctg-length.py summary-count-mapped.tsv 151 TABLE_normalised.txt

# In this step, an additional normalization is performed on the previously normalized table to adjust coverage values based on a predefined scale by columns. 
# The script normalize_table_by_columns.py takes the previously generated file (TABLE_normalised.txt) as input, selects a specific column (in this case, column 2), 
# and applies a normalization factor (1,000,000) to scale the values per sample. The output is saved in a file named TABLE_normalised_per_sample.txt.

# Download normalize_table_by_columns.py script
gdrive_download 1c_fD520xtrCNlUIq9VqqqSvY2OryMXTU normalize_table_by_columns.py

python2.7 normalize_table_by_columns.py TABLE_normalised.txt 2 1000000 TABLE_normalised_per_sample.txt

########
## 12 ##
########

################
## ANNOTATION ## 
################

# Functional and taxonomic annotation are processes used to characterize genetic sequences by assigning biological information and classification. 
# - Functional annotation involves identifying the roles or functions of genes and proteins, such as their involvement in specific pathways, cellular processes, 
# or molecular interactions. 
# - Taxonomic annotation, on the other hand, assigns sequences to their corresponding organisms or taxonomic groups, 
# providing insights into the evolutionary and ecological context of the data. 
# Together, these annotations allow researchers to understand both the biological role and the origin of the sequences, 
# which is critical in fields such as genomics, metagenomics, and molecular biology.

# In this step, gene annotation tasks are performed by integrating alignment results with fungal protein and NCBI databases. 
# First, sample information and genomic sequence data are prepared and organized. 
# Then, these sequences are aligned with fungal proteins and NCBI proteins to obtain the best matches. 
# Subsequently, taxonomic information is added to the alignment results through the download and processing of taxonomy files. 
# The results are formatted, taxonomy tables are combined, and the best matches are selected based on bit score among the annotations. 
# Finally, a table is generated containing taxonomic data and KOG functions of the annotated proteins, completing the annotation and classification process.

cd ..
mkdir 7_ANNOTATION
cp ./6_NORMALISE_MAPPING/TABLE_normalised_per_sample.txt ./7_ANNOTATION/
cd ./7_ANNOTATION

# Download contig_mapping_to_genecall_mapping.py script
gdrive_download 1DakO7roc9C2GJ-SkZuy8AZKTV3QHan14 contig_mapping_to_genecall_mapping.py

python2.7 contig_mapping_to_genecall_mapping.py ~/metagenomic_analysis/4_FGS/belen_MG_Megahit_genecalling_fgs.faa TABLE_normalised_per_sample.txt

# Add “#” to the name of the samples:
head -1 TABLE_normalised_per_sample.txt_genecall.txt| awk -F'\t' '{printf $1"\t"$2 ;for(i=3; i<=NF; ++i) printf "\t%s", "#"$i }' |  awk -F '\t' '{print $0}' > header.txt

tail -n +2  TABLE_normalised_per_sample.txt_genecall.txt > table.txt

cat header.txt table.txt > TABLE_NORM_SAMPLES_GENECALL.txt

#####################################
## JGI FUNGAL PROTEINS - BIOCEV PC ##  
#####################################

cd /mnt/DATA/DATABASES/FUNGAL_PROTEINS_JGI/
cp JGI_FUNGAL_PROTEINS_ANNOTATED_20210312.faa.zip /mnt/DATA/belen/7_ANNOTATION/
cd /mnt/DATA/belen/7_ANNOTATION/
unzip JGI_FUNGAL_PROTEINS_ANNOTATED_20210312.faa.zip

diamond blastp -d ~/metagenomic_analysis/7_ANNOTATION/JGI_FUNGAL_PROTEINS_ANNOTATED_20210312.faa -q ~/metagenomic_analysis/4_FGS/belen_MG_Megahit_genecalling_fgs.faa -e 1E-5 -o genecalling_JGI_FUN_20210312.txt -f 6 -p 120 -b12 -c1

###########################################################
## Total time = 1406.2s                                  ##
## Reported 72585274 pairwise alignments, 72585274 HSPs  ##
## 4541138 queries aligned.                              ##
###########################################################

export LANG=en_US.UTF-8
export LC_ALL=en_US.UTF-8

sort -t$'\t' -k1,1 -k12,12gr -k11,11g -k3,3gr genecalling_JGI_FUN_20210312.txt | sort -u -k1,1 --merge > genecalling_JGI_FUN_20210312_best.txt

# GENERA DEFINED
diamond blastp -d /mnt/DATA/DATABASES/NCBI_nr_DIAMOND/NCBI_nr_20210225_diamond_GENERA -q /mnt/DATA/belen/metagenomic_analysis/4_FGS/belen_MG_Megahit_genecalling_fgs.faa -e 1E-5 -o belen_MG_genecalling_NCBI_nr_PROTEINS_GENERA.txt -f 6 -p 120 -b12 -c1

#############################################################
## Total time = 57412.9s                                   ##
## Reported 280564387 pairwise alignments, 280564387 HSPs. ##
## 12605946 queries aligned.                               ##
#############################################################

export LANG=en_US.UTF-8
export LC_ALL=en_US.UTF-8
sort -t$'\t' -k1,1 -k12,12gr -k11,11g -k3,3gr belen_MG_genecalling_NCBI_nr_PROTEINS_GENERA.txt | sort -u -k1,1 --merge > genecalling_NCBI_nr_PROTEINS_best.txt

#  ADD TAXONOMY TO BLAST RESULTS

function gdrive_download () {
 CONFIRM=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate "https://docs.google.com/uc?export=download&id=$1" -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')
 wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$CONFIRM&id=$1" -O $2
 rm -rf /tmp/cookies.txt
}

# Download jgi_abr_org_list.txt
gdrive_download 12c28kgIw4mPBIhQutNGladdAXwNLtvlR jgi_abr_org_list.txt

# Download replace_fungal_annot_by_taxname.py script
gdrive_download 1XBTtiC1JYl2rzeV7idN2WrveEZknmnQi replace_fungal_annot_by_taxname.py

python2.7 replace_fungal_annot_by_taxname.py genecalling_JGI_FUN_20210312_best.txt jgi_abr_org_list.txt genecalling_JGI_FUNGAL_PROTEINS_best_reformate.txt

# FUNGAL
awk -F'\t' '{print $2}' genecalling_JGI_FUNGAL_PROTEINS_best_reformate.txt | sort | uniq > FUNGAL_NAMES.txt

# NCBI
awk -F'\t' '{print $2}' genecalling_NCBI_nr_PROTEINS_best.txt | sort | uniq > ALL_ACCESSIONS.txt

# Download get_taxonomy_offline.py script
gdrive_download 1o8KmSbwzOsjjeouK3dR0RNmWkMjdfFow get_taxonomy_offline.py

python2.7 get_taxonomy_offline.py ALL_ACCESSIONS.txt /mnt/DATA/DATABASES/NCBI_nr_DIAMOND/ACC2TAXID_nr_current.txt /mnt/DATA/DATABASES/NCBI_nr_DIAMOND/TAXONOMY_TAXID_ALL_fixed.txt taxa_all_accessions.txt 

########################################################################
## accession list loaded...(4106547)                                  ##
## taxonomy loaded...(907158)                                         ##
## taxonomy retrieved... 4106547 vs acc (4106547) - should be equal!  ##
## Done :]                                                            ##
########################################################################

# REFORMAT

# Download replace_acc_by_sp_from_taxonomy.py script
gdrive_download 1jQ3F3ZuA0sBJVy3eJiRhwAaxh31LKqSF replace_acc_by_sp_from_taxonomy.py

python2.7 replace_acc_by_sp_from_taxonomy.py genecalling_NCBI_nr_PROTEINS_best.txt taxa_all_accessions.txt genecalling_NCBI_nr_PROTEINS_best_reformat.txt

######################################################
## number of taxa: 4106547 (4106547)                ##
## DONE :) Processed blast: 12605946 - NOT FOUND 0  ##
######################################################

# COMBINE TAXONOMY TABLES

# Download JGI_TAXA_TAB_2021.txt
gdrive_download 1VtSyy7OutKZ6fAZMTYY2HkZUDNcpfY6V JGI_TAXA_TAB_2021.txt

# Download combine_taxonomy_tables.py script
gdrive_download 1F5p28LpaHrSYWwNI_82V9eHoKIb63y8G combine_taxonomy_tables.py

python2.7 combine_taxonomy_tables.py FUNGAL_NAMES.txt JGI_TAXA_TAB_2021.txt taxa_all_accessions.txt TAX_TAB.tab 

########################################################################
## names loaded...                                                    ##              
## FUNGAL TAXONOMY PROCESSED - NOT FOUND 0 vs. FOUND 1498             ##
## OTHER TAXONOMY PROCESSED - REDUCING TO 34416 vs. ORIGINAL 4106548  ##
## DONE :)                                                            ##
########################################################################

# Download get_best_hit_by_bitscore_multi.py script
gdrive_download 1-3XE5Le8I1_HzQdWbaAHev4ZlrSs6lUi get_best_hit_by_bitscore_multi.py

python2.7 get_best_hit_by_bitscore_multi.py genecalling_NCBI_nr_PROTEINS_best_reformat.txt genecalling_JGI_FUNGAL_PROTEINS_best_reformate.txt

###############################################################################
## FILE: genecalling_NCBI_nr_PROTEINS_best_reformat.txt - HITS: 12605946     ##
## NEW ANNOTATIONS: 12605946 - REPLACED: 0 - CURRENT BEST HITS: 12605946     ##
##                                                                           ##
## FILE: genecalling_JGI_FUNGAL_PROTEINS_best_reformate.txt - HITS: 4541138  ##
## NEW ANNOTATIONS: 1400 - REPLACED: 7597 - CURRENT BEST HITS: 12607346      ##
##                                                                           ##
## done :)                                                                   ##
###############################################################################

awk -F'\t' '{print $2}' best_of_the_blast.txt | sort | uniq > ALL_TAXA_NAMES.txt

# Download get_taxonomy_basedonnames.py script
gdrive_download 1XruvN2qGN2-dUHZNSn0jXOYmUxoJ3Uz0 get_taxonomy_basedonnames.py

python2.7 get_taxonomy_basedonnames.py ALL_TAXA_NAMES.txt TAX_TAB.tab TAX_TAB_FINAL.tab 

#######################################################
## names loaded...                                   ##
## TAXONOMY PROCESSED - NOT FOUND 0 vs. FOUND 35600  ##
## DONE :)                                           ##
#######################################################

awk -F'\t' '{print $1"\t"$12"\t"$2""}' best_of_the_blast.txt > TAXONOMY_BEST_OF_SIMPLE.txt

# KOGG FROM JGI-MYCO-GENOMES

awk -F'[|\t]' '{print $1"\t"$15"\t["$5"]"}' genecalling_JGI_FUN_20210312_best.txt > FUNCTION_JGI_KOG_SIMPLE.txt

########
## 13 ##
########

##############
## SPLIT IT ##
##############

# The annotated genomic sequences in the FASTA file are divided into smaller groups. 
# This is performed using a script that fragments the file belen_MG_Megahit_genecalling_fgs.faa into defined-sized parts (in this case, 83,000 sequences per group). 
# This segmentation facilitates the management and processing of large volumes of genomic data.

# Download split_fasta_by_group_size.py script 
gdrive_download 1mGbdx3OBumymosW24WaYfZT9nq_a1z1z split_fasta_by_group_size.py

python2.7 split_fasta_by_group_size.py /mnt/DATA/belen/metagenomic_analysis_chapter_4/4_FGS/belen_MG_Megahit_genecalling_fgs.faa 83000

cd ..
mkdir 8_SPLIT 
cd ./7_ANNOTATION
mv *.fas ../8_SPLIT
cd ../8_SPLIT/

########
## 14 ## 
########

######################
## dbCAN ANNOTATION ##
######################

# In this step, the annotation of CAZy (Carbohydrate-Active Enzymes) is performed using the local dbCAN database. 
# First, the FASTA files are processed with the script run_dbcan.py, which searches for Hidden Markov Model (HMM) profiles within the local dbCAN database. 
# The analysis is executed in parallel to optimize processing. The results are consolidated into a single file (all_dbCAN.txt), 
# from which the best matches are selected based on e-value to generate a filtered file (all_dbCAN_best.txt). 
# Finally, unique gene names are extracted, and a simplified table (CAZy_BEST_SIMPLE.txt) is created, 
# containing the best annotations and identifying carbohydrate-active enzymes present in the samples.

# dbCAN local database
conda activate run_dbcan

for file in *.fas
do
 sample=${file%%.fas}
 mkdir ${sample}
done

for file in *.fas
do
  sample=${file%%.fas}
  echo "run_dbcan.py ${file} protein --db_dir /mnt/DATA/DATABASES/run_dbcan_master/db/ -t hmmer --out_dir ${sample} --hmm_cpu 1 --dia_cpu 1"
done > dbcan.sh

cat dbcan.sh | parallel

echo "" > all_dbCAN.txt
for file in *.fas
do
 sample=${file%%.fas}
 wc -l ${sample}/hmmer.out
 cat ${sample}/hmmer.out >> all_dbCAN.txt
done

# dbCAN annotation
export LC_ALL=en_US.UTF-8
export LANG=en_US.UTF-8

sort -t$'\t' -k3,3 -k5,5g all_dbCAN.txt | sort -u -k3,3 --merge > all_dbCAN_best.txt

awk -F'[.\t]' '{print $1}' all_dbCAN_best.txt | sort | uniq > hmm_names_uniq.txt

awk -F'[.\t]' '{print $1}' all_dbCAN_best.txt > hmm_names.txt

awk -F'\t' '{print $3"\t"$5}' all_dbCAN_best.txt > all_dbCAN_best_gene_eval.txt

paste -d"\t" all_dbCAN_best_gene_eval.txt hmm_names.txt > CAZy_BEST_SIMPLE.txt

########
## 15 ##
########

#################
## KOFAM - KOs ##
#################

# In this step, functional annotation is performed using the KOfam database, which assigns KEGG Orthology (KO) functions to genes based on HMM profiles. 
# The workflow begins by organizing directories for results (ko_tbl) and temporary files (tmp). 
# The hmmsearch tool is then used to compare KOfam HMM profiles against genomic FASTA sequences, with tasks executed in parallel for efficiency. 
# All results are merged into a single file (KOFAM_all.out.txt).
# A Python script is then used to filter the results based on predefined thresholds and e-values, 
# producing a final table (hmmsearch_KOFAM_multi_best.txt) containing the most reliable functional annotations, linking genes to their corresponding biological roles.

mkdir ko_tbl
mkdir tmp

for i in /mnt/DATA1/priscila/kofamKOALA/db/profiles/*.hmm
do
  file=${i##*/}
  ko=${file%%.hmm}
  echo "hmmsearch --tblout ko_tbl/${ko}.out.txt --noali --cpu 1 -E 1e-5 ${i} ~/metagenomic_analysis/4_FGS/belen_MG_Megahit_genecalling_fgs.faa >/dev/null 2>&1"
done > hmmsearch_kofam.sh

cat hmmsearch_kofam.sh | parallel -j 70 --tmpdir tmp

cat ko_tbl/*.out.txt > KOFAM_all.out.txt

python2.7 kegg_multi_from_kofamkoala_raw_filterby_thresholds_evalues.py KOFAM_all.out.txt /mnt/DATA1/priscila/kofamKOALA/db/ko_list hmmsearch_KOFAM_multi_best.txt

########
## 16 ## 
########

#########################
## KEGG AND dbCAN tree ##
#########################

# In this step, the KEGG ontology tree is generated and filtered for unique KOs (KEGG Orthologies) based on the functional annotations from the previous step. 
# First, a script is used to extract unique KOs from the KEGG annotations (KO_UNIQUE_from_KO_simple.py). 
# Then, the KEGG ontology table (kegg_tab.txt) is processed to retain only the KOs present in the data, creating a subtable with relevant KOs (KOFAM_KOs_tree.tab).

# Similarly, a tree for CAZy (Carbohydrate-Active Enzymes) is generated. 
# Unique CAZy identifiers are extracted from the annotations and a script (get_CAZy_tree.py) is used to create a CAZy-specific ontology tree (CAZy_tree.tab),
# which provides a structured representation of the identified carbohydrate-active enzymes in the data.

# Download kegg_tab.txt
gdrive_download 11CVgwqy6O2mJ5rc4vxevYQVp04oJrQsl kegg_tab.txt

# Download GET_KEGG_ontology_subtable.py script
gdrive_download 1AmOiMqLHE8nberbvYEDPT12JShAfSW_w GET_KEGG_ontology_subtable.py

# Download KO_UNIQUE_from_KO_simple.py
gdrive_download 1aqENUDPh3dCoxnO8YfBZaUjwwxkmlMM1 KO_UNIQUE_from_KO_simple.py

python2.7 KO_UNIQUE_from_KO_simple.py hmmsearch_KOFAM_multi_best.txt

python2.7 GET_KEGG_ontology_subtable.py kegg_tab.txt hmmsearch_KOFAM_multi_best.txt.unique.txt KOFAM_KOs_tree.tab

awk -F'\t' '{print $3}' ~/metagenomic_analysis/8_SPLIT/CAZy_BEST_SIMPLE.txt | sort | uniq > CAZy_BEST_unique.txt

# Download get_CAZy_tree.py script
gdrive_download 1SGVK2cqWCLozEPGNLvPRs0YG-ckrF_CZ get_CAZy_tree.py

python2.7 get_CAZy_tree.py CAZy_BEST_unique.txt CAZy_tree.tab

########
## 17 ##
########

##############################
## LINK ANNOTATION TO TABLE ##
##############################

# In this step, the annotation results are linked to the sample table, integrating multiple sources of functional and taxonomic data. 
# We obtain the final tables of the metagenomic analysis, where we can appreciate the abundance of each sample, the taxonomy and associated functionality. 

# Download link_simple_table_to_mapping_table.py script
gdrive_download 198TDGsV1cBfLEZorb5znFHysG47XEj5t link_simple_table_to_mapping_table.py

cd ../7_ANNOTATION/
mv TABLE_NORM_SAMPLES_GENECALL.txt ~/metagenomic_analysis/8_SPLIT/
cp TAXONOMY_BEST_OF_SIMPLE.txt ~/metagenomic_analysis/belen/8_SPLIT/
cd ../8_SPLIT

TABLE="~/metagenomic_analysis/8_SPLIT/TABLE_NORM_SAMPLES_GENECALL.txt"

echo "${TABLE}"

TABLE_BASE=${TABLE%%.${TABLE##*.}}

echo "${TABLE_BASE}"

python2.7 link_simple_table_to_mapping_table.py ${TABLE} TAXONOMY_BEST_OF_SIMPLE.txt TAX_BEST bitscore ${TABLE_BASE}_TAX.txt

python2.7 link_simple_table_to_mapping_table.py ${TABLE_BASE}_TAX.txt CAZy_BEST_SIMPLE.txt CAZy e-val ${TABLE_BASE}_TAX_CAZy.txt

python2.7 link_simple_table_to_mapping_table.py ${TABLE_BASE}_TAX_CAZy.txt ../7_ANNOTATION/FUNCTION_JGI_KOG_SIMPLE.txt KOG e-val ${TABLE_BASE}_TAX_CAZy_KOG.tab

python2.7 link_simple_table_to_mapping_table.py ${TABLE_BASE}_TAX_CAZy_KOG.tab hmmsearch_KOFAM_multi_best.txt KEGG e-val ${TABLE_BASE}_TAX_CAZy_KOG_KEGG.tab

cd ../7_ANNOTATION/
cp TAX_TAB_FINAL.tab ../8_SPLIT/
cd ../8_SPLIT

# Download add_higher_taxonomy.py script
gdrive_download 1AWuqqPaP2rUMpF_uMHOGEs8aE3Iy7JS2 add_higher_taxonomy.py

python2.7 add_higher_taxonomy.py ${TABLE_BASE}_TAX_CAZy_KOG_KEGG.tab TAX_TAB_FINAL.tab TAX_BEST ${TABLE_BASE}_TAX2_CAZy_KOG_KEGG.tab TAX_tree_genus.tab

cd ..
mkdir 9_FINAL_TABLES
cp ./8_SPLIT/CAZy_tree.tab ./9_FINAL_TABLES
cp ./8_SPLIT/TABLE_NORM_SAMPLES_GENECALL_TAX_CAZy_KOG_KEGG.tab ./9_FINAL_TABLES
cp ./8_SPLIT/TAX_TAB_FINAL.tab ./9_FINAL_TABLES
cp ./8_SPLIT/TABLE_NORM_SAMPLES_GENECALL_TAX2_CAZy_KOG_KEGG.tab ./9_FINAL_TABLES

########################
## ANALYSIS COMPLETED ##
########################
