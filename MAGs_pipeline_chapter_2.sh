###################
## MAGs ANALYSIS ##
###################

# The analysis of MAGs (Metagenome-Assembled Genomes) is an approach used in metagenomics to reconstruct complete genomes of microorganisms present in an environmental 
# sample, without the need for prior isolation in culture. Using raw metagenomic sequences, MAGs are obtained by an assembly and binning process, 
# in which contigs (DNA fragments) are grouped into bins representing individual genomes. 
# These genomes can come from bacteria, archaea or other microbes present in the sample. 

#######
## 1 ##
#######

#############
## BINNING ##
#############

# Binning is a key step in metagenomic analysis. The main objective of binning is to group contigs (assembled DNA fragments) into bins, 
# where each bin represents a possible individual genome. 

conda activate metawrap

metawrap binning -a /mnt/DATA/belen/MAGS_assembly/final.contigs.fa -o binning_metawrap -t 120 -m 1000 --metabat2 --maxbin2 --concoct --universal --run-checkm --interleaved /mnt/DATA/belen/MAGS_assembly/*.pe.qc.fq.gz

#######
## 2 ##
#######

####################
## BIN REFINEMENT ##
####################

# In this step you pick the best version of each bin. You can be more or less stringent in this step by lowering the completeness a bit.

metawrap bin_refinement -o /mnt/DATA/belen/MAGS_assembly/metawrap_refined_bins/ -A /mnt/DATA/belen/MAGS_assembly/binning_metawrap/metabat2_bins/ -B /mnt/DATA/belen/MAGS_assembly/binning_metawrap_concot/concoct_bins/ -C /mnt/DATA/belen/MAGS_assembly/binning_metawrap/maxbin2_bins/ -m 1000 -t 120 -c 50 -x 10

#######
## 3 ##
#######

##################
## CheckM2 STEP ##
##################

# CheckM2 is used to evaluate the quality of the refined bins (MAGs) obtained. CheckM2 is a tool that estimates the completeness and contamination of MAGs.

conda activate checkm2

checkm2 predict -i /mnt/DATA/belen/MAGS_assembly/metawrap_refined_bins/metawrap_50_10_bins/ -x fa --output-directory refinded_checkm2 --database_path /mnt/DATA1/priscila/checkm2/database/CheckM2_database/uniref100.KO.1.dmnd --tmpdir ./ --threads 240

awk '$2 >= 50 && $3 <=10' refinded_checkm2/quality_report.tsv > good_bins_checkm2.tsv

awk '$2 >= 50 && $3 <=10' refinded_checkm2/quality_report.tsv | cut -f1 > bins_list

for i in $(cat bins_list); do cp metawrap_50_10_bins/$i.fa selected_bins/ ;done 

#######
## 4 ##
#######

##########################
## TAXONOMIC ANNOTATION ##
##########################

# GTDB-Tk (Genome Taxonomy Database Toolkit) is used to assign taxonomy to refined MAGs using the GTDB database version 2.4.0 (v220). 
# This tool classifies microbial genomes from complete genomic data and provides a standardized taxonomy based on the phylogenetic tree proposed by GTDB.

conda activate /mnt/DATA1/priscila/condaenvs/gtdbtk220

gtdbtk classify_wf --genome_dir  /mnt/DATA/belen/MAGS_assembly/metawrap_refined_bins/metawrap_50_10_bins/ --out_dir refined_checkm2_gtdb220 --cpus 240 --pplacer_cpus 60 -x .fa --tmpdir ./ --skip_ani_screen

#######
## 5 ##
#######

##########
## GUNC ##
##########

# GUNC (Genomic UNcertainty Calculator), a tool designed to assess the taxonomic contamination and consistency of MAGs, is used. 
# This analysis is crucial to verify the quality of the refined MAGs and ensure that they represent unique and consistent genomes rather than 
# mixtures of genetic material from different organisms.

conda activate gunc

export TMPDIR="/mnt/DATA/projects/priscila/tmp/"
echo $TMPDIR

mkdir selected_bins_gunc

gunc run --input_dir /mnt/DATA/belen/MAGS_assembly/metawrap_refined_bins/metawrap_50_10_bins/ --detailed_output --contig_taxonomy_output --use_species_level --out_dir selected_bins_gunc --threads 120 --db_file /mnt/DATA1/priscila/database/gunc_db_progenomes2.1.dmnd --file_suffix .fa

#######
## 6 ##
#######

#####################################
## RELATIVE QUANTIFICATION OF MAGs ##
#####################################

# The relative quantification of MAGs is performed using the Minimap2 and CoverM tools. 
# In this step, the relative abundance of each MAG in the microbial community is determined based on the mapping of MAGs.

mkdir mags_bams

conda activate coverm

export TMPDIR="/mnt/DATA/projects/priscila/tmp/"
export TMPDIR="/mnt/DATA/belen/MAGS_assembly"

echo $TMPDIR

coverm genome --mapper minimap2-sr --methods relative_abundance -o coverm_relative_abundance_selected.txt --bam-file-cache-directory mags_bams --interleaved /mnt/DATA/belen/MAGS_assembly/*.pe.qc.fq.gz --genome-fasta-directory /mnt/DATA/belen/MAGS_assembly/metawrap_refined_bins/metawrap_50_10_bins/ -x fa --threads 240

#######
## 7 ##
#######

###########################
## FUNCTIONAL ANNOTATION ##
###########################

conda activate DRAM

# We have
lrwxrwxrwx. 1 belen belen   32 Dec 11 17:39 gtdbtk.ar53.summary.tsv -> classify/gtdbtk.ar53.summary.tsv
lrwxrwxrwx. 1 belen belen   34 Dec 11 18:29 gtdbtk.bac120.summary.tsv -> classify/gtdbtk.bac120.summary.tsv

# Combine both files
head -n 1 /mnt/DATA/belen/MAGS_assembly/metawrap_refined_bins/refined_checkm2_gtdb220/gtdbtk.bac120.summary.tsv > gtdbtk_summary_bac_arch.tsv
tail -n +2 /mnt/DATA/belen/MAGS_assembly/metawrap_refined_bins/refined_checkm2_gtdb220/gtdbtk.bac120.summary.tsv >> gtdbtk_summary_bac_arch.tsv
tail -n +2 /mnt/DATA/belen/MAGS_assembly/metawrap_refined_bins/refined_checkm2_gtdb220/gtdbtk.ar53.summary.tsv >> gtdbtk_summary_bac_arch.tsv

# Functional annotation
DRAM.py annotate -i '/mnt/DATA/belen/MAGS_assembly/metawrap_refined_bins/metawrap_50_10_bins/*.fa' \
-o dram_annotation/ \
--min_contig_size 2000 \
--gtdb_taxonomy /mnt/DATA/belen/MAGS_assembly/metawrap_refined_bins/refined_checkm2_gtdb220/gtdbtk_summary_bac_arch.tsv \
--checkm_quality /mnt/DATA/belen/MAGS_assembly/refinded_checkm2/quality_report.tsv \
--threads 512 \
--verbose \
--kofam_use_dbcan2_thresholds \
--keep_tmp_dir

cd /mnt/DATA/belen/MAGS_assembly/dram_annotation

DRAM.py distill -i /mnt/DATA/belen/MAGS_assembly/dram_annotation/annotations.tsv -o genome_summaries --trna_path /mnt/DATA/belen/MAGS_assembly/dram_annotation/trnas.tsv --rrna_path /mnt/DATA/belen/MAGS_assembly/dram_annotation/rrnas.tsv

#######
## 8 ## 
#######

################################
## RELATIVE ABUBDANCE OF MAGs ##
################################

# In this step, quantification of the relative abundance of MAGs in each sample is performed using the Salmon tool. 
# The aim is to determine how many reads from each sample map to the different MAGs, which provides information on the relative abundance of the assembled genomes 
# in the different samples.

conda activate metawrap

metawrap quant_bins2 -a /mnt/DATA/belen/MAGS_assembly/final.contigs.fa -o quantified_bins -t 240 -b /mnt/DATA/belen/MAGS_assembly/metawrap_refined_bins/metawrap_50_10_bins/ /mnt/DATA/belen/MAGS_assembly/binning_metawrap/work_files/*.bam

#############################
## MAGs ANALYSIS COMPLETED ##
#############################
