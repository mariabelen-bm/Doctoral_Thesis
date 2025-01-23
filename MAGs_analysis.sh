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
