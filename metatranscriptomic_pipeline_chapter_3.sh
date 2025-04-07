#################################
## METATRANSCRIPTOMIC ANALYSIS ##
#################################

conda activate fastp


for file in *_R1.fq.gz
do
   sample=${file%%_R1.fq.gz}
   fastp --detect_adapter_for_pe --adapter_sequence=AGATCGGAAGAG --adapter_sequence_r2=AGATCGGAAGAG -W 1 -M 3 -5 -3 -g -q 30 -u 50 -l 50 -h ${sample}.html --thread=16 --dont_eval_duplication  -i ${sample}_R1.fq.gz -I ${sample}_R2.fq.gz --unpaired1=filtered/${sample}.se.fq --unpaired2=filtered/${sample}.se.fq  --stdout > filtered/${sample}.pe.trim.qc.fq
done

# cat se reads into same file as pe

for file in *.pe.trim.qc.fq
do
   sample=${file%%.pe.trim.qc.fq}
cat ${sample}.se.fq >> ${file}
done


##REMOVE rRNAs WITH - bbduk.sh

bbdir='/home/kdanielmorais/bioinformatics/tools/BBtools/'
echo ${bbdir}                                            
for file in *.pe.trim.qc.fq
do
    sample=${file%%.pe.trim.qc.fq}	
    echo "${bbdir}bbmap/bbduk.sh ordered k=31 ref=${bbdir}ribokmers.fa.gz ow=true in=${file} out=rRNA_remove/${sample}_NO_rRNA.pe.fq outm=rRNA_remove/${sample}_rRNA.pe.fq"
done > remove_rrna.sh

cat remove_rrna.sh | parallel

# fix paired R1 and R2 files 

for file in *_NO_rRNA.pe.fq
do
   sample=${file%%_NO_rRNA.pe.fq}
   echo "/home/kdanielmorais/bioinformatics/tools/BBtools/bbmap/repair.sh in=${file} out1=for_assembly/${sample}.pe.fq.1 out2=for_assembly/${sample}.pe.fq.2 outsingle=for_assembly/${sample}.se.fq repair"
done > extract_command.sh

cat extract_command.sh | parallel

#rm -rf *.tr.qc.fq.cut

# trinity assembly 
## obs about rnaSeq libs: the kit used normaly here are stranded rna preps, this generates strand-specific transcripts and should be assembled a bit different. Trinity has an option for this --SS_lib_type (can be RF or FR) the kit "TruSeq" uses dUTP method (FR according to Trinity documents). Should try to compare this data assebled normally and considering strand-specifi option as well????

## put all reads together

cat *.1 > all.qc.fq.1
cat *.2 > all.qc.fq.2
cat *.se.fq > all.se.qc.fq
#put unpaired reads into file .1
cat all.se.qc.fq >> all.qc.fq.1

### TRINITY ASSEMBLY  
mkdir trinity_assembly
#trinity can't use more than 200G of ram
conda activate TrinityEnv
Trinity --NO_SEQTK  --seqType fq  --left all.qc.fq.1 --right all.qc.fq.2 --CPU 120 --max_memory 200G --output trinity_assembly/


## SAMPLE MAPPING 
#find a way for this loop to work across different folders!!!!!

REF=../trinity_assembly.Trinity.fasta
reference=${REF%%.fasta}
echo "reference is" ${reference}
mkdir ${reference}_build
bowtie2-build ${REF} ${reference}_build/${reference}.build


ref=trinity_assembly.Trinity_build
reference=${ref%%_build}

for file in *_NO_rRNA.pe.fq
 do
 sample=${file%%_NO_rRNA.pe.fq}
 echo "processing ${sample}... reference ${reference}"

 bowtie2 -p 70 -x for_assembly/trinity_assembly.Trinity_build/trinity_assembly.Trinity.build -q ${file} -S ${sample}.sam
 echo "sam file is done..."
 
 #rm -rf ${sample}.pe.fq 

 samtools view -Sb ${sample}.sam > ${sample}.bam
 echo "bam file is done..."

 rm -rf ${sample}.sam

 samtools view -c -f 4 ${sample}.bam > ${sample}.reads-unmapped.count.txt
 echo "unmapped reads info done..."

 samtools view -c -F 4 ${sample}.bam > ${sample}.reads-mapped.count.txt
 echo "mapped reads info done..."

 samtools sort -o ${sample}.sorted.bam ${sample}.bam
 echo "bam file was sorted..."

 rm -rf ${sample}.bam

 samtools index ${sample}.sorted.bam
 echo "soerted bam file was indexed..."
 samtools idxstats ${sample}.sorted.bam > ${sample}.reads.by.contigs.txt

 echo "sample ${sample} is done..."
done

######GETTING COUNT MAPPING TABLE CONTIG/SAMPLE
gdrive_download 1HDB2EF-pq-EJxQxsI1uv1-iVjTl6tAlo count-up-mapped-from-results-txt-with-ctg-length.py

python2.7 count-up-mapped-from-results-txt-with-ctg-length.py *.reads.by.contigs.txt

#####################
# get coverage file #
#####################

https://drive.google.com/file/d/1S2AQHd2YIjnxZz2kIa2avo1RSAuj-pWT/view?usp=sharing

gdrive_download 1S2AQHd2YIjnxZz2kIa2avo1RSAuj-pWT get_assembly_coverage.py

# adjusted to 145bp because used the trimmed and filtered reads _NO_rRNA files

python2.7 get_assembly_coverage.py summary-count-mapped.tsv 145 Ruben_substrateMT_assembly_coverage.txt


###GENE CALLING - FragGeneScan

#####################
# GENECALLING - FGS #
#####################

# Link creation
ln -s /home/kdanielmorais/bioinformatics/tools/fraggenescan/FragGeneScan1.31/train/ ./

# Command
FragGeneScan -s  ../for_assembly/trinity_assembly.Trinity.fasta -w 1 -o Ruben_substrateMT_trinity_genecalling -t complete -p 120


###FALTA FAZER A PARTE FINAL DO MAPPING E NORMALIZACAO
######################################
# 2a NORMALISE MAPPING TABLE PER BASE
######################################
gdrive_download 1w0bfttjXFZ64NHD8bP7UDaQcS1yd20qR normalize-mapping-table-by-read-length-and-ctg-length.py
 
python2.7 normalize-mapping-table-by-read-length-and-ctg-length.py summary-count-mapped.tsv 145 TABLE_normalised_145.txt

WARNING: gene length is 0 bp - *
done...


########################################
# 2b NORMALISE MAPPING TABLE PER SAMPLE
########################################
gdrive_download 1c_fD520xtrCNlUIq9VqqqSvY2OryMXTU normalize_table_by_columns.py

python2.7 normalize_table_by_columns.py TABLE_normalised_145.txt 2 1000000 TABLE_normalised_per_sample.txt



##########################
# 3 MULTIPLY BY GENECALL # -->         
##########################

#k141_43039
gdrive_download 1DakO7roc9C2GJ-SkZuy8AZKTV3QHan14 contig_mapping_to_genecall_mapping.py

python2.7 contig_mapping_to_genecall_mapping.py genecalling/Ruben_substrateMT_trinity_genecalling.faa TABLE_normalised_per_sample.txt

ls
############################
# ADD "#" TO SAMPLE NAMES  #
############################  

head -1 TABLE_normalised_per_sample.txt_genecall.txt | awk -F'\t' '{printf $1"\t"$2 ;for(i=3; i<=NF; ++i) printf "\t%s", "#"$i }' |  awk -F '\t' '{print $0}' > header.txt

tail -n +2  TABLE_normalised_per_sample.txt_genecall.txt > table.txt

cat header.txt table.txt > TABLE_NORM_SAMPLES_GENECALL.txt


#####################
## DRAM ANNOTATION ##
#####################

#error in description can be fixed with 
DRAM-setup.py update_description_db
# runs for a few hours and takes a few hundred Gb disk space
conda activate DRAM
DRAM.py annotate_genes -i genecalling/Ruben_substrateMT_trinity_genecalling.faa -o annotation_DRAM --threads 240 --verbose --use_uniref 

###obs for this step
https://github.com/WrightonLabCSU/DRAM/issues/62
dbcan-CAZy uses filtering indicated by them - dbCAN2 suggestions for thresholds" suggested by dbcan2:
(see http://bcb.unl.edu/dbCAN2/blast.php)
E-Value < 1e-15, coverage > 0.35



#summarize results
DRAM.py distill -i annotation_DRAM2/annotations.tsv -o annotation_DRAM2/distilled


# next time run with FOAM-hmm_rel1a.hmm database 



#######################----TAXONOMY NCBI-JGI----################################


###################################
# JGI FUNGAL PROTEINS - BIOCEV PC #
###################################


/home/kdanielmorais/bioinformatics/tools/diamond blastp -d /mnt/DATA/DATABASES/FUNGAL_PROTEINS_JGI/JGI_FUNGAL_PROTEINS_ANNOTATED_20210312 -q genecalling/Ruben_substrateMT_trinity_genecalling.faa -e 1E-5 -o taxonomy/genecalling_JGI_FUN_20210312.txt -f 6 -p 256 -b12 -c1

export LANG=en_US.UTF-8
export LC_ALL=en_US.UTF-8
sort -t$'\t' -k1,1 -k12,12gr -k11,11g -k3,3gr genecalling_JGI_FUN_20210312.txt | sort -u -k1,1 --merge > genecalling_JGI_FUN_20210312_best.txt


# GENERA DEFINED   ##########RUNNING THIS -- started on 20.07.2022 at
/home/kdanielmorais/bioinformatics/tools/diamond blastp -d /mnt/DATA/DATABASES/NCBI_nr_DIAMOND/NCBI_nr_20210225_diamond_GENERA -q CLEMENTINE_MT_Trinity_genecalling_fgs.faa -e 1E-5 -o genecalling_NCBI_nr_PROTEINS_GENERA.txt -f 6 -p 512 -b12 -c1
##
Total time = 7076.34s
Reported 27886590 pairwise alignments, 27886590 HSPs.
1229412 queries aligned.
##

export LANG=en_US.UTF-8
export LC_ALL=en_US.UTF-8
sort -t$'\t' -k1,1 -k12,12gr -k11,11g -k3,3gr genecalling_NCBI_nr_PROTEINS_GENERA.txt | sort -u -k1,1 --merge > genecalling_NCBI_nr_GENERA_PROTEINS_best.txt



#################################
# ADD TAXONOMY TO BLAST RESULTS #
#################################

# JGI
# latest 20210406
gdrive_download 12c28kgIw4mPBIhQutNGladdAXwNLtvlR jgi_abr_org_list.txt

gdrive_download 1XBTtiC1JYl2rzeV7idN2WrveEZknmnQi replace_fungal_annot_by_taxname.py

python2.7 replace_fungal_annot_by_taxname.py genecalling_JGI_FUN_20210312_best.txt jgi_abr_org_list.txt genecalling_JGI_FUNGAL_PROTEINS_best_reformate.txt

awk -F'\t' '{print $2}' genecalling_JGI_FUNGAL_PROTEINS_best_reformate.txt | sort | uniq > FUNGAL_NAMES.txt




# NCBI
awk -F'\t' '{print $2}' genecalling_NCBI_nr_GENERA_PROTEINS_best.txt | sort | uniq > ALL_ACCESSIONS.txt

#gdrive_download 1FQdQ2Oh3sgyc2IKymYz_B6mKQwlklia_ retrieve_taxonomy_by_accession_with_taxid_library.py
#python2.7 retrieve_taxonomy_by_accession_with_taxid_library.py ALL_ACCESSIONS.txt /mnt/DATA/DATABASES/ACC2TAXID/ACC2TAXID_nr_current.txt taxa_all_accessions.txt

https://drive.google.com/file/d/1o8KmSbwzOsjjeouK3dR0RNmWkMjdfFow/view?usp=sharing

gdrive_download 1o8KmSbwzOsjjeouK3dR0RNmWkMjdfFow get_taxonomy_offline.py

python2.7 get_taxonomy_offline.py ALL_ACCESSIONS.txt /mnt/DATA/DATABASES/NCBI_nr_DIAMOND/ACC2TAXID_nr_current.txt /mnt/DATA/DATABASES/NCBI_nr_DIAMOND/TAXONOMY_TAXID_ALL_fixed.txt taxa_all_accessions.txt 

accession list loaded...(779978)
taxonomy loaded...(907158)
taxonomy retrieved... 779978 vs acc (779978) - should be equal!
Done :]


# REFORMAT
https://drive.google.com/file/d/1jQ3F3ZuA0sBJVy3eJiRhwAaxh31LKqSF/view?usp=sharing

gdrive_download 1jQ3F3ZuA0sBJVy3eJiRhwAaxh31LKqSF replace_acc_by_sp_from_taxonomy.py

python2.7 replace_acc_by_sp_from_taxonomy.py genecalling_NCBI_nr_GENERA_PROTEINS_best.txt TAXONOMY_ALL.txt genecalling_NCBI_GENERA_PROTEINS_best_reformat.txt
#
number of taxa: 779978 (779978)
DONE :) Processed blast: 1229412 - NOT FOUND 0
#

# COMBINE TAXONOMY TABLES
# GET TAXONOMY FOR ALL
gdrive_download 1VtSyy7OutKZ6fAZMTYY2HkZUDNcpfY6V JGI_TAXA_TAB_2021.txt
#https://drive.google.com/file/d/1F5p28LpaHrSYWwNI_82V9eHoKIb63y8G/view?usp=sharing

gdrive_download 1F5p28LpaHrSYWwNI_82V9eHoKIb63y8G combine_taxonomy_tables.py

python2.7 combine_taxonomy_tables.py FUNGAL_NAMES.txt JGI_TAXA_TAB_2021.txt TAXONOMY_ALL.txt TAX_TAB.tab 
##
names loaded...
FUNGAL TAXONOMY PROCESSED - NOT FOUND 0 vs. FOUND 1498
OTHER TAXONOMY PROCESSED - REDUCING TO 23625 vs. ORIGINAL 779979
DONE :)
###

####
gdrive_download 1-3XE5Le8I1_HzQdWbaAHev4ZlrSs6lUi get_best_hit_by_bitscore_multi.py

python2.7 get_best_hit_by_bitscore_multi.py genecalling_NCBI_GENERA_PROTEINS_best_reformat.txt genecalling_JGI_FUNGAL_PROTEINS_best_reformate.txt
#
FILE: genecalling_NCBI_GENERA_PROTEINS_best_reformat.txt - HITS: 1229412
NEW ANNOTATIONS: 1229412 - REPLACED: 0 - CURRENT BEST HITS: 1229412

FILE: genecalling_JGI_FUNGAL_PROTEINS_best_reformate.txt - HITS: 714548
NEW ANNOTATIONS: 12375 - REPLACED: 308679 - CURRENT BEST HITS: 1241787
#

awk -F'\t' '{print $2}' best_of_the_blast.txt | sort | uniq > ALL_TAXA_NAMES.txt

gdrive_download 1XruvN2qGN2-dUHZNSn0jXOYmUxoJ3Uz0 get_taxonomy_basedonnames.py

python2.7 get_taxonomy_basedonnames.py ALL_TAXA_NAMES.txt TAX_TAB.tab TAX_TAB_FINAL.tab 
#names loaded...
TAXONOMY PROCESSED - NOT FOUND 0 vs. FOUND 24671
DONE :)
#
awk -F'\t' '{print $1"\t"$12"\t"$2""}' best_of_the_blast.txt > TAXONOMY_BEST_OF_SIMPLE.txt

##############################
# KOGG FROM JGI-MYCO-GENOMES #
##############################

# name e-val KOG
awk -F'[|\t]' '{print $1"\t"$15"\t["$5"]"}' genecalling_JGI_FUN_20210312_best.txt > FUNCTION_JGI_KOG_SIMPLE.txt

##########---- SPLITING FOR FOAM -----#########	


############
# split it #
############

gdrive_download 1mGbdx3OBumymosW24WaYfZT9nq_a1z1z split_fasta_by_group_size.py

python2.7 split_fasta_by_group_size.py /mnt/DATA1/priscila/rubenMT/filtered/rRNA_remove/genecalling/Ruben_substrateMT_trinity_genecalling.faa 83000

mkdir SPLIT 
mv *.fas SPLIT

########
# FOAM #
########
cd SPLIT
for file in *.fas
do
   output=${file%%.fas}
   echo "/home/kdanielmorais/bioinformatics/tools/hmmer-3.0/src/hmmsearch --tblout ${output}.txt --noali --cpu 1 -E 1e-5 /mnt/DATA/DATABASES/FOAM_db/FOAM-hmm_rel1.hmm ${file} >/dev/null 2>&1"
done > foam.sh

mkdir tmp
cat foam.sh | parallel --tmpdir tmp

##### PROCESS OUTPUT ######

for file in *.txt
do
 sample=${file%%.txt}
 grep -v '#' ${file} | awk -F' ' '{print $1"\t"$3"\t"$5"\t"$6}' > ${sample}.for_sort.txt
done


export LC_ALL=en_US.UTF-8
export LANG=en_US.UTF-8

for file in *.for_sort.txt
do
 sample=${file%%.for_sort.txt}
 echo "sort -t$'\t' -k1,1 -k4,4gr -k3,3g ${file} | sort -u -k1,1 --merge > ${sample}.sorted_best.txt"
done > sort.sh

cat sort.sh | parallel

cat *.sorted_best.txt > FOAM_BEST.txt

##########################################################
# FOAM ANNOTATION
##########################################################

gdrive_download 1ckuIqWVVarFgtcEUQOne2z-TdkDO3aXU FOAM_simple_multi_from_raw.py

python2.7 FOAM_simple_multi_from_raw.py FOAM_BEST.txt FUNCTION_FOAM_KO_SIMPLE_MULTI.txt


#####################################
#  USE dbCAN LOCAL DATABASE - CONDA # 
#####################################

#First activate the dbcan environment with
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

#####################################################
# dbCAN ANNOTATION
#####################################################

export LC_ALL=en_US.UTF-8
export LANG=en_US.UTF-8

sort -t$'\t' -k3,3 -k5,5g all_dbCAN.txt | sort -u -k3,3 --merge > all_dbCAN_best.txt

awk -F'[.\t]' '{print $1}' all_dbCAN_best.txt | sort | uniq > hmm_names_uniq.txt

awk -F'[.\t]' '{print $1}' all_dbCAN_best.txt > hmm_names.txt

awk -F'\t' '{print $3"\t"$5}' all_dbCAN_best.txt > all_dbCAN_best_gene_eval.txt

paste -d"\t" all_dbCAN_best_gene_eval.txt hmm_names.txt > CAZy_BEST_SIMPLE.txt


###################################
# LINK ANNOTATION TO TABLE
###################################

gdrive_download 198TDGsV1cBfLEZorb5znFHysG47XEj5t link_simple_table_to_mapping_table.py

#python2.7 link_simple_table_to_mapping_table.py mapping_table_normalised_per_sample_genecall.txt best_of_the_blast_simple.txt BESTTAX bitscore MAPTAB_NORMPERSAMPLE_GENES_BESTTAX.ta
#mv MG_CLEMENTINE_normalised_per_sample.txt_genecall.txt TABLE_NORM_SAMPLES_GENECALL.txt

TABLE="TABLE_NORM_SAMPLES_GENECALL.txt"

echo "${TABLE}"

TABLE_BASE=${TABLE%%.${TABLE##*.}}

echo "${TABLE_BASE}"

python2.7 link_simple_table_to_mapping_table.py ../${TABLE} TAXONOMY_BEST_OF_SIMPLE.txt TAX_BEST bitscore ${TABLE_BASE}_TAX.txt
Functions processed... 1241787
Done... used functions: 1241787/1241787

python2.7 link_simple_table_to_mapping_table.py ${TABLE_BASE}_TAX.txt CAZy_BEST_SIMPLE.txt CAZy e-val ${TABLE_BASE}_TAX_CAZy.txt
Functions processed... 10697
Done... used functions: 10696/10697

python2.7 link_simple_table_to_mapping_table.py ${TABLE_BASE}_TAX_CAZy.txt FUNCTION_FOAM_KO_SIMPLE_MULTI.txt FOAM e-val ${TABLE_BASE}_TAX_CAZy_FOAM.tab
Functions processed... 295860
Done... used functions: 295860/295860

python2.7 link_simple_table_to_mapping_table.py ${TABLE_BASE}_TAX_CAZy_FOAM.tab FUNCTION_JGI_KOG_SIMPLE.txt KOG e-val ${TABLE_BASE}_TAX_CAZy_FOAM_KOG.tab
Functions processed... 714548
Done... used functions: 714548/714548
# genus taxonomy

https://drive.google.com/file/d/1AWuqqPaP2rUMpF_uMHOGEs8aE3Iy7JS2/view?usp=sharing

gdrive_download 1AWuqqPaP2rUMpF_uMHOGEs8aE3Iy7JS2 add_higher_taxonomy.py

#python2.7 add_higher_taxonomy.py TABLE_NORM_SAMPLES_GENECALL_TAX_CAZy_FOAM_KOG.tab TAX_tree.tab TAX_BEST TABLE_NORM_SAMPLES_GENECALL_TAX2_CAZy_FOAM_KOG.tab TAX_tree_genus.tab

python2.7 add_higher_taxonomy.py TABLE_NORM_SAMPLES_GENECALL_TAX_CAZy_FOAM_KOG.tab TAX_TAB_FINAL.tab TAX_BEST TABLE_NORM_SAMPLES_GENECALL_TAX2_CAZy_FOAM_KOG.tab TAX_tree_genus.tab

>>>duplicate<<<
new: Eukaryota  Nematoda        Chromadorea     Rhabditida      Setariidae      [Setaria]
old: Eukaryota  Streptophyta    Magnoliopsida   Poales  Poaceae [Setaria]
>>>duplicate<<<
new: Eukaryota  Mucoromycota    Mucoromycetes   Mucorales       Lichtheimiaceae [Fennellomyces]
old: Eukaryota  Mucoromycota    Mucoromycetes   Mucorales       Syncephalastraceae      [Fennellomyces]
>>>duplicate<<<
new: Eukaryota  Basidiomycota   Agaricomycetes  Polyporales     Meruliaceae     [Phlebia]
old: Eukaryota  Basidiomycota   Agaricomycetes  Polyporales     Steccherinaceae [Phlebia]
>>>duplicate<<<
new: Eukaryota  Basidiomycota   Agaricomycetes  Agaricales      Tricholomataceae        [Infundibulicybe]
old: Eukaryota  Basidiomycota   Agaricomycetes  Agaricales      undefined Agaricales    [Infundibulicybe]
>>>duplicate<<<
new: Eukaryota  Ascomycota      Saccharomycetes Saccharomycetales       Trichomonascaceae       [Blastobotrys]
old: Eukaryota  Ascomycota      Saccharomycetes Saccharomycetales       Trigonopsidaceae        [Blastobotrys]
>>>duplicate<<<
new: Eukaryota  Ascomycota      Saccharomycetes Saccharomycetales       Debaryomycetaceae       [Candida]
old: Eukaryota  Ascomycota      Saccharomycetes Saccharomycetales       undefined Saccharomycetales     [Candida]
>>>duplicate<<<
new: Eukaryota  Ascomycota      Eurotiomycetes  Eurotiales      Thermoascaceae  [Paecilomyces]
old: Eukaryota  Ascomycota      Sordariomycetes Hypocreales     Clavicipitaceae [Paecilomyces]
>>>duplicate<<<
new: Eukaryota  Ascomycota      Dothideomycetes Pleosporales    Cucurbitariaceae        [Pyrenochaeta]
old: Eukaryota  Ascomycota      Dothideomycetes Pleosporales    Neopyrenochaetaceae     [Pyrenochaeta]
>>>duplicate<<<
new: Bacteria   Calditrichaeota Calditrichae    Calditrichales  Calditrichaceae [Caldithrix]
old: Bacteria   Calditrichaeota Calditrichia    Calditrichales  Calditrichaceae [Caldithrix]
done :]


######################
#unique dbCAN models #

awk -F'\t' '{print $3}' CAZy_BEST_SIMPLE.txt | sort | uniq > CAZy_BEST_unique.txt

https://drive.google.com/file/d/1SGVK2cqWCLozEPGNLvPRs0YG-ckrF_CZ/view?usp=sharing

gdrive_download 1SGVK2cqWCLozEPGNLvPRs0YG-ckrF_CZ get_CAZy_tree.py

python2.7 get_CAZy_tree.py CAZy_BEST_unique.txt CAZy_tree.tab

#############
# KOG TREE #
#############
gdrive_download 1uoUDoD5El-Gdlv9VNQ6BGATtacGsq2MF KOG_TAB_2021_03_03.txt


##################
### KofamKOALA ###  use 0.00001 for e-val
##################
in "/mnt/DATA1/priscila/rubenMT/filtered/rRNA_remove/kofam_KOs/kofam_koala_KOs_rubenMT.txt" 
/mnt/DATA1/priscila/kofamKOALA/bin/kofam_scan-1.3.0/exec_annotation ../genecalling/Ruben_substrateMT_trinity_genecalling.faa -o kofam_koala_KOs_rubenMT.txt -p /mnt/DATA1/priscila/kofamKOALA/db/profiles/ -k /mnt/DATA1/priscila/kofamKOALA/db/ko_list --cpu=120 --tmp-dir=tmp -E 1e-5 -f detail-tsv

###extract KOs and evals 
python2.7 KO_simple_tab_from_koala.py kofam_koala_KOs_rubenMT.txt KO_simple_table.txt

#add this to the bigTable

python2.7 link_simple_table_to_mapping_table.py TABLE_NORM_SAMPLES_GENECALL_TAX2_CAZy_FOAM_KOG.tab kofam_KOs/KO_simple_table.txt  KEGG e-val TABLE_NORM_SAMPLES_GENECALL_TAX2_CAZy_FOAM_KOG_KEGG.tab

#############
# KEGG TREE #
#############
gdrive_download 11CVgwqy6O2mJ5rc4vxevYQVp04oJrQsl kegg_tab.txt
gdrive_download 1AmOiMqLHE8nberbvYEDPT12JShAfSW_w GET_KEGG_ontology_subtable.py

gdrive_download 1aqENUDPh3dCoxnO8YfBZaUjwwxkmlMM1 KO_UNIQUE_from_KO_simple.py
cat kofam_KOs/KO_simple_table.txt FUNCTION_FOAM_KO_SIMPLE_MULTI.txt > KEGG_ALL_MULTI.txt

python2.7 KO_UNIQUE_from_KO_simple.py KEGG_ALL_MULTI.txt
python2.7 GET_KEGG_ontology_subtable.py kegg_tab.txt KEGG_ALL_MULTI.txt.unique.txt FOAM_KEGG_tree.tab
