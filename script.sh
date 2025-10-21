#! /bin/bash

####################### END-TO-END WORKFLOW FOR ANALYZING TARGETED EXOME SEQ DATA ##########################################################
################################################ Following GATK Best Practices Workflow #######################################


#--------------------------------------------------------------------------------------------------------------------------------------
######################################### SETTING UP DIRECTORIES AND DEFINE PATHS #######################################################
mkdir reads ref_files aligned_reads bed_files results #### make the required directories

########################### Defining variables for directory paths ########################################## 
ref=~/tar_es_analysis/ref_files/b37/human_g1k_v37.fasta 
known_sites_snps=~/tar_es_analysis/ref_files/b37/dbsnp_138.b37.vcf
known_sites_indels=~/tar_es_analysis/ref_files/b37/Mills_and_1000G_gold_standard.indels.b37.vcf
aligned_reads=~/tar_es_analysis/aligned_reads
reads=~/tar_es_analysis/reads
results=~/tar_es_analysis/results
bed_files=~/tar_es_analysis/bed_files
final_results=~/tar_es_analysis/final_results

#----------------------------------------------------------------------------------------------------------------------------------------

################################################## DOWNLOAD REFERENCE FILES ##############################################################

################################ Download the refrence genome ##################################


#-----------------------------------------------------------------------------------------------------------------------------------

############################################### HG19 Genome Files ####################################################################
################## Downloading hg19 reference genome (Latest Version) ######################
################## Alternatively, can use newer versions of the genome (hg38 or T2T consortium, available on UCSC) #################

###################### If you want to download the reference from UCSC genome browser ###############################
wget -P ${ref}/ https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/latest/hg19.fa.gz ### download the ref file

gunzip ${ref}/hg19.fa.gz ### unzip the fasta file


########################### Need to index and make dictionary file #####################################
samtools faidx ${ref}/hg19.fa
gatk CreateSequenceDictionary R=$${ref}/hg19.fa O=${ref}/hg19.dict


################## However, the UCSC HG19 genome does not have known sites in GATK Resource bundle #####################

#-----------------------------------------------------------------------------------------------------------------------------------


################ Hence downloading b37 (HG19) genome from GATK resource bundle ########################

#################################### b37 (HG19) Genome Files from GATK Resources ############################
################### Download the reference sequence file ########################

wget -P ${ref}/b37 ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz

gunzip ${ref}/b37/human_g1k_v37.fasta.gz #### unzip the file
 

##################### Index the reference file and create sequence dictionary files ######################

samtools faidx ${ref}/b37/human_g1k_v37.fasta ### index

gatk CreateSequenceDictionary R=${ref}/b37/human_g1k_v37.fasta O=${ref}/b37/human_g1k_v37.dict ### create seq dict


################## Alternatively can download the index and dictionary files from resource bundle #######################

wget -P ${ref}/b37 https://console.cloud.google.com/storage/browser/_details/broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai

wget -P ${ref}/b37 https://console.cloud.google.com/storage/browser/_details/broad-references/hg19/v0/Homo_sapiens_assembly19.dict


############################# Download known sites for BQSR from GATK Resource bundle ################################

wget -P ${ref}/b37 ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/dbsnp_138.b37.vcf.gz

wget -P ${ref}/b37 ftp://gsapubftp-anonymous@ftp.broadinstitute.org/bundle/b37/Mills_and_1000G_gold_standard.indels.b37.vcf.gz

gunzip ${ref}/b37/dbsnp_138.b37.vcf.gz

gunzip ${ref}/b37/Mills_and_1000G_gold_standard.indels.b37.vcf.gz


#-------------------------------------------------------------------------------------------------------------------------------------


########################################## HG38 Genome Files ##########################################################
######################### Download and unzip hg38 genome files ####################################

wget -P ${ref} https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz

gunzip ${ref}/hg38.fa.gz


########################### Index reference file  ###########################

samtools faidx ${ref}/hg38.fa

######################## Create sequence dictionary ###########################

gatk CreateSequenceDictionary R=${ref}/hg38.fa O=${ref}/hg38.dict


########################### Download known sites files for BQSR from GATK resource bundle ################################

wget -P ${ref} https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf

wget -P ${ref} https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf.idx



#--------------------------------------------------------------------------------------------------------------------------------------

################################################ VARIANT CALLING ###############################################

##################### Make sure your raw reads are in the reads directory ################################


#---------------------------------------------------------------------------------------------------------------------------------------

############################################# STEP-1: QUALITY CHECK ################################################################

#################### Checking quality of fastq files through FASTQC tool 

fastqc ${reads}/sample_R1.fastq.gz -o ${reads}/
fastqc ${reads}/sample_R2.fastq.gz -o ${reads}/


############## The overall sequence quality is good
############### No trimming required 
############## However, we will run fastp for an overall QC

#---------------------------------------------------------------------------------------------------------------------------------------



################################# STEP-2: QUALITY CONTROL #############################

#### Performing an overall all-in-one-QC using FASTP (for paired-end data)

fastp -i ${reads}/sample_R1.fastq.gz -I ${reads}/sample_R2.fastq.gz -o ${reads}/sample_R1_filt.fastq.gz -O ${reads}/sample_R2_filt.fastq.gz


################# After fastp, we can re-run fastqc to check if the pre-processing has taken effect. 
############### Not re-running FASTQC as the original reads were also of good quality. 

#----------------------------------------------------------------------------------------------------------------------------------------



########################################## STEP-3: ALIGNMENT ######################################


############## Aliging the reads with reference genome using BWA-MEM
############# Before alignment, index the reference genome

bwa index ${ref}

############### Before moving on to alignment, we would need read group information from fastq files
############### Extract the first four lines from fastq.gz files as follows

zcat ${reads}/sample_R1_filt.fastq.gz | head -n 4

########### The header shows read group info
########### Header for the above files is:
######## @MN01877:55:000H7CWWN:1:11102:23610:1050 1:N:0:55


######### If you want to find out if the reads are from different lanes

zcat ${reads}/sample_R1_filt.fastq.gz | grep "^@" | cut -d ":" -f 4 | sort | uniq -c

###### All reads are arising from same lane (1), so its not necessary to specify every field. Only sample ID, instrument name and sample name is enough
### @RG\tID:sample\tPL:ILLUMINA\tSM:sample


############################## Run BWA-MEM #################################

bwa mem -t 4 -R "@RG\tID:sample\tPL:ILLUMINA\tSM:sample" ${ref} ${reads}/sample_R1_filt.fastq.gz ${reads}/sample_R2_filt.fastq.gz > ${aligned_reads}/sample.paired.sam

############ To check the number of threads on system

nproc

################### Check first 100 lines of sam file ######################

samtools view ${aligned_reads}/sample.paired.sam | less

###################### Check the basic statistics of alignment ####################

samtools flagstat ${aligned_reads}/sample.paired.sam


#----------------------------------------------------------------------------------------------------------------------------------------



############################# STEP-4: MARK DUPLICATES AND SORT SAM FILES #######################################
###################### Flag duplicate reads, remove them and convert sam to bam ############################

gatk MarkDuplicatesSpark -I ${aligned_reads}/sample.paired.sam -O ${aligned_reads}/sample_sorted_dedup_reads.bam


#---------------------------------------------------------------------------------------------------------------------------------------



################################### STEP-5: PERFORM BASE QUALITY SCORE RECALIBRATION ##############################
#################### Since we have SNPs and Indels in different VCFs, we will pass both vcf files in the known-sites function ##########
################# Index the VCF files, otherwise gatk will generate an error ####################

gatk IndexFeatureFile -I ${known_sites_snps} ### index the snps

gatk IndexFeatureFile -I ${known_sites_indels} ### index the indels

##### BQSR 

gatk BaseRecalibrator -I ${aligned_reads}/sample_sorted_dedup_reads.bam -R ${ref} --known-sites ${known_sites_snps} --known-sites ${known_sites_indels} -O ${results}/recal_data.table

## Apply the model to recalibrate base quality scores

gatk ApplyBQSR -I ${aligned_reads}/sample_sorted_dedup_reads.bam -R ${ref} --bqsr-recal-file ${results}/recal_data.table -O ${aligned_reads}/sample_sorted_dedup_bqsr_reads.bam 


#-----------------------------------------------------------------------------------------------------------------------------------------



################################## STEP-6: VARIANT CALLING USING GATK HAPLOTYPE CALLER #####################################
##################### Run HaplotypeCaller with bed file #############################

gatk HaplotypeCaller -R ${ref} -I ${aligned_reads}/sample_sorted_dedup_bqsr_reads.bam -L ${bed_files}/capturekit_modified.bed -O ${results}/raw_variants_panel.vcf


################## After running Haplotype Caller, view VCF files ######################

bcftools view ${results}/raw_variants_panel.vcf

#----------------------------------------------------------------------------------------------------------------------------------------



############################################# STEP:7 EXTRACT SNPS AND INDELS ################################################
#################### Extract SNPs and Indels from VCF file generated with gene panel #########################

gatk SelectVariants -R ${ref} -V ${results}/raw_variants_panel.vcf --select-type SNP -O ${results}/raw_snps_panel.vcf

gatk SelectVariants -R ${ref} -V ${results}/raw_variants_panel.vcf --select-type INDEL -O ${results}/raw_indels_panel.vcf


#----------------------------------------------------------------------------------------------------------------------------------------



#################################### STEP-8: FILTERING THE VARIANTS ####################################
#################### Hard Filtering using gatk thresholds ###################################

################ Filter SNPs ###################

gatk VariantFiltration \
	-R ${ref} \
	-V ${results}/raw_snps_panel.vcf \
	-O ${results}/filtered_snps.vcf \
	-filter-name "QD" -filter "QD < 2.0" \
    -filter-name "QUAL" -filter "QUAL < 30.0" \
	-filter-name "FS" -filter "FS > 60.0" \
	-filter-name "MQ" -filter "MQ < 40.0" \
	-filter-name "SOR" -filter "SOR > 3.0" \
	-filter-name "MQRankSum" -filter "MQRankSum < -12.5" \
	-filter-name "ReadPosRankSum" -filter "ReadPosRankSum < -8.0" \
	-genotype-filter-expression "GQ < 20" \
	-genotype-filter-name "GQ"
    
    

########### Filtering Indels ########################

gatk VariantFiltration \
	-R ${ref} \
	-V ${results}/raw_indels_panel.vcf \
	-O ${results}/filtered_indels.vcf \
	-filter-name "QD" -filter "QD < 2.0" \
    -filter-name "QUAL" -filter "QUAL < 30.0" \
	-filter-name "FS" -filter "FS > 200.0" \
	-filter-name "ReadPosRankSum" -filter "ReadPosRankSum < -20.0" \
	-genotype-filter-expression "GQ < 20" \
	-genotype-filter-name "GQ"




###################### Filter the variants that pass the thresholds ###################################

####### SNPs ###############

gatk SelectVariants \
	--exclude-filtered \
	-V ${results}/filtered_snps.vcf \
	-O ${results}/analysis-ready-snps.vcf

############ Indels ###############

gatk SelectVariants \
	--exclude-filtered \
	-V ${results}/filtered_indels.vcf \
	-O ${results}/analysis-ready-indels.vcf

bcftools view ${results}/analysis-ready-indels.vcf ### view the VCFs

bcftools view ${results}/analysis-ready-snps.vcf


############## The analysis-ready VCFs contain variants that passed the quality filters ########################

#--------------------------------------------------------------------------------------------------------------------------------------



#################################### STEP-9: ANNOTATING THE VARIANTS #######################################
########## Using VariantAnnotater tool of gatk ######################
############# Adds dbsnp ID to the known sites ####################

######### Indels ###############

gatk VariantAnnotator \
   -R ${ref} \
   -V ${results}/analysis-ready-indels.vcf \
   -O ${results}/annotated-indels.vcf \
   -L ${bed_files}/capturekit_modified.bed \
   --dbsnp ${known_sites_indels}

################ SNPs #####################

gatk VariantAnnotator \
   -R ${ref} \
   -V ${results}/analysis-ready-snps.vcf \
   -O ${results}/annotated-snps.vcf \
   -L ${bed_files}/capturekit_modified.bed \
   --dbsnp ${known_sites_snps}

########### The VCF doesn't contain gene names ############################
######## Add gene names from BED Files through bcftools annotate ################################

############ Indels ####################

bcftools annotate -a ${bed_files}/capturekit_modified.bed -h <(echo '##INFO=<ID=GENE,Number=1,Type=String,Description="Genes">') -c CHROM,FROM,TO,GENE -O z -o ${results}/annotated_indels_final.vcf.gz ${results}/annotated-indels.vcf

gunzip ${results}/annotated_indels_final.vcf.gz


############ SNPs ################

bcftools annotate -a ${bed_files}/capturekit_modified.bed -h <(echo '##INFO=<ID=GENE,Number=1,Type=String,Description="Genes">') -c CHROM,FROM,TO,GENE -O z -o ${results}/annotated_snps_final.vcf.gz ${results}/annotated-snps.vcf

gunzip ${results}/annotated_snps_final.vcf.gz


#----------------------------------------------------------------------------------------------------------------------------------------


########################################## STEP-10: EXTRACT VARIANTS TO A TABLE ####################################################

######### Extract the known variants to a table for easy visualization ###############################

########## SNPs ##############

gatk VariantsToTable \
     -V ${results}/annotated_snps_final.vcf \
     -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F FILTER -F GENE \
     -O ${results}/snps.table

########### Indels #############

gatk VariantsToTable \
     -V ${results}/annotated_indels_final.vcf \
     -F CHROM -F POS -F ID -F REF -F ALT -F QUAL -F FILTER -F GENE \
     -O ${results}/indels.table


########### Merge snps and indels to a single table ####################

(cat ${results}/indels.table; tail -n +2 ${results}/snps.table) > ${results}/final_variants.table

#---------------------------------------------------------------------------------------------------------------------------------------


########################################################### END OF SCRIPT ##############################################################