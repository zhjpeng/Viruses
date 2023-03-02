#!/bin/bash

###################################################################
############################# 开工大吉#############################
###################################################################


###################################################################
######################### meta-genomics ###########################
###################################################################

MAIN=$(pwd)
BACU=${MAIN}/00-BaculoviridaeGenome
DATA=${MAIN}/01-reference
RAW=${MAIN}/02-rawFastQ
QUAL=${MAIN}/03-qualityAssesment
TRIM=${MAIN}/04-cleanFastQ
HOST=${MAIN}/05-removeHost
CONTIG=${MAIN}/06-primaryAssembly
VIRUSE=${MAIN}/07-DatabaseSearch
NAME=${MAIN}/06-taxonClassification

FLAG="[MetaTrans "$(date "+%Y-%m-%d %H:%M:%S")"] "

CONDA=/home/zhangJP/miniconda3/envs/MetaTranscriptome
ROUTE=${CONDA}/bin
ADAPTER=${CONDA}/share/trimmomatic-0.39-2/adapters
DB=/home/zhangJP/database

## format the name of raw sequencing data
:<<H
for K in $(ls)
do
    cp ${K}/*_1.fq.gz ${K//[-_]/}_1.fq.gz && cp ${K}/*_2.fq.gz ${K//[-_]/}_2.fq.gz
done

for K in $(ls ${RAW} | awk -F '_' '/gz$/{print $1}' | sort | uniq)
do
    mv ${RAW}/${K}_1.fq.gz ${RAW}/${K,,}_R1.fq.gz && mv ${RAW}/${K}_2.fq.gz ${RAW}/${K,,}_R2.fq.gz
done

for K in $(ls ${RAW} | awk -F '_' '/gz$/{print $1}' | sort | uniq)
do
    mv ${RAW}/${K}_R1.fq.gz ${RAW}/${K^}_R1.fq.gz && mv ${RAW}/${K}_R2.fq.gz ${RAW}/${K^}_R2.fq.gz
done

FILE1=$(ls ${RAW} | grep "fq.gz$" | wc -l)

echo -e "\n" ${FLAG} "A total of" ${FILE1} "files were detected...\n"

echo -e ${FLAG} "Start to asses the quality of these data...\n"

############################ Quality Control ######################

if [[ ! -d ${QUAL} ]]; then
    echo -e ${FLAG} ${QUAL} "is not existing, create it now...\n"
    mkdir -p ${QUAL}
else
    echo -e ${FLAG} ${QUAL} "is existing, delete and create it now...\n"
    rm -r ${QUAL} && mkdir -p ${QUAL}
fi

## FastQC
## https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

cd ${QUAL}

${ROUTE}/fastqc -t 30 \
                -o ${QUAL} \
                ${RAW}/*.gz

${ROUTE}/multiqc ${QUAL}/

cd -

echo -e ${FLAG} "Complete to asses the quality of these data...\n"

if [[ ! -d ${TRIM} ]]; then
    echo -e ${FLAG} ${TRIM} "is not existing, create it now...\n"
    mkdir -p ${TRIM}
else
    echo -e ${FLAG} ${TRIM} "is existing, delete and create it now...\n"
    rm -r ${TRIM} && mkdir -p ${TRIM}
fi

for K in $(ls ${RAW} | grep "gz$" | awk -F '_' '{print $1}'| sort | uniq)
do

    echo -e ${FLAG} "Start to trim the low-quality reads of" ${K} "...\n"

    ############################ Trim LowQuality Bases ######################
    ${ROUTE}/trimmomatic PE ${RAW}/${K}_R1.fq.gz ${RAW}/${K}_R2.fq.gz \
          ${TRIM}/${K}-R1.clean.fq.gz ${TRIM}/${K}-R1.unpair.fq.gz \
          ${TRIM}/${K}-R2.clean.fq.gz ${TRIM}/${K}-R2.unpair.fq.gz \
          ILLUMINACLIP:${ADAPTER}/TruSeq3-PE-2.fa:2:20:10 \
          TOPHRED33 HEADCROP:10 \
          SLIDINGWINDOW:4:20 MINLEN:36 \
          -threads 4

    echo -e ${FLAG} "Complete to trim the low-quality reads of" ${K} "...\n"
done

############################ Identification of Baculoviruse--positive samples ######################

cd 00-BaculoviridaeGenome

bwa index ${BACU}/Baculoviridae.fa \
          -p ${BACU}/Baculoviridae

samtools faidx ${BACU}/Baculoviridae.fa

cd -

if [[ ! -d ${HOST} ]]; then
    echo -e ${FLAG} ${HOST} "is not existing, create it now...\n"
    mkdir -p ${HOST}
else
    echo -e ${FLAG} ${HOST} "is existing, delete and create it now...\n"
    rm -r ${HOST} && mkdir -p ${HOST}
fi

for K in $(ls ${TRIM} | grep ".fq.gz$" | awk -F '-' '{print $1}'| sort | uniq)
do

    bwa mem -R "@RG\tID:foo\tSM:${K}\tLB:library1\tPL:illumina" \
            -t 30 ${BACU}/Baculoviridae \
            ${TRIM}/${K}-R1.clean.fq.gz \
            ${TRIM}/${K}-R2.clean.fq.gz |
    samtools sort -@ 2 \
                  -m 10G \
                  -o ${HOST}/${K}-BAC-sort.bam

    samtools index -@ 2 ${HOST}/${K}-BAC-sort.bam

    cd ${HOST}

    mosdepth -t 10 ${K} ${K}-BAC-sort.bam && rm ${K}-BAC-sort.bam

    cat ${K}.mosdepth.summary.txt | \
                         cut -f 4 | sed '1d'| \
                         awk -v sample=${K} '{if ($1>2) print sample}' |\
                         sort | uniq >> ${MAIN}/bac-positive.txt

    cd -
done

############################ Remove Host Sequences ######################
cd ${DATA}

bwa index ${DATA}/Helicoverpa_armigera.fa \
          -p ${DATA}/Helicoverpa_armigera

samtools faidx ${DATA}/Helicoverpa_armigera.fa

cd -

if [[ ! -d ${CONTIG} ]]; then
    echo -e ${FLAG} ${CONTIG} "is not existing, create it now...\n"
    mkdir -p ${CONTIG}
else
    echo -e ${FLAG} ${CONTIG} "is existing, delete and create it now...\n"
    rm -r ${CONTIG} && mkdir -p ${CONTIG}
fi

for K in $(cat ${MAIN}/bac-positive.txt | sort | uniq)
do
    bwa mem -R "@RG\tID:foo\tSM:${K}\tLB:library1\tPL:illumina" \
            -t 30 ${DATA}/Helicoverpa_armigera \
            ${TRIM}/${K}-R1.clean.fq.gz \
            ${TRIM}/${K}-R2.clean.fq.gz |
        samtools sort -@ 2 -m 10G -o ${HOST}/${K}-sort.bam

    ## samtools sort -n represents sort reads by their name
    samtools view -b -h -f 4 ${HOST}/${K}-sort.bam > ${HOST}/${K}_Un.bam

    samtools sort -n -@ 2 -m 10G -o ${HOST}/${K}_Un.sort.bam ${HOST}/${K}_Un.bam

    bedtools bamtofastq -i ${HOST}/${K}_Un.sort.bam \
                        -fq ${HOST}/${K}-1.Un.fq \
                        -fq2 ${HOST}/${K}-2.Un.fq

    gzip ${HOST}/${K}-1.Un.fq && gzip ${HOST}/${K}-2.Un.fq && rm ${HOST}/${K}-sort.bam ${HOST}/${K}_Un.bam ${HOST}/${K}_Un.sort.bam

    cd ${CONTIG}

    ${ROUTE}/megahit -1 ${HOST}/${K}-1.Un.fq.gz \
                     -2 ${HOST}/${K}-2.Un.fq.gz \
                     --num-cpu-threads 30 \
                     --out-dir ${CONTIG}/${K}-Un-megahit \
                     --out-prefix ${K} && cp ${CONTIG}/${K}-Un-megahit/*.fa ${CONTIG}/${K}-Un.fa

    cd -
done
H

if [[ ! -d ${VIRUSE} ]]; then
    echo -e ${FLAG} ${VIRUSE} "is not existing, create it now...\n"
    mkdir -p ${VIRUSE}
else
    echo -e ${FLAG} ${VIRUSE} "is existing, delete and create it now...\n"
    rm -r ${VIRUSE} && mkdir -p ${VIRUSE}
fi

makeblastdb -in ${BACU}/Baculoviridae.fa -dbtype nucl -out ${BACU}/Baculoviridae

for K in $(cat ${MAIN}/bac-positive.txt | sort | uniq)
do
    blastn -db ${BACU}/Baculoviridae \
           -query ${CONTIG}/${K}-Un.fa \
           -out ${VIRUSE}/${K} \
           -outfmt 6 \
           -num_threads 30 \
           -evalue 1e-8 \
           -num_alignments 1

    cat ${VIRUSE}/${K} | cut -f 1 | sort | uniq > ${VIRUSE}/${K}.hits

    seqkit grep -f ${VIRUSE}/${K}.hits ${CONTIG}/${K}-Un.fa >> ${VIRUSE}/${K}.hits.fa

    diamond blastx --db ${DB}/nr \
                   --query ${VIRUSE}/${K}.hits.fa \
                   --out ${VIRUSE}/${K}.nr \
                   --evalue 1e-5 \
                   --outfmt 6 qseqid qlen sseqid slen pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames salltitles \
                   --sensitive \
                   --threads 30 \
                   --max-target-seqs 1 \
                   --max-hsps 1
done

if [[ ! -d ${HOST} ]]; then
    echo ${HOST} "is not existing, create it now"
    mkdir -p ${HOST}
else
    echo ${HOST} "is existing, delete and create it now"
    rm -r ${HOST} && mkdir -p ${HOST}
fi

if [[ ! -d ${CONTIG} ]]; then
    echo ${CONTIG} "is not existing, create it now"
    mkdir -p ${CONTIG}
else
    echo ${CONTIG} "is existing, delete and create it now"
    rm -r ${CONTIG} && mkdir -p ${CONTIG}
fi

if [[ ! -d ${RNAViruse} ]]; then
    echo ${RNAViruse} "is not existing, create it now"
    mkdir -p ${RNAViruse}
else
    echo ${RNAViruse} "is existing, delete and create it now"
    rm -r ${RNAViruse} && mkdir -p ${RNAViruse}
fi

if [[ ! -d ${NAME} ]]; then
    echo ${NAME} "is not existing, create it now"
    mkdir -p ${NAME}
else
    echo ${NAME} "is existing, delete and create it now"
    rm -r ${NAME} && mkdir -p ${NAME}
fi

source ~/miniconda3/bin/activate MetaTranscriptome

## softwares
## Quality Control: fastQC; MultiQC; Trimmomatic;
## Remove Host Reads: Hisat2; samtools; bedtools;
## Assemble Genome: megahit;
## classify virus:
## Taxon assign: Diamond; blast;
##
## Mafft; IQtree;

## wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/902/850/365/GCA_902850365.1_PGI_CHILSU_v2/GCA_902850365.1_PGI_CHILSU_v2_genomic.fna.gz

DB=/home/zhangJP/database
RdRp=/home/zhangJP/biosoft/NeoRdRp

CONDA=/home/zhangJP/miniconda3/envs/MetaTranscriptome
ROUTE=${CONDA}/bin
ADAPTER=${CONDA}/share/trimmomatic-0.39-2/adapters
T=26

############################   Anticipation  ######################

#hisat2-build -p ${T} ${DATA}/GCA_902850365.1_PGI_CHILSU_v2_genomic.fa ${DATA}/HostGenome

############################ Quality Control ######################

## FastQC
## https://www.bioinformatics.babraham.ac.uk/projects/fastqc/

#cd ${RAW}

#${ROUTE}/fastqc -t ${T} *.gz

#${ROUTE}/multiqc ./

#cd -

for K in $(ls ${RAW} | grep "gz$" | awk -F '-' '{print $1}'| sort | uniq)
do

    ############################ Trim LowQuality Bases ######################
    ${ROUTE}/trimmomatic PE ${RAW}/${K}-1.fq.gz ${RAW}/${K}-2.fq.gz \
          ${TRIM}/${K}-1.clean.fq.gz ${TRIM}/${K}-1.unpair.fq.gz \
          ${TRIM}/${K}-2.clean.fq.gz ${TRIM}/${K}-2.unpair.fq.gz \
          ILLUMINACLIP:${ADAPTER}/TruSeq3-PE-2.fa:2:20:10 \
          TOPHRED33 HEADCROP:10 \
          SLIDINGWINDOW:4:20 MINLEN:36 \
          -threads ${T}

    ############################ Remove Host Sequences ######################

    ## samtools sort -n represents sort reads by their name
    hisat2 -q --dta -x ${DATA}/HostGenome \
                    -p ${T} \
                    -1 ${TRIM}/${K}-1.clean.fq.gz \
                    -2 ${TRIM}/${K}-2.clean.fq.gz |\
    samtools sort -@ 2 -m 10G > ${TEMP}/${K}_sort.bam

    samtools view -b -h -f 4 ${TEMP}/${K}_sort.bam > ${TEMP}/${K}_Un.bam

    samtools sort -n -@ 2 -m 10G -o ${TEMP}/${K}_Un.sort.bam ${TEMP}/${K}_Un.bam

    bedtools bamtofastq -i ${TEMP}/${K}_Un.sort.bam \
                        -fq ${HOST}/${K}-1.Un.fq \
                        -fq2 ${HOST}/${K}-2.Un.fq

    gzip ${HOST}/${K}-1.Un.fq && gzip ${HOST}/${K}-2.Un.fq

    cd ${CONTIG}

    ${ROUTE}/megahit -1 ${TRIM}/${K}-1.clean.fq.gz \
                     -2 ${TRIM}/${K}-2.clean.fq.gz \
                     --num-cpu-threads ${T} \
                     --out-dir ${CONTIG}/${K}-megahit \
                     --out-prefix ${K} && cp ${CONTIG}/${K}-megahit/*.fa ${CONTIG}/${K}.fa

    ${ROUTE}/megahit -1 ${HOST}/${K}-1.Un.fq.gz \
                     -2 ${HOST}/${K}-2.Un.fq.gz \
                     --num-cpu-threads ${T} \
                     --out-dir ${CONTIG}/${K}-Un-megahit \
                     --out-prefix ${K}-Un && cp ${CONTIG}/${K}-Un-megahit/*.fa ${CONTIG}/${K}-Un.fa

    cd ${MAIN}
done

## Taxon Assigned
for K in $(ls ${CONTIG} | grep "fa$" | awk -F '.' '{print $1}')
do
    ## remove contigs shorter than 1kb
    seqkit seq -m 1000 ${CONTIG}/${K}.fa -o ${RNAViruse}/${K}_filter.fa -g

    ## RdRp protein sequences
    diamond blastx --query ${RNAViruse}/${K}_filter.fa \
                   --db ${RdRp}/NeoRdRp-seq \
                   --out ${RNAViruse}/${K}_RdRp.txt \
                   --evalue 1e-5 \
                   --outfmt 6 \
                   --sensitive \
                   --threads 6 \
                   --max-target-seqs 1

    cat ${RNAViruse}/${K}_RdRp.txt | cut -f 1 > ${RNAViruse}/${K}_RdRp_Contig.txt

    seqkit grep -f ${RNAViruse}/${K}_RdRp_Contig.txt --threads 2 ${RNAViruse}/${K}_filter.fa > ${RNAViruse}/${K}_RdRp_Contig.fa

    diamond blastx --db ${DB}/nr \
                   --query ${RNAViruse}/${K}_RdRp_Contig.fa \
                   --out ${RNAViruse}/${K}_Nr.txt \
                   --evalue 1e-5 \
                   --outfmt 6 qseqid qlen sseqid slen pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames salltitles \
                   --sensitive \
                   --threads ${T} \
                   --max-target-seqs 1 \
                   --max-hsps 1

    ## 获取与Nr数据库二次比对上的记录
    cat ${RNAViruse}/${K}_Nr.txt | csvtk -t grep -f 15 -P /home/zhangJP/database/Sub-Nr/Virus.list | cut -f 15 > ${RNAViruse}/${K}_Nr_reciprocal.txid

    cat ${RNAViruse}/${K}_Nr.txt | csvtk -t grep -f 15 -P /home/zhangJP/database/Sub-Nr/Virus.list | cut -f 1 > ${RNAViruse}/${K}_Nr_reciprocal.txt

    ## 获取命中病毒的完整分类信息
    ## taxonkit reformat -F 补充不完整的分类信息
    taxonkit lineage -c ${RNAViruse}/${K}_Nr_reciprocal.txid --delimiter "\t" | taxonkit reformat -F > ${RNAViruse}/${K}_Nr_reciprocal.lineage

    ## 根据比对记录，提取对应的序列
    seqkit grep -f ${RNAViruse}/${K}_Nr_reciprocal.txt --threads 2 ${RNAViruse}/${K}_filter.fa > ${RNAViruse}/${K}_Nr_reciprocal.fa

    ## viral genome
    tblastx -query ${RNAViruse}/${K}_Nr_reciprocal.fa \
           -db ${DB}/viral/viralGenomeMerge \
           -out ${RNAViruse}/${K}_Genome_reciprocal \
           -evalue 1e-5 \
           -outfmt 6 \
           -num_threads ${T}

    cut -f 1 2 ${RNAViruse}/${K}_Genome_reciprocal > ${RNAViruse}/${K}_BestGenome_reciprocal.txt
done

cat ${RNAViruse}/*_Nr_reciprocal.lineage | sort | uniq > ${RNAViruse}/merge_Nr_reciprocal.lineage

for K in $(ls ${CONTIG} | grep "fa$" | awk -F '.' '{print $1}')
do
    ## viral genome
    tblastx -query ${CONTIG}/${K}.fa \
           -db ${DB}/viral/viralGenomeMerge \
           -out ${NAME}/${K}_viralGenome \
           -evalue 1e-5 \
           -outfmt 6 \
           -num_threads ${T}

    ## viral protein
    diamond blastx --query ${CONTIG}/${K}.fa \
                   --db ${DB}/viral/viralProteinMerge \
                   --out ${NAME}/${K}_viralProtein.txt \
                   --evalue 1e-5 \
                   --outfmt 6 \
                   --sensitive \
                   --threads ${T} \
                   --max-target-seqs 1 \
                   --max-hsps 1

    ## Nr
    diamond blastx --db ${DB}/nr \
                   --query ${CONTIG}/${K}.fa \
                   --out ${NAME}/${K}_nr.txt \
                   --evalue 1e-5 \
                   --outfmt 6 qseqid qlen sseqid slen pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames salltitles \
                   --sensitive \
                   --threads ${T} \
                   --max-target-seqs 1 \
                   --max-hsps 1
done
