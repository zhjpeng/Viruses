#!/bin/bash

###################################################################
############################# 开工大吉#############################
###################################################################


###################################################################
######################### meta-genomics ###########################
###################################################################

MAIN=$(pwd)
GENOME=${MAIN}/10-reference
RAW=${MAIN}/20-rawData
TRIM=${MAIN}/30-cleanData
CONTIG=${MAIN}/40-primaryAssembly
VIRUSE=${MAIN}/50-DatabaseSearch
NAME=${MAIN}/60-taxonClassification

FLAG="[MetaGenomics "$(date "+%Y-%m-%d %H:%M:%S")"] "

THREAD=32

ViralFlye=/home/zhangJP/biosoft/viralFlye
CONDA=/home/zhangJP/miniconda3/envs/MetaTranscriptome
ROUTE=${CONDA}/bin

############################ Download data  ######################
:<<K
source ~/miniconda3/bin/activate aspera
#ascp -QT -k 1 -l 300m -P33001 -i /home/zhangJP/miniconda3/envs/aspera/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/run/ERR694/ERR6942233/XRKUK_20190708_S64018_PL100132152-1_B01.subreads.bam 20-rawData/

##Cnaphalocrocis medinalis
#ascp -QT -k 1 -l 300m -P33001 -i /home/zhangJP/miniconda3/envs/aspera/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR120/031/SRR12014231/SRR12014231_subreads.fastq.gz 20-rawData/Cnaphalocrocis_medinalis/
#ascp -QT -k 1 -l 300m -P33001 -i /home/zhangJP/miniconda3/envs/aspera/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR120/032/SRR12014232/SRR12014232_subreads.fastq.gz 20-rawData/Cnaphalocrocis_medinalis/
#ascp -QT -k 1 -l 300m -P33001 -i /home/zhangJP/miniconda3/envs/aspera/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR120/033/SRR12014233/SRR12014233_subreads.fastq.gz 20-rawData/Cnaphalocrocis_medinalis/
#ascp -QT -k 1 -l 300m -P33001 -i /home/zhangJP/miniconda3/envs/aspera/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR120/034/SRR12014234/SRR12014234_subreads.fastq.gz 20-rawData/Cnaphalocrocis_medinalis/
#ascp -QT -k 1 -l 300m -P33001 -i /home/zhangJP/miniconda3/envs/aspera/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR120/035/SRR12014235/SRR12014235_subreads.fastq.gz 20-rawData/Cnaphalocrocis_medinalis/
#ascp -QT -k 1 -l 300m -P33001 -i /home/zhangJP/miniconda3/envs/aspera/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR120/036/SRR12014236/SRR12014236_subreads.fastq.gz 20-rawData/Cnaphalocrocis_medinalis/
#ascp -QT -k 1 -l 300m -P33001 -i /home/zhangJP/miniconda3/envs/aspera/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR120/036/SRR12014236/SRR12014236_subreads.fastq.gz 20-rawData/Cnaphalocrocis_medinalis/
#ascp -QT -k 1 -l 300m -P33001 -i /home/zhangJP/miniconda3/envs/aspera/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR120/037/SRR12014237/SRR12014237_subreads.fastq.gz 20-rawData/Cnaphalocrocis_medinalis/

## Helicoverpa zea
#ascp -QT -k 1 -l 300m -P33001 -i /home/zhangJP/miniconda3/envs/aspera/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR161/029/SRR16121629/SRR16121629_subreads.fastq.gz 20-rawData/Helicoverpa_zea

## Spodoptera frugiperda
#ascp -QT -k 1 -l 300m -P33001 -i /home/zhangJP/miniconda3/envs/aspera/etc/asperaweb_id_dsa.openssh era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR111/042/SRR11149142/SRR11149142_subreads.fastq.gz 20-rawData/Spodoptera_frugiperda
K

############################ Index reference genome  ######################
source ~/miniconda3/bin/activate viralFlye

## install dependency
#conda install -y hisat2 gffread
:<<U
cd ${GENOME}

for K in $(ls ${GENOME} | grep "fa$" | awk -F '.' '{print $1}'| sort | uniq)
do

    gffread ${GENOME}/${K}.gff3 -T -o ${GENOME}/${K}.gtf
    extract_splice_sites.py ${GENOME}/${K}.gtf > ${GENOME}/${K}.ss
    extract_exons.py ${GENOME}/${K}.gtf > ${GENOME}/${K}.exon

    hisat2-build --ss ${GENOME}/${K}.ss --exon ${GENOME}/${K}.exon ${GENOME}/${K}.fa ${GENOME}/${K}_trans
done

cd ${MAIN}

############################ Quality Control ######################
if [[ ! -d ${TRIM} ]]; then
    echo -e ${FLAG} ${TRIM} "is not existing, create it now...\n"
    mkdir -p ${TRIM}
else
    echo -e ${FLAG} ${TRIM} "is existing, delete and create it now...\n"
    rm -r ${TRIM} && mkdir -p ${TRIM}
fi

if [[ ! -d ${CONTIG} ]]; then
    echo -e ${FLAG} ${CONTIG} "is not existing, create it now...\n"
    mkdir -p ${CONTIG}
else
    echo -e ${FLAG} ${CONTIG} "is existing, delete and create it now...\n"
    rm -r ${CONTIG} && mkdir -p ${CONTIG}
fi

if [[ ! -d ${VIRUSE} ]]; then
    echo -e ${FLAG} ${VIRUSE} "is not existing, create it now...\n"
    mkdir -p ${VIRUSE}
else
    echo -e ${FLAG} ${VIRUSE} "is existing, delete and create it now...\n"
    rm -r ${VIRUSE} && mkdir -p ${VIRUSE}
fi
U

for K in $(cat ${MAIN}/mapping.txt | grep -v "#" | awk -F ':' '{print $1}' | sort | uniq)
do
:<<H
    if [[ ! -d ${TRIM}/${K} ]]; then
        echo -e ${FLAG} ${K} "is not existing, create it now...\n"
        mkdir -p ${TRIM}/${K}
    else
        echo -e ${FLAG} ${K} "is existing, delete and create it now...\n"
        rm -r ${TRIM}/${K} && mkdir -p ${TRIM}/${K}
    fi
H
    ############################ Remove Host Reads ######################

    for J in $(cat ${MAIN}/mapping.txt | grep ${K} | awk -F ':' '{print $2}' | tr ',' '\n' | sort | uniq)
    do
:<<Y
       echo -e ${FLAG} "Start to remove the host reads (" ${K} ") from the sample of" ${J} "...\n"

        hisat2 -p 20 --dta \
               -x ${GENOME}/${K}_trans \
               -1 ${RAW}/${J}_1.fq.gz \
               -2 ${RAW}/${J}_2.fq.gz | \
        samtools sort -@ 2 -m 10G > ${TRIM}/${K}_${J}.sort.bam

        samtools view -b -h -f 4 ${TRIM}/${K}_${J}.sort.bam > ${TRIM}/${K}_${J}.un.bam

        samtools sort -n -@ 2 -m 10G -o ${TRIM}/${K}_${J}.un.sort.bam ${TRIM}/${K}_${J}.un.bam

        bedtools bamtofastq -i ${TRIM}/${K}_${J}.un.sort.bam \
                        -fq ${TRIM}/${K}_${J}.un-1.fq \
                        -fq2 ${TRIM}/${K}_${J}.un-2.fq

        gzip ${TRIM}/${K}_${J}.un-1.fq && gzip ${TRIM}/${K}_${J}.un-2.fq

        echo -e ${FLAG} "Complete to remove the host reads (" ${K} ") from the sample of" ${J} "...\n"

        cd ${TRIM}

        ${ROUTE}/megahit -1 ${TRIM}/${K}_${J}.un-1.fq.gz \
                         -2 ${TRIM}/${K}_${J}.un-2.fq.gz \
                         --num-cpu-threads 30 \
                         --out-dir ${CONTIG}/${K}_${J} \
                         --out-prefix ${K}_${J}.megahit && cp ${CONTIG}/${K}_${J}/*.fa ${CONTIG}/${K}_${J}.megahit.fa

        diamond blastx --query ${CONTIG}/${K}_${J}.megahit.fa \
                       --db  ~/database/viral/viralProteinMerge.dmnd \
                       --out ${VIRUSE}/viral_${K}_${J}.blast \
                       --evalue 1e-5 \
                       --outfmt 6 \
                       --sensitive \
                       --threads 30 \
                       --max-target-seqs 1

        echo -e ${FLAG} "Complete to detect potential viruses in" ${K} "...\n"
Y
        cat ${VIRUSE}/viral_${K}_${J}.blast | cut -f 1 > ${VIRUSE}/viral_${K}_${J}.hits.id

        seqkit grep -f ${VIRUSE}/viral_${K}_${J}.hits.id --threads 2 ${CONTIG}/${K}_${J}.megahit.fa > ${VIRUSE}/viral_${K}_${J}.hits.contig

        diamond blastx --query ${VIRUSE}/viral_${K}_${J}.hits.contig \
                       --db ~/database/nr.dmnd \
                       --out ${VIRUSE}/nr_${K}_${J}.blast \
                       --evalue 1e-5 \
                       --outfmt 6 qseqid qlen sseqid slen pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames salltitles \
                       --sensitive \
                       --threads 32 \
                       --max-target-seqs 1 \
                       --max-hsps 1

        ## 获取与Nr数据库二次比对上的记录
        cat ${VIRUSE}/nr_${K}_${J}.blast | csvtk -t grep -f 15 -P /home/zhangJP/database/Sub-Nr/Virus.list | cut -f 15 | uniq > ${VIRUSE}/nr_${K}_${J}.name

        cat ${VIRUSE}/nr_${K}_${J}.blast | csvtk -t grep -f 15 -P /home/zhangJP/database/Sub-Nr/Virus.list | cut -f 1 | uniq > ${VIRUSE}/nr_${K}_${J}.id    
    done
done

## viral identification

if [[ ! -d ${NAME} ]]; then
    echo -e ${FLAG} ${NAME} "is not existing, create it now...\n"
    mkdir -p ${NAME}
else
    echo -e ${FLAG} ${NAME} "is existing, delete and create it now...\n"
    rm -r ${NAME} && mkdir -p ${NAME}
fi

for K in $(ls ${CONTIG} | grep "fa$" | awk -F '.' '{print $1}')
do
    ## 获取命中病毒的完整分类信息
    ## taxonkit reformat -F 补充不完整的分类信息
    taxonkit lineage -c ${VIRUSE}/nr_${K}_${J}.name --delimiter "  " > ${VIRUSE}/${K}-Nr-hits.lineage

    ## 根据比对记录，提取对应的序列
    seqkit grep -f ${RNAViruse}/${K}-Nr-hits.id --threads 2 ${RNAViruse}/${K}-cdhit.fa > ${NAME}/${K}-Nr-hits.fa
done
