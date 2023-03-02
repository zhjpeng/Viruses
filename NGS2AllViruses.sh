#!/bin/bash

###################################################################
######################### meta-genomics ###########################
###################################################################

MAIN=$(pwd)
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
CDHIT=/home/zhangJP/biosoft/cd-hit-v4.8.1-2019-0228
ROUTE=${CONDA}/bin
ADAPTER=${CONDA}/share/trimmomatic-0.39-2/adapters
DB=/home/zhangJP/database

:<<L

if [[ ! -d ${RAW} ]]; then
    echo -e ${FLAG} ${RAW} "is not existing, create it now...\n"
    mkdir -p ${RAW}
else
    echo -e ${FLAG} ${RAW} "is existing, delete and create it now...\n"
    rm -r ${RAW} && mkdir -p ${RAW}
fi

cp ${MAIN}/*_tsv.txt ${RAW}

source ~/miniconda3/bin/activate aspera

for K in $(ls ${RAW} | grep "_tsv.txt" | awk -F "_" '{print $4}')
do
    cat ${RAW}/filereport_read_run_${K}_tsv.txt | awk -F "\t" '{print $2}' | sed -e '1d' -e 's/;/\n/' > ${RAW}/${K}_tsv.txt

    cat ${RAW}/${K}_tsv.txt | while read id; do ascp -QT -k 1 -l 300m -P33001 -i /home/zhangJP/miniconda3/envs/aspera/etc/asperaweb_id_dsa.openssh era-fasp@$id ${RAW};done
done

## download reference genome and make index
if [[ ! -d ${DATA} ]]; then
    echo -e ${FLAG} ${DATA} "is not existing, create it now...\n"
    mkdir -p ${DATA}
else
    echo -e ${FLAG} ${DATA} "is existing, delete and create it now...\n"
    rm -r ${DATA} && mkdir -p ${DATA}
fi

cd ${DATA}

### A. cerana
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/442/555/GCF_001442555.1_ACSNU-2.0/GCF_001442555.1_ACSNU-2.0_genomic.fna.gz && mv GCF_001442555.1_ACSNU-2.0_genomic.fna.gz Acerana.fa.gz

bwa index Acerana.fa.gz -p Acerana

### A. mellifera
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/254/395/GCF_003254395.2_Amel_HAv3.1/GCF_003254395.2_Amel_HAv3.1_genomic.fna.gz && mv GCF_003254395.2_Amel_HAv3.1_genomic.fna.gz Amellifera.fa.gz

bwa index Amellifera.fa.gz -p Amellifera

cd -

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
                ${RAW}/PRJNA787435/*.gz

${ROUTE}/multiqc ${QUAL}

cd -

echo -e ${FLAG} "Complete to asses the quality of these data...\n"
L

if [[ ! -d ${TRIM} ]]; then
    echo -e ${FLAG} ${TRIM} "is not existing, create it now...\n"
    mkdir -p ${TRIM}
else
    echo -e ${FLAG} ${TRIM} "is existing, delete and create it now...\n"
    rm -r ${TRIM} && mkdir -p ${TRIM}
fi

if [[ ! -d ${HOST} ]]; then
    echo -e ${FLAG} ${HOST} "is not existing, create it now...\n"
    mkdir -p ${HOST}
else
    echo -e ${FLAG} ${HOST} "is existing, delete and create it now...\n"
    rm -r ${HOST} && mkdir -p ${HOST}
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

for K in $(ls ${RAW}/PRJNA787435 | grep "gz$" | awk -F '_' '{print $1}'| sort | uniq)
do
    echo -e ${FLAG} "Start to trim the low-quality reads of" ${K} "...\n"

    ############################ Trim LowQuality Bases ######################
    ${ROUTE}/trimmomatic PE ${RAW}/PRJNA787435/${K}_1.fastq.gz ${RAW}/PRJNA787435/${K}_2.fastq.gz \
          ${TRIM}/${K}-R1.clean.fq.gz ${TRIM}/${K}-R1.unpair.fq.gz \
          ${TRIM}/${K}-R2.clean.fq.gz ${TRIM}/${K}-R2.unpair.fq.gz \
          ILLUMINACLIP:${ADAPTER}/TruSeq3-PE-2.fa:2:20:10 \
          TOPHRED33 HEADCROP:10 \
          SLIDINGWINDOW:4:20 MINLEN:36 \
          -threads 4

    echo -e ${FLAG} "Complete to trim the low-quality reads of" ${K} "...\n"

    bwa mem -t 30 \
            -R "@RG\tID:foo\tSM:${K}\tLB:library1" \
            ${DATA}/Amellifera \
            ${TRIM}/${K}-R1.clean.fq.gz \
            ${TRIM}/${K}-R2.clean.fq.gz | \
    samtools sort -@ 2 -O bam -o ${HOST}/${K}.bam

    samtools view -b -h -f 4 ${HOST}/${K}.bam > ${HOST}/${K}-Un.bam

    ## samtools sort -n represents sort reads by their name
    samtools sort -n -@ 2 -m 10G -o ${HOST}/${K}_Un.sort.bam ${HOST}/${K}_Un.bam

    bedtools bamtofastq -i ${HOST}/${K}-Un.bam -fq ${TRIM}/${K}-1.Un.fq -fq2 ${TRIM}/${K}-2.Un.fq

    gzip ${TRIM}/${K}-1.Un.fq && gzip ${TRIM}/${K}-2.Un.fq

    cd ${CONTIG}

    ${ROUTE}/megahit -1 ${TRIM}/${K}-1.Un.fq.gz \
                     -2 ${TRIM}/${K}-2.Un.fq.gz \
                     --num-cpu-threads 30 \
                     --out-dir ./${K}-megahit \
                     --out-prefix ${K}

    ${CDHIT}/cd-hit-est -i ${K}-megahit/${K}.contigs.fa -o ${CONTIG}/${K}.cdhit.fa -c 0.95 -n 8 -M 28000 -T 30

    seqkit seq -m 1000 ${CONTIG}/${K}.cdhit.fa -o ${CONTIG}/${K}.1000.fa -g

    cd -

    diamond blastx --query ${CONTIG}/${K}.1000.fa \
                   --db  ~/database/viral/viralProteinMerge.dmnd \
                   --out ${VIRUSE}/viral_${K}.blast \
                   --evalue 1e-5 \
                   --outfmt 6 \
                   --sensitive \
                   --threads 30 \
                   --max-target-seqs 1

    echo -e ${FLAG} "Complete to detect potential viruses in" ${K} "...\n"

    cat ${VIRUSE}/viral_${K}.blast | cut -f 1 > ${VIRUSE}/viral_${K}.hits.id

    seqkit grep -f ${VIRUSE}/viral_${K}.hits.id --threads 2 ${CONTIG}/${K}.1000.fa > ${VIRUSE}/viral_${K}.hits.contig

    diamond blastx --query ${VIRUSE}/viral_${K}.hits.contig \
                   --db ~/database/nr.dmnd \
                   --out ${VIRUSE}/nr_${K}.blast \
                   --evalue 1e-5 \
                   --outfmt 6 qseqid qlen sseqid slen pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames salltitles \
                   --sensitive \
                   --threads 30 \
                   --max-target-seqs 1 \
                   --max-hsps 1

    ## 获取与Nr数据库二次比对上的记录
    cat ${VIRUSE}/nr_${K}.blast | csvtk -t grep -f 15 -P /home/zhangJP/database/Sub-Nr/Virus.list | cut -f 15 | uniq > ${VIRUSE}/nr_${K}.name

    cat ${VIRUSE}/nr_${K}.blast | csvtk -t grep -f 15 -P /home/zhangJP/database/Sub-Nr/Virus.list | cut -f 1 | uniq > ${VIRUSE}/nr_${K}.id
done
