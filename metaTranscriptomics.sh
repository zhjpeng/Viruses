#!/bin/bash

###################################################################
############################# 开工大吉#############################
###################################################################


###################################################################
####################### meta-transcriptome ########################
###################################################################

MAIN=$(pwd)
TEMP=${MAIN}/00-tempFile
RAW=${MAIN}/01-rawFastQ
QUAL=${MAIN}/02-qualityAssesment
TRIM=${MAIN}/03-cleanFastQ
CONTIG=${MAIN}/04-primaryAssembly
RNAViruse=${MAIN}/05-RdRpBlast
NAME=${MAIN}/06-taxonClassification
ABUND=${MAIN}/07-abundenceMeasure

DB=/home/zhangJP/database
RdRp=/home/zhangJP/biosoft/NeoRdRp

CONDA=/home/zhangJP/miniconda3/envs/MetaTranscriptome
ROUTE=${CONDA}/bin
ADAPTER=${CONDA}/share/trimmomatic-0.39-2/adapters
CDHIT=/home/zhangJP/biosoft/cd-hit-v4.8.1-2019-0228
T=26

FLAG="[MetaTrans "$(date "+%Y-%m-%d %H:%M:%S")"] "

if [[ ! -d ${TEMP} ]]; then
    echo -e ${FLAG} ${TEMP} "is not existing, create it now...\n"
    mkdir -p ${TEMP}
else
    echo -e ${FLAG} ${TEMP} "is existing, delete and create it now...\n"
    rm -r ${TEMP} && mkdir -p ${TEMP}
fi

FILE1=$(ls ${RAW} | grep "fq.gz$" | wc -l)
echo -e ${FLAG} "A total of" ${FILE1} "files were detected...\n"

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

${ROUTE}/fastqc -t ${T} ${RAW}/*.gz -o ${QUAL}

${ROUTE}/multiqc ${RAW}/

cd -

echo -e ${FLAG} "Complete to asses the quality of these data...\n"

if [[ ! -d ${TRIM} ]]; then
    echo -e ${FLAG} ${TRIM} "is not existing, create it now...\n"
    mkdir -p ${TRIM}
else
    echo -e ${FLAG} ${TRIM} "is existing, delete and create it now...\n"
    rm -r ${TRIM} && mkdir -p ${TRIM}
fi

for K in $(ls ${RAW} | grep "gz$" | awk -F '_' '{print $1"_"$2}'| sort | uniq)
do

    echo -e ${FLAG} "Start to trim the low-quality reads of"${K}"...\n"

    ############################ Trim LowQuality Bases ######################
    ${ROUTE}/trimmomatic PE ${RAW}/${K}_1.fq.gz ${RAW}/${K}_2.fq.gz \
          ${TRIM}/${K}-1.clean.fq.gz ${TRIM}/${K}-1.unpair.fq.gz \
          ${TRIM}/${K}-2.clean.fq.gz ${TRIM}/${K}-2.unpair.fq.gz \
          ILLUMINACLIP:${ADAPTER}/TruSeq3-PE-2.fa:2:20:10 \
          TOPHRED33 HEADCROP:10 \
          SLIDINGWINDOW:4:20 MINLEN:36 \
          -threads ${T}

    echo -e ${FLAG} "Complete to trim the low-quality reads of" ${K} "...\n"

done

############################ Sequence Assembly ######################

if [[ ! -d ${CONTIG} ]]; then
    echo -e ${FLAG} ${CONTIG} "is not existing, create it now...\n"
    mkdir -p ${CONTIG}
else
    echo -e ${FLAG} ${CONTIG} "is existing, delete and create it now...\n"
    rm -r ${CONTIG} && mkdir -p ${CONTIG}
fi

cd ${CONTIG}
for K in $(ls ${TRIM} | grep "gz$" | awk -F '-' '{print $1}' | sort | uniq)
do
    echo -e ${FLAG} "Start to assemble the reads of" ${K} "...\n"

    ${ROUTE}/megahit -1 ${TRIM}/${K}-1.clean.fq.gz \
                     -2 ${TRIM}/${K}-2.clean.fq.gz \
                     --num-cpu-threads ${T} \
                     --out-dir ${CONTIG}/${K}-megahit \
                     --out-prefix ${K} && cp ${CONTIG}/${K}-megahit/*.fa ${CONTIG}/${K}.fa

    echo -e ${FLAG} "Complete to assemble the reads of" ${K} "...\n"
done

############################ Viral Identification ######################

if [[ ! -d ${RNAViruse} ]]; then
    echo -e ${FLAG} ${RNAViruse} "is not existing, create it now...\n"
    mkdir -p ${RNAViruse}
else
    echo -e ${FLAG} ${RNAViruse} "is existing, delete and create it now...\n"
    rm -r ${RNAViruse} && mkdir -p ${RNAViruse}
fi

## viral taxonomy
for K in $(ls ${CONTIG} | grep "fa$" | awk -F '.' '{print $1}')
do
    ## remove contigs shorter than 1kb
    echo -e ${FLAG} "Start to filter out contigs shorter than 1 kb in" ${K} "...\n"

    seqkit seq -m 1000 ${CONTIG}/${K}.fa -o ${RNAViruse}/${K}-filter-1k.fa -g

    echo -e ${FLAG} "Complete to filter out contigs shorter than 1 kb in" ${K} "...\n"

    ## collaspe contigs redundancy with 99% identity
    echo -e ${FLAG} "Start to remove redundant contigs with 99% identity in" ${K} "...\n"

    ${CDHIT}/cd-hit-est -i ${RNAViruse}/${K}-filter-1k.fa \
               -o ${RNAViruse}/${K}-cdhit.fa \
               -c 0.99 \
               -n 5 \
               -M 28000 \
               -T 32

    echo -e ${FLAG} "Complete to remove redundant contigs with 99% identity in" ${K} "...\n"

    ## RdRp protein sequences

    echo -e ${FLAG} "Start to detect potential RNA viruses in" ${K} "...\n"

    diamond blastx --query ${RNAViruse}/${K}-cdhit.fa \
                   --db ${RdRp}/NeoRdRp-seq \
                   --out ${RNAViruse}/${K}-RdRp-hits.txt \
                   --evalue 1e-5 \
                   --outfmt 6 \
                   --sensitive \
                   --threads 30 \
                   --max-target-seqs 1

    echo -e ${FLAG} "Complete to detect potential RNA viruses in" ${K} "...\n"

    cat ${RNAViruse}/${K}-RdRp-hits.txt | cut -f 1 > ${RNAViruse}/${K}-RdRp-hits.id

    seqkit grep -f ${RNAViruse}/${K}-RdRp-hits.id --threads 2 ${RNAViruse}/${K}-cdhit.fa > ${RNAViruse}/${K}-RdRp-hits.contig

    diamond blastx --db ${DB}/nr \
                   --query ${RNAViruse}/${K}-RdRp-hits.contig \
                   --out ${RNAViruse}/${K}-Nr-hits.txt \
                   --evalue 1e-5 \
                   --outfmt 6 qseqid qlen sseqid slen pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames salltitles \
                   --sensitive \
                   --threads ${T} \
                   --max-target-seqs 1 \
                   --max-hsps 1

    ## 获取与Nr数据库二次比对上的记录
    cat ${RNAViruse}/${K}-Nr-hits.txt | csvtk -t grep -f 15 -P /home/zhangJP/database/Sub-Nr/Virus.list | cut -f 15 | uniq > ${RNAViruse}/${K}-Nr-hits.name

    cat ${RNAViruse}/${K}-Nr-hits.txt | csvtk -t grep -f 15 -P /home/zhangJP/database/Sub-Nr/Virus.list | cut -f 1 | uniq > ${RNAViruse}/${K}-Nr-hits.id

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
    taxonkit lineage -c ${RNAViruse}/${K}-Nr-hits.name --delimiter "  " > ${NAME}/${K}-Nr-hits.lineage

    ## 根据比对记录，提取对应的序列
    seqkit grep -f ${RNAViruse}/${K}-Nr-hits.id --threads 2 ${RNAViruse}/${K}-cdhit.fa > ${NAME}/${K}-Nr-hits.fa
done

## the abundence of RNA viruses

if [[ ! -d ${ABUND} ]]; then
    echo -e ${FLAG} ${ABUND} "is not existing, create it now...\n"
    mkdir -p ${ABUND}
else
    echo -e ${FLAG} ${ABUND} "is existing, delete and create it now...\n"
    rm -r ${ABUND} && mkdir -p ${ABUND}
fi

for K in $(ls ${CONTIG} | grep "fa$" | awk -F '.' '{print $1}')
do
    bwa index ${NAME}/${K}-Nr-hits.fa -p ${ABUND}/${K}

    bwa mem -t 30 ${ABUND}/${K} \
            ${TRIM}/${K}-1.clean.fq.gz \
            ${TRIM}/${K}-2.clean.fq.gz |\
    samtools sort -m 1g --threads 5 -o ${ABUND}/${K}.bam

    samtools index ${ABUND}/${K}.bam -@ 10

    samtools flagstat -@ 10 ${ABUND}/${K}.bam > ${ABUND}/${K}.flagstat

    mosdepth -t 30 ${ABUND}/${K}  ${ABUND}/${K}.bam
done

############################ phylogenetic relationship ######################
## download the reference genome
mkdir Amalgaviridae Articulavirales Astroviridae Birnaviridae Bunyavirales Cystoviridae Flaviviridae Hepe-Virga Hypoviridae Luteo-Sobemo Mono-Chu Narna-Levi Nidovirales Nodaviridae Partiti-Picobirna Picornavirales Potyviridae Reovirales Tolivirales Toti-Chryso

## 1511860 Amalgaviridae family
taxonkit list -j 30 --ids 1511860 --indent "" --data-dir ~/database/Sub-Nr/taxdmp/ > Amalgaviridae/Amalgaviridae.list

ncbi-genome-download --section genbank --formats fasta --taxids Amalgaviridae/Amalgaviridae.list viral

cp ./genbank/viral/*/*fna.gz Amalgaviridae && rm -r ./genbank/

## 2499411 Articulavirales order
taxonkit list -j 30 --ids 2499411 --indent "" --data-dir ~/database/Sub-Nr/taxdmp/ > Articulavirales/Articulavirales.list

ncbi-genome-download --section genbank --formats fasta --taxids Articulavirales/Articulavirales.list viral

cp ./genbank/viral/*/*fna.gz Articulavirales && rm -r ./genbank/

## 39733 Astroviridae family
taxonkit list -j 30 --ids 39733 --indent "" --data-dir ~/database/Sub-Nr/taxdmp/ > Astroviridae/Astroviridae.list

ncbi-genome-download --section genbank --formats fasta --taxids Astroviridae/Astroviridae.list viral

cp ./genbank/viral/*/*fna.gz Astroviridae && rm -r ./genbank/

## 10993 Birnaviridae family, dsRNA
taxonkit list -j 30 --ids 10993 --indent "" --data-dir ~/database/Sub-Nr/taxdmp/ > Birnaviridae/Birnaviridae.list

ncbi-genome-download --section genbank --formats fasta --taxids Birnaviridae/Birnaviridae.list viral

cp ./genbank/viral/*/*fna.gz Birnaviridae && rm -r ./genbank/

## 1980410 Bunyavirales order, RNA virus
taxonkit list -j 30 --ids 1980410 --indent "" --data-dir ~/database/Sub-Nr/taxdmp/ > Bunyavirales/Bunyavirales.list

ncbi-genome-download --section genbank --formats fasta --taxids Bunyavirales/Bunyavirales.list viral

cp ./genbank/viral/*/*fna.gz Bunyavirales && rm -r ./genbank/

## 10877 Cystoviridae family, dsRNA
taxonkit list -j 30 --ids 10877 --indent "" --data-dir ~/database/Sub-Nr/taxdmp/ > Cystoviridae/Cystoviridae.list

ncbi-genome-download --section genbank --formats fasta --taxids Cystoviridae/Cystoviridae.list viral

cp ./genbank/viral/*/*fna.gz Cystoviridae && rm -r ./genbank/

## 11050 Flaviviridae family, ssRNA(+)
taxonkit list -j 30 --ids 11050 --indent "" --data-dir ~/database/Sub-Nr/taxdmp/ > Flaviviridae/Flaviviridae.list

ncbi-genome-download --section genbank --formats fasta --taxids Flaviviridae/Flaviviridae.list viral

cp ./genbank/viral/*/*fna.gz Flaviviridae && rm -r ./genbank/

## (Hepe-Virga) Virgaviridae, Bromoviridae, Closteroviridae, Togaviridae, Endornaviridae, Tymovirales, Hepeviridae, Alphatetraviridae
## 675071 Virgaviridae family, ssRNA(+)
## 39740 Bromoviridae family, ssRNA(+)
## 69973 Closteroviridae family, ssRNA(+)
## 11018 Togaviridae family, ssRNA(+)
## 564644 Endornaviridae family, dsRNA
## 675063 Tymovirales order, ssRNA(+)
## 291484 Hepeviridae family, ssRNA(+)
## 1283207 Alphatetraviridae family, ssRNA(+)
for K in 675071 39740 69973 11018 564644 675063 291484 1283207
do
    taxonkit list -j 30 --ids ${K} --indent "" --data-dir ~/database/Sub-Nr/taxdmp/ > Hepe-Virga/${K}.list

    ncbi-genome-download --section genbank --formats fasta --taxids Hepe-Virga/${K}.list viral

    cp ./genbank/viral/*/*fna.gz Hepe-Virga/ && rm -r ./genbank/
done

## 39748 Hypoviridae family, ssRNA(+)
taxonkit list -j 30 --ids 39748 --indent "" --data-dir ~/database/Sub-Nr/taxdmp/ > Hypoviridae/Hypoviridae.list

ncbi-genome-download --section genbank --formats fasta --taxids Hypoviridae/Hypoviridae.list viral

cp ./genbank/viral/*/*fna.gz Hypoviridae/ && rm -r ./genbank/

## (Luteo-Sobemo) Luteovirus, Solemoviridae, Barnaviridae
## 12036 Luteovirus genus, ssRNA(+)
## 2169577 Solemoviridae family, ssRNA(+)
## 39741 Barnaviridae family, ssRNA(+)
for K in 12036 2169577 39741
do
    taxonkit list -j 30 --ids ${K} --indent "" --data-dir ~/database/Sub-Nr/taxdmp/ > Luteo-Sobemo/${K}.list

    ncbi-genome-download --section genbank --formats fasta --taxids Luteo-Sobemo/${K}.list viral

    cp ./genbank/viral/*/*fna.gz Luteo-Sobemo/ && rm -r ./genbank/
done

## (Mono-Chu) Chuviridae, Nyamiviridae, Mymonaviridae, Rhabdoviridae, Paramyxoviridae
## 2501952 Chuviridae family, ssRNA(-)
## 1513294 Nyamiviridae family, ssRNA(-)
## 1883067 Mymonaviridae family, ssRNA(-)
## 11270 Rhabdoviridae family, ssRNA(-)
## 11158 Paramyxoviridae family, ssRNA(-)
for K in 2501952 1513294 1883067 11270 11158
do
    taxonkit list -j 30 --ids ${K} --indent "" --data-dir ~/database/Sub-Nr/taxdmp/ > Luteo-Sobemo/${K}.list

    ncbi-genome-download --section genbank --formats fasta --taxids Luteo-Sobemo/${K}.list viral

    cp ./genbank/viral/*/*fna.gz Luteo-Sobemo/ && rm -r ./genbank/
done

## （Narna-Levi）Narnaviridae, Mitoviridae, Botourmiaviridae, Leviviridae
## Narnaviridae 186766 family, ssRNA(+)
## Mitoviridae 2732892 family, ssRNA(+)
## Botourmiaviridae 2560063 family, ssRNA(+)
## Leviviridae 2842319 family, ssRNA(+)
for K in 186766 2732892 2560063 2842319
do
    taxonkit list -j 30 --ids ${K} --indent "" --data-dir ~/database/Sub-Nr/taxdmp/ > Narna-Levi/${K}.list

    ncbi-genome-download --section genbank --formats fasta --taxids Narna-Levi/${K}.list viral

    cp ./genbank/viral/*/*fna.gz Narna-Levi/ && rm -r ./genbank/
done

## Nidovirales 76804 order, ssRNA(+)
taxonkit list -j 30 --ids 76804 --indent "" --data-dir ~/database/Sub-Nr/taxdmp/ > Nidovirales/Nidovirales.list

ncbi-genome-download --section genbank --formats fasta --taxids Nidovirales/Nidovirales.list viral

cp ./genbank/viral/*/*fna.gz Nidovirales/ && rm -r ./genbank/

## Nodaviridae 12283 family, ssRNA(+)
taxonkit list -j 30 --ids 12283 --indent "" --data-dir ~/database/Sub-Nr/taxdmp/ > Nodaviridae/Nodaviridae.list

ncbi-genome-download --section genbank --formats fasta --taxids Nodaviridae/Nodaviridae.list viral

cp ./genbank/viral/*/*fna.gz Nodaviridae/ && rm -r ./genbank/

## (Partiti-Picobirna) Partitiviridae, Picobirnaviridae
## Partitiviridae 11012 family, dsRNA
## Picobirnaviridae 585893 family, dsRNA
for K in 11012 585893
do
    taxonkit list -j 30 --ids ${K} --indent "" --data-dir ~/database/Sub-Nr/taxdmp/ > Partiti-Picobirna/${K}.list

    ncbi-genome-download --section genbank --formats fasta --taxids Partiti-Picobirna/${K}.list viral

    cp ./genbank/viral/*/*fna.gz Partiti-Picobirna/ && rm -r ./genbank/
done

## Picornavirales 464095 order, ssRNA(+)
taxonkit list -j 30 --ids 464095 --indent "" --data-dir ~/database/Sub-Nr/taxdmp/ > Picornavirales/Picornavirales.list

ncbi-genome-download --section genbank --formats fasta --taxids Picornavirales/Picornavirales.list viral

cp ./genbank/viral/*/*fna.gz Picornavirales/ && rm -r ./genbank/

## Potyviridae 39729 family, ssRNA(+)
taxonkit list -j 30 --ids 39729 --indent "" --data-dir ~/database/Sub-Nr/taxdmp/ > Potyviridae/Potyviridae.list

ncbi-genome-download --section genbank --formats fasta --taxids Potyviridae/Potyviridae.list viral

cp ./genbank/viral/*/*fna.gz Potyviridae/ && rm -r ./genbank/

## Reovirales 2732541 order, dsRNA
taxonkit list -j 30 --ids 39729 --indent "" --data-dir ~/database/Sub-Nr/taxdmp/ > Reovirales/Reovirales.list

ncbi-genome-download --section genbank --formats fasta --taxids Reovirales/Reovirales.list viral

cp ./genbank/viral/*/*fna.gz Reovirales/ && rm -r ./genbank/

## Tolivirales 2732548 order, ssRNA(+)
taxonkit list -j 30 --ids 2732548 --indent "" --data-dir ~/database/Sub-Nr/taxdmp/ > Tolivirales/Tolivirales.list

ncbi-genome-download --section genbank --formats fasta --taxids Tolivirales/Tolivirales.list viral

cp ./genbank/viral/*/*fna.gz Tolivirales/ && rm -r ./genbank/

## (Toti-Chryso) Totiviridae, Chrysoviridae, Quadriviridae
## Totiviridae 11006 family, dsRNA
## Chrysoviridae 249310 family, dsRNA
## Quadriviridae 1299296 family, dsRNA
for K in 11006 249310 1299296
do
    taxonkit list -j 30 --ids ${K} --indent "" --data-dir ~/database/Sub-Nr/taxdmp/ > Toti-Chryso/${K}.list

    ncbi-genome-download --section genbank --formats fasta --taxids Toti-Chryso/${K}.list viral

    cp ./genbank/viral/*/*fna.gz Toti-Chryso/ && rm -r ./genbank/
done
