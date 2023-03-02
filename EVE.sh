Current=$(pwd)

GENOME=/home/zhangJP/LepPublicGenomes/00-genomesDB
ViralProtein=/home/zhangJP/database/viral
Nr=/home/zhangJP/database

#for K in Chrysodeixis_includens Mamestra_brassicae Autographa_gamma Xestia_c-nigrum Mythimna_impura Phlogophora_meticulosa Hydraecia_micacea Xestia_xanthographa Diachrysia_chrysitis Agrotis_puta Xestia_sexstrigata Ochropleura_plecta Autographa_pulchrina Tholera_decimalis Apamea_monoglypha Mesoligia_furuncula Noctua_pronuba Cosmia_trapezina Griposia_aprilina Noctua_fimbriata Mythimna_ferrago Atethmia_centrago Mythimna_albipuncta Noctua_janthe Hecatera_dysodea Luperina_testacea Euplexia_lucipara Abrostola_tripartita Amphipyra_tragopoginis Dryobotodes_eremita Agrochola_macilenta Omphaloscelis_lunosa Diarsia_rubi Brachylomia_viminalis Aporophyla_lueneburgensis Amphipyra_berbera Agrochola_circellaris Xylocampa_areola
for K in Spodoptera_exigua Spodoptera_frugiperda Spodoptera_littoralis Trichoplusia_ni Helicoverpa_armigera Helicoverpa_zea Craniophora_ligustri Caradrina_clavipalpis Allophyes_oxyacanthae Acronicta_aceris Eupsilia_transversa Protodeltote_pygarga Mythimna_separata Chrysodeixis_includens Mamestra_brassicae Autographa_gamma Xestia_c-nigrum Mythimna_impura Phlogophora_meticulosa Hydraecia_micacea Xestia_xanthographa Diachrysia_chrysitis Agrotis_puta Xestia_sexstrigata Ochropleura_plecta Autographa_pulchrina Tholera_decimalis Apamea_monoglypha Mesoligia_furuncula Noctua_pronuba Cosmia_trapezina Griposia_aprilina Noctua_fimbriata Mythimna_ferrago Atethmia_centrago Mythimna_albipuncta Noctua_janthe Hecatera_dysodea Luperina_testacea Euplexia_lucipara Abrostola_tripartita Amphipyra_tragopoginis Dryobotodes_eremita Agrochola_macilenta Omphaloscelis_lunosa Diarsia_rubi Brachylomia_viminalis Aporophyla_lueneburgensis Amphipyra_berbera Agrochola_circellaris Xylocampa_areola
do
    cd ${Current}/10-HighQualityGenome/
    makeblastdb -in ${K}.fa -dbtype nucl -out ${K}

    ## similarity search
    cd ${Current}/20-tblastnSearch
    tblastn -query ${ViralProtein}/viralProteinMerge.fa \
            -db ${Current}/10-HighQualityGenome/${K} \
            -out ${K}-tblastn \
            -evalue 1e-20 \
            -outfmt 6 \
            -num_threads 30
done

## Phylogenetic Tree
for file in $(find 10-HighQualityGenome -name "full_table.tsv")
do
    grep -v "^#" ${file} | awk '$2=="Complete" {print $1}' >> 30-PhylogeneticTree/complete_busco_ids.txt
done

### 获取所有基因组都有的单拷贝基因
sort complete_busco_ids.txt | uniq -c | grep -P "211 " | awk '{print $2}' > 30-PhylogeneticTree/final_busco_ids.txt

### 修改每个文件的名称和其包含序列的名称，即加上物种名称，如Blastobasis_lacticolella
mkdir -p 30-PhylogeneticTree/BUSCOs_Single_Gene
for dir in $(find . -type d -name "single_copy_busco_sequences")
do
    sppname=$(dirname $dir | cut -d "/" -f 3 )
    for file in ${dir}/*.faa
    do
        filename=$(basename ${file})
        cat ${file} | sed 's/^>/>'${sppname}'|/g' > ./30-PhylogeneticTree/BUSCOs_Single_Gene/${sppname}_${filename}
    done
done

### 将不同物种的同一个蛋白质序列整合到单独的文件中
mkdir -p 30-PhylogeneticTree/BUSCOs_Multiple_Genes
for ID in $(ls ./30-PhylogeneticTree/BUSCOs_Single_Gene | grep "faa$" | awk -F "_" '{print $3}'| sort | uniq)
do
    cat ./30-PhylogeneticTree/BUSCOs_Single_Gene/*${ID}*.faa > ./30-PhylogeneticTree/BUSCOs_Multiple_Genes/${ID}
done

mkdir 10-HighQualityGenome
mkdir 20-tblastnSearch

for K in $(ls ${GENOME} | grep "fa.gz$" | awk -F "." '{print $1}')
do
    cd ${Current}/10-HighQualityGenome/

    source ~/miniconda3/bin/activate base
    ## cut off sequences shorter than 10k
    seqkit seq -m 10000 ${GENOME}/${K}.fa.gz -o ${K}.fa -g

    ## Quality Assessment --BUSCOs

    source ~/miniconda3/bin/activate BUSCO5

    busco --mode genome \
          --in ${K}.fa \
          --lineage_dataset /home/zhangJP/biosoft/BUSCO5/lepidoptera_odb10/ \
          --out ${K} \
          --cpu 30 \
          --offline

    ## Phylogenetic Tree
    for file in $(find 10-HighQualityGenome -name "full_table.tsv")
    do
        grep -v "^#" ${file} | awk '$2=="Complete" {print $1}' >> 30-PhylogeneticTree/complete_busco_ids.txt
    done

    ### 获取所有基因组都有的单拷贝基因
    sort complete_busco_ids.txt |uniq -c | grep -P "211 " | awk '{print $2}' > final_busco_ids.txt

    ### 修改每个文件的名称和其包含序列的名称，即加上物种名称，如Blastobasis_lacticolella

    mkdir BUSCOs_Single_Gene
    for dir in $(find . -type d -name "single_copy_busco_sequences")
    do
        sppname=$(dirname $dir | cut -d "/" -f 3 )
        for file in ${dir}/*.faa
        do
            filename=$(basename ${file})
            cat ${file} | sed 's/^>/>'${sppname}'|/g' > ./30-PhylogeneticTree/BUSCOs_Single_Gene/${sppname}_${filename}
        done
    done

    ### 将不同物种的同一个蛋白质序列整合到单独的文件中
    mkdir BUSCOs_Multiple_Genes
    for ID in $(ls ./30-PhylogeneticTree/BUSCOs_Genes | awk -F "_" '{print $3}'| sort | uniq)
    do
        cat *${ID}* > ./30-PhylogeneticTree/BUSCOs_Multiple_Genes/${ID}.faa
    done

    ## Quality Assessment --LAI

    ## git clone https://github.com/oushujun/EDTA.git
    ## cd EDTA/ && conda env create -f EDTA.yml

    source ~/miniconda3/bin/activate EDTA

    perl ~/biosoft/EDTA/EDTA.pl --genome ${K}.fa \
                      --species others \
                      --step all \
                      --threads 30 \
                      --anno 1 \
                      --overwrite 1

    ## build index for different genomes

    source ~/miniconda3/bin/activate base

    makeblastdb -in ${K}.fa -dbtype nucl -out ${K}

    ## similarity search
    cd ${Current}/20-tblastnSearch
    tblastn -query ${ViralProtein}/viralProteinMerge.fa \
            -db ${Current}/10-HighQualityGenome/${K} \
            -out ${K}-tblastn \
            -evalue 1e-20 \
            -outfmt 6 \
            -num_threads 30

    ## extract target nuclear sequences
    ### convert blast result to bed file
    #cut -f 2,9,10 ${K}-tblastn > ${K}_hit.bed

    #bedtools merge -i ${K}_hit.bed -d 10

    #seqkit subseq --bed ${K}_hit.bed -w 0 ${K}.fa -o ${K}_hit.fa

    ## reciprocal similarity search
    #diamond blastx --db ${Nr}/nr \
    #               --query ${K}_hit.fa \
    #               --out ${K}-blastx \
    #               --evalue 1e-5 \
    #               --outfmt 6 qseqid qlen sseqid slen pident length mismatch gapopen qstart qend sstart send evalue bitscore staxids sscinames salltitles \
    #               --sensitive \
    #               --threads 30 \
    #               --max-target-seqs 1 \
    #               --max-hsps 1
#done
