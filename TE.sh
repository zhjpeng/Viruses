RMER=/home/zhangJP/biosoft/RepeatModeler-2.0.3
RMK=/home/zhangJP/biosoft/RepeatMasker
MELT=/home/zhangJP/biosoft/MELTv2.2.2

#GENOME=/home/zhangJP/LepPublicGenomes/00-Chilo_suppressalis
GENOME=/home/zhangJP/LepPublicGenomes/00-genomesDB
BUSCO_DB=/home/zhangJP/biosoft/BUSCO5/lepidoptera_odb10/

ASSES=/home/zhangJP/LepPublicGenomes/30-TE/10-HighQualityGenome

FLAG="[TE-Finder "$(date "+%Y-%m-%d %H:%M:%S")"] "

:<<P

if [[ ! -d ${ASSES} ]]; then
    echo -e ${FLAG} ${ASSES} "is not existing, create it now...\n"
    mkdir -p ${ASSES}
else
    echo -e ${FLAG} ${ASSES} "is existing, delete and create it now...\n"
    rm -r ${ASSES} && mkdir -p ${ASSES}
fi

P

source ~/miniconda3/bin/activate EDTA

${RMK}/famdb.py -i ${RMK}/Libraries/RepeatMaskerLib.h5 families -f embl -ad 50557 --include-class-in-name > ${ASSES}/insecta.embl

${RMK}/util/buildRMLibFromEMBL.pl ${ASSES}/insecta.embl > ${ASSES}/insecta.fasta

#for K in $(ls ${GENOME} | grep "fa.gz$" | awk -F "." '{print $1}')
for K in Spodoptera_exigua Spodoptera_frugiperda Spodoptera_littoralis Trichoplusia_ni Helicoverpa_armigera Helicoverpa_zea Craniophora_ligustri Caradrina_clavipalpis Allophyes_oxyacanthae Acronicta_aceris Eupsilia_transversa Protodeltote_pygarga Mythimna_separata Chrysodeixis_includens Mamestra_brassicae Autographa_gamma Xestia_c-nigrum Mythimna_impura Phlogophora_meticulosa Hydraecia_micacea Xestia_xanthographa Diachrysia_chrysitis Agrotis_puta Xestia_sexstrigata Ochropleura_plecta Autographa_pulchrina Tholera_decimalis Apamea_monoglypha Mesoligia_furuncula Noctua_pronuba Cosmia_trapezina Griposia_aprilina Noctua_fimbriata Mythimna_ferrago Atethmia_centrago Mythimna_albipuncta Noctua_janthe Hecatera_dysodea Luperina_testacea Euplexia_lucipara Abrostola_tripartita Amphipyra_tragopoginis Dryobotodes_eremita Agrochola_macilenta Omphaloscelis_lunosa Diarsia_rubi Brachylomia_viminalis Aporophyla_lueneburgensis Amphipyra_berbera Agrochola_circellaris Xylocampa_areola
do
    cd ${ASSES}

    if [[ ! -d ${ASSES}/${K} ]]; then
        echo -e ${FLAG} ${K} "is not existing, create it now...\n"
        mkdir -p ${K}
    else
        echo -e ${FLAG} ${K} "is existing, create a new file to replace it...\n"
        rm -r ${K} && mkdir -p ${K}
    fi

    cd ${ASSES}/${K}

    ## filter out sequences shorter than 10k
    ## BUSCO doesn't support gzip file

    zcat ${GENOME}/${K}.fa.gz | seqkit sort -lr | seqkit replace -p ".+" -r "Seq{nr}" --nr-width 5 | seqkit seq -m 10000 -g -o ${ASSES}/${K}/${K}.fa

    ## rename these sequences

    ## Quality Assessment --BUSCOs

    source ~/miniconda3/bin/activate BUSCO5

    busco --mode genome \
          --in ${ASSES}/${K}/${K}.fa \
          --lineage_dataset ${BUSCO_DB} \
          --out ${K} \
          --cpu 30 \
          --offline

    ## Quality Assessment --LAI

    ## install EDTA
    ## git clone https://github.com/oushujun/EDTA.git
    ## cd EDTA/ && conda env create -f EDTA.yml

    source ~/miniconda3/bin/activate EDTA

    perl ~/biosoft/EDTA/EDTA.pl --genome ${ASSES}/${K}/${K}.fa \
                      --species others \
                      --step all \
                      --threads 30 \
                      --anno 1

    LAI -genome ${ASSES}/${K}/${K}.fa \
        -t 30 \
        -intact ${ASSES}/${K}/${K}.fa.mod.EDTA.raw/LTR/${K}.fa.mod.pass.list \
        -all ${ASSES}/${K}/${K}.fa.mod.EDTA.anno/${K}.fa.mod.out

    ## install repeatmodeler2 and repeatmasker
    # conda install repeatmodeler repeatmasker

    ## prediction TE with repeatmodeler2
    ${RMER}/BuildDatabase -name ${K}-TE-DB \
                  -engine ncbi \
                  ${ASSES}/${K}/${K}.fa

    # xxx-families.fa: Consensus sequences
    ${RMER}/RepeatModeler -database ${K}-TE-DB -engine ncbi -pa 20 >& run.out

    # repeatmasker
    cat ../insecta.fasta ${K}-TE-DB-families.fa > ${K}.repeatlib.fa

    ${RMK}/RepeatMasker -e ncbi -pa 8 -xsmall -lib ${K}.repeatlib.fa -gff -dir ./ ${ASSES}/${K}/${K}.fa

done

## estimate the frequency of TE InDel
## download url: https://sourceforge.net/projects/popoolation-te2/files/latest/download
### tar -zvxf MELTv2.2.2.tar.gz
### conda install bowtie2

#java -jar ${MELT}/MELT.jar BuildTransposonZIP
