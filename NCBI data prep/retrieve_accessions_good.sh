#BIOSAMPLES=$(esearch -db nuccore -query "rep_origin [FEATURE] AND bacteria[Organism]" | elink -target biosample | efetch -format docsum | xtract.Linux -pattern DocumentSummary -block Accession -element Accession)
BIOSAMPLES=$(esearch -db nuccore -query "rep_origin [FEATURE] AND bacteria[Organism]" | elink -target assembly | efetch -format docsum | xtract.Linux -pattern DocumentSummary -block Synonym -element RefSeq)
MDATA="mdata.tab"
echo -e "Assembly\tChromosomeAccessions\tfromRefSeq" >> ${MDATA} #Put a header in the file
for BIOSAMPLE in ${BIOSAMPLES}
do
    RESULT=$(esearch -db nuccore -query ${BIOSAMPLE} | efetch -format docsum | xtract.Linux -pattern DocumentSummary -element Genome,AccessionVersion)
    declare -a CHROMOSOMES=()
    if [[ -z "$RESULT" ]]; then
        REFSEQ="False"
    else
        ACC_LIST=(${RESULT//\s+/\s})
        REFSEQ="True"
        flag=false
        for i in ${ACC_LIST[@]}
        do
            if $flag; then
                CHROMOSOMES+=("$i")
                echo "$i"
                flag=false
            elif [[ "$i" == "chromosome" ]]; then
                flag=true
            fi
        done
    fi
    echo "chromies: ${CHROMOSOMES[@]}"
    echo -e ${BIOSAMPLE}'\t'${CHROMOSOMES[@]}'\t'${REFSEQ} >> ${MDATA}
done
