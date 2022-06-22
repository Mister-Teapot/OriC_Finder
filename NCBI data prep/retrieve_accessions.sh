# Make sure to have EntrezDirect: https://www.ncbi.nlm.nih.gov/books/NBK179288/
# This script is based in part on the README from: https://github.com/schultzm/entrez_direct_tut
# Some entries will return no chromosomes, but still show as fromRefSeq=true. AFAIK, these entries are either shotgun sequences or sequences that have been replaced by different entries.
BIOSAMPLES=$(esearch -db nuccore -query "rep_origin [FEATURE] AND bacteria[Organism]" | elink -target assembly | efetch -format docsum | xtract.Linux -pattern DocumentSummary -block Synonym -element RefSeq)
MDATA="mdata.tab"
echo -e "Assembly\tChromosomeAccession\tfromRefSeq" >> ${MDATA}
for BIOSAMPLE in ${BIOSAMPLES}
do
    RESULT=$(esearch -db nuccore -query ${BIOSAMPLE} | efetch -format docsum | xtract.Linux -pattern DocumentSummary -element Genome,AccessionVersion)
    declare -a CHROMOSOMES=()
    if [[ -z "$RESULT" ]]; then
        REFSEQ="False"
    else
        # This script does not work for shotgun sequences, since they usually have 50+ partial sequences and 1 main one that is not annotated as a chromosome
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
    echo -e ${BIOSAMPLE}'\t'${CHROMOSOMES[@]}'\t'${REFSEQ} >> ${MDATA}
done
