BIOSAMPLES=$(esearch -db nuccore -query "rep_origin [FEATURE] AND bacteria[Organism]" | elink -target biosample | efetch -format docsum | xtract.Linux -pattern DocumentSummary -block Accession -element Accession)
MDATA="mdata2.tab"
echo -e "BioSample\tAccession\tRefSeq" >> ${MDATA} #Put a header in the file
for BIOSAMPLE in ${BIOSAMPLES[@]}
do
    QUERY=$(esearch -db assembly -query ${BIOSAMPLE} | efetch -format docsum)
    ACCESSION=$(echo ${QUERY} | xtract.Linux -pattern DocumentSummary -block Synonym -element RefSeq)
    REFSEQ="True"
    if [[ -z "$ACCESSION" ]]; then
        ACCESSION=$(echo ${QUERY} | xtract.Linux -pattern DocumentSummary -block Synonym -element Genbank)
        REFSEQ="False"
    fi
    echo -e ${BIOSAMPLE}'\t'${ACCESSION}'\t'${REFSEQ} >> ${MDATA}
done
