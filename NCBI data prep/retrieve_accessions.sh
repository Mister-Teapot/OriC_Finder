BIOSAMPLES=$(esearch -db nuccore -query "rep_origin [FEATURE] AND bacteria[Organism]" | elink -target biosample | efetch -format docsum | xtract.Linux -pattern DocumentSummary -block Accession -element Accession)
MDATA="mdata.tab"
#echo -e "BioSample\tAssembly" >> ${MDATA} #Put a header in the file
for BIOSAMPLE in ${BIOSAMPLES[@]}
do
    DOCSUM=$(esearch -db assembly -query ${BIOSAMPLE} | efetch -format docsum)
    ASSEMBLY=$(echo ${DOCSUM} | xtract.Linux -pattern DocumentSummary -block Synonym -element RefSeq)
    echo -e ${BIOSAMPLE}'\t'${ASSEMBLY} >> ${MDATA}
done
