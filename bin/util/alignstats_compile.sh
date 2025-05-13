#!/bin/bash

# Variables - replace MDIR with your directory path
MDIR="/path/to/your/directory/"
LOGFILE="${MDIR}logs/alignstats_summary_compile.log"

# Ensure logs directory exists
mkdir -p ${MDIR}logs/as

echo "STARTcompileAstats" > $LOGFILE

# Find first file and extract header
FIRST_FILE=$(find ${MDIR}*/align/*/alignqc/alignstats/*alignstats.tsv | head -n 1)
head -n 1 "$FIRST_FILE" > ${MDIR}other_reports/alignstats_bsummary.tsv
echo "a_${FIRST_FILE}" >> $LOGFILE

# Find all files and concatenate their last lines
find ${MDIR}*/align/*/alignqc/alignstats/*alignstats.tsv | while read file; do
    tail -n 1 "$file" >> ${MDIR}other_reports/alignstats_bsummary.tsv
    echo "b_${file}" >> $LOGFILE
done

# Copy file to create other output TSV files
cp ${MDIR}other_reports/alignstats_bsummary.tsv ${MDIR}other_reports/alignstats_csummary.tsv
cp ${MDIR}other_reports/alignstats_bsummary.tsv ${MDIR}other_reports/alignstats_combo_mqc.tsv

# Convert spaces to tabs in alignstats_combo_mqc.tsv
perl -pi -e 's/ /\t/g;' ${MDIR}other_reports/alignstats_combo_mqc.tsv

# Create final output by copying alignstats_combo_mqc.tsv
cp ${MDIR}other_reports/alignstats_combo_mqc.tsv ${MDIR}other_reports/alignstats_gs_mqc.tsv

# Error handling
if [ $? -ne 0 ]; then
    touch ${MDIR}logs/ALIGNSTATSCOMPIEFAILEDw_$?
fi
