#!/bin/bash

# Check if output file is provided
if [ $# -ne 1 ]; then
    echo "Usage: $0 <output_file>"
    exit 1
fi

OUTFILE="$1"

# Extract directory from output file
OUTDIR=$(dirname "${OUTFILE}")
mkdir -p "${OUTDIR}"  # Create directory if it doesn't exist

# Define paths - should probs not hardcode these change later
SNP_LIST=../resources/ldsr/ldsr_hg19_refs//hm3_no_MHC.list.txt
HG19_DIR=../resources/ldsr/ldsr_hg19_refs/plink_files
HG38_DIR=../resources/ldsr/ldsr_hg38_refs/plink_files

# Step A: Get hg37 (hg19) coordinates
echo "Extracting hg19 coordinates..."

# Concatenate all hg19 bim files
cat ${HG19_DIR}/1000G.EUR.QC.*.bim > "${OUTDIR}/all_hg19.bim"
cat ${HG38_DIR}/1000G.EUR.hg38.*.bim > "${OUTDIR}/all_hg38.bim"

for BIM in "${OUTDIR}/all_hg19.bim" "${OUTDIR}/all_hg38.bim"; do
    echo "Checking for rsID, or chr and base pos, duplicates: $BIM"

    PREFIX=$(basename "$BIM" .bim)

    # Check for duplicate rsIDs
    dup_rsid_count=$(awk '{print $2}' "$BIM" | sort | uniq -d | wc -l)
    if [ "$dup_rsid_count" -gt 0 ]; then
        echo "  Found $dup_rsid_count duplicate rsIDs"
        awk '{print $2}' "$BIM" | sort | uniq -d > "${PREFIX}_dup_rsid.txt"
        awk 'NR==FNR {rsid[$1]; next} $2 in rsid' "${PREFIX}_dup_rsid.txt" "$BIM" > "${OUTDIR}/${PREFIX}_lines_with_dup_rsid.txt"
    else
        echo "  No duplicate rsIDs found"
    fi

    # Check for duplicate chr+pos
    dup_pos_count=$(awk '{print $1"\t"$4}' "$BIM" | sort | uniq -d | wc -l)
    if [ "$dup_pos_count" -gt 0 ]; then
        echo "  Found $dup_pos_count duplicate chr:pos pairs"
        awk '{print $1"\t"$4}' "$BIM" | sort | uniq -d > "${PREFIX}_dup_pos.txt"
        awk 'NR==FNR {pos[$1"\t"$2]; next} ($1"\t"$4) in pos' "${PREFIX}_dup_pos.txt" "$BIM" > "${OUTDIR}/${PREFIX}_lines_with_dup_pos.txt"
    else
        echo "  No duplicate chr:pos pairs found"
    fi

    echo
done

# Extract SNPs matching the rsID list
awk 'NR==FNR {s[$1]; next} $2 in s' ${SNP_LIST} "${OUTDIR}/all_hg19.bim" > "${OUTDIR}/selected_snps_hg19.txt"

# Create BED file for liftOver (chr, start, end, rsID)
# - bim format: chr, rsID, cm, pos, a1, a2
# - BED needs "chr" prefix and 0-based start (pos-1), 1-based end (pos)
awk '{print "chr"$1 "\t" ($4-1) "\t" $4 "\t" $2}' "${OUTDIR}/selected_snps_hg19.txt" > "${OUTDIR}/hg19_snps.bed"

# Step B: Lift over to hg38
echo "Lifting over to hg38..."

# Download chain file if not present (uncomment if needed)
#wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/liftOver
#wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz
#wget http://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz

# Run liftOver
../resources/liftover/liftOver "${OUTDIR}/hg19_snps.bed" ../resources/liftover/hg19ToHg38.over.chain.gz "${OUTDIR}/hg38_snps.bed" "${OUTDIR}/unmapped.txt"

# Extract lifted rsIDs
cut -f4 "${OUTDIR}/hg38_snps.bed" > "${OUTDIR}/lifted_hg38_rsids.txt"

# Track SNPs lost during liftOver
UNMAPPED_COUNT=$(grep -v "^#" "${OUTDIR}/unmapped.txt" | wc -l)
echo "Number of rsIDs not lifted: ${UNMAPPED_COUNT}"

# Verify rsIDs in hg38 reference files
echo "Verifying rsIDs in hg38 reference..."

# Concatenate all hg38 bim rsIDs
cat ${HG38_DIR}/1000G.EUR.hg38.*.bim | cut -f2 > "${OUTDIR}/all_hg38_rsids.txt"

# Find common rsIDs between lifted list and hg38 reference
comm -12 <(sort "${OUTDIR}/lifted_hg38_rsids.txt") <(sort "${OUTDIR}/all_hg38_rsids.txt") > "${OUTFILE}"

# Count common rsIDs
COMMON_COUNT=$(wc -l < "${OUTFILE}")
echo "Number of lifted rsIDs present in hg38 reference: ${COMMON_COUNT}"

# Output final list
echo "Final rsID list for LDSC saved as ${OUTFILE}"

