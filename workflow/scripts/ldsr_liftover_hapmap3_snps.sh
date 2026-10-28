#!/bin/bash
#--------------------------------------------------------------------------------------
#
#    Lift HapMap3 SNPs from hg19 to hg38 for LDSC
#
#--------------------------------------------------------------------------------------
#
# Pipeline:  09SLDSR | Rule: lift_hapmap3_snps
#            Upstream:   get_hg19_refs, get_hg38_refs
#            Downstream: ldsr_ld_scores_hg38 (--print-snps)
#
# Purpose:   - Report duplicate rsIDs and chr:pos pairs in the hg19 and hg38 .bim files
#            - Select HapMap3 no-MHC SNPs from the hg19 .bim files
#            - Lift them to hg38 with UCSC liftOver
#            - Keep lifted rsIDs that are present in the hg38 .bim files
#
# Usage:     ldsr_liftover_hapmap3_snps.sh <output_file> [snp_list] [hg19_dir]
#                                          [hg38_dir] [liftover] [chain]
#            (run from workflow/; optional arguments override the defaults below)
#
# Inputs:    snp_list  HapMap3 no-MHC rsID list
#            hg19_dir  Directory of 1000G.EUR.QC.<chr>.bim files
#            hg38_dir  Directory of 1000G.EUR.hg38.<chr>.bim files
#            liftover  UCSC liftOver binary
#            chain     hg19ToHg38 chain file
#
# Outputs:   output_file  Lifted rsIDs present in the hg38 reference, one per line
#            Intermediate files (concatenated .bim, BED files, unmapped SNPs,
#            duplicate reports) are written next to output_file.
#
# Notes:     liftOver and the chain file are downloaded once by hand from
#            http://hgdownload.soe.ucsc.edu (admin/exe/linux.x86_64/liftOver and
#            goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz).
#
#--------------------------------------------------------------------------------------

set -euo pipefail

if [[ $# -lt 1 || $# -gt 6 ]]; then
    echo "Usage: $0 <output_file> [snp_list] [hg19_dir] [hg38_dir] [liftover] [chain]"
    exit 1
fi

OUTFILE="$1"
SNP_LIST="${2:-../resources/ldsr/ldsr_hg19_refs/hm3_no_MHC.list.txt}"
HG19_DIR="${3:-../resources/ldsr/ldsr_hg19_refs/plink_files}"
HG38_DIR="${4:-../resources/ldsr/ldsr_hg38_refs/plink_files}"
LIFTOVER="${5:-../resources/liftover/liftOver}"
CHAIN="${6:-../resources/liftover/hg19ToHg38.over.chain.gz}"

OUTDIR=$(dirname "${OUTFILE}")
mkdir -p "${OUTDIR}"

echo "Variables"
echo "============================"
printf "%-10s %s\n" OUTFILE "${OUTFILE}" SNP_LIST "${SNP_LIST}" HG19_DIR "${HG19_DIR}" \
    HG38_DIR "${HG38_DIR}" LIFTOVER "${LIFTOVER}" CHAIN "${CHAIN}"
echo "============================"
echo

##  Concatenate .bim files and check for duplicates  ----------------------------------
echo "Extracting hg19 coordinates..."
cat "${HG19_DIR}"/1000G.EUR.QC.*.bim > "${OUTDIR}/all_hg19.bim"
cat "${HG38_DIR}"/1000G.EUR.hg38.*.bim > "${OUTDIR}/all_hg38.bim"

for BIM in "${OUTDIR}/all_hg19.bim" "${OUTDIR}/all_hg38.bim"; do
    echo "Checking for rsID, or chr and base pos, duplicates: ${BIM}"

    PREFIX="${OUTDIR}/$(basename "${BIM}" .bim)"

    dup_rsid_count=$(awk '{print $2}' "${BIM}" | sort | uniq -d | wc -l)
    if [[ "${dup_rsid_count}" -gt 0 ]]; then
        echo "  Found ${dup_rsid_count} duplicate rsIDs"
        awk '{print $2}' "${BIM}" | sort | uniq -d > "${PREFIX}_dup_rsid.txt"
        awk 'NR==FNR {rsid[$1]; next} $2 in rsid' "${PREFIX}_dup_rsid.txt" "${BIM}" \
            > "${PREFIX}_lines_with_dup_rsid.txt"
    else
        echo "  No duplicate rsIDs found"
    fi

    dup_pos_count=$(awk '{print $1"\t"$4}' "${BIM}" | sort | uniq -d | wc -l)
    if [[ "${dup_pos_count}" -gt 0 ]]; then
        echo "  Found ${dup_pos_count} duplicate chr:pos pairs"
        awk '{print $1"\t"$4}' "${BIM}" | sort | uniq -d > "${PREFIX}_dup_pos.txt"
        awk 'NR==FNR {pos[$1"\t"$2]; next} ($1"\t"$4) in pos' "${PREFIX}_dup_pos.txt" "${BIM}" \
            > "${PREFIX}_lines_with_dup_pos.txt"
    else
        echo "  No duplicate chr:pos pairs found"
    fi

    echo
done

##  Select HapMap3 SNPs and make hg19 BED  --------------------------------------------
awk 'NR==FNR {s[$1]; next} $2 in s' "${SNP_LIST}" "${OUTDIR}/all_hg19.bim" \
    > "${OUTDIR}/selected_snps_hg19.txt"

# .bim: chr, rsID, cM, pos, A1, A2 -> BED: chr-prefixed, 0-based start, 1-based end
awk '{print "chr"$1 "\t" ($4-1) "\t" $4 "\t" $2}' "${OUTDIR}/selected_snps_hg19.txt" \
    > "${OUTDIR}/hg19_snps.bed"

##  Lift over to hg38  ----------------------------------------------------------------
echo "Lifting over to hg38..."
"${LIFTOVER}" "${OUTDIR}/hg19_snps.bed" "${CHAIN}" \
    "${OUTDIR}/hg38_snps.bed" "${OUTDIR}/unmapped.txt"

cut -f4 "${OUTDIR}/hg38_snps.bed" > "${OUTDIR}/lifted_hg38_rsids.txt"

# grep exits 1 when every SNP lifted, which would stop the script under pipefail
UNMAPPED_COUNT=$({ grep -v "^#" "${OUTDIR}/unmapped.txt" || true; } | wc -l)
echo "Number of rsIDs not lifted: ${UNMAPPED_COUNT}"

##  Keep lifted rsIDs present in hg38 reference  --------------------------------------
echo "Verifying rsIDs in hg38 reference..."
cut -f2 "${OUTDIR}/all_hg38.bim" > "${OUTDIR}/all_hg38_rsids.txt"

comm -12 <(sort "${OUTDIR}/lifted_hg38_rsids.txt") <(sort "${OUTDIR}/all_hg38_rsids.txt") \
    > "${OUTFILE}"

COMMON_COUNT=$(wc -l < "${OUTFILE}")
echo "Number of lifted rsIDs present in hg38 reference: ${COMMON_COUNT}"
echo "Final rsID list for LDSC saved as ${OUTFILE}"

#--------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------
