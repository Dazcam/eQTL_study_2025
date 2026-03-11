set -euo pipefail

# Set variables -----
OUTDIR=$1
OUTFILE_NAME=$2
CELL_TYPES_STR=$3
THREADS=$4

BASE_RESULTS_DIR="../results/12CTWAS"
SOURCE_DIR_NAME="weights"
ARCHIVE_ROOT_NAME="fusion_weights"
COMPRESSION_LEVEL=6

SOURCE_DIR="${BASE_RESULTS_DIR}/${SOURCE_DIR_NAME}"
OUTFILE="${OUTDIR}/${OUTFILE_NAME}"

mkdir -p "${OUTDIR}"

if [[ -z "$THREADS" || ! "$THREADS" =~ ^[0-9]+$ ]]; then
    echo "Error: THREADS must be a positive integer, got: '$THREADS'" >&2
    exit 2
fi

if [[ -z "$CELL_TYPES_STR" ]]; then
    echo "Error: No cell types provided" >&2
    exit 2
fi

read -r -a CELL_TYPES <<< "$CELL_TYPES_STR"
echo "Parsed cell types (${#CELL_TYPES[@]}): ${CELL_TYPES[*]}" >&2

# Build tar args -----
tar_args=()

for ct in "${CELL_TYPES[@]}"; do
  tar_args+=( "${SOURCE_DIR_NAME}/${ct}" )
done

# Create tarball -----
echo "Creating archive: ${OUTFILE}"
echo "Including cell types: ${CELL_TYPES[*]}"
echo "Source directory: ${SOURCE_DIR}"
echo "Archive root name: ${ARCHIVE_ROOT_NAME}"

tar \
  -C "${BASE_RESULTS_DIR}" \
  --transform="s|^${SOURCE_DIR_NAME}|${ARCHIVE_ROOT_NAME}|" \
  -cf - \
  "${tar_args[@]}" \
| pigz -"${COMPRESSION_LEVEL}" -p "${THREADS}" \
> "${OUTFILE}"

# Verify -----
echo "Verifying gzip integrity..."
pigz -t "${OUTFILE}"

echo "Verifying tar structure..."
tar -tzf "${OUTFILE}" > /dev/null

echo "Archive created and verified successfully."
