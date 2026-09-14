#!/bin/bash

# Runs MEMOTE quality verification on the final SBML models, producing one
# HTML snapshot report per organism.
#
# Usage: bash scripts/verification/memote_verification.sh [--cores N]
# Example: bash scripts/verification/memote_verification.sh --cores 4

source "$(dirname "$0")/../config.sh"

cores=$(( $(sysctl -n hw.logicalcpu 2>/dev/null || nproc 2>/dev/null || echo 4) - 1 ))
[ "$cores" -lt 1 ] && cores=1

while [ $# -gt 0 ]; do
    case "$1" in
        --cores) cores="$2"; shift 2 ;;
        *) echo "Unknown argument: $1" >&2; exit 1 ;;
    esac
done

sbml_dir="$pipeline_root/data/models/final/xml"
reports_dir="$pipeline_root/data/models/final/reports"
mkdir -p "$reports_dir"

run_memote() {
    local f="$1"
    local name
    name=$(basename "$f" .xml)
    DISABLE_PANDERA_IMPORT_WARNING=True "$memote_bin" report snapshot \
        --filename "$reports_dir/$name.html" "$f" \
        > /dev/null 2>&1 \
        || echo "memote failed on $name" >&2
}
export -f run_memote
export memote_bin reports_dir

files=("$sbml_dir"/*.xml)
if [ ! -e "${files[0]}" ]; then
    echo "No .xml files found in $sbml_dir -- nothing to process." >&2
    exit 1
fi

if [ ${#files[@]} -eq 1 ]; then
    run_memote "${files[0]}"
elif command -v parallel > /dev/null 2>&1; then
    printf '%s\n' "${files[@]}" | parallel -j "$cores" --bar run_memote {}
else
    for f in "${files[@]}"; do
        run_memote "$f"
    done
fi

echo "Reports generated: $(ls "$reports_dir"/*.html 2>/dev/null | wc -l | tr -d ' ')"
