#!/bin/bash

# Runs genomes -> draft models: CarveMe (draft reconstruction from
# protein_sequences/*.faa and genomic_sequences/*.fna) followed by a bare
# CobraPy conversion (SBML -> clean SBML + .mat).
#
# Moves any Diamond hit tables to diamond_hits/.
#
# Output filenames are the input filename up to the first '.' -- for
# RefSeq-style names (GCF_001421785.1_Leaf33_protein.faa) this keeps only
# the bare accession (GCF_001421785).
#
# Usage: bash scripts/draft_model_generation/genome_to_draftmodels.sh [--cores N] [--solver NAME]
# Example: bash scripts/draft_model_generation/genome_to_draftmodels.sh --cores 4 --solver gurobi

source "$(dirname "$0")/../config.sh"

cores=$(( $(sysctl -n hw.logicalcpu 2>/dev/null || nproc 2>/dev/null || echo 4) - 1 ))
[ "$cores" -lt 1 ] && cores=1
solver="gurobi"

while [ $# -gt 0 ]; do
    case "$1" in
        --cores) cores="$2"; shift 2 ;;
        --solver) solver="$2"; shift 2 ;;
        *) echo "Unknown argument: $1" >&2; exit 1 ;;
    esac
done

mkdir -p "$diamond_hits_dir" "$carveme_dir" "$cobrapy_xml_dir" "$cobrapy_mat_dir"

## --- Step 1: CarveMe ---

run_carve() {
    local flag=$1
    local files=("${@:2}")
    [ ${#files[@]} -eq 0 ] && return

    if [ ${#files[@]} -eq 1 ]; then
        outfile=$(basename "${files[0]}" | cut -d. -f1)
        err=$("$carve_bin" $flag "${files[0]}" -o "$carveme_dir/$outfile.xml" --solver "$solver" 2>&1 >/dev/null) \
            || echo "carve failed on $(basename "${files[0]}"): $err" >&2

    elif command -v parallel > /dev/null 2>&1; then
        printf '%s\n' "${files[@]}" | parallel -j "$cores" --bar \
            'outfile=$(basename {} | cut -d. -f1); err=$("'"$carve_bin"'" '"$flag"' "{}" -o "'"$carveme_dir"'/$outfile.xml" --solver '"$solver"' 2>&1 >/dev/null) || echo "carve failed on $(basename {}): $err" >&2'

    else
        for f in "${files[@]}"; do
            outfile=$(basename "$f" | cut -d. -f1)
            err=$("$carve_bin" $flag "$f" -o "$carveme_dir/$outfile.xml" --solver "$solver" 2>&1 >/dev/null) \
                || echo "carve failed on $(basename "$f"): $err" >&2
        done
    fi
}

protein_files=("$protein_sequences_dir"/*.faa)
[ -e "${protein_files[0]}" ] && run_carve "" "${protein_files[@]}"

dna_files=("$genomic_sequences_dir"/*.fna)
[ -e "${dna_files[0]}" ] && run_carve "--dna" "${dna_files[@]}"

find "$protein_sequences_dir" "$genomic_sequences_dir" -maxdepth 1 -name "*.tsv" -exec mv {} "$diamond_hits_dir/" \; 2>/dev/null

## --- Step 2: CobraPy conversion (bare -- SBML -> clean SBML + .mat, no mismatch detection) ---

tmppy=$(mktemp /tmp/cobrapy_convert.XXXXXX.py)
trap 'rm -f "$tmppy"' EXIT

cat > "$tmppy" << 'PYEOF'
import sys, warnings
warnings.filterwarnings('ignore')
import cobra.io

xml_in, xml_out, mat_out = sys.argv[1], sys.argv[2], sys.argv[3]
model = cobra.io.read_sbml_model(xml_in)
cobra.io.write_sbml_model(model, xml_out)
cobra.io.save_matlab_model(model, mat_out, varname='model')
PYEOF

draft_files=("$carveme_dir"/*.xml)
if [ -e "${draft_files[0]}" ]; then
    if [ ${#draft_files[@]} -eq 1 ]; then
        name=$(basename "${draft_files[0]}" .xml)
        err=$("$python_env" -W ignore "$tmppy" "${draft_files[0]}" "$cobrapy_xml_dir/$name.xml" "$cobrapy_mat_dir/$name.mat" 2>&1 >/dev/null) \
            || echo "cobrapy conversion failed on $name: $err" >&2

    elif command -v parallel > /dev/null 2>&1; then
        printf '%s\n' "${draft_files[@]}" | parallel -j "$cores" --bar \
            'name=$(basename {} .xml); err=$("'"$python_env"'" -W ignore "'"$tmppy"'" {} "'"$cobrapy_xml_dir"'/$name.xml" "'"$cobrapy_mat_dir"'/$name.mat" 2>&1 >/dev/null) || echo "cobrapy conversion failed on $name: $err" >&2'

    else
        for f in "${draft_files[@]}"; do
            name=$(basename "$f" .xml)
            err=$("$python_env" -W ignore "$tmppy" "$f" "$cobrapy_xml_dir/$name.xml" "$cobrapy_mat_dir/$name.mat" 2>&1 >/dev/null) \
                || echo "cobrapy conversion failed on $name: $err" >&2
        done
    fi
fi
