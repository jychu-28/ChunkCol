#!/usr/bin/env bash
# ============================================================================
# Genome collinearity comparison -> SyRI structural variation detection
# -> plotsr visualization.
#
# Input:
#   Chromosome-level genome FASTA files and .fai indexes.
#
# Output:
#   syri_full/<comparison>/<comparison>.png
# ============================================================================

set -euo pipefail

# ============================================
# 0. Configuration
# ============================================

MINIMAP2=${MINIMAP2:-/home/jychu/miniconda3/bin/minimap2}
SYRI=${SYRI:-/tmp/syri_pandas2/bin/syri}
PLOTSR=${PLOTSR:-/home/jychu/miniconda3/envs/syri/bin/plotsr}
WORKDIR=${WORKDIR:-$(pwd)}

CHUNK_SIZE=${CHUNK_SIZE:-100000000}
OVERLAP=${OVERLAP:-10000}
THREADS=${THREADS:-10}
PARALLEL_CHUNKS=${PARALLEL_CHUNKS:-7}

GENOME_LUDLOWII=${GENOME_LUDLOWII:-"$WORKDIR/Paeonia.ludlowii.fna.chr1-5.fa"}
GENOME_POLISHED=${GENOME_POLISHED:-"$WORKDIR/final.polish.merged.all.fasta.chr1-5.fa"}
GENOME_FENGDAN=${GENOME_FENGDAN:-"$WORKDIR/Paeoniaostii.Fengdan.genome.fa.chr1-5.fa"}

CHR_LUDLOWII=(Chr1 Chr2 Chr3 Chr4 Chr5)
CHR_POLISHED=(chr01 chr02 chr03 chr04 chr05)
CHR_FENGDAN=(Chr01 Chr02 Chr03 Chr04 Chr05)

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

check_file() {
    local path=$1
    if [[ ! -s "$path" ]]; then
        echo "Missing required file: $path" >&2
        exit 1
    fi
}

check_command() {
    local path=$1
    local name=$2
    if [[ ! -x "$path" ]]; then
        echo "Missing executable $name: $path" >&2
        exit 1
    fi
}

# ============================================
# 1. Extract chromosome FASTA and split chunks
# ============================================

extract_and_split() {
    local genome_fa=$1
    local chr_name=$2
    local out_prefix=$3

    local chr_fa="${out_prefix}.fa"
    if [[ ! -s "$chr_fa" ]]; then
        echo "  Extracting $chr_name from $(basename "$genome_fa") ..."
        python3 - "$genome_fa" "$chr_name" "$chr_fa" <<'PY'
import sys
from Bio import SeqIO

genome_fa, chr_name, chr_fa = sys.argv[1:4]
for rec in SeqIO.parse(genome_fa, "fasta"):
    if rec.id == chr_name:
        SeqIO.write(rec, chr_fa, "fasta")
        print(f"  -> {rec.id}: {len(rec.seq):,} bp")
        break
else:
    raise SystemExit(f"chromosome not found: {chr_name}")
PY
    fi

    local first_chunk="${out_prefix}_chunk0000.fa"
    if [[ -s "$first_chunk" ]]; then
        echo "  [SKIP] chunks exist"
        return
    fi

    echo "  Splitting into ${CHUNK_SIZE} bp chunks ..."
    python3 - "$chr_fa" "$out_prefix" "$CHUNK_SIZE" "$OVERLAP" <<'PY'
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

chr_fa, out_prefix, chunk_size, overlap = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
for rec in SeqIO.parse(chr_fa, "fasta"):
    seq = str(rec.seq)
    n_chunks = 0
    for i, start in enumerate(range(0, len(seq), chunk_size)):
        end = min(start + chunk_size + overlap, len(seq))
        chunk = SeqRecord(
            Seq(seq[start:end]),
            id=f"{rec.id}_chunk{i:04d}",
            description=f"start={start}",
        )
        with open(f"{out_prefix}_chunk{i:04d}.fa", "w") as handle:
            SeqIO.write(chunk, handle, "fasta")
        n_chunks += 1
    print(f"  -> {rec.id}: {n_chunks} chunks")
PY
}

# ============================================
# 2. minimap2 chunk alignment
# ============================================

run_minimap2_chunks() {
    local ref_prefix=$1
    local qry_prefix=$2
    local out_dir=$3

    mkdir -p "$out_dir"

    local n_ref
    local n_qry
    n_ref=$(find "$(dirname "$ref_prefix")" -maxdepth 1 -name "$(basename "$ref_prefix")_chunk*.fa" | wc -l)
    n_qry=$(find "$(dirname "$qry_prefix")" -maxdepth 1 -name "$(basename "$qry_prefix")_chunk*.fa" | wc -l)

    local taskfile="$out_dir/tasks.txt"
    : > "$taskfile"
    for ((r = 0; r < n_ref; r++)); do
        for ((q = 0; q < n_qry; q++)); do
            local expected=$((r * n_qry / n_ref))
            if [[ $q -ge $((expected - 1)) && $q -le $((expected + 1)) ]]; then
                printf -v r_pad "%04d" "$r"
                printf -v q_pad "%04d" "$q"
                printf "%s\t%s\t%s\n" \
                    "${ref_prefix}_chunk${r_pad}.fa" \
                    "${qry_prefix}_chunk${q_pad}.fa" \
                    "${out_dir}/ref${r}_qry${q}.paf" >> "$taskfile"
            fi
        done
    done

    echo "  Aligning $(wc -l < "$taskfile") chunk pairs ..."
    parallel --colsep '\t' -j "$PARALLEL_CHUNKS" \
        "$MINIMAP2 -x asm20 --cs --eqx -t 2 {1} {2} > {3} 2>/dev/null" \
        :::: "$taskfile"

    cat "$out_dir"/ref*_qry*.paf > "$out_dir/merged.paf"
    echo "  -> merged.paf: $(wc -l < "$out_dir/merged.paf") lines"
}

# ============================================
# 3. Restore chromosome-level PAF coordinates
# ============================================

fix_coords() {
    local merged_paf=$1
    local ref_fai=$2
    local qry_fai=$3
    local fixed_paf=$4

    if [[ -s "$fixed_paf" ]]; then
        echo "  [SKIP] $(basename "$fixed_paf") exists"
        return
    fi

    echo "  Fixing coordinates ..."
    python3 "$SCRIPT_DIR/fix_paf_coords.py" \
        --input "$merged_paf" \
        --output "$fixed_paf" \
        --fai "$ref_fai" \
        --fai "$qry_fai" \
        --chunk-size "$CHUNK_SIZE"
}

# ============================================
# 4. Convert cs:Z: tags to cg:Z: tags
# ============================================

add_cigar() {
    local fixed_paf=$1
    local cg_paf=$2

    if [[ -s "$cg_paf" && "$(wc -l < "$cg_paf")" -eq "$(wc -l < "$fixed_paf")" ]]; then
        echo "  [SKIP] $(basename "$cg_paf") exists"
        return
    fi

    echo "  Adding cg:Z: tags ..."
    python3 "$SCRIPT_DIR/add_cigar_to_paf.py" \
        --input "$fixed_paf" \
        --output "$cg_paf"
}

# ============================================
# 5. SyRI structural variation detection
# ============================================

run_syri() {
    local cg_paf=$1
    local ref_fa=$2
    local qry_fa=$3
    local out_dir=$4
    local prefix=$5
    local syri_out="$out_dir/${prefix}syri.out"

    if [[ -s "$syri_out" ]]; then
        echo "  [SKIP] $syri_out exists ($(du -h "$syri_out" | cut -f1))"
        return
    fi

    echo "  Running SyRI: $prefix ($(wc -l < "$cg_paf") alignments) ..."
    "$SYRI" -c "$cg_paf" -r "$ref_fa" -q "$qry_fa" \
        -F P --nc 1 --nosnp \
        --prefix "$prefix" --dir "$out_dir"

    if [[ -s "$syri_out" ]]; then
        echo "  -> $(basename "$syri_out"): $(du -h "$syri_out" | cut -f1), $(wc -l < "$syri_out") lines"
    else
        echo "  FAILED: $prefix" >&2
        return 1
    fi
}

# ============================================
# 6. plotsr visualization
# ============================================

run_plotsr() {
    local out_dir=$1
    local plot_name=$2
    local ref_genome=$3
    local qry_genome=$4
    local ref_label=$5
    local qry_label=$6
    shift 6

    local merged_out="$out_dir/${plot_name}.syri.out"
    cat "$@" > "$merged_out"
    echo "  Merged syri.out: $(wc -l < "$merged_out") lines"

    local genomes_txt="$out_dir/genomes.txt"
    cat > "$genomes_txt" <<EOF
#file	name	tags
$ref_genome	$ref_label	lw:1.5
$qry_genome	$qry_label	lw:1.5
EOF

    echo "  Running plotsr ..."
    "$PLOTSR" \
        --sr "$merged_out" \
        --genomes "$genomes_txt" \
        -o "$out_dir/${plot_name}.png" \
        -H 10 -W 14 -f 10

    echo "  -> $out_dir/${plot_name}.png ($(du -h "$out_dir/${plot_name}.png" | cut -f1))"
}

# ============================================
# Process one genome comparison
# ============================================

run_comparison() {
    local cmp_name=$1
    local ref_genome_fa=$2
    local qry_genome_fa=$3
    local -n ref_chrs=$4
    local -n qry_chrs=$5
    local ref_label=$6
    local qry_label=$7

    local ref_fai="${ref_genome_fa}.fai"
    local qry_fai="${qry_genome_fa}.fai"

    check_file "$ref_genome_fa"
    check_file "$qry_genome_fa"
    check_file "$ref_fai"
    check_file "$qry_fai"

    echo ""
    echo "=============================================="
    echo "  $cmp_name"
    echo "=============================================="

    local split_dir="$WORKDIR/split_100M/$cmp_name"
    local paf_dir="$WORKDIR/01_paf/$cmp_name"
    local syri_dir="$WORKDIR/syri_full/$cmp_name"
    mkdir -p "$split_dir" "$paf_dir" "$syri_dir"

    local syri_outputs=()

    for i in 0 1 2 3 4; do
        local ref_chr="${ref_chrs[$i]}"
        local qry_chr="${qry_chrs[$i]}"
        local pair="${ref_chr}_vs_${qry_chr}"
        local pair_dir="$paf_dir/$pair"

        echo ""
        echo "--- $pair ---"

        extract_and_split "$ref_genome_fa" "$ref_chr" "$split_dir/${ref_chr}"
        extract_and_split "$qry_genome_fa" "$qry_chr" "$split_dir/${qry_chr}"

        local merged_paf="$pair_dir/merged.paf"
        if [[ ! -s "$merged_paf" ]]; then
            mkdir -p "$pair_dir"
            run_minimap2_chunks "$split_dir/${ref_chr}" "$split_dir/${qry_chr}" "$pair_dir"
        else
            echo "  [SKIP] minimap2: $pair/merged.paf exists ($(wc -l < "$merged_paf") lines)"
        fi

        local fixed_paf="$pair_dir/merged.fixed.paf"
        fix_coords "$merged_paf" "$ref_fai" "$qry_fai" "$fixed_paf"

        local cg_paf="$pair_dir/merged.fixed.cg.paf"
        add_cigar "$fixed_paf" "$cg_paf"

        local ref_chr_fa="$split_dir/${ref_chr}.fa"
        local qry_chr_fa="$split_dir/${qry_chr}.fa"
        run_syri "$cg_paf" "$ref_chr_fa" "$qry_chr_fa" "$syri_dir" "$pair"

        syri_outputs+=("$syri_dir/${pair}syri.out")
    done

    echo ""
    echo "=== Plotting $cmp_name ==="
    run_plotsr "$syri_dir" "$cmp_name" \
        "$ref_genome_fa" "$qry_genome_fa" \
        "$ref_label" "$qry_label" \
        "${syri_outputs[@]}"

    echo ""
    echo "=== $cmp_name DONE ==="
}

# ============================================
# Main
# ============================================

check_command "$MINIMAP2" "minimap2"
check_command "$SYRI" "syri"
check_command "$PLOTSR" "plotsr"

echo "Pipeline start: $(date)"
echo "Minimap2: $("$MINIMAP2" --version 2>&1 | head -1)"
echo "Chunk size: ${CHUNK_SIZE} bp"
echo "Overlap: ${OVERLAP} bp"
echo ""

run_comparison \
    "ludlowii_vs_polished" \
    "$GENOME_LUDLOWII" "$GENOME_POLISHED" \
    CHR_LUDLOWII CHR_POLISHED \
    "ludlowii" "polished"

run_comparison \
    "polished_vs_fengdan" \
    "$GENOME_POLISHED" "$GENOME_FENGDAN" \
    CHR_POLISHED CHR_FENGDAN \
    "polished" "fengdan"

echo ""
echo "=============================================="
echo "  ALL DONE at $(date)"
echo "=============================================="
echo "Output:"
echo "  syri_full/ludlowii_vs_polished/ludlowii_vs_polished.png"
echo "  syri_full/polished_vs_fengdan/polished_vs_fengdan.png"

