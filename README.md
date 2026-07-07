# ChunkCol: Chunk-based collinearity plotting for ultra-large genomes

A chunk-based strategy designed **specifically for ultra-large genomes, using peony (*Paeonia*) as an example**, especially when a single chromosome can exceed **2 Gb** and direct whole-chromosome alignment becomes difficult to run efficiently.

The core idea is to split each chromosome into 100 Mb overlapping chunks, align only neighboring chunk pairs, and then restore chunk-level PAF coordinates back to whole-chromosome coordinates before running **SyRI** and **plotsr**. This keeps the alignment workload manageable while still producing chromosome-scale synteny and structural variation visualizations.

This project was developed for chromosome-level comparative analysis among *Paeonia* genomes. It uses **minimap2** for chunk alignment, converts minimap2 `cs:Z:` tags into `cg:Z:` CIGAR tags required by SyRI, detects structural rearrangements with **SyRI**, and finally visualizes genome collinearity with **plotsr**.

## Workflow

```text
Whole-genome FASTA files (chr1-5.fa)
   |
   |-- [1] Extract chromosomes and split into 100 Mb chunks
   |
   |-- [2] Align adjacent chunks with minimap2 -x asm20
   |        -> ref*_qry*.paf
   |        -> merged.paf
   |
   |-- [3] Restore chunk-level coordinates to chromosome-level coordinates
   |        -> merged.fixed.paf
   |
   |-- [4] Convert cs:Z: tags to cg:Z: tags
   |        -> merged.fixed.cg.paf
   |
   |-- [5] Run SyRI with PAF input
   |        -> syri.out, syri.summary, syri.vcf
   |
   `-- [6] Plot collinearity with plotsr
            -> PNG synteny plot
```

## Repository Structure

```text
.
├── README.md
├── docs/
│   └── workflow.md
└── scripts/
    ├── pipeline_full.sh
    ├── fix_paf_coords.py
    └── add_cigar_to_paf.py
```

## Dependencies

Required command-line tools:

- minimap2
- SyRI
- plotsr
- GNU parallel
- samtools, for creating FASTA `.fai` indexes if needed

Required Python packages:

- biopython

Example installation with conda:

```bash
conda create -n syri-pipeline -c bioconda -c conda-forge \
  minimap2 syri plotsr parallel samtools biopython
conda activate syri-pipeline
```

Depending on your SyRI installation, you may need to use a separate environment or a specific pandas version.

## Input Files

The pipeline expects chromosome-level FASTA files and corresponding `.fai` indexes:

```text
Paeonia.ludlowii.fna.chr1-5.fa
Paeonia.ludlowii.fna.chr1-5.fa.fai

final.polish.merged.all.fasta.chr1-5.fa
final.polish.merged.all.fasta.chr1-5.fa.fai

Paeoniaostii.Fengdan.genome.fa.chr1-5.fa
Paeoniaostii.Fengdan.genome.fa.chr1-5.fa.fai
```

If indexes are missing, create them with:

```bash
samtools faidx genome.fa
```

## Quick Start

Edit the configuration section in `scripts/pipeline_full.sh`:

```bash
MINIMAP2=/path/to/minimap2
SYRI=/path/to/syri
PLOTSR=/path/to/plotsr

GENOME_LUDLOWII="$WORKDIR/Paeonia.ludlowii.fna.chr1-5.fa"
GENOME_POLISHED="$WORKDIR/final.polish.merged.all.fasta.chr1-5.fa"
GENOME_FENGDAN="$WORKDIR/Paeoniaostii.Fengdan.genome.fa.chr1-5.fa"
```

Then run:

```bash
bash scripts/pipeline_full.sh
```

## Output

The main outputs are:

```text
syri_full/ludlowii_vs_polished/ludlowii_vs_polished.png
syri_full/polished_vs_fengdan/polished_vs_fengdan.png
```

Intermediate outputs:

```text
split_100M/        chromosome FASTA chunks
01_paf/            chunk-level and corrected PAF files
syri_full/         SyRI results and plotsr figures
```

## Notes

- The default chunk size is 100 Mb with 10 kb overlap.
- Only neighboring chunks are aligned to reduce computational cost.
- `fix_paf_coords.py` converts chunk-local coordinates back to whole-chromosome coordinates.
- `add_cigar_to_paf.py` converts minimap2 `cs:Z:` tags to `cg:Z:` tags because SyRI requires CIGAR information for PAF input.
- Large FASTA, PAF, SyRI, and plot output files are intentionally ignored by `.gitignore`.

## Citation

If you use this workflow, please cite the underlying tools:

- minimap2
- SyRI
- plotsr



