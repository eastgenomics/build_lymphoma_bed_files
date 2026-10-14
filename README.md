# build_lymphoma_bed_files

## What are typical use cases for this script?

Building the athena exons file and panel BED for the lymphoma (gaea) assay
from approved inputs, using fixed rules plus a reviewed overrides file.

---

## What inputs are required for this script to run?

### Required

- `--tsv` (`file`): approved coding exons +-`--pad` bp
  (`GCF_000001405.39_GRCh38.p13_genomic_20211119.symbols.exon_5bp.tsv`)
- `--gff` (`file`): gzipped RefSeq GFF the TSV was made from
  (`GCF_000001405.39_GRCh38.p13_genomic_20211119.gff.gz`)
- `--probes` (`file`): probe BED (`Probes_GRCh38_Lymphoma_v1.0.bed`)
- `--panel-transcripts` (`file`): transcripts in the panel, one per line
  (`panel_transcripts.txt`)
- `--utr3-transcripts` (`file`): gene<TAB>transcript pairs whose 3'UTR goes in
  the panel (`utr3_transcripts.tsv`)
- `--overrides` (`file`): overrides TSV with a header row: action, chrom, start,
  end, gene, transcript, label, reason (`overrides.tsv`)
- `--curated-transcripts` (`file`): gene<TAB>transcript<TAB>requested rows whose
  targeted non-coding sequence goes in the exons file (`curated_transcripts.tsv`)
- `--pad` (`int`): padding (bp) either side of the CDS in the `--tsv` rows (`5`)
- `--exons-out` (`str`): exons file to write
- `--panel-out` (`str`): panel BED to write

_n.b. none of the input files are kept in git: put them in `data/`, which is
ignored, along with `output/`._

## How does this script work?

The exons file and panel are built from four rules:

- **Panel** = coding exons (with replaces) + 3'UTRs + override adds, for the
  transcripts in `--panel-transcripts`, minus excludes
- **Exons file** = coding exons (with replaces) + 3'UTRs + override adds +
  targeted non-coding sequence

### Coding exons

- The approved `--tsv` rows, unchanged

### 3'UTRs

- For the transcripts in `--utr3-transcripts`
- Taken from the RefSeq GFF, starting at the padded CDS edge, trimmed to
  probe-covered sequence
- Labelled `utr3`, or `utr3a`/`utr3b`/... if the trim splits it

### Overrides

For anything the rules can't produce. Each override must have a `reason`.

- `add`: add a row to the exons file (and so the panel); must be fully
  probe-covered
- `replace`: swap the TSV row with the same transcript and label for these
  coordinates, in the exons file and the panel; must contain the TSV row, and
  the extra sequence must be fully probe-covered
  - _n.b. for exceptions only - currently MGA exon 14, extended to the MANE
    NM_001400225.1 exon start so it stays one row_
- `exclude`: drop a row from the panel; must have no probe coverage

### Targeted non-coding sequence

- For the transcripts in `--curated-transcripts` (the genes manually curated in
  v1, DI-1335), derived from the probes
- Keeps probe-covered sequence around the transcript that isn't already covered
  by the rules above, if it's longer than the longest probe overhang past a
  coding exon on any other panel transcript (measured from the data, ~110bp)
- Also keeps any intron the probes cover completely
- Labelled from the GFF as `flank5`/`flank3`, `utr5_N`/`utr3_N` (UTR part of
  exon N) or `intron_N`
- Exons file only, never the panel

## Setup

Requires Python 3.9+ and `bedtools` on PATH.

```
python3 -m venv .venv
source .venv/bin/activate
pip install -e .
```

Install bedtools (e.g. `sudo apt install bedtools`),
or download a release binary and add its folder to PATH in `~/.bashrc`:

```
export PATH="$HOME/tools/bedtools2/bin:$PATH"
```

Check with `bedtools --version`.

---

## Example commands

### Build

```
python3 build_gaea_beds.py \
    --tsv data/GCF_000001405.39_GRCh38.p13_genomic_20211119.symbols.exon_5bp.tsv \
    --gff data/GCF_000001405.39_GRCh38.p13_genomic_20211119.gff.gz \
    --probes data/Probes_GRCh38_Lymphoma_v1.0.bed \
    --panel-transcripts data/panel_transcripts.txt \
    --utr3-transcripts data/utr3_transcripts.tsv \
    --overrides data/overrides.tsv \
    --curated-transcripts data/curated_transcripts.tsv \
    --pad 5 \
    --exons-out output/GCF_000001405.39_GRCh38.p13_genomic_20211119.gaea_exons_5bp_athena_v2.0.0.bed \
    --panel-out output/GCF_000001405.39_GRCh38.p13_genomic_20211119.gaea_panel_5bp_athena_coding_3utr_v2.0.0.bed
```

## What does this script output

- `--exons-out`: the athena exons file, headerless tab-separated (chrom, start,
  end, gene, transcript, exon), sorted by position
  (`GCF_000001405.39_GRCh38.p13_genomic_20211119.gaea_exons_5bp_athena_v2.0.0.bed`)
- `--panel-out`: the panel BED, headerless tab-separated (chrom, start, end,
  transcript), sorted by position
  (`GCF_000001405.39_GRCh38.p13_genomic_20211119.gaea_panel_5bp_athena_coding_3utr_v2.0.0.bed`)

The script also prints the probe overhang cut-off it used and the row count of
each file.
