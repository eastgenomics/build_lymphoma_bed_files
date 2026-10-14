#!/usr/bin/env python3
"""
Build the athena exons file and panel BED for the lymphoma (gaea) assay.

Rules
  1. Coding exons: the approved symbols.exon_5bp.tsv, unchanged.
  2. 3'UTRs for the transcripts in utr3_transcripts.tsv: taken from the
     RefSeq GFF, starting at the padded CDS edge, trimmed to probe-covered
     sequence. Labelled utr3, or utr3a/utr3b/... if the trim splits it.
  3. Overrides (overrides.tsv) for anything the rules can't produce:
       add      add a row to the exons file (and so the panel)
       replace  swap the TSV row with the same transcript and label for
                these coordinates, in the exons file and the panel
       exclude  drop a row from the panel
  4. Targeted non-coding sequence for the transcripts in
     curated_transcripts.tsv, derived from the probes: probe-covered
     sequence around the transcript that isn't already covered by 1-3,
     keeping pieces longer than the longest probe overhang past a coding
     exon on any other panel transcript, plus any intron the probes cover
     completely. Labelled from the GFF as flank5/flank3, utr5_N/utr3_N (UTR
     part of exon N) or intron_N. Exons file only, never the panel.
Panel      = 1 (with replaces) + 2 + override adds, for panel_transcripts.txt,
             minus excludes.
Exons file = 1 (with replaces) + 2 + override adds + 4.

Requires bedtools on PATH. See README.md for an example command.
"""
import argparse
import gzip
import subprocess
import tempfile
from io import StringIO
from pathlib import Path

import pandas as pd

COLS = ["chrom", "start", "end", "gene", "transcript", "exon"]
CHROM_ORDER = {**{str(i): i for i in range(1, 23)}, "X": 23, "Y": 24, "MT": 25}


def bedtools(*args):
    """
    Run a bedtools command and return its output as a DataFrame.

    Parameters
    ----------
    *args : str, int, Path or pandas.DataFrame
        Arguments to pass to ``bedtools``, e.g. ``"intersect", "-a", df1,
        "-b", df2``. DataFrames are written to temporary headerless BED
        files and replaced by their paths.

    Returns
    -------
    pandas.DataFrame
        The command's output, with integer column names and column 0
        (chromosome) read as str. Empty if the command printed nothing.

    Raises
    ------
    subprocess.CalledProcessError
        If bedtools exits with a non-zero status.
    """
    with tempfile.TemporaryDirectory() as tmp:
        cmd = ["bedtools"]
        for i, arg in enumerate(args):
            if isinstance(arg, pd.DataFrame):
                path = f"{tmp}/{i}.bed"
                arg.to_csv(path, sep="\t", header=False, index=False)
                arg = path
            cmd.append(str(arg))
        out = subprocess.run(
            cmd, check=True, capture_output=True, text=True
        ).stdout
    if not out:
        return pd.DataFrame()
    return pd.read_csv(StringIO(out), sep="\t", header=None, dtype={0: str})


def read_bed(path, names):
    """
    Read a headerless tab-separated file, skipping ``#`` comment lines.

    Parameters
    ----------
    path : str or Path
        File to read.
    names : list of str
        Column names, one per column in the file.

    Returns
    -------
    pandas.DataFrame
        The file's rows, with ``chrom`` and ``exon`` columns (if present)
        read as str.
    """
    return pd.read_csv(path, sep="\t", header=None, names=names, comment="#",
                       dtype={"chrom": str, "exon": str})


def anti_join(df, other, on):
    """
    Keep the rows of ``df`` that have no match in ``other``.

    Parameters
    ----------
    df, other : pandas.DataFrame
        Tables to compare.
    on : list of str
        Columns to match on.

    Returns
    -------
    pandas.DataFrame
        Rows of ``df`` whose ``on`` values aren't in ``other``.
    """
    merged = df.merge(
        other[on].drop_duplicates(), on=on, how="left", indicator=True
    )
    return merged[merged._merge == "left_only"].drop(columns="_merge")


def read_overrides(path):
    """
    Read the overrides file.

    Parameters
    ----------
    path : str or Path
        Overrides TSV with a header row: action, chrom, start, end, gene,
        transcript, label, reason.

    Returns
    -------
    pandas.DataFrame
        One row per override, with ``label`` renamed to ``exon`` to match
        the exons-file columns (``COLS``).
    """
    overrides = pd.read_csv(path, sep="\t", dtype={"chrom": str, "label": str})
    return overrides.rename(columns={"label": "exon"})


def apply_replaces(coding, replaces):
    """
    Swap TSV rows for the coordinates given by ``replace`` overrides.

    Each replace row matches the TSV row with the same transcript and exon
    label, and must contain it (a replace can only extend a row).

    Parameters
    ----------
    coding : pandas.DataFrame
        Coding exon rows from the TSV, in ``COLS`` format.
    replaces : pandas.DataFrame
        ``replace`` override rows, in ``COLS`` format.

    Returns
    -------
    pandas.DataFrame
        ``coding`` with each matched row's chrom, start and end replaced.

    Raises
    ------
    ValueError
        If a replace row doesn't match exactly one TSV row, or doesn't
        contain the row it matches.
    """
    coding = coding.copy()
    for r in replaces.itertuples(index=False):
        hit = (coding.transcript == r.transcript) & (coding.exon == r.exon)
        if hit.sum() != 1:
            raise ValueError(
                f"replace {r.transcript} exon {r.exon} "
                f"matches {hit.sum()} TSV rows, not 1"
            )
        old = coding[hit].iloc[0]
        if not (
            old.chrom == r.chrom and r.start <= old.start and r.end >= old.end
        ):
            raise ValueError(
                f"replace {r.transcript} exon {r.exon} "
                f"doesn't contain the TSV row "
                f"{old.chrom}:{old.start}-{old.end}"
            )
        coding.loc[hit, ["chrom", "start", "end"]] = [r.chrom, r.start, r.end]
    return coding


def merged_probes(path):
    """
    Read the probe BED and merge it into contiguous targeted blocks.

    Abutting probes in the design are often separated by a 1bp seam, so
    probes up to 1bp apart are merged (``bedtools merge -d 1``) rather
    than leaving the seam as untargeted sequence.

    Parameters
    ----------
    path : str or Path
        Probe BED (chrom, start, end).

    Returns
    -------
    pandas.DataFrame
        Merged probe intervals, columns 0-2 (chrom, start, end).
    """
    probes = read_bed(path, ["chrom", "start", "end"]).sort_values(
        ["chrom", "start"]
    )
    return bedtools("merge", "-d", 1, "-i", probes)


def gff_features(path, transcripts):
    """
    Extract exon and CDS features for the given transcripts from the GFF.

    Only primary chromosomes (``NC_0000*`` accessions) are read, so PAR /
    alt-scaffold copies are ignored. Coordinates are converted from GFF
    (1-based, closed) to BED (0-based, half-open), and accessions to the
    chromosome names used in the TSV and probes (e.g. ``NC_000016.10`` ->
    ``16``, ``NC_000023.11`` -> ``X``).

    Parameters
    ----------
    path : str or Path
        Gzipped RefSeq GFF.
    transcripts : iterable of str
        Transcript IDs, e.g. ``NM_003745.2``.

    Returns
    -------
    pandas.DataFrame
        Columns chrom, start, end, feature (``exon`` or ``CDS``),
        transcript, strand.

    Raises
    ------
    ValueError
        If any transcript has no features in the GFF.
    """
    keys = {f"Parent=rna-{tx};": tx for tx in transcripts}
    rows = []
    with gzip.open(path, "rt") as fh:
        for line in fh:
            f = line.split("\t")
            if (
                len(f) < 9
                or f[2] not in ("exon", "CDS")
                or not f[0].startswith("NC_0000")
            ):
                continue
            tx = next((tx for key, tx in keys.items() if key in f[8]), None)
            if tx:
                # NC_000016.10 -> 16, NC_000023.11 -> X;
                # GFF 1-based -> BED 0-based
                accession = str(int(f[0][3:9]))
                chrom = {"23": "X", "24": "Y"}.get(accession, accession)
                rows.append((chrom, int(f[3]) - 1, int(f[4]), f[2], tx, f[6]))
    feats = pd.DataFrame(
        rows,
        columns=["chrom", "start", "end", "feature", "transcript", "strand"],
    )
    missing = set(transcripts) - set(feats.transcript)
    if missing:
        raise ValueError(f"not found in the GFF: {sorted(missing)}")
    return feats


def transcript_parts(feats, tx):
    """
    Unpack one transcript's features from ``gff_features``.

    Parameters
    ----------
    feats : pandas.DataFrame
        Features from ``gff_features`` including this transcript.
    tx : str
        Transcript ID.

    Returns
    -------
    chrom : str
        Chromosome.
    strand : {"+", "-"}
        Strand.
    cds_lo, cds_hi : int
        Lowest CDS start and highest CDS end (unpadded).
    exons : list of tuple of int
        Exon (start, end) intervals, sorted by position.
    """
    f = feats[feats.transcript == tx]
    cds = f[f.feature == "CDS"]
    exons = [
        tuple(e)
        for e in f[f.feature == "exon"]
        .sort_values("start")[["start", "end"]]
        .values.tolist()
    ]
    return (
        f.chrom.iloc[0],
        f.strand.iloc[0],
        cds.start.min(),
        cds.end.max(),
        exons,
    )


def labelled_rows(intervals, gene, tx, labels):
    """
    Make exons-file rows, adding ``a``, ``b``, ... to repeated labels.

    Parameters
    ----------
    intervals : pandas.DataFrame
        Intervals in transcript order, columns 0-2 (chrom, start, end).
    gene : str
        Gene symbol.
    tx : str
        Transcript ID.
    labels : list of str
        One label per interval. A label used more than once gets ``a``,
        ``b``, ... suffixes in order, e.g. ``utr3a``, ``utr3b``.

    Returns
    -------
    pandas.DataFrame
        Rows in ``COLS`` format.
    """
    counts = pd.Series(labels).value_counts()
    seen = {}
    exon = []
    for label in labels:
        seen[label] = seen.get(label, -1) + 1
        exon.append(
            label
            if counts[label] == 1
            else f"{label}{chr(ord('a') + seen[label])}"
        )
    return pd.DataFrame(
        {
            "chrom": intervals[0].values,
            "start": intervals[1].values,
            "end": intervals[2].values,
            "gene": gene,
            "transcript": tx,
            "exon": exon,
        }
    )


def utr3_rows(gene, tx, feats, probes, pad):
    """
    Build the panel 3'UTR rows for one transcript (rule 2).

    The 3'UTR is the exon sequence on the 3' side of the CDS, starting at
    the padded CDS edge (CDS end + ``pad`` on the + strand, CDS start -
    ``pad`` on the - strand) so it abuts the last coding row without
    overlapping it. Only probe-covered sequence is kept.

    Parameters
    ----------
    gene : str
        Gene symbol for the output rows.
    tx : str
        Transcript ID.
    feats : pandas.DataFrame
        Features from ``gff_features`` including this transcript.
    probes : pandas.DataFrame
        Merged probe intervals from ``merged_probes``.
    pad : int
        Padding (bp) either side of the CDS in the TSV coding rows.

    Returns
    -------
    pandas.DataFrame
        Rows in ``COLS`` format, in transcript order, labelled ``utr3``,
        or ``utr3a``, ``utr3b``, ... if probe gaps split the 3'UTR.

    Raises
    ------
    ValueError
        If none of the 3'UTR is probe-covered.
    """
    chrom, strand, cds_lo, cds_hi, exons = transcript_parts(feats, tx)
    lo, hi = cds_lo - pad, cds_hi + pad  # padded CDS envelope

    # exon sequence outside the padded CDS, on the 3' side only
    exons = pd.DataFrame([(chrom, s, e) for s, e in exons])
    utr = bedtools(
        "subtract", "-a", exons, "-b", pd.DataFrame([[chrom, lo, hi]])
    )
    utr = utr[utr[1] >= hi] if strand == "+" else utr[utr[2] <= lo]

    covered = bedtools("intersect", "-a", utr, "-b", probes)
    if covered.empty:
        raise ValueError(f"{gene} {tx}: 3'UTR has no probe coverage")
    covered = covered.sort_values(1, ascending=(strand == "+"))
    return labelled_rows(covered, gene, tx, ["utr3"] * len(covered))


def overhang_limit(tsv, transcripts, probes):
    """
    Measure how far probes extend past coding exons by design.

    The capture design targets coding exons plus a small margin, so probe
    blocks touching a coding row usually run a little way past it. The
    longest such overhang on ordinary (non-curated) transcripts is the
    cut-off above which non-coding probe coverage counts as deliberately
    targeted.

    Parameters
    ----------
    tsv : pandas.DataFrame
        Coding exon rows from the TSV, in ``COLS`` format.
    transcripts : iterable of str
        Transcripts to measure, i.e. the panel transcripts that weren't
        manually curated.
    probes : pandas.DataFrame
        Merged probe intervals from ``merged_probes``.

    Returns
    -------
    int
        Longest stretch (bp) of probe-covered sequence outside any TSV
        coding row, in probe blocks touching these transcripts' coding rows.
    """
    coding = tsv[tsv.transcript.isin(transcripts)][["chrom", "start", "end"]]
    blocks = bedtools("intersect", "-u", "-a", probes, "-b", coding)
    beyond = bedtools(
        "subtract", "-a", blocks, "-b", tsv[["chrom", "start", "end"]]
    )
    return int((beyond[2] - beyond[1]).max())


def targeted_noncoding(tx, feats, existing, probes, limit):
    """
    Find the deliberately targeted non-coding sequence around one transcript.

    Takes the probe blocks overlapping the transcript's span, removes
    anything already covered by ``existing`` rows, and keeps the pieces
    longer than ``limit`` (so ordinary probe overhang is ignored). Introns
    the probes cover completely are always kept, since short targeted
    introns can be under ``limit``.

    Parameters
    ----------
    tx : str
        Transcript ID.
    feats : pandas.DataFrame
        Features from ``gff_features`` including this transcript.
    existing : pandas.DataFrame
        This transcript's rows already in the exons file (coding, 3'UTR
        and override adds), with chrom, start and end columns.
    probes : pandas.DataFrame
        Merged probe intervals from ``merged_probes``.
    limit : int
        Overhang cut-off from ``overhang_limit``.

    Returns
    -------
    pandas.DataFrame
        Merged targeted intervals, columns 0-2 (chrom, start, end).
        Empty if nothing is targeted beyond ``existing``.
    """
    chrom, _, _, _, exons = transcript_parts(feats, tx)
    span = pd.DataFrame([[chrom, exons[0][0], exons[-1][1]]])

    blocks = bedtools("intersect", "-u", "-a", probes, "-b", span)
    pieces = bedtools(
        "subtract", "-a", blocks, "-b", existing[["chrom", "start", "end"]]
    )
    if not pieces.empty:
        pieces = pieces[pieces[2] - pieces[1] > limit]

    introns = pd.DataFrame(
        [[chrom, e1, s2] for (_, e1), (s2, _) in zip(exons, exons[1:])]
    )
    if not introns.empty:
        cov = bedtools("coverage", "-a", introns, "-b", probes)
        introns = cov[cov.iloc[:, -1] == 1].iloc[:, :3]

    found = pd.concat(
        [pieces.iloc[:, :3] if not pieces.empty else pieces, introns]
    )
    if found.empty:
        return found
    return bedtools("merge", "-i", found.sort_values([1, 2]))


def targeted_rows(gene, tx, targeted, feats, existing):
    """
    Turn one transcript's targeted intervals into labelled exons-file rows.

    The transcript and its flanks are tiled by feature from the GFF -
    5'/3' flank, the UTR part of each exon, and each intron - and the
    tiles are cut to the targeted intervals, minus anything already
    covered by ``existing`` rows so the new rows abut them. Exons and
    introns are numbered in transcript order, matching the GFF and TSV.

    Parameters
    ----------
    gene : str
        Gene symbol for the output rows.
    tx : str
        Transcript ID.
    targeted : pandas.DataFrame
        Intervals from ``targeted_noncoding``, columns 0-2 (chrom, start,
        end).
    feats : pandas.DataFrame
        Features from ``gff_features`` including this transcript.
    existing : pandas.DataFrame
        This transcript's rows already in the exons file, with chrom,
        start and end columns.

    Returns
    -------
    pandas.DataFrame
        Rows in ``COLS`` format, in transcript order, labelled ``flank5``,
        ``flank3``, ``utr5_N``, ``utr3_N`` (UTR part of exon N) or
        ``intron_N``. A feature split into several pieces gets ``a``,
        ``b``, ... suffixes. Empty if nothing is left to add.
    """
    if targeted.empty:
        return pd.DataFrame(columns=COLS)
    chrom, strand, cds_lo, cds_hi, exons = transcript_parts(feats, tx)
    # exon numbers in transcript order, as in the GFF and TSV
    numbers = list(range(1, len(exons) + 1))
    if strand == "-":
        numbers.reverse()
    low, high = ("5", "3") if strand == "+" else ("3", "5")

    # tile the transcript and its flanks by feature
    lo_bound = min(targeted[1].min(), exons[0][0])
    hi_bound = max(targeted[2].max(), exons[-1][1])
    tiles = [
        (lo_bound, exons[0][0], f"flank{low}"),
        (exons[-1][1], hi_bound, f"flank{high}"),
    ]
    for (s, e), num in zip(exons, numbers):
        tiles.append((s, min(e, cds_lo), f"utr{low}_{num}"))
        tiles.append((max(s, cds_hi), e, f"utr{high}_{num}"))
    for (_, e1), (s2, _), n1, n2 in zip(
        exons, exons[1:], numbers, numbers[1:]
    ):
        tiles.append((e1, s2, f"intron_{min(n1, n2)}"))
    tiles = pd.DataFrame(
        [(chrom, s, e, label) for s, e, label in tiles if e > s]
    )

    rows = bedtools("intersect", "-a", tiles, "-b", targeted)
    if not rows.empty:
        rows = bedtools(
            "subtract", "-a", rows, "-b", existing[["chrom", "start", "end"]]
        )
    if rows.empty:
        return pd.DataFrame(columns=COLS)
    rows = rows.sort_values(1, ascending=(strand == "+"))  # transcript order
    return labelled_rows(rows, gene, tx, list(rows[3]))


def build(
    tsv,
    gff,
    probes,
    panel_transcripts,
    utr3_transcripts,
    overrides,
    curated_transcripts,
    pad,
):
    """
    Build the exons file and panel from the input files, by the rules above.

    Parameters
    ----------
    tsv : str or Path
        Approved coding exons +-``pad`` bp (e.g. symbols.exon_5bp.tsv).
    gff : str or Path
        Gzipped RefSeq GFF the TSV was made from.
    probes : str or Path
        Probe BED.
    panel_transcripts : str or Path
        Transcripts in the panel, one per line.
    utr3_transcripts : str or Path
        gene<TAB>transcript pairs whose 3'UTR goes in the panel.
    overrides : str or Path
        Overrides TSV (action, chrom, start, end, gene, transcript, label,
        reason).
    curated_transcripts : str or Path
        gene<TAB>transcript<TAB>requested rows whose targeted non-coding
        sequence goes in the exons file.
    pad : int
        Padding (bp) either side of the CDS in the TSV coding rows.

    Returns
    -------
    exons : pandas.DataFrame
        Exons file rows in ``COLS`` format, sorted by position.
    panel : pandas.DataFrame
        Panel rows (chrom, start, end, transcript), sorted by position.
    limit : int
        Probe overhang cut-off used for rule 4, from ``overhang_limit``.
    """
    coding = read_bed(tsv, COLS)
    panel_tx = set(read_bed(panel_transcripts, ["transcript"]).transcript)
    utr3_tx = read_bed(utr3_transcripts, ["gene", "transcript"])
    overrides = read_overrides(overrides)
    probes = merged_probes(probes)
    curated = read_bed(
        curated_transcripts, ["gene", "transcript", "requested"]
    )

    feats = gff_features(
        gff, set(utr3_tx.transcript) | set(curated.transcript)
    )
    utr3 = pd.concat([utr3_rows(gene, tx, feats, probes, pad)
                      for gene, tx in utr3_tx.itertuples(index=False)])

    coding = apply_replaces(
        coding, overrides[overrides.action == "replace"][COLS]
    )
    adds = overrides[overrides.action == "add"][COLS]
    excludes = overrides[overrides.action == "exclude"][COLS]

    scored = pd.concat([coding, utr3, adds])
    panel = anti_join(scored[scored.transcript.isin(panel_tx)], excludes, COLS)

    limit = overhang_limit(coding, panel_tx - set(curated.transcript), probes)
    extra = []
    for gene, tx, _ in curated.itertuples(index=False):
        existing = scored[scored.transcript == tx]
        targeted = targeted_noncoding(tx, feats, existing, probes, limit)
        extra.append(targeted_rows(gene, tx, targeted, feats, existing))

    order = ["chrom", "start", "end", "transcript", "exon"]

    def key(column):
        return column.map(CHROM_ORDER) if column.name == "chrom" else column

    exons = pd.concat([scored, *extra]).sort_values(order, key=key)
    panel = panel.sort_values(order, key=key)
    return exons, panel[["chrom", "start", "end", "transcript"]], limit


def parse_args(argv=None):
    """
    Parse the command-line arguments.

    Parameters
    ----------
    argv : list of str, optional
        Arguments to parse. Defaults to ``sys.argv[1:]``.

    Returns
    -------
    argparse.Namespace
        tsv, gff, probes, panel_transcripts, utr3_transcripts, overrides,
        curated_transcripts, exons_out and panel_out as ``Path``, and pad
        as ``int``.
    """
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("--tsv", type=Path, required=True, metavar="TSV",
                    help="approved coding exons +-PAD bp "
                    "(symbols.exon_5bp.tsv)")
    ap.add_argument("--gff", type=Path, required=True, metavar="GFF",
                    help="gzipped RefSeq GFF the TSV was made from")
    ap.add_argument("--probes", type=Path, required=True, metavar="BED",
                    help="probe BED, e.g. Probes_GRCh38_Lymphoma_v1.0.bed")
    ap.add_argument(
        "--panel-transcripts", type=Path, required=True, metavar="TXT",
        help="transcripts in the panel, one per line",
    )
    ap.add_argument(
        "--utr3-transcripts", type=Path, required=True, metavar="TSV",
        help="gene<TAB>transcript pairs whose 3'UTR goes in the panel",
    )
    ap.add_argument(
        "--overrides", type=Path, required=True, metavar="TSV",
        help="overrides TSV (action, chrom, start, end, gene, transcript, "
        "label, reason)",
    )
    ap.add_argument(
        "--curated-transcripts", type=Path, required=True, metavar="TSV",
        help="gene<TAB>transcript<TAB>requested rows whose targeted "
        "non-coding sequence goes in the exons file",
    )
    ap.add_argument("--pad", type=int, required=True, metavar="BP",
                    help="padding (bp) either side of the CDS in the TSV "
                    "coding rows")
    ap.add_argument("--exons-out", type=Path, required=True, metavar="BED",
                    help="exons file to write")
    ap.add_argument("--panel-out", type=Path, required=True, metavar="BED",
                    help="panel BED to write")
    return ap.parse_args(argv)


def main():
    """
    Parse the command-line inputs, build both files and write them out.

    Writes ``--exons-out`` and ``--panel-out`` as headerless tab-separated
    files, and prints the overhang cut-off and row counts.
    """
    args = parse_args()
    exons, panel, limit = build(
        tsv=args.tsv,
        gff=args.gff,
        probes=args.probes,
        panel_transcripts=args.panel_transcripts,
        utr3_transcripts=args.utr3_transcripts,
        overrides=args.overrides,
        curated_transcripts=args.curated_transcripts,
        pad=args.pad,
    )

    for df, path in [(exons, args.exons_out), (panel, args.panel_out)]:
        path.parent.mkdir(parents=True, exist_ok=True)
        df.to_csv(path, sep="\t", header=False, index=False)
    print(
        "Probe overhang past coding exons on non-curated panel transcripts: "
        f"up to {limit}bp"
    )
    print(f"Wrote {args.exons_out} ({len(exons)} rows)")
    print(f"Wrote {args.panel_out} ({len(panel)} rows)")


if __name__ == "__main__":
    main()
