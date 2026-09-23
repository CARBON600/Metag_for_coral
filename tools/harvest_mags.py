#!/usr/bin/env python3
"""Harvest a per-combination status matrix and a MAG quality table.

Read-only with respect to the pipeline ``output/`` tree; it only writes into
``--outdir`` (default ``harvest/``).  The (binner, treat, method, group, sample)
grid is enumerated from disk, so combinations that were never executed are still
reported (``status == "not_run"``) rather than silently missing.

Usage:
    python3 tools/harvest_mags.py --output-dir output --outdir harvest
"""
import argparse
import ast
import csv
import glob
import json
import os
import sys

DEFAULT_BINNERS = ["unitem", "comebin", "metadecoder", "semibin2_single", "metawrap"]

# binner -> (raw bin dir relative to the combo dir, extension, gzipped)
# The metawrap raw dir MUST equal rules/common.smk::binner_raw_genome_dir.
BINNER_LAYOUT = {
    "unitem": ("consensus/bins", "fna", True),
    "comebin": ("comebin_res/comebin_res_bins", "fa", False),
    "metadecoder": (".", "fasta", False),
    "semibin2_single": ("output_bins", "fa", True),
    "metawrap": ("04_FINAL_BINS_FOR_GTDB", "fa", False),
}
BINNER_DONE = {
    "unitem": "consensus.done",
    "comebin": "comebin.done",
    "metadecoder": "cluster.done",
    "semibin2_single": "semibin.done",
    "metawrap": "mw.done",
}

# quality column -> candidate CheckM bin_stats_ext keys (first match wins).
# `marker lineage` is lowercase in CheckM 1.0.x AND 1.2.x (the title-case
# "Marker lineage" only appears in the qa table header, handled separately by
# _load_checkm_qa). The alias therefore fixes marker_lineage for ALL arms.
CHECKM_KEYS = [
    ("completeness", ("Completeness",)),
    ("contamination", ("Contamination",)),
    ("strain_heterogeneity", ("Strain heterogeneity",)),
    ("genome_size", ("Genome size",)),
    ("gc", ("GC",)),
    ("contigs", ("# contigs", "contigs")),
    ("marker_lineage", ("Marker lineage", "marker lineage")),
]

# columns that only exist in the `checkm qa -o 2` table (results.tsv), not in
# storage/bin_stats_ext.tsv; keyed by the harvest quality column name.
CHECKM_QA_KEYS = {"strain_heterogeneity", "marker_lineage"}

COMBO_HEADER = [
    "binner", "treat", "method", "group", "sample",
    "contigs_r2000", "sorted_bam", "raw_bins", "prepared_bins",
    "refinem", "checkm", "gtdb", "status",
]
QUALITY_HEADER = [
    "binner", "treat", "method", "group", "sample", "bin_id",
    "completeness", "contamination", "strain_heterogeneity",
    "genome_size", "gc", "contigs", "marker_lineage",
    "gtdb_domain", "gtdb_taxonomy",
]

# Columns of the MIMAG tables. mimag_bins.tsv holds the passing MAGs and
# mimag_excluded.tsv the rest, so the two reconcile to the full bin set.
MIMAG_HEADER = [
    "binner", "treat", "method", "group", "sample", "bin_id", "fasta_relpath",
    "completeness", "contamination", "strain_heterogeneity",
    "genome_size", "gc", "contigs", "marker_lineage",
    "gtdb_domain", "gtdb_taxonomy", "mimag_pass", "exclude_reason",
]
MIMAG_SUMMARY_HEADER = [
    "binner", "treat", "method", "group", "sample",
    "n_total", "n_mimag", "n_excluded", "status",
]


def _rel_parts(path, root):
    return os.path.relpath(path, root).split(os.sep)


def _glob(*parts):
    # Normalise separators so the glob works on both POSIX and Windows.
    return glob.glob(os.path.join(*parts).replace("\\", "/"))


def discover_grid(output_dir):
    """Return sorted dimension lists discovered from the output tree."""
    treats, methods, groups, samples, binners = set(), set(), set(), set(), set()

    for path in _glob(output_dir, "assemble", "*", "*", "*", "*", "contigs.fasta"):
        p = _rel_parts(path, os.path.join(output_dir, "assemble"))
        if len(p) >= 4:
            treats.add(p[0]); methods.add(p[1]); groups.add(p[2]); samples.add(p[3])

    for path in _glob(output_dir, "assemble", "*", "*", "*", "*", "*_sorted.bam"):
        p = _rel_parts(path, os.path.join(output_dir, "assemble"))
        if len(p) >= 4:
            treats.add(p[0]); methods.add(p[1]); groups.add(p[2]); samples.add(p[3])

    for sub in ("binning", "bins_prepared", "refinem", "checkm", "gtdb"):
        for path in _glob(output_dir, sub, "*", "*", "*", "*", "*"):
            p = _rel_parts(path, os.path.join(output_dir, sub))
            if len(p) >= 5:
                binners.add(p[0]); treats.add(p[1]); methods.add(p[2]); groups.add(p[3]); samples.add(p[4])

    for path in _glob(output_dir, "fq4dep", "*", "*", "*_1.fq.gz"):
        p = _rel_parts(path, os.path.join(output_dir, "fq4dep"))
        if len(p) >= 3:
            methods.add(p[0]); groups.add(p[1]); samples.add(p[2][: -len("_1.fq.gz")])

    for path in _glob(output_dir, "fastqs", "*", "*"):
        if os.path.isdir(path):
            p = _rel_parts(path, os.path.join(output_dir, "fastqs"))
            if len(p) >= 2:
                groups.add(p[0]); samples.add(p[1])

    for path in _glob(output_dir, "bt2", "*", "*.bam"):
        p = _rel_parts(path, os.path.join(output_dir, "bt2"))
        if len(p) >= 2:
            groups.add(p[0]); samples.add(p[1][: -len(".bam")])

    if not binners:
        binners = set(DEFAULT_BINNERS)

    return (sorted(treats), sorted(methods), sorted(groups), sorted(samples), sorted(binners))


def _count_files(directory, ext, gzipped):
    if not os.path.isdir(directory):
        return 0
    suffix = "." + ext + (".gz" if gzipped else "")
    n = 0
    for name in os.listdir(directory):
        if os.path.isfile(os.path.join(directory, name)) and name.endswith(suffix):
            n += 1
    return n


def _count_nonempty_lines(path):
    if not os.path.isfile(path):
        return 0
    n = 0
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.strip():
                n += 1
    return n


def _count_gtdb_rows(path):
    # The GTDB-Tk summary's first non-empty line is the header, so count data rows.
    n = _count_nonempty_lines(path)
    return n - 1 if n > 0 else 0


def combo_status(output_dir, binner, treat, method, group, sample):
    base = os.path.join(output_dir, "assemble", treat, method, group, sample)
    has_assembly = os.path.isfile(os.path.join(base, "contigs.fasta"))
    has_contigs = os.path.isfile(os.path.join(base, "contigs_r2000bp.fasta"))
    has_bam = os.path.isfile(os.path.join(base, sample + "_sorted.bam"))

    bin_base = os.path.join(output_dir, "binning", binner, treat, method, group, sample)
    raw_rel, ext, gz = BINNER_LAYOUT[binner]
    raw_dir = bin_base if raw_rel == "." else os.path.join(bin_base, raw_rel)
    raw_bins = _count_files(raw_dir, ext, gz)
    bin_done = os.path.isfile(os.path.join(bin_base, BINNER_DONE[binner]))

    prep_dir = os.path.join(output_dir, "bins_prepared", binner, treat, method, group, sample)
    prepared_bins = 0
    if os.path.isdir(prep_dir):
        prepared_bins = sum(1 for n in os.listdir(prep_dir) if n != "prepare.done")

    refinem = os.path.isfile(os.path.join(output_dir, "refinem", binner, treat, method, group, sample, "refinem.done"))
    checkm = os.path.isfile(os.path.join(output_dir, "checkm", binner, treat, method, group, sample, "checkm.done"))
    gtdb = os.path.isfile(os.path.join(output_dir, "gtdb", binner, treat, method, group, sample, "gtdb.done"))

    # Actual QC numbers, not just ".done" markers: an empty CheckM bin_stats_ext
    # and empty GTDB summaries mean QC produced nothing (e.g. bins were filtered
    # out), which must not be reported as "ok".
    checkm_table = os.path.join(output_dir, "checkm", binner, treat, method, group, sample, "storage", "bin_stats_ext.tsv")
    gtdb_dir = os.path.join(output_dir, "gtdb", binner, treat, method, group, sample)
    qc_rows = _count_nonempty_lines(checkm_table)
    qc_rows += _count_gtdb_rows(os.path.join(gtdb_dir, "{0}_{1}.bac120.summary.tsv".format(sample, binner)))
    qc_rows += _count_gtdb_rows(os.path.join(gtdb_dir, "{0}_{1}.ar122.summary.tsv".format(sample, binner)))

    if not (has_assembly or has_contigs or has_bam):
        status = "not_run"
    elif not bin_done:
        status = "partial"
    elif raw_bins == 0:
        status = "empty"
    elif prepared_bins == 0 or qc_rows == 0:
        status = "qc_empty"
    else:
        status = "ok"

    return {
        "contigs_r2000": int(has_contigs),
        "sorted_bam": int(has_bam),
        "raw_bins": raw_bins,
        "prepared_bins": prepared_bins,
        "refinem": int(refinem),
        "checkm": int(checkm),
        "gtdb": int(gtdb),
        "status": status,
    }


def _load_checkm(path):
    """Parse a CheckM storage/bin_stats_ext.tsv into {bin_id: metrics dict}."""
    out = {}
    if not os.path.isfile(path):
        return out
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            parts = line.split("\t", 1)
            if len(parts) != 2:
                continue
            bin_id, blob = parts[0].strip(), parts[1].strip()
            rec = None
            try:
                rec = json.loads(blob)
            except ValueError:
                try:
                    rec = ast.literal_eval(blob)
                except (ValueError, SyntaxError):
                    rec = None
            if isinstance(rec, dict):
                out[bin_id] = rec
    return out


def _load_gtdb(path):
    """Return {user_genome: classification} from a GTDB-Tk summary TSV."""
    out = {}
    if not os.path.isfile(path):
        return out
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        rows = [ln.rstrip("\n") for ln in fh if ln.strip()]
    if not rows:
        return out
    header = rows[0].split("\t")
    if "user_genome" not in header:
        return out
    i_genome = header.index("user_genome")
    i_class = header.index("classification") if "classification" in header else -1
    for line in rows[1:]:
        fields = line.split("\t")
        if i_genome >= len(fields):
            continue
        cls = fields[i_class] if 0 <= i_class < len(fields) else ""
        out[fields[i_genome]] = cls
    return out


def _load_checkm_qa(path):
    """Return {bin_id: {strain_heterogeneity, marker_lineage}} from a
    ``checkm qa -o 2 --tab_table`` results.tsv.

    That table is the only place CheckM writes 'Strain heterogeneity'; it is
    produced by the ``checkm_qa`` rule. Columns are looked up by name, so a
    build that reorders them is still handled. Empty/missing files yield {}.
    """
    out = {}
    if not os.path.isfile(path):
        return out
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        rows = [ln.rstrip("\n") for ln in fh if ln.strip()]
    if not rows:
        return out
    header = rows[0].split("\t")
    if "Bin Id" not in header:
        return out
    i_bin = header.index("Bin Id")
    i_strh = header.index("Strain heterogeneity") if "Strain heterogeneity" in header else -1
    i_mark = header.index("Marker lineage") if "Marker lineage" in header else -1
    for line in rows[1:]:
        fields = line.split("\t")
        if i_bin >= len(fields):
            continue
        out[fields[i_bin]] = {
            "strain_heterogeneity": fields[i_strh] if 0 <= i_strh < len(fields) else "",
            "marker_lineage": fields[i_mark] if 0 <= i_mark < len(fields) else "",
        }
    return out


def _domain(classification):
    if classification:
        first = classification.split(";")[0].strip()
        if first.startswith("d__"):
            return first[len("d__"):]
    return ""


def _metric(rec, keys):
    """Return the first present key's value from a CheckM record, else ''."""
    for key in keys:
        if key in rec:
            return rec[key]
    return ""


def _qa_or(rec, keys, qa_row, qa_key):
    """`_metric(rec, keys)`, falling back to the qa table when the ext file has
    no value for this column (e.g. `Strain heterogeneity`)."""
    val = _metric(rec, keys)
    if val == "" and qa_row:
        val = qa_row.get(qa_key, "")
    return val


def _refinem_fasta_relpath(binner, treat, method, group, sample, bin_id):
    """Path (relative to output/) of the RefineM bin that CheckM/GTDB scored.

    RefineM ``filter_bins`` writes ``<basename>.filtered.<ext>`` (see
    refinem/main.py), so the CheckM/GTDB ``bin_id`` already carries the
    ``.filtered`` suffix and the extension is the binner extension.  The
    extension is taken from ``BINNER_LAYOUT`` (the same mapping the raw-bin
    counter uses); ``emit_mimag`` verifies the constructed path exists so a
    drift between this mapping and ``common.smk``'s ``binner_extension`` is
    caught by the existence check below.
    """
    ext = BINNER_LAYOUT[binner][1]
    return "refinem/{0}/{1}/{2}/{3}/{4}/{5}.{6}".format(
        binner, treat, method, group, sample, bin_id, ext)


def emit_mimag(output_dir, outdir, min_completeness, max_contamination,
               binners=None, treats=None, methods=None, groups=None, samples=None):
    """Write the MIMAG tables into ``outdir``.

    Thresholds are inclusive (``>=`` / ``<=``) to match how dRep compares them.
    Only passing MAGs land in ``mimag_bins.tsv``; the rest go to
    ``mimag_excluded.tsv`` with an ``exclude_reason``.
    """
    d_treats, d_methods, d_groups, d_samples, d_binners = discover_grid(output_dir)
    treats = sorted(treats) if treats else d_treats
    methods = sorted(methods) if methods else d_methods
    groups = sorted(groups) if groups else d_groups
    samples = sorted(samples) if samples else d_samples
    binners = sorted(binners) if binners else (d_binners or list(DEFAULT_BINNERS))

    min_comp = float(min_completeness)
    max_con = float(max_contamination)

    pass_rows, excl_rows, summ_rows = [], [], []
    missing = []
    for binner in binners:
        for treat in treats:
            for method in methods:
                for group in groups:
                    for sample in samples:
                        cbase = os.path.join(output_dir, "checkm", binner, treat, method, group, sample)
                        gbase = os.path.join(output_dir, "gtdb", binner, treat, method, group, sample)
                        checkm = _load_checkm(os.path.join(cbase, "storage", "bin_stats_ext.tsv"))
                        qa = _load_checkm_qa(os.path.join(cbase, "results.tsv"))
                        gtdb = {}
                        gtdb.update(_load_gtdb(os.path.join(gbase, "{0}_{1}.bac120.summary.tsv".format(sample, binner))))
                        gtdb.update(_load_gtdb(os.path.join(gbase, "{0}_{1}.ar122.summary.tsv".format(sample, binner))))

                        c_total = c_mimag = c_excl = 0
                        for bid in sorted(set(checkm) | set(gtdb)):
                            c_total += 1
                            rec = checkm.get(bid, {})
                            comp = _metric(rec, ("Completeness",))
                            con = _metric(rec, ("Contamination",))
                            cls = gtdb.get(bid, "")

                            reason = ""
                            passed = False
                            if not rec:
                                reason = "no_checkm"
                            else:
                                try:
                                    low = float(comp) < min_comp
                                    high = float(con) > max_con
                                    if low and high:
                                        reason = "both"
                                    elif low:
                                        reason = "low_completeness"
                                    elif high:
                                        reason = "high_contamination"
                                    else:
                                        passed = True
                                except (TypeError, ValueError):
                                    reason = "no_checkm"

                            relpath = ""
                            if passed:
                                relpath = _refinem_fasta_relpath(binner, treat, method, group, sample, bid)
                                if not os.path.isfile(os.path.join(output_dir, relpath)):
                                    missing.append(os.path.join(output_dir, relpath))

                            qa_row = qa.get(bid, {})
                            row = [
                                binner, treat, method, group, sample, bid, relpath,
                                comp, con, _qa_or(rec, ("Strain heterogeneity",), qa_row, "strain_heterogeneity"),
                                _metric(rec, ("Genome size",)), _metric(rec, ("GC",)),
                                _metric(rec, ("# contigs", "contigs")),
                                _qa_or(rec, ("Marker lineage", "marker lineage"), qa_row, "marker_lineage"),
                                _domain(cls), cls,
                                "True" if passed else "False", reason,
                            ]
                            if passed:
                                pass_rows.append(row)
                                c_mimag += 1
                            else:
                                excl_rows.append(row)
                                c_excl += 1
                        # The combo status distinguishes a never-run branch from
                        # a binner that legitimately produced 0 bins (status
                        # values come from combo_status, the same parser used by
                        # the standalone harvest workflow).
                        st = combo_status(output_dir, binner, treat, method, group, sample)
                        summ_rows.append([binner, treat, method, group, sample,
                                          c_total, c_mimag, c_excl, st["status"]])

    if missing:
        raise RuntimeError(
            "mimag: {0} passing MAG(s) have no RefineM FASTA (e.g. {1}); "
            "check BINNER_LAYOUT vs common.smk:binner_extension".format(len(missing), missing[0]))

    bins_path = os.path.join(outdir, "mimag_bins.tsv")
    excl_path = os.path.join(outdir, "mimag_excluded.tsv")
    summ_path = os.path.join(outdir, "mimag_summary.tsv")
    with open(bins_path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(MIMAG_HEADER)
        writer.writerows(pass_rows)
    with open(excl_path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(MIMAG_HEADER)
        writer.writerows(excl_rows)
    with open(summ_path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(MIMAG_SUMMARY_HEADER)
        writer.writerows(summ_rows)

    sys.stderr.write(
        "harvest: mimag wrote {0} pass / {1} excluded over {2} combos -> {3}\n".format(
            len(pass_rows), len(excl_rows), len(summ_rows), bins_path))
    return {"n_pass": len(pass_rows), "n_excluded": len(excl_rows), "n_combo": len(summ_rows)}


def main(argv=None):
    ap = argparse.ArgumentParser(description="Harvest combo status and MAG quality.")
    ap.add_argument("--output-dir", default="output", help="pipeline output directory")
    ap.add_argument("--outdir", default="harvest", help="directory for the summary tables")
    ap.add_argument("--emit-mimag", action="store_true",
                    help="write mimag_{bins,summary,excluded}.tsv instead of the harvest tables")
    ap.add_argument("--min-completeness", type=float, default=50.0,
                    help="MIMAG completeness threshold (gate 1)")
    ap.add_argument("--max-contamination", type=float, default=10.0,
                    help="MIMAG contamination threshold (gate 1)")
    ap.add_argument("--binners", nargs="*", default=None)
    ap.add_argument("--treats", nargs="*", default=None)
    ap.add_argument("--methods", nargs="*", default=None)
    ap.add_argument("--groups", nargs="*", default=None)
    ap.add_argument("--samples", nargs="*", default=None)
    args = ap.parse_args(argv)

    if not os.path.isdir(args.output_dir):
        sys.stderr.write("ERROR: output dir not found: {0}\n".format(args.output_dir))
        return 1
    os.makedirs(args.outdir, exist_ok=True)

    if args.emit_mimag:
        emit_mimag(
            args.output_dir, args.outdir, args.min_completeness, args.max_contamination,
            binners=args.binners, treats=args.treats, methods=args.methods,
            groups=args.groups, samples=args.samples)
        return 0

    treats, methods, groups, samples, binners = discover_grid(args.output_dir)
    sys.stderr.write(
        "harvest: grid binners={0} treats={1} methods={2} groups={3} samples={4}\n".format(
            len(binners), len(treats), len(methods), len(groups), len(samples)
        )
    )

    combo_path = os.path.join(args.outdir, "combo_status.tsv")
    quality_path = os.path.join(args.outdir, "mag_quality.tsv")

    n_combo = 0
    with open(combo_path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(COMBO_HEADER)
        for binner in binners:
            for treat in treats:
                for method in methods:
                    for group in groups:
                        for sample in samples:
                            st = combo_status(args.output_dir, binner, treat, method, group, sample)
                            writer.writerow(
                                [binner, treat, method, group, sample]
                                + [st["contigs_r2000"], st["sorted_bam"], st["raw_bins"],
                                   st["prepared_bins"], st["refinem"], st["checkm"], st["gtdb"],
                                   st["status"]]
                            )
                            n_combo += 1

    n_mag = 0
    with open(quality_path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(QUALITY_HEADER)
        for binner in binners:
            for treat in treats:
                for method in methods:
                    for group in groups:
                        for sample in samples:
                            cbase = os.path.join(args.output_dir, "checkm", binner, treat, method, group, sample)
                            gbase = os.path.join(args.output_dir, "gtdb", binner, treat, method, group, sample)
                            checkm = _load_checkm(os.path.join(cbase, "storage", "bin_stats_ext.tsv"))
                            qa = _load_checkm_qa(os.path.join(cbase, "results.tsv"))
                            gtdb = {}
                            gtdb.update(_load_gtdb(os.path.join(gbase, "{0}_{1}.bac120.summary.tsv".format(sample, binner))))
                            gtdb.update(_load_gtdb(os.path.join(gbase, "{0}_{1}.ar122.summary.tsv".format(sample, binner))))
                            for bid in sorted(set(checkm) | set(gtdb)):
                                rec = checkm.get(bid, {})
                                qa_row = qa.get(bid, {})
                                row = [binner, treat, method, group, sample, bid]
                                for _col, keys in CHECKM_KEYS:
                                    if _col in CHECKM_QA_KEYS:
                                        val = _qa_or(rec, keys, qa_row, _col)
                                    else:
                                        val = _metric(rec, keys)
                                    row.append(val)
                                cls = gtdb.get(bid, "")
                                row.append(_domain(cls))
                                row.append(cls)
                                writer.writerow(row)
                                n_mag += 1

    sys.stderr.write("harvest: wrote {0} combo rows -> {1}\n".format(n_combo, combo_path))
    sys.stderr.write("harvest: wrote {0} MAG rows -> {1}\n".format(n_mag, quality_path))
    return 0


if __name__ == "__main__":
    sys.exit(main())
