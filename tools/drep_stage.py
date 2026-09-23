#!/usr/bin/env python3
"""Staging, verification and summary helpers for the dRep layer.

Reads the MIMAG table (output/summary/mimag_bins.tsv) written by
tools/harvest_mags.py rather than re-parsing CheckM/GTDB here. Stdlib only, so
it runs in any of the workflow's python environments.

Subcommands
-----------
stage      collect MAGs for one pool (per_workflow) or for all pools
           (cross_workflow) into a private, freshly wiped work directory, and
           write source/, manifest.tsv, genomeInfo.csv, g_source.list, status.txt.
verify     check the dRep log reports 100.00% of genomes passed checkM filtering.
summary    write a per-pool summary.tsv (n_total/n_mimag/n_rep/n_merged/status).
collect    concatenate every per-pool summary.tsv into per_workflow.tsv.
taxonomy   L9: cluster membership (Cdb) + winners (Wdb) + taxon propagation.
"""
import argparse
import csv
import glob
import os
import shutil
import sys

MANIFEST_HEADER = [
    "staged_basename", "binner", "treat", "method", "group", "sample",
    "bin_id", "orig_relpath",
]
# dRep's --genomeInfo requires the columns genome/completeness/contamination;
# strain_heterogeneity is carried as well so the documented winner score
# (which contains con*(strh/100)) is actually evaluated.
GENOMEINFO_HEADER = ["genome", "completeness", "contamination", "strain_heterogeneity"]
SUMMARY_HEADER = [
    "binner", "treat", "method", "group", "sample",
    "n_total", "n_mimag", "n_rep", "n_merged", "status",
]
FINAL_HEADER = [
    "binner", "treat", "method", "group", "sample", "bin_id",
    "cluster", "secondary_cluster", "is_representative", "score",
    "own_taxonomy", "representative_taxonomy", "final_taxonomy", "taxonomy_agreement",
]
FIVE_TUPLE = ["binner", "treat", "method", "group", "sample"]
CHECKM_PASS_NEEDLE = "100.00% of genomes passed checkM filtering"


def _read_tsv(path):
    if not os.path.isfile(path):
        return []
    with open(path, "r", encoding="utf-8", errors="replace", newline="") as fh:
        return [dict(row) for row in csv.DictReader(fh, delimiter="\t")]


def _read_csv(path):
    if not os.path.isfile(path):
        return []
    with open(path, "r", encoding="utf-8", errors="replace", newline="") as fh:
        return [dict(row) for row in csv.DictReader(fh)]


def _count_rows(path):
    if not os.path.isfile(path):
        return 0
    n = 0
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for i, line in enumerate(fh):
            if i == 0:
                continue
            if line.strip():
                n += 1
    return n


def _staged_name(row, mode):
    if mode == "cross_workflow":
        # Double underscores: binner can be "semibin2_single" and bin_id can be
        # "bin_1.filtered", so single underscores cannot be parsed back.
        return "{0}__{1}__{2}__{3}__{4}__{5}.fasta".format(
            row["binner"], row["treat"], row["method"], row["group"], row["sample"], row["bin_id"])
    return os.path.basename(row["fasta_relpath"])


def _select(rows, args):
    if args.mode == "cross_workflow":
        return list(rows)
    missing = [k for k in FIVE_TUPLE if getattr(args, k, None) in (None, "")]
    if missing:
        raise RuntimeError("per_workflow staging needs --{0}".format(", --".join(missing)))
    keep = []
    for row in rows:
        if all(row.get(k) == getattr(args, k) for k in FIVE_TUPLE):
            keep.append(row)
    return keep


def cmd_stage(args):
    rows = _read_tsv(args.mimag_bins)
    if not rows:
        sys.stderr.write("[drep_stage] WARNING: no rows in {0}\n".format(args.mimag_bins))
    rows = _select(rows, args)
    # mimag_bins.tsv only holds passing MAGs; keep the check anyway.
    rows = [r for r in rows if r.get("mimag_pass", "True") == "True"]

    workdir = args.workdir
    if os.path.isdir(workdir):
        shutil.rmtree(workdir)
    elif os.path.exists(workdir):
        os.remove(workdir)
    source = os.path.join(workdir, "source")
    os.makedirs(source)

    manifests, geninfo, gpaths = [], [], []
    for row in rows:
        src = os.path.join(args.output_dir, row["fasta_relpath"])
        if not os.path.isfile(src):
            raise RuntimeError("missing RefineM FASTA: {0}".format(src))
        staged = _staged_name(row, args.mode)
        dst = os.path.join(source, staged)
        if os.path.exists(dst):
            raise RuntimeError("duplicate staged basename: {0}".format(staged))
        shutil.copy2(src, dst)
        manifests.append({
            "staged_basename": staged, "binner": row["binner"], "treat": row["treat"],
            "method": row["method"], "group": row["group"], "sample": row["sample"],
            "bin_id": row["bin_id"], "orig_relpath": row["fasta_relpath"],
        })
        # A missing strain-heterogeneity value must not reach dRep as an empty
        # field: _validate_genomeInfo calls astype(float) on it and would crash.
        # 0 reproduces dRep's own "column absent" behaviour for that genome.
        geninfo.append({
            "genome": staged, "completeness": row["completeness"],
            "contamination": row["contamination"],
            "strain_heterogeneity": row.get("strain_heterogeneity") or "0",
        })
        gpaths.append(os.path.abspath(dst))

    with open(os.path.join(workdir, "manifest.tsv"), "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=MANIFEST_HEADER, delimiter="\t")
        writer.writeheader()
        for m in manifests:
            writer.writerow(m)

    with open(os.path.join(workdir, "genomeInfo.csv"), "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=GENOMEINFO_HEADER)
        writer.writeheader()
        for g in geninfo:
            writer.writerow(g)

    # One absolute path per line, no empty lines: dRep's load_genomes treats a
    # single positional argument as a path list and asserts each line is a file.
    with open(os.path.join(workdir, "g_source.list"), "w") as fh:
        for p in gpaths:
            fh.write(p + "\n")

    n = len(manifests)
    n_files = len([x for x in os.listdir(source)
                   if os.path.isfile(os.path.join(source, x))])
    # The four views of the pool must all have the same size.
    if not (n == len(geninfo) == len(gpaths) == n_files):
        raise RuntimeError(
            "staged pool mismatch: manifest={0} genomeInfo={1} list={2} files={3}".format(
                n, len(geninfo), len(gpaths), n_files))
    # dRep reads -g as a path list, so every line has to be an existing file.
    if n >= 2 and any(not os.path.isfile(p) for p in gpaths):
        raise RuntimeError("g_source.list references a path that is not a file")

    if n == 0:
        status = "skipped_empty"
    elif n == 1:
        status = "passthrough_single"
    else:
        status = "dereplicated"

    # dRep is not invoked for 0 or 1 MAG, so the declared outputs are written
    # here. For the dereplicated branch data_tables/ is left absent on purpose:
    # an empty *.csv would break dRep's WorkDirectory load and its Bdb check.
    if status in ("skipped_empty", "passthrough_single"):
        data_tables = os.path.join(workdir, "data_tables")
        os.makedirs(data_tables, exist_ok=True)
        for name in ("Bdb", "Cdb", "Wdb"):
            open(os.path.join(data_tables, name + ".csv"), "w").close()
        derep = os.path.join(workdir, "dereplicated_genomes")
        os.makedirs(derep, exist_ok=True)
        if status == "passthrough_single":
            shutil.copy2(gpaths[0], os.path.join(derep, os.path.basename(gpaths[0])))

    with open(os.path.join(workdir, "status.txt"), "w") as fh:
        fh.write(status + "\n")

    sys.stderr.write("[drep_stage] stage {0}: n_mimag={1} status={2}\n".format(
        workdir, n, status))
    return 0


def cmd_verify(args):
    text = ""
    if os.path.isfile(args.log):
        with open(args.log, "r", encoding="utf-8", errors="replace") as fh:
            text = fh.read()
    if CHECKM_PASS_NEEDLE not in text:
        raise RuntimeError(
            "checkM pass-rate line '{0}' missing from {1}".format(CHECKM_PASS_NEEDLE, args.log))
    # A dereplicated run has to leave at least one representative; otherwise an
    # empty dereplicated_genomes/ is reported as a legitimate n_rep=0 row.
    if args.derep_dir is not None:
        n_rep = 0
        if os.path.isdir(args.derep_dir):
            n_rep = len([x for x in os.listdir(args.derep_dir)
                         if os.path.isfile(os.path.join(args.derep_dir, x))])
        if n_rep < 1:
            raise RuntimeError(
                "no representative MAG in {0}".format(args.derep_dir))
    sys.stderr.write("[drep_stage] {0} reports 100% passed checkM filtering\n".format(args.log))
    return 0


def cmd_summary(args):
    n_mimag = _count_rows(os.path.join(args.workdir, "manifest.tsv"))
    derep = os.path.join(args.workdir, "dereplicated_genomes")
    n_rep = 0
    if os.path.isdir(derep):
        n_rep = len([x for x in os.listdir(derep)
                     if os.path.isfile(os.path.join(derep, x))])
    status = ""
    status_file = os.path.join(args.workdir, "status.txt")
    if os.path.isfile(status_file):
        with open(status_file, "r", encoding="utf-8", errors="replace") as fh:
            status = fh.read().strip()

    n_total = ""
    for row in _read_tsv(args.mimag_summary):
        if all(row.get(k) == getattr(args, k) for k in FIVE_TUPLE):
            n_total = row.get("n_total", "")
            break

    row = [args.binner, args.treat, args.method, args.group, args.sample,
           n_total, n_mimag, n_rep, int(n_mimag) - int(n_rep), status]
    with open(os.path.join(args.workdir, "summary.tsv"), "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(SUMMARY_HEADER)
        writer.writerow(row)
    sys.stderr.write("[drep_stage] summary {0}: n_mimag={1} n_rep={2} status={3}\n".format(
        args.workdir, n_mimag, n_rep, status))
    return 0


def cmd_collect(args):
    rows = []
    pattern = os.path.join(args.per_workflow_root, "*", "*", "*", "*", "*", "summary.tsv")
    for path in sorted(glob.glob(pattern)):
        rows.extend(_read_tsv(path))

    n_bins = _count_rows(args.mimag_bins)
    total_mimag = 0
    for row in rows:
        try:
            total_mimag += int(row.get("n_mimag", 0))
        except (TypeError, ValueError):
            pass
    # Both sides of this check come from the same file, so they have to agree.
    if total_mimag != n_bins:
        raise RuntimeError(
            "funnel mismatch: sum(n_mimag)={0} != rows(mimag_bins)={1}".format(total_mimag, n_bins))

    with open(args.out, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(SUMMARY_HEADER)
        for row in rows:
            writer.writerow([row.get(col, "") for col in SUMMARY_HEADER])
    sys.stderr.write("[drep_stage] collect: {0} pools, sum(n_mimag)={1} -> {2}\n".format(
        len(rows), total_mimag, args.out))
    return 0


def cmd_taxonomy(args):
    mimag = _read_tsv(args.mimag_bins)
    manifest = _read_tsv(args.manifest)
    cdb = _read_csv(args.cdb)
    wdb = _read_csv(args.wdb)

    mimag_by_key = {}
    for row in mimag:
        key = tuple(row.get(k, "") for k in FIVE_TUPLE) + (row.get("bin_id", ""),)
        mimag_by_key[key] = row
    m_by_staged = {m.get("staged_basename", ""): m for m in manifest}

    # Membership comes from Cdb's secondary_cluster column; winners from Wdb
    # (genome/cluster/score, where 'cluster' holds the set id).
    g2cluster = {}
    for row in cdb:
        genome = row.get("genome", "")
        if genome:
            g2cluster[genome] = row.get("secondary_cluster", "")
    representatives = {}
    for row in wdb:
        cluster = row.get("cluster", "")
        if cluster and cluster not in representatives:
            representatives[cluster] = (row.get("genome", ""), row.get("score", ""))

    out_rows = []
    sets = {}
    for genome, cluster in g2cluster.items():
        member = m_by_staged.get(genome)
        if member is None:
            continue
        mkey = tuple(member.get(k, "") for k in FIVE_TUPLE) + (member.get("bin_id", ""),)
        own = mimag_by_key.get(mkey, {}).get("gtdb_taxonomy", "")

        rep_genome, rep_score = representatives.get(cluster, ("", ""))
        rep = m_by_staged.get(rep_genome)
        rep_tax = ""
        if rep is not None:
            rkey = tuple(rep.get(k, "") for k in FIVE_TUPLE) + (rep.get("bin_id", ""),)
            rep_tax = mimag_by_key.get(rkey, {}).get("gtdb_taxonomy", "")

        out_rows.append([
            member.get("binner", ""), member.get("treat", ""), member.get("method", ""),
            member.get("group", ""), member.get("sample", ""), member.get("bin_id", ""),
            cluster, cluster,
            "True" if genome == rep_genome else "False",
            rep_score if genome == rep_genome else "",
            own, rep_tax, rep_tax,
            "True" if own == rep_tax else "False",
        ])
        sets.setdefault(cluster, set()).add(own)

    discordant = sum(1 for taxes in sets.values() if len(taxes) > 1)
    with open(args.out, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(FINAL_HEADER)
        writer.writerows(out_rows)
    sys.stderr.write(
        "[drep_stage] taxonomy: rows={0} sets={1} discordant_sets={2} -> {3}\n".format(
            len(out_rows), len(sets), discordant, args.out))
    return 0


def _add_common_pool_args(p):
    for key in FIVE_TUPLE:
        p.add_argument("--" + key, default=None)


def main(argv=None):
    ap = argparse.ArgumentParser(description="dRep layer staging and summarisation.")
    sub = ap.add_subparsers(dest="command")

    p = sub.add_parser("stage", help="stage MAGs into a private work directory")
    p.add_argument("--mode", choices=["per_workflow", "cross_workflow"], required=True)
    p.add_argument("--mimag-bins", default="output/summary/mimag_bins.tsv")
    p.add_argument("--output-dir", default="output")
    p.add_argument("--workdir", required=True)
    _add_common_pool_args(p)
    p.set_defaults(func=cmd_stage)

    p = sub.add_parser("verify", help="assert the dRep checkM-pass log line")
    p.add_argument("--log", required=True)
    p.add_argument("--derep-dir", default=None,
                   help="if set, require at least one representative in this directory")
    p.set_defaults(func=cmd_verify)

    p = sub.add_parser("summary", help="write a per-pool summary.tsv")
    p.add_argument("--workdir", required=True)
    p.add_argument("--mimag-summary", default="output/summary/mimag_summary.tsv")
    _add_common_pool_args(p)
    p.set_defaults(func=cmd_summary)

    p = sub.add_parser("collect", help="concatenate per-pool summaries")
    p.add_argument("--per-workflow-root", default="output/drep_per_workflow")
    p.add_argument("--mimag-bins", default="output/summary/mimag_bins.tsv")
    p.add_argument("--out", default="output/drep_per_workflow/per_workflow.tsv")
    p.set_defaults(func=cmd_collect)

    p = sub.add_parser("taxonomy", help="L9 membership + taxon propagation")
    p.add_argument("--mimag-bins", default="output/summary/mimag_bins.tsv")
    p.add_argument("--cdb", required=True)
    p.add_argument("--wdb", required=True)
    p.add_argument("--manifest", required=True)
    p.add_argument("--out", required=True)
    p.set_defaults(func=cmd_taxonomy)

    args = ap.parse_args(argv)
    if not getattr(args, "func", None):
        ap.error("a subcommand is required")
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
