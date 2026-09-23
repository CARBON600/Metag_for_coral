#!/usr/bin/env python3
"""Assert the outputs of a Metag_for_coral-main validation run.

Two modes:

  outputs      walk <clone>/output for the given (narrow) matrix and print a
               PASS/WARN/FAIL table; exit 1 if anything FAILs. This is the part
               that catches the failures that are otherwise silent: missing
               stage markers, bin-name/CheckM-key drift, RefineM accidentally
               running for metawrap, a CheckM table that was copied from the
               guard instead of re-scored, an empty Strain heterogeneity column.
  fingerprint  parse a `snakemake -n` log, and either write the per-rule job
               counts as the expected fingerprint (--bless) or compare against
               it (--check). This is the regression part: after a code change,
               the DAG shape must still be what you expect.

Examples:
  python3 check_outputs.py outputs --output-dir <clone>/output --sample MINI1
  python3 tools/mini/check_outputs.py fingerprint --log <clone>/logs/dry.log --expect <clone>/tools/mini/fixture/expected_jobs.json --bless
  python3 tools/mini/check_outputs.py fingerprint --log <clone>/logs/dry.log --expect <clone>/tools/mini/fixture/expected_jobs.json --check
"""
import argparse
import glob
import hashlib
import json
import os
import re
import sys

# binner -> bin extension used by that arm downstream of prepare/refinem
EXT = {"unitem": "fna", "comebin": "fa", "metadecoder": "fasta",
       "semibin2_single": "fa", "metawrap": "fa"}

RESULTS = []


def add(level, name, detail=""):
    RESULTS.append((level, name, detail))


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def nonempty_lines(path):
    if not os.path.isfile(path):
        return 0
    n = 0
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.strip():
                n += 1
    return n


def read_header(path):
    if not os.path.isfile(path):
        return []
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.strip():
                return line.rstrip("\n").split("\t")
    return []


def checkm_keys(path):
    keys = []
    if not os.path.isfile(path):
        return keys
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.strip():
                keys.append(line.split("\t", 1)[0].strip())
    return keys


def tsv_rows(path):
    rows = []
    if not os.path.isfile(path):
        return rows
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        head = None
        for line in fh:
            if not line.strip():
                continue
            f = line.rstrip("\n").split("\t")
            if head is None:
                head = f
                continue
            rows.append(dict(zip(head, f)))
    return rows


def matrix_from_config(path):
    """Derive the (binner, treat, method, group, sample) matrix from a config file,
    mirroring the Snakefile: other binners run on assembly_treats, metawrap on
    metawrap_treats (default = assembly_treats)."""
    import yaml
    with open(path, "r", encoding="utf-8") as fh:
        cfg = yaml.safe_load(fh) or {}
    data_dir = cfg["data_dir"]
    samples = sorted(os.path.basename(p)[: -len("_1.fq.gz")]
                     for p in glob.glob(os.path.join(data_dir, "*_1.fq.gz")))
    if not samples:
        raise SystemExit("ERROR: no <sample>_1.fq.gz under data_dir={0}".format(data_dir))
    treats = list(cfg["assembly_treats"])
    methods = list(cfg["filter_methods"])
    groups = list(cfg["groups"])
    binners = list(cfg["binners"])
    mw_treats = list(cfg.get("metawrap_treats") or treats)
    pairs = ([(b, t) for b in binners if b != "metawrap" for t in treats]
             + [(b, t) for b in binners if b == "metawrap" for t in mw_treats])
    return [(b, t, m, g, s) for (b, t) in pairs for m in methods
            for g in groups for s in samples]


def combo_dirs(out, prefix, b, t, m, g, s):
    return os.path.join(out, prefix, b, t, m, g, s)


def check_outputs(args):
    out = args.output_dir
    if args.config:
        combos = matrix_from_config(args.config)
        print("[check_outputs] matrix from {0}: {1} combination(s)".format(
            args.config, len(combos)))
        for c in combos:
            print("    {0}/{1}/{2}/{3}/{4}".format(*c))
    else:
        binners = [x.strip() for x in args.binners.split(",") if x.strip()]
        treats = [x.strip() for x in args.treats.split(",") if x.strip()]
        methods = [x.strip() for x in args.methods.split(",") if x.strip()]
        groups = [x.strip() for x in args.groups.split(",") if x.strip()]
        combos = [(b, t, m, g, args.sample) for b in binners for t in treats
                  for m in methods for g in groups]

    # Arms that are allowed to come back empty: their own rules print "(allowed)"
    # (rules/binning.smk:31,66,105,208,245; rules/metawrap.smk:180).  Override with
    # FIXTURE_EMPTY_OK="a,b,c"; set it to "" to demand bins from every arm again.
    EMPTY_OK = {x.strip() for x in os.environ.get(
        "FIXTURE_EMPTY_OK", "unitem,comebin,metadecoder").split(",") if x.strip()}

    def verdict(n, arm, what):
        """PASS when n>0, PASS with a note when the arm is allowed to be empty,
        FAIL otherwise.  Returns (level, detail)."""
        if n:
            return "PASS", "{0} {1}".format(n, what)
        if arm in EMPTY_OK:
            return "PASS", "0 {0} (allowed empty for arm '{1}')".format(what, arm)
        return "FAIL", "0 {0} (needs >0)".format(what)

    # ---- 1. per-combination stage markers + bin/QC consistency -------------
    for (b, t, m, g, s) in combos:
        tag = "{0}/{1}/{2}/{3}/{4}".format(b, t, m, g, s)
        ext = EXT.get(b, "fa")

        prep = combo_dirs(out, "bins_prepared", b, t, m, g, s)
        n_prep = 0
        prep_ran = os.path.isfile(os.path.join(prep, "prepare.done"))
        if prep_ran:
            n_prep = len([x for x in os.listdir(prep) if x != "prepare.done"])
            level, detail = verdict(n_prep, b, "bin file(s)")
            add(level, "prepare_bins " + tag, detail)
        else:
            add("FAIL", "prepare_bins " + tag, "missing prepare.done")
        # "this arm ran and legitimately produced nothing" - distinct from "blocked"
        arm_empty = prep_ran and n_prep == 0 and b in EMPTY_OK

        rd = combo_dirs(out, "refinem", b, t, m, g, s)
        if os.path.isfile(os.path.join(rd, "refinem.done")):
            bins = sorted(os.path.basename(p)[: -(len(ext) + 1)]
                          for p in os.listdir(rd) if p.endswith("." + ext))
            level, detail = verdict(len(bins), b, "bin(s)")
            add(level, "refinem " + tag, detail)
        else:
            bins = []
            add("FAIL", "refinem " + tag, "missing refinem.done")

        cd = combo_dirs(out, "checkm", b, t, m, g, s)
        keys = checkm_keys(os.path.join(cd, "storage", "bin_stats_ext.tsv"))
        if os.path.isfile(os.path.join(cd, "checkm.done")):
            level, detail = verdict(len(keys), b, "row(s) in storage/bin_stats_ext.tsv")
            add(level, "checkm " + tag, detail)
        else:
            add("FAIL", "checkm " + tag, "missing checkm.done")

        qa = os.path.join(cd, "results.tsv")
        hdr = read_header(qa)
        need = ["Bin Id", "Strain heterogeneity", "Marker lineage"]
        if os.path.isfile(qa):
            missing = [c for c in need if c not in hdr]
            if missing and arm_empty:
                add("PASS", "checkm_qa " + tag,
                    "no scored bins (arm '{0}' is allowed to be empty)".format(b))
            else:
                add("PASS" if not missing else "FAIL", "checkm_qa " + tag,
                    "header ok" if not missing else "columns missing: {0}".format(missing))
        else:
            add("FAIL", "checkm_qa " + tag, "missing results.tsv (rule checkm_qa)")

        gd = combo_dirs(out, "gtdb", b, t, m, g, s)
        if os.path.isfile(os.path.join(gd, "gtdb.done")):
            n = nonempty_lines(os.path.join(gd, "{0}_{1}.bac120.summary.tsv".format(s, b)))
            n += nonempty_lines(os.path.join(gd, "{0}_{1}.ar122.summary.tsv".format(s, b)))
            add("PASS", "gtdbtk " + tag, "{0} summary row(s) (0 is legitimate)".format(
                max(0, n - 2)))
        else:
            add("FAIL", "gtdbtk " + tag, "missing gtdb.done")

        dd = combo_dirs(out, "drep_per_workflow", b, t, m, g, s)
        if os.path.isfile(os.path.join(dd, "derep.done")):
            add("PASS", "drep_per_workflow " + tag, "")
        elif arm_empty:
            add("PASS", "drep_per_workflow " + tag,
                "no dRep input (arm '{0}' is allowed to be empty)".format(b))
        else:
            add("FAIL", "drep_per_workflow " + tag, "missing derep.done")

        # bin basenames == CheckM keys <- the join the MIMAG gate relies on
        if bins and keys:
            if sorted(bins) == sorted(keys):
                add("PASS", "key==filename " + tag, "{0} id(s)".format(len(bins)))
            else:
                add("FAIL", "key==filename " + tag,
                    "refinem={0} checkm={1}".format(sorted(bins)[:3], sorted(keys)[:3]))

        if b == "metawrap":
            gdir = os.path.join(out, "binning", "metawrap", t, m, g, s)
            add("PASS" if os.path.isfile(os.path.join(gdir, "mw.done")) else "FAIL",
                "metawrap_bins " + tag, "")
            gb = os.path.join(gdir, "04_FINAL_BINS_FOR_GTDB")
            n = len([p for p in os.listdir(gb)]) if os.path.isdir(gb) else 0
            add("PASS" if n else "FAIL", "metawrap final bins " + tag,
                "{0} file(s) in 04_FINAL_BINS_FOR_GTDB".format(n))
            # RefineM must NOT have run for this arm
            if os.path.isfile(os.path.join(rd, "scaffold_stats.tsv")):
                add("FAIL", "refinem skipped " + tag,
                    "scaffold_stats.tsv present -> RefineM ran for metawrap")
            else:
                add("PASS", "refinem skipped " + tag, "no RefineM output (bypass ok)")
            # the reported table must be CheckM 1.2.3 output, not the guard's copy
            gtsv = os.path.join(gdir, "03_BIN_REASSEMBLY", "reassembled_bins.checkm",
                                "storage", "bin_stats_ext.tsv")
            ptsv = os.path.join(cd, "storage", "bin_stats_ext.tsv")
            if os.path.isfile(gtsv) and os.path.isfile(ptsv):
                if md5(gtsv) == md5(ptsv):
                    add("FAIL", "checkm rescored " + tag,
                        "pipeline table is byte-identical to the guard's 1.0.12 table")
                else:
                    add("PASS", "checkm rescored " + tag, "1.2.3 table differs from guard's")

    # ---- 2. summary tables --------------------------------------------------
    summ = os.path.join(out, "summary", "mimag_summary.tsv")
    rows = tsv_rows(summ)
    want = set(combos)
    got = set((r.get("binner"), r.get("treat"), r.get("method"), r.get("group"),
               r.get("sample")) for r in rows)
    miss = want - got
    add("PASS" if not miss else "FAIL", "mimag_summary rows",
        "{0}/{1} expected combos present{2}".format(len(got & want), len(want),
        "" if not miss else " missing={0}".format(sorted(miss)[:3])))
    if rows:
        bad = [r for r in rows if r.get("status") not in ("ok", "empty", "qc_empty")]
        add("PASS" if not bad else "WARN", "mimag_summary status",
            "statuses={0}".format(sorted(set(r.get("status") for r in rows))))
        for r in rows:
            if (r.get("binner"), r.get("treat"), r.get("method"), r.get("group"),
                    r.get("sample")) in want and r.get("status") != "ok":
                add("WARN", "status not ok: {0}/{1}/{2}/{3}/{4}".format(
                    r.get("binner"), r.get("treat"), r.get("method"), r.get("group"),
                    r.get("sample")), "status={0}".format(r.get("status")))

    for name in ("mimag_bins.tsv", "mimag_excluded.tsv"):
        p = os.path.join(out, "summary", name)
        add("PASS" if nonempty_lines(p) >= 1 else "FAIL", name + " exists", "")

    # ---- 3. Strain heterogeneity actually populated (the results.tsv path) --
    n_checked = n_strh = n_mark = 0
    for name in ("mimag_bins.tsv", "mimag_excluded.tsv"):
        for r in tsv_rows(os.path.join(out, "summary", name)):
            if (r.get("binner"), r.get("treat"), r.get("method"), r.get("group"),
                    r.get("sample")) not in want:
                continue
            if not r.get("completeness"):
                continue
            n_checked += 1
            if r.get("strain_heterogeneity", ""):
                n_strh += 1
            if r.get("marker_lineage", ""):
                n_mark += 1
    if n_checked == 0:
        add("WARN", "strain_heterogeneity populated", "no scored bin in this matrix")
    else:
        add("PASS" if n_strh == n_checked else "FAIL",
            "strain_heterogeneity populated",
            "{0}/{1} scored bins non-empty".format(n_strh, n_checked))
        add("PASS" if n_mark == n_checked else "WARN",
            "marker_lineage populated",
            "{0}/{1} scored bins non-empty".format(n_mark, n_checked))

    n_total = sum(int(r.get("n_total") or 0) for r in rows)
    n_mimag = sum(int(r.get("n_mimag") or 0) for r in rows)
    if n_total and n_mimag == 0:
        add("WARN", "MIMAG outcome",
            "0/{0} bins passed: expected for --synthesize fixtures, investigate if you "
            "used real genomes".format(n_total))
    elif n_mimag:
        add("PASS", "MIMAG outcome", "{0}/{1} bins passed".format(n_mimag, n_total))

    # ---- 4. downstream ------------------------------------------------------
    pw = tsv_rows(os.path.join(out, "drep_per_workflow", "per_workflow.tsv"))
    mw_rows = [r for r in pw if r.get("binner") == "metawrap"]
    add("PASS" if pw else "FAIL", "per_workflow.tsv", "{0} row(s)".format(len(pw)))
    add("PASS" if mw_rows else "WARN", "per_workflow has metawrap",
        "{0} metawrap row(s)".format(len(mw_rows)))
    for rel in ("summary/final_mags.tsv", "drep_cross_workflow/data_tables/Bdb.csv",
                "drep_cross_workflow/data_tables/Cdb.csv",
                "drep_cross_workflow/data_tables/Wdb.csv"):
        p = os.path.join(out, rel)
        add("PASS" if os.path.isfile(p) else "FAIL", os.path.basename(rel), "")


def parse_jobs(log_path):
    counts, total = {}, None
    if not os.path.isfile(log_path):
        return None, None
    with open(log_path, "r", encoding="utf-8", errors="replace") as fh:
        txt = fh.read()
    m = re.search(r"Job stats:\n(.*?)\n\n", txt, re.S)
    if not m:
        return None, None
    for line in m.group(1).splitlines():
        mm = re.match(r"^(\S+)\s+(\d+)\s*$", line.strip())
        if not mm:
            continue
        if mm.group(1) == "total":
            total = int(mm.group(2))
        elif mm.group(1) != "job":
            counts[mm.group(1)] = int(mm.group(2))
    return counts, total


def check_fingerprint(args):
    counts, total = parse_jobs(args.log)
    if counts is None:
        print("WARN: no 'Job stats' block in {0} (dry run may have failed earlier)".format(
            args.log))
        if os.path.isfile(args.log):
            with open(args.log, "r", encoding="utf-8", errors="replace") as fh:
                print(fh.read()[-2000:])
        return 2
    cur = {"total": total, "rules": counts}
    if args.bless:
        with open(args.expect, "w", encoding="utf-8") as fh:
            json.dump(cur, fh, indent=2, sort_keys=True)
        print("[fingerprint] blessed {0} rule(s), total={1} -> {2}".format(
            len(counts), total, args.expect))
        for k in sorted(counts):
            print("  {0:34s} {1}".format(k, counts[k]))
        return 0
    if not os.path.isfile(args.expect):
        print("WARN: no expected fingerprint at {0}; run with --bless first".format(args.expect))
        return 2
    with open(args.expect, "r", encoding="utf-8") as fh:
        exp = json.load(fh)
    ok = True
    for k in sorted(set(list(exp["rules"]) + list(counts))):
        a, b = exp["rules"].get(k, 0), counts.get(k, 0)
        if a != b:
            ok = False
            print("  DIFF {0:32s} expected={1} got={2}".format(k, a, b))
    if exp.get("total") != total:
        print("  DIFF {0:32s} expected={1} got={2}".format("TOTAL", exp.get("total"), total))
        ok = False
    print("[fingerprint] {0} (expected total={1}, got {2})".format(
        "MATCH" if ok else "MISMATCH", exp.get("total"), total))
    return 0 if ok else 1


def main(argv=None):
    ap = argparse.ArgumentParser(description="Assert validation-run outputs.")
    sub = ap.add_subparsers(dest="mode", required=True)

    o = sub.add_parser("outputs")
    o.add_argument("--output-dir", required=True)
    o.add_argument("--config", default=None,
                   help="mini config: derive the matrix (preferred)")
    o.add_argument("--sample", default="MINI1")
    o.add_argument("--binners", default="semibin2_single,metawrap")
    o.add_argument("--treats", default="control_assemble")
    o.add_argument("--methods", default="bt2")
    o.add_argument("--groups", default="exp")

    f = sub.add_parser("fingerprint")
    f.add_argument("--log", required=True)
    f.add_argument("--expect", required=True)
    f.add_argument("--bless", action="store_true")
    f.add_argument("--check", action="store_true")

    args = ap.parse_args(argv)
    if args.mode == "outputs":
        check_outputs(args)
        width = max(len(n) for _, n, _ in RESULTS)
        n_fail = 0
        for level, name, detail in RESULTS:
            if level == "FAIL":
                n_fail += 1
            print("{0:5s} {1:{2}s} {3}".format(level, name, width, detail))
        print("\n[check_outputs] {0} PASS / {1} WARN / {2} FAIL".format(
            sum(1 for r in RESULTS if r[0] == "PASS"),
            sum(1 for r in RESULTS if r[0] == "WARN"), n_fail))
        return 1 if n_fail else 0
    return check_fingerprint(args)


if __name__ == "__main__":
    sys.exit(main())
