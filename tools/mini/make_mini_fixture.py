#!/usr/bin/env python3
"""Build a tiny, format-correct fixture for Metag_for_coral-main.

It writes:

  <out>/data/<SAMPLE>_1.fq.gz, <SAMPLE>_2.fq.gz   paired 2x150 bp reads, seeded
  <out>/genomes/<name>.fna                        the genomes the reads came from
  <out>/fixture_manifest.tsv                      genome / length / depth / pairs
  <out>/config.mini.yaml                          partial Snakemake config override
  <out>/fixture_selfcheck.txt                     what this script verified itself

Genome source (one of):

  --genome A.fna [--genome B.fna ...]   real FASTAs. Best: one or two of your own
                                        MAGs, e.g.
                                          output/drep_cross_workflow/dereplicated_genomes/*.fna
                                        Complete-ish genomes so CheckM scores them and
                                        GTDB-Tk can classify them.
  --synthesize                          no external input: two random 2 Mb genomes with
                                        different GC (helps binning separate them).
                                        Format-correct, but CheckM scores ~0%
                                        completeness, so no MAG can pass the MIMAG
                                        gate: plumbing-only test.

Usage (on the cluster, snakemake env active):

  python3 tools/mini/make_mini_fixture.py --out <clone>/tools/mini/fixture --repo <repo> \
      --genome MAG_A.fna --genome MAG_B.fna \
      --guard <repo>/tools/mw_guard.sh --local-base /tmp/$USER/mw
"""
import argparse
import gzip
import io
import os
import random
import sys

COMP = str.maketrans("ACGTacgtN", "TGCAtgcaN")


def revcomp(s):
    return s.translate(COMP)[::-1]


def read_fasta(path):
    """Read a FASTA into [(name, SEQ)]."""
    out, name, chunks = [], None, []
    with open(path, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    out.append((name, "".join(chunks).upper()))
                name = line[1:].strip().split()[0] or "seq"
                chunks = []
            elif line:
                chunks.append(line.strip())
    if name is not None:
        out.append((name, "".join(chunks).upper()))
    if not out:
        raise SystemExit("ERROR: no FASTA records in {0}".format(path))
    return out


def collapse(records):
    return "".join(seq for _, seq in records)


def open_gz_text(path):
    # mtime=0 so identical seeds give byte-identical files.
    return io.TextIOWrapper(gzip.GzipFile(path, "wb", compresslevel=6, mtime=0),
                            encoding="ascii", newline="\n")


def synthesize(name, length, gc, seed):
    rnd = random.Random(seed)
    w = [(1 - gc) / 2, gc / 2, gc / 2, (1 - gc) / 2]
    return name, "".join(rnd.choices("ACGT", weights=w, k=length))


def qual_string(rnd, read_len, mean_q=36.0):
    """A realistic, deterministic Phred+33 string for ONE read.

    The fixture must NOT emit a flat quality string.  SPAdes' BayesHammer -- the
    error corrector that metaWRAP's per-bin reassembly always runs (it calls
    spades with --careful) -- infers the Phred offset from the observed quality
    alphabet and aborts with "Failed to determine offset!" (spades-hammer, err
    code 255) when that alphabet is too narrow.  metaWRAP then reports
    "was not successfully reassembled", leaves reassembled_bins/ empty and exits
    1.  The pipeline's own assembly is immune only because rules/assembly.smk
    passes --only-assembler, so it never runs BayesHammer.
    So: vary the qualities per base AND force a low-quality 3' tail, so that
    Phred+33 is unambiguous (chars below 59 appear in every read).
    """
    decay = 8.0 / max(1, read_len)              # mild 3' drop-off
    q = [int(max(2.0, min(40.0, rnd.gauss(mean_q, 3.0) - decay * i)))
         for i in range(read_len)]
    for k, tail_q in enumerate((2, 7, 15, 25)):  # force a low-quality tail
        if len(q) > k:
            q[len(q) - 1 - k] = tail_q
    return "".join(chr(33 + v) for v in q)


def simulate_reads(sample_id, genomes, depths, args, outdir):
    rnd = random.Random(args.seed)
    rndq = random.Random(args.seed + 977)   # separate stream: read layout unchanged
    flat_qual = args.qual_char * args.read_len if args.qual_mode == "flat" else None
    frags = []
    for gi, ((gname, seq), depth) in enumerate(zip(genomes, depths)):
        n = int(round(depth * len(seq) / (2.0 * args.read_len)))
        frags.extend([gi] * n)
    rnd.shuffle(frags)

    r1_path = os.path.join(outdir, "data", "{0}_1.fq.gz".format(sample_id))
    r2_path = os.path.join(outdir, "data", "{0}_2.fq.gz".format(sample_id))
    os.makedirs(os.path.dirname(r1_path), exist_ok=True)

    written = 0
    with open_gz_text(r1_path) as f1, open_gz_text(r2_path) as f2:
        for i, gi in enumerate(frags):
            gname, seq = genomes[gi]
            ins = int(rnd.gauss(args.insert, args.insert_sd))
            ins = max(args.read_len + 10, min(ins, args.max_insert))
            if ins > len(seq):
                ins = len(seq)
            if len(seq) - ins <= 0:
                continue
            start = rnd.randrange(0, len(seq) - ins)
            frag = seq[start:start + ins]
            r1 = frag[:args.read_len]
            r2 = revcomp(frag[-args.read_len:])
            if len(r1) != args.read_len or len(r2) != args.read_len:
                continue
            tag = "{0}:{1}:{2}".format(sample_id, gname, i)
            q = flat_qual if flat_qual is not None else qual_string(rndq, args.read_len)
            f1.write("@{0}/1\n{1}\n+\n{2}\n".format(tag, r1, q))
            f2.write("@{0}/2\n{1}\n+\n{2}\n".format(tag, r2, q))
            written += 1
    return r1_path, r2_path, written


def selfcheck(r1_path, r2_path, genomes, n_sample=200):
    """Re-open both FASTQs and verify structure, provenance and pairing geometry."""
    rnd = random.Random(1234)
    gmap = dict(genomes)
    lines = []
    for path in (r1_path, r2_path):
        with gzip.open(path, "rt", encoding="ascii") as fh:
            recs = [ln.rstrip("\n") for ln in fh]
        if len(recs) % 4 != 0:
            raise SystemExit("SELFCHECK FAIL: {0} line count {1} not divisible by 4".format(
                path, len(recs)))
        n = len(recs) // 4
        for k in range(0, min(n, n_sample)):
            h, s, p, q = recs[4 * k:4 * k + 4]
            if not h.startswith("@"):
                raise SystemExit("SELFCHECK FAIL: header {0!r}".format(h[:40]))
            if p != "+":
                raise SystemExit("SELFCHECK FAIL: separator {0!r}".format(p))
            if len(s) != len(q):
                raise SystemExit("SELFCHECK FAIL: seq/qual length {0} != {1}".format(len(s), len(q)))
            if set(s) - set("ACGTN"):
                raise SystemExit("SELFCHECK FAIL: non-ACGTN bases in {0}".format(h))
            qs = set(q)
            if len(qs) < 3:
                raise SystemExit("SELFCHECK FAIL: quality alphabet has only {0} distinct "
                                 "value(s); BayesHammer cannot infer the Phred offset "
                                 "(metaWRAP reassembly runs spades --careful)".format(len(qs)))
            if min(ord(c) for c in qs) >= 59:
                raise SystemExit("SELFCHECK FAIL: no quality char below 59 in {0}; Phred+33 "
                                 "is not unambiguous for BayesHammer".format(h[:40]))
            if k % 50 == 0:
                gname = h[1:].split(":")[1].split("/")[0]
                gseq = gmap[gname]
                # R1 is forward, R2 is the reverse complement of its genome segment,
                # so accept either orientation here (pairing is checked below).
                if gseq.find(s[:31]) < 0 and gseq.find(revcomp(s)[:31]) < 0:
                    raise SystemExit("SELFCHECK FAIL: {0} read absent from genome {1}".format(
                        path, gname))
        lines.append("{0}\tn_records={1}".format(os.path.basename(path), n))

    with gzip.open(r1_path, "rt", encoding="ascii") as f1, \
            gzip.open(r2_path, "rt", encoding="ascii") as f2:
        r1 = [ln.rstrip("\n") for ln in f1]
        r2 = [ln.rstrip("\n") for ln in f2]
    if len(r1) != len(r2):
        raise SystemExit("SELFCHECK FAIL: R1/R2 record count differ")
    n = len(r1) // 4
    ok = 0
    for k in rnd.sample(range(n), min(200, n)):
        h1, s1 = r1[4 * k], r1[4 * k + 1]
        h2, s2 = r2[4 * k], r2[4 * k + 1]
        if h1[:-2] != h2[:-2]:
            raise SystemExit("SELFCHECK FAIL: pair ids differ: {0} vs {1}".format(h1, h2))
        if not h1.endswith("/1") or not h2.endswith("/2"):
            raise SystemExit("SELFCHECK FAIL: /1 /2 suffixes missing")
        gname = h1[1:].split(":")[1]
        seq = gmap[gname]
        p1 = seq.find(s1[:31])
        p2 = seq.find(revcomp(s2)[:31])
        if p1 < 0 or p2 < 0 or not (1 <= (p2 - p1) <= 1500):
            raise SystemExit("SELFCHECK FAIL: pairing geometry wrong for {0} (p1={1} p2={2})".format(
                h1, p1, p2))
        ok += 1
    lines.append("paired-geometry samples ok: {0}".format(ok))
    return lines


def write_config_mini(path, args, fixture_data_dir, real_cfg):
    threads = {}
    for k, v in (real_cfg.get("threads") or {}).items():
        try:
            threads[k] = min(int(v), args.thread_cap)
        except (TypeError, ValueError):
            threads[k] = v
    threads["metawrap"] = args.thread_cap
    res = dict(real_cfg.get("resources") or {})
    res["metawrap_mem_mb"] = args.metawrap_mem_gb * 1024
    try:
        res["drep_pw_mem_mb"] = min(int(res.get("drep_pw_mem_mb", 16000)), 8000)
        res["drep_xw_mem_mb"] = min(int(res.get("drep_xw_mem_mb", 64000)), 16000)
    except (TypeError, ValueError):
        pass

    def ymap(d, indent=2):
        pad = " " * indent
        return "".join("{0}{1}: {2}\n".format(pad, k, v) for k, v in d.items())

    def ylist(items):
        return "\n".join("  - {0}".format(x) for x in items)

    binners = [b.strip() for b in args.binners.split(",") if b.strip()]
    text = """# Partial config OVERRIDE for the tiny end-to-end validation run.
# Merged over config.yaml by Snakemake (CLI configfile wins; dicts merge deeply).
#   bash run.sh --configfile tools/mini/fixture/config.mini.yaml -F --cores {cap}
#
# The matrix is deliberately minimal: every rule runs for at least one
# combination, and the expensive steps (SPAdes, metaWRAP reassembly, GTDB-Tk
# pplacer) see ~4 Mbp of sequence instead of a whole metagenome.
data_dir: "{data_dir}"

groups:
{groups}

filter_methods:
{methods}

assembly_treats:
{treats}

binners:
{binners}

metawrap_treats:
{treats}

metawrap:
  env: {mw_env}
  guard: "{guard}"
  local_base: "{local_base}"
  refine_mem_gb: {mem_gb}
  reassemble_mem_gb: {mem_gb}
  keep_short_work: 0
  archive_existing: 0

# Same keys as config.yaml (scaled down) so nothing disappears when merged.
threads:
{threads}
resources:
{resources}
""".format(
        cap=args.thread_cap,
        data_dir=fixture_data_dir,
        groups=ylist(args.groups.split(",")),
        methods=ylist(args.methods.split(",")),
        treats=ylist(args.treats.split(",")),
        binners=ylist(binners),
        mw_env=args.metawrap_env,
        guard=args.guard,
        local_base=args.local_base,
        mem_gb=args.metawrap_mem_gb,
        threads=ymap(threads),
        resources=ymap(res),
    )
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(text)
    return binners


def main(argv=None):
    ap = argparse.ArgumentParser(description="Build a tiny fixture for Metag_for_coral-main.")
    ap.add_argument("--out", required=True, help="fixture dir to create")
    ap.add_argument("--repo", required=True, help="Metag_for_coral-main checkout (for config.yaml)")
    ap.add_argument("--sample", default="MINI1")
    ap.add_argument("--genome", action="append", default=[],
                    help="genome FASTA (repeatable); omit to use --synthesize")
    ap.add_argument("--synthesize", action="store_true")
    ap.add_argument("--depths", default="30,12", help="per-genome depth, comma separated")
    ap.add_argument("--read-len", type=int, default=150)
    ap.add_argument("--insert", type=int, default=350)
    ap.add_argument("--insert-sd", type=int, default=50)
    ap.add_argument("--max-insert", type=int, default=700)
    ap.add_argument("--qual-char", default="I",
                    help="Phred char for --qual-mode flat only (I = Q40)")
    ap.add_argument("--qual-mode", choices=("illumina", "flat"), default="illumina",
                    help="illumina = varied per-base Phred with a low 3' tail (default; required by "
                         "metaWRAP's spades --careful / BayesHammer); flat = one repeated --qual-char")
    ap.add_argument("--seed", type=int, default=20260921)
    ap.add_argument("--groups", default="exp")
    ap.add_argument("--methods", default="bt2")
    ap.add_argument("--treats", default="control_assemble")
    ap.add_argument("--binners", default="semibin2_single,metawrap")
    ap.add_argument("--thread-cap", type=int, default=16)
    ap.add_argument("--metawrap-mem-gb", type=int, default=16)
    ap.add_argument("--metawrap-env", default="metawrap-env")
    ap.add_argument("--guard", default="/path/to/mw_guard.sh")
    ap.add_argument("--local-base", default="/tmp/$USER/mw")
    args = ap.parse_args(argv)

    try:
        import yaml
    except ImportError:
        raise SystemExit("ERROR: PyYAML missing; run this with the snakemake env python3 "
                         "(conda activate snakemake).")

    cfg_path = os.path.join(args.repo, "config.yaml")
    if not os.path.isfile(cfg_path):
        raise SystemExit("ERROR: {0} not found (--repo)".format(cfg_path))
    with open(cfg_path, "r", encoding="utf-8") as fh:
        real_cfg = yaml.safe_load(fh) or {}

    os.makedirs(os.path.join(args.out, "data"), exist_ok=True)
    os.makedirs(os.path.join(args.out, "genomes"), exist_ok=True)

    depths = [float(x) for x in args.depths.split(",")]
    genomes = []
    if args.synthesize:
        specs = [("MINIGEN_A", 2000000, 0.66, args.seed + 1),
                 ("MINIGEN_B", 2000000, 0.34, args.seed + 2),
                 ("MINIGEN_C", 2000000, 0.50, args.seed + 3)]
        for name, length, gc, sd in specs[:max(1, len(depths))]:
            genomes.append(synthesize(name, length, gc, sd))
    else:
        if not args.genome:
            raise SystemExit("ERROR: give --genome (repeatable) or --synthesize")
        for i, p in enumerate(args.genome):
            recs = read_fasta(p)
            name = "GENOME{0}_{1}".format(i + 1, os.path.basename(p).split(".")[0])[:40]
            genomes.append((name, collapse(recs)))
    if len(depths) < len(genomes):
        depths = depths + [depths[-1]] * (len(genomes) - len(depths))
    depths = depths[:len(genomes)]

    for name, seq in genomes:
        with open(os.path.join(args.out, "genomes", name + ".fna"), "w", encoding="utf-8") as fh:
            for i in range(0, len(seq), 80):
                fh.write(">{0}\n{1}\n".format(name, seq[i:i + 80]))

    r1, r2, n_pairs = simulate_reads(args.sample, genomes, depths, args, args.out)
    lines = selfcheck(r1, r2, genomes)

    with open(os.path.join(args.out, "fixture_manifest.tsv"), "w", encoding="utf-8") as fh:
        fh.write("genome\tlength_bp\ttarget_depth\tgc\n")
        for (name, seq), d in zip(genomes, depths):
            gc = (seq.count("G") + seq.count("C")) / float(len(seq))
            fh.write("{0}\t{1}\t{2}\t{3:.3f}\n".format(name, len(seq), d, gc))
        fh.write("# read_len\t{0}\n# insert\t{1}+-{2}\n# pairs\t{3}\n# seed\t{4}\n".format(
            args.read_len, args.insert, args.insert_sd, n_pairs, args.seed))

    mini_cfg = os.path.join(args.out, "config.mini.yaml")
    binners = write_config_mini(mini_cfg, args, os.path.join(os.path.abspath(args.out), "data"),
                               real_cfg)

    with open(os.path.join(args.out, "fixture_selfcheck.txt"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")

    print("[fixture] genomes   : {0} ({1} bp total)".format(
        len(genomes), sum(len(s) for _, s in genomes)))
    print("[fixture] read pairs: {0} (2 x {1} bp)".format(n_pairs, args.read_len))
    print("[fixture] data      : {0}".format(os.path.join(args.out, "data")))
    print("[fixture] config    : {0}".format(mini_cfg))
    print("[fixture] binners   : {0}".format(", ".join(binners)))
    print("[fixture] selfcheck : {0}".format(" | ".join(lines)))
    for tag, val in (("guard", args.guard), ("local-base", args.local_base)):
        if val.startswith("/path/to"):
            print("[fixture] WARNING: --{0} is still a placeholder: {1}".format(tag, val))
    return 0


if __name__ == "__main__":
    sys.exit(main())
