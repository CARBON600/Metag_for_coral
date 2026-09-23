# metaWRAP arm (the 5th binner): an external guard wrapper (config
# `metawrap.guard`) runs metaWRAP's binning (metaBAT2 + MaxBin2 + CONCOCT),
# bin_refinement and reassemble_bins, and stages the final bins under
# <BASE_FINAL>/<sample>/04_FINAL_BINS_FOR_GTDB/.
# Two non-obvious constraints: BASE_FINAL is set per (treat, method, group), so the
# guard's FINAL_ROOT equals this combination's directory
# (output/binning/metawrap/{treat}/{method}/{group}/{sample}); and the guard's
# MIN_COMPLETENESS/MAX_CONTAMINATION come from `mimag`, so its internal
# preselection cannot drift from the downstream MIMAG gate. The full design
# (CheckM 1.0.12 archive vs the shared 1.2.3 rescoring, RefineM bypass, metaWRAP
# filters) is in the README "metaWRAP arm" section.

import hashlib
import os

MW = config.get("metawrap", {})
MW_ENV = MW.get("env", "metawrap-env")
MW_MEM_MB = config.get("resources", {}).get("metawrap_mem_mb", 64000)

if MW_BINNERS and not MW.get("guard"):
    raise ValueError(
        "config.yaml: 'metawrap.guard' must be set when 'metawrap' is in 'binners'.")
if MW_BINNERS and not MW.get("local_base"):
    raise ValueError(
        "config.yaml: 'metawrap.local_base' must be set when 'metawrap' is in 'binners'.")

# Resolve a relative guard against the repository root (the directory of the
# main Snakefile), so `tools/mw_guard.sh` stays valid wherever the checkout lives.
MW_GUARD = MW.get("guard", "")
if MW_GUARD and not os.path.isabs(MW_GUARD):
    MW_GUARD = os.path.normpath(os.path.join(workflow.basedir, MW_GUARD))

if MW_BINNERS:
    # The guard is invoked as `bash <guard>`, so +x is not required; existence and
    # byte-identity (md5) are the real gate. The md5 pins the exact frozen script,
    # which is a stronger identity check than a version string (guard_version in
    # config.yaml is documentation, kept in sync with it).
    if not os.path.isfile(MW_GUARD):
        raise ValueError(
            "config.yaml: metawrap.guard not found: {0}".format(MW_GUARD))
    if not os.access(MW_GUARD, os.R_OK):
        raise ValueError(
            "config.yaml: metawrap.guard is not readable: {0}".format(MW_GUARD))
    expected_md5 = MW.get("guard_md5")
    if not expected_md5:
        raise ValueError(
            "config.yaml: 'metawrap.guard_md5' must be set when 'metawrap' is in 'binners'.")
    with open(MW_GUARD, "rb") as _guard_fh:
        actual_md5 = hashlib.md5(_guard_fh.read()).hexdigest()
    if actual_md5 != expected_md5:
        raise ValueError(
            "config.yaml: metawrap.guard md5 mismatch: expected {0}, got {1} ({2})".format(
                expected_md5, actual_md5, MW_GUARD))


rule metawrap_bins:
    input:
        r1=assembly_reads_r1,
        r2=assembly_reads_r2,
        contigs="output/assemble/{treat}/{method}/{group}/{sample}/contigs_r2000bp.fasta"
    output:
        done="output/binning/metawrap/{treat}/{method}/{group}/{sample}/mw.done"
    threads: config["threads"]["metawrap"]
    conda:
        MW_ENV
    resources:
        mem_mb=MW_MEM_MB
    log:
        "logs/metawrap_bins/{treat}/{method}/{group}/{sample}.log"
    params:
        guard=MW_GUARD,
        env=MW_ENV,
        # group-level base; the guard appends /<sample> to form FINAL_ROOT, which
        # is exactly this combination's directory.
        base_final=lambda wc: f"output/binning/metawrap/{wc.treat}/{wc.method}/{wc.group}",
        local_base=MW["local_base"],
        min_comp=config["mimag"]["min_completeness"],
        max_con=config["mimag"]["max_contamination"],
        refine_mem_gb=MW.get("refine_mem_gb", 40),
        reassemble_mem_gb=MW.get("reassemble_mem_gb", 40),
        keep_short_work=MW.get("keep_short_work", 0),
        archive_existing=MW.get("archive_existing", 0),
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {log})
        exec > {log} 2>&1

        final_root="{params.base_final}/{wildcards.sample}"
        # The guard writes into <BASE_FINAL>/<SAMPLE_ID>; wipe only THIS
        # sample's subtree, never the shared group directory other samples use.
        rm -rf "$final_root"
        mkdir -p "{params.base_final}"

        METAWRAP_ENV={params.env} \
        THREADS={threads} \
        MIN_COMPLETENESS={params.min_comp} \
        MAX_CONTAMINATION={params.max_con} \
        REFINE_MEM_GB={params.refine_mem_gb} \
        REASSEMBLE_MEM_GB={params.reassemble_mem_gb} \
        KEEP_SHORT_WORK={params.keep_short_work} \
        ARCHIVE_EXISTING={params.archive_existing} \
        RESUME=0 \
        MW_LOCAL_BASE={params.local_base} \
        BASE_FINAL={params.base_final} \
        bash {params.guard} run \
          --samples {wildcards.sample} \
          --r1 {input.r1} \
          --r2 {input.r2} \
          --assembly {input.contigs}

        # The guard exits non-zero when metaWRAP fails; metaWRAP itself
        # hard-errors when no bin survives, so a 0-MAG run surfaces as a failed
        # job here (not as an `empty` status). Verify the payload this rule
        # promises downstream.
        bins_dir="$final_root/04_FINAL_BINS_FOR_GTDB"
        if [ ! -d "$bins_dir" ]; then
          echo "ERROR: metaWRAP did not produce $bins_dir"
          exit 1
        fi
        nbins=$(ls "$bins_dir"/*.fa 2>/dev/null | wc -l || true)
        if [ "$nbins" -eq 0 ]; then
          echo "ERROR: metaWRAP produced 0 final bins under $bins_dir"
          exit 1
        fi
        echo "metawrap_bins: $nbins final bins under $bins_dir"

        # Validate the guard's own summary when it is present (soft: the guard's
        # exit code already gates success).
        vsum=$(find "$final_root" -maxdepth 1 -name 'VALIDATION_SUMMARY.txt' | head -n 1 || true)
        if [ -n "$vsum" ] && ! grep -q "status: SUCCESS" "$vsum"; then
          echo "ERROR: $vsum does not report 'status: SUCCESS'"
          exit 1
        fi

        touch {output.done}
        """


# Pass the prepared metaWRAP bins through to the text pipeline's `refinem`
# location WITHOUT running RefineM: metaWRAP already refined/consolidated and
# reassembled the bins, and its reassembled contigs are not in
# contigs_r2000bp.fasta / {sample}_sorted.bam, which RefineM requires. This rule
# exists only so that checkm_lineage_wf / gtdbtk_classify (which both read
# output/refinem/{binner}/.../) have their input.
rule metawrap_ingest:
    input:
        checked="output/check_qc_inputs.done",
        prepared="output/bins_prepared/metawrap/{treat}/{method}/{group}/{sample}/prepare.done"
    output:
        done="output/refinem/metawrap/{treat}/{method}/{group}/{sample}/refinem.done"
    resources:
        mem_mb=2000
    log:
        "logs/metawrap_ingest/{treat}/{method}/{group}/{sample}.log"
    params:
        genomes=lambda wc: f"output/bins_prepared/metawrap/{wc.treat}/{wc.method}/{wc.group}/{wc.sample}",
        outdir=lambda wc: f"output/refinem/metawrap/{wc.treat}/{wc.method}/{wc.group}/{wc.sample}",
        # guard's archived CheckM 1.0.12 table (from the guard's 03_BIN_REASSEMBLY
        # tree); used only for the key cross-check below, never reported.
        guard_checkm=lambda wc: (
            "output/binning/metawrap/{0}/{1}/{2}/{3}/03_BIN_REASSEMBLY/"
            "reassembled_bins.checkm/storage/bin_stats_ext.tsv".format(
                wc.treat, wc.method, wc.group, wc.sample)),
        ext="fa",
    shell:
        r"""
        set -euo pipefail
        mkdir -p $(dirname {log})
        exec > {log} 2>&1
        rm -rf {params.outdir}
        mkdir -p {params.outdir}
        n=0
        for f in {params.genomes}/*.{params.ext}; do
          [ -e "$f" ] || continue
          cp "$f" {params.outdir}/
          n=$((n + 1))
        done
        echo "metawrap_ingest: $n bins passed through to {params.outdir} (RefineM skipped)"
        [ "$n" -eq 0 ] && echo "WARNING: 0 bins after ingest (allowed)"

        # Cross-check the final bin basenames against the guard's CheckM 1.0.12
        # key column when that archive is present. They must be identical; a
        # guard that renamed or dropped a bin is caught here instead of being
        # silently mis-joined by the downstream harvest/dRep layers.
        gc={params.guard_checkm}
        if [ -s "$gc" ]; then
          tmp_a=$(mktemp)
          tmp_b=$(mktemp)
          ls {params.outdir}/*.fa | sed 's#.*/##; s/\.fa$//' | sort > "$tmp_a"
          cut -f1 "$gc" | sort > "$tmp_b"
          if ! diff -u "$tmp_a" "$tmp_b"; then
            rm -f "$tmp_a" "$tmp_b"
            echo "ERROR: bin basenames != metaWRAP CheckM(1.0.12) keys ($gc)"
            exit 1
          fi
          rm -f "$tmp_a" "$tmp_b"
        else
          echo "WARNING: no metaWRAP CheckM(1.0.12) archive at $gc; skipping key cross-check"
        fi

        touch {output.done}
        """
