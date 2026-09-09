#!/usr/bin/env python3
"""
Generate a minimal but *realistic* test dataset for the preprocessing stage.

The FASTQ files shipped for Phase 2 contained two synthetic reads and were only
ever meant to exercise samplesheet parsing. They cannot validate preprocessing:
fastp has nothing to trim and Bowtie2 has nothing to align. This script builds a
dataset small enough to run in seconds, yet structured so that every step has an
observable, predictable effect:

  * a mini "host" genome and a distinct "microbial" genome
  * reads drawn from both, in a known ratio, so host removal has a target
    retention rate to check against
  * adapter read-through on a known fraction of reads, so fastp reports a
    non-zero adapter-trimming count
  * degraded quality tails, so fastp's quality trimming does measurable work

Because the generated FASTQ/FASTA/index files are excluded by .gitignore, this
script is the reproducible source of truth — regenerate with:

    python3 test/data/make_test_data.py --outdir test/data

Determinism: a fixed RNG seed means the expected counts printed at the end stay
stable across runs and machines.
"""

import gzip
import random
import argparse
from pathlib import Path

# Illumina TruSeq adapter — the sequence fastp should detect and clip.
ADAPTER_R1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
ADAPTER_R2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"

READ_LEN = 100
FRAGMENT_LEN = 300

# Per-sample composition: how many read pairs, and what fraction comes from the
# host genome. Different ratios between samples make it obvious in the summary
# table that host removal is sample-specific rather than a fixed filter.
SAMPLES = {
    "S01": {"pairs": 2000, "host_fraction": 0.40},
    "S02": {"pairs": 2000, "host_fraction": 0.15},
}

# Fraction of pairs that get adapter read-through and severe quality decay.
ADAPTER_FRACTION = 0.20
BAD_TAIL_FRACTION = 0.15


def revcomp(seq):
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def random_genome(rng, length, gc_content):
    """
    Build a random sequence at a target GC content.

    Host and microbial genomes are given clearly different GC contents so that
    reads from one do not spuriously align to the other, and so the GC shift
    after host removal is visible in the fastp/QC metrics.
    """
    at, gc = (1.0 - gc_content) / 2, gc_content / 2
    return "".join(
        rng.choices("ACGT", weights=[at, gc, gc, at], k=length)
    )


def quality_string(rng, length, degrade_tail):
    """
    Phred+33 qualities: high quality up front, optionally decaying at the 3' end.

    Real Illumina reads lose quality along the read; reproducing that is what
    gives fastp's --cut_tail something meaningful to act on.
    """
    quals = []
    for i in range(length):
        if degrade_tail and i > length * 0.6:
            # Ramp down toward Q10-ish over the final 40% of the read
            progress = (i - length * 0.6) / (length * 0.4)
            mean_q = 36 - int(26 * progress)
        else:
            mean_q = 36
        q = max(2, min(40, int(rng.gauss(mean_q, 2))))
        quals.append(chr(q + 33))
    return "".join(quals)


def make_read_pair(rng, genome, add_adapter, degrade_tail):
    """Draw one fragment from `genome` and derive an R1/R2 pair from its ends."""
    start = rng.randrange(0, len(genome) - FRAGMENT_LEN)
    fragment = genome[start:start + FRAGMENT_LEN]

    r1 = fragment[:READ_LEN]
    r2 = revcomp(fragment[-READ_LEN:])

    if add_adapter:
        # Short insert => sequencer reads through into the adapter. Truncate the
        # insert and append adapter so the read stays READ_LEN long.
        insert = rng.randrange(40, 70)
        r1 = (fragment[:insert] + ADAPTER_R1)[:READ_LEN]
        r2 = (revcomp(fragment[:insert]) + ADAPTER_R2)[:READ_LEN]
        # Pad if the adapter did not fill the read
        r1 = r1 + "".join(rng.choices("ACGT", k=READ_LEN - len(r1)))
        r2 = r2 + "".join(rng.choices("ACGT", k=READ_LEN - len(r2)))

    q1 = quality_string(rng, len(r1), degrade_tail)
    q2 = quality_string(rng, len(r2), degrade_tail)
    return r1, q1, r2, q2


def write_sample(rng, sample, config, host_genome, microbe_genome, outdir):
    """Write <sample>_R1.fastq.gz / _R2.fastq.gz and report its composition."""
    pairs = config["pairs"]
    n_host = int(pairs * config["host_fraction"])

    r1_path = outdir / f"{sample}_R1.fastq.gz"
    r2_path = outdir / f"{sample}_R2.fastq.gz"

    with gzip.open(r1_path, "wt") as out1, gzip.open(r2_path, "wt") as out2:
        for i in range(pairs):
            source = "host" if i < n_host else "microbe"
            genome = host_genome if source == "host" else microbe_genome

            add_adapter = rng.random() < ADAPTER_FRACTION
            degrade_tail = rng.random() < BAD_TAIL_FRACTION

            r1, q1, r2, q2 = make_read_pair(rng, genome, add_adapter, degrade_tail)

            # Encode the true origin in the read name so failures are diagnosable
            name = f"{sample}:{source}:{i}"
            out1.write(f"@{name} 1:N:0:1\n{r1}\n+\n{q1}\n")
            out2.write(f"@{name} 2:N:0:1\n{r2}\n+\n{q2}\n")

    return {"pairs": pairs, "host": n_host, "microbe": pairs - n_host}


def make_dummy_dbs(outdir):
    """
    Stub-run DB gate directories (Phase 16).

    Kraken2/Bracken/HUMAnN guards pass the DB path through
    `file(..., checkIfExists: true)`, so a stub-run that must cover their
    branch logic needs *existing* paths — contents are never read (stub mode
    runs no real tool). Empty directories are the honest minimum: nothing
    pretends to be a real database. Guard-exercising stub runs point at:

        --kraken2_db      <outdir>/dummy_dbs/kraken2
        --bracken_db      <outdir>/dummy_dbs/kraken2
        --humann_db       <outdir>/dummy_dbs/humann
        --metaphlan_db    <outdir>/dummy_dbs/humann/metaphlan
    """
    root = outdir / "dummy_dbs"
    for sub in ("kraken2", "humann/chocophlan", "humann/uniref", "humann/metaphlan"):
        (root / sub).mkdir(parents=True, exist_ok=True)
    print(f"\nStub-run DB gates: {root}/ (kraken2/ + humann/{{chocophlan,uniref,metaphlan}}/)"
          " — empty dirs, only used by -stub-run branch guards")


def write_pathway_fixture(outdir):
    """
    KO→pathway mapping fixture (Phase 16).

    Two-column, header-less table in the format consumed by
    bin/integrate_functional.py (KO, pathway per line; a "ko:" prefix on the
    KO column is stripped on load, so one row carries it on purpose). Rows are
    synthetic and clearly labeled test_pathway_* — this exercises the Phase 14
    Pathway-join wiring inside the standard test profile, not biology.
    """
    path = outdir / "pathway_test.tsv"
    rows = [
        "ko:K00001\tko00010 test_pathway_alpha",
        "K00001\tko02010 test_pathway_beta",
        "K00002\tko00020 test_pathway_alpha",
    ]
    with open(path, "w") as handle:
        handle.write("\n".join(rows) + "\n")
    print(f"Pathway fixture: {path} (3 synthetic rows; "
          "wired via test.config pathway_db)")


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, default=Path("test/data"))
    parser.add_argument("--seed", type=int, default=42)
    return parser.parse_args()


def main():
    args = parse_args()
    rng = random.Random(args.seed)

    args.outdir.mkdir(parents=True, exist_ok=True)
    ref_dir = args.outdir / "host_reference"
    ref_dir.mkdir(exist_ok=True)

    # Two well-separated genomes: 45% vs 60% GC keeps cross-alignment negligible
    host_genome = random_genome(rng, 60_000, gc_content=0.45)
    microbe_genome = random_genome(rng, 60_000, gc_content=0.60)

    host_fasta = ref_dir / "mini_host.fa"
    with open(host_fasta, "w") as handle:
        handle.write(">mini_host_chr1 synthetic host genome for pipeline testing\n")
        for i in range(0, len(host_genome), 60):
            handle.write(host_genome[i:i + 60] + "\n")

    print(f"Host reference: {host_fasta} ({len(host_genome):,} bp)")

    summary = {}
    for sample, config in SAMPLES.items():
        summary[sample] = write_sample(
            rng, sample, config, host_genome, microbe_genome, args.outdir
        )

    print("\nGenerated samples (expected composition):")
    print(f"  {'sample':<8} {'pairs':>8} {'host':>8} {'microbe':>8}  expected retention")
    for sample, counts in summary.items():
        retention = 100.0 * counts["microbe"] / counts["pairs"]
        print(
            f"  {sample:<8} {counts['pairs']:>8,} {counts['host']:>8,} "
            f"{counts['microbe']:>8,}  ~{retention:.0f}% after host removal"
        )

    make_dummy_dbs(args.outdir)
    write_pathway_fixture(args.outdir)

    print(f"\nNext: build the Bowtie2 index\n"
          f"  bowtie2-build {host_fasta} {ref_dir / 'mini_host'}")


if __name__ == "__main__":
    main()
