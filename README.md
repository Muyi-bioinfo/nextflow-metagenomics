# nextflow-metagenomics

**English** | [中文](README_CN.md)

[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)
[![Nextflow 26.04.4](https://img.shields.io/badge/Nextflow-26.04.4-0dc09d?logo=nextflow&logoColor=white)](https://www.nextflow.io/)
![conda](https://img.shields.io/badge/conda-supported-3EB049?logo=anaconda&logoColor=white)
![Docker](https://img.shields.io/badge/Docker-supported-2496ED?logo=docker&logoColor=white)
![Singularity](https://img.shields.io/badge/Singularity-supported-5E97F6)

A modular, reproducible shotgun-metagenomics workflow built on
[Nextflow DSL2](https://www.nextflow.io/): taxonomic classification, functional
profiling, cross-sample integration, MAG reconstruction, and result
visualization.

**Project status:** The full pipeline is implemented (input → preprocessing →
read-based ∥ assembly → MAG → QC → dereplication → classification → gene
prediction → annotation → abundance → integration → reporting → visualization).

> **Verification scope** (what has actually been executed, versus what is
> pending): preprocessing, assembly (MEGAHIT), mapping, binning, gene
> prediction, abundance, integration and MultiQC have been verified with real
> runs on synthetic test data. Kraken2, Bracken, HUMAnN, CheckM2, dRep,
> GTDB-Tk, DIAMOND, eggNOG-mapper and RGI depend on multi-GB databases that are
> not available on the development machine — these have been verified via
> **stub-run** only (channel topology, parameter guards, output structure).
> Cross-sample merge logic and plotting (except the
> database-free CoverM abundance heatmap) are likewise database-dependent and
> verified with stub/synthetic tables. **Execution against real databases and
> container runs on real hardware have not been verified.** See
> [Limitations](#limitations).

---

## Table of contents

- [Overview](#overview)
- [Workflow](#workflow)
- [Repository structure](#repository-structure)
- [Features](#features)
- [Quick start](#quick-start)
- [Input & output](#input--output)
- [Configuration & databases](#configuration--databases)
- [Testing & reproducibility](#testing--reproducibility)
- [Limitations](#limitations)
- [Tools & citations](#tools--citations)
- [Documentation](#documentation)
- [License](#license)

## Overview

`nextflow-metagenomics` is a Nextflow DSL2 workflow for shotgun metagenomics.
It consumes raw FASTQ, performs QC and optional host removal, then runs
**read-based analysis** (taxonomic classification, abundance estimation,
functional profiling) and the **assembly/MAG path** (assembly, MAG
reconstruction, QC, dereplication, classification, gene prediction, functional
annotation, abundance) **in parallel** from the same clean reads. Per-sample
read-based results are merged into cross-sample matrices, all results are joined
into core tables, and a MultiQC report plus a set of figures are generated.

The workflow is built as a layered, single-responsibility Nextflow DSL2 project:
a thin entry point, one orchestrating workflow, per-stage subworkflows, and
local process modules — with `tuple(meta, …)` metadata threading every channel
from input to report.

## Workflow

| Stage | Tools |
| ----- | ----- |
| Quality control & host removal | FastQC → fastp → Bowtie2 + samtools |
| Read-based profiling | Kraken2 → Bracken ∥ HUMAnN |
| Assembly | MEGAHIT ∥ metaSPAdes → QUAST |
| Mapping & coverage | Bowtie2 → samtools → depth |
| Binning | MetaBAT2 |
| MAG QC | CheckM2 |
| Dereplication | dRep |
| Taxonomy | GTDB-Tk |
| Gene prediction | Prodigal |
| Functional annotation | DIAMOND ∥ eggNOG-mapper ∥ RGI (CARD) |
| MAG abundance | CoverM |
| Integration | custom Python — two core tables |
| Reporting | MultiQC + 7 figures (custom Python) |

![nextflow-metagenomics workflow overview](docs/images/workflow_compact.png)

- **Parallel:** read-based analysis and the assembly/MAG chain fork from clean
  reads and never depend on each other.
- **Serial:** assembly → mapping → binning, and MAG QC → dereplication →
  classification (dRep scoring consumes QC output, and GTDB-Tk classifies only
  the dereplicated representatives).
- **Fan-out:** gene prediction / annotation / abundance consume the
  representative set in parallel; an integration stage joins their summaries.
- A full vertical version is available at
  [docs/images/workflow.png](docs/images/workflow.png).

## Repository structure

```text
nextflow-metagenomics/
├── main.nf                    # thin DSL2 entry → workflows/mag.nf
├── workflows/mag.nf           # main workflow: all stages wired here
├── nextflow.config            # parameter definitions + profiles + db_dir derivation
├── conf/                      # base / conda / docker / singularity / slurm / test
├── subworkflows/local/        # 14 subworkflows (one per stage)
├── modules/local/             # 45 local modules / 47 processes (single responsibility)
├── bin/                       # 14 Python parsing/orchestration scripts + 1 diagram utility
├── assets/                    # sentinel files (e.g. empty.tsv)
├── test/                      # synthetic data, test config, unit tests
├── docs/                      # architecture / channels / workflow / input / output / database / parameters
├── environment.yml            # pinned development environment (nf-meta)
├── setup_env.sh               # environment creation + post-install fixes
└── LICENSE
```

## Features

- **Nextflow DSL2, layered architecture** — a thin `main.nf` → `workflows/mag.nf`
  → 14 subworkflows → 45 local modules (47 process definitions). Each process
  has a single responsibility (`FASTQC` / `FASTP`, not `RUN_ALL_ANALYSIS`) and
  declares its own conda/container directives, publishDir, stub block and
  `versions.yml`.

- **Channel & metadata contract** — every channel is `tuple(meta, …)`;
  `meta = [id, group, batch, host, single_end]` is carried and extended
  throughout the workflow (`assembler`, `assembly_mode`, `samples` appended
  after assembly).

- **Explicit data-flow design (`join` / `combine` / `collect` / `groupTuple`):**
  - `combine` scatters Bracken over taxonomic levels (sample × level).
  - `combine(by: 0)` pairs contigs with the reads that assembled them, keyed on
    sample identity (correct under `--assembler both`).
  - `groupTuple` + `groupKey` aggregate sorted BAMs per assembly unit.
  - `collectFile` / `toSortedList` aggregate batch-level inputs with deterministic
    ordering so `-resume` hits its cache.
  - `flatMap` expands a bin directory into one record per MAG; empty channels are
    a legal state (0-bin samples).

- **Parameterization & skip strategy:**
  - All ~80 parameters are declared in `nextflow.config`; database/index paths
    are passed via `params.*` only (no hardcoded paths).
  - Per-tool `--skip_*` switches plus a master switch per branch; a missing
    read-based database **warns and skips** that parallel branch, while a missing
    MAG-level database (CheckM2 / GTDB-Tk / annotation) **errors explicitly** —
    no silent degradation and no fabricated output.

- **Resource management** — four label tiers
  (`process_single/low/medium/high`) in `conf/base.config`, capped by
  `resourceLimits`; per-tool overrides for memory-hungry tools (KRAKEN2, GTDBTK
  at 64 GB) and retry-with-backoff on resource exit codes.

- **Reproducibility & run profiles** — version pinning (`environment.yml` +
  per-process conda/container directives), a fixed MetaBAT2 seed, `versions.yml`
  per process, and `-stub-run` support. Profiles for `conda`, `docker`,
  `singularity`, `slurm` and `test`, plus a centralized `--db_dir` layout that
  derives eight database paths from one tree.

- **Per-sample → cross-sample integration** — per-sample Bracken/HUMAnN tables
  are merged into sample×taxa and sample×pathway matrices (plus a Bray–Curtis
  distance matrix); QC/classification/annotation/abundance are joined into two
  core tables, and 7 figures are rendered as a read-only consumer layer.

## Quick start

Installation — use `setup_env.sh` (it applies post-install fixes that
`mamba env create` alone cannot express; see the notes at the top of
`environment.yml`):

```bash
bash setup_env.sh          # create + repair + per-tool startup verification
conda activate nf-meta
```

Run:

```bash
# Standard run (read-based and assembly branches; branches without a database
# are skipped automatically with a warning)
nextflow run main.nf --input samplesheet.csv

# Specify a batch id (output goes to results/<batch_id>/)
nextflow run main.nf --input samplesheet.csv --batch_id batch_001

# Test run on built-in synthetic data (requires the nf-meta environment)
nextflow run main.nf -profile test \
  --skip_read_based --skip_mag_qc --skip_dereplication \
  --skip_taxonomy --skip_annotation

# Stub run: verify channel topology + output structure with no data or databases
# (the dummy/fake DB flags exercise every branch without real tools)
nextflow run main.nf -profile test -stub-run \
  --kraken2_db test/data/dummy_dbs/kraken2 --bracken_db test/data/dummy_dbs/kraken2 \
  --humann_db test/data/dummy_dbs/humann \
  --metaphlan_db test/data/dummy_dbs/humann/metaphlan \
  --checkm2_db /fake --gtdbtk_db /fake --diamond_db /fake \
  --eggnog_db /fake --card_db /fake

# Alternative runtimes / scheduler (add databases as needed)
nextflow run main.nf --input samplesheet.csv -profile docker
nextflow run main.nf --input samplesheet.csv -profile singularity     # Apptainer/Singularity
nextflow run main.nf --input samplesheet.csv -profile slurm
nextflow run main.nf --input samplesheet.csv -profile conda

# Centralized database layout: one parameter instead of eight
nextflow run main.nf --input samplesheet.csv --db_dir /path/to/databases
```

Profiles can be combined (e.g. `-profile test,docker -stub-run`).

## Input & output

### Input

`--input` points to a samplesheet CSV (one sample per line; leave `fastq_2`
empty for single-end):

```csv
sample,fastq_1,fastq_2,group,batch,host
S01,/data/reads/S01_R1.fastq.gz,/data/reads/S01_R2.fastq.gz,case,batch01,human
S02,/data/reads/S02_R1.fastq.gz,/data/reads/S02_R2.fastq.gz,control,batch01,human
```

| Column | Meaning |
| ------ | ------- |
| `sample` | unique sample id (carried as `meta.id`) |
| `fastq_1` | path to forward reads (R1), required |
| `fastq_2` | path to reverse reads (R2); leave empty for single-end |
| `group` | sample group label (e.g. `case` / `control`) |
| `batch` | batch label, carried in metadata |
| `host` | host genome label for host removal (e.g. `human`) |

`single_end` is inferred automatically (empty `fastq_2` → single-end) and added
to the validated output. The `CHECK_SAMPLESHEET` process validates the file,
completes relative paths to absolute, and emits
`00_metadata/validated_samplesheet.csv`. See [docs/input.md](docs/input.md).

### Output

Results are organized per batch under `results/<batch_id>/` in fixed numbered
directories (`00_metadata/` … `14_integrated/`, `99_multiqc/`). The two core
deliverables are:

- `mag_metadata.tsv` — QC / classification / genome statistics / abundance,
  one row per representative MAG.
- `mag_functional_annotation.tsv` — KO / COG / GO / Pathway / ARG, one row per
  gene.

Cross-sample integration adds matrices under `03_taxonomy/combined/` and
`04_function/combined/`; plotting adds a `figures/` subdirectory under the
relevant stage directories (7 PNGs). See
[docs/output.md](docs/output.md) for the per-file inventory.

## Configuration & databases

### Configuration & skip strategy

All parameters live in `nextflow.config`; the complete reference (grouped by
stage, with defaults) is in [docs/parameters.md](docs/parameters.md). Key points:

- **publishDir closures** defer `${params.outdir}/${params.batch_id}/...`
  evaluation so profile/CLI overrides take effect.
- **skip semantics** are transitive: skipping an upstream stage leaves its
  downstream with empty channels (which terminate cleanly). `--skip_mag_qc`
  requires `--skip_dereplication` (dRep scoring consumes QC output).
- **databases** are passed via `--*_db` (or derived from a single `--db_dir`
  tree); see [docs/database.md](docs/database.md).

### Databases

| Branch / stage | Database param | Missing-database behavior |
| --- | --- | --- |
| Kraken2 / Bracken | `--kraken2_db` (`--bracken_db` defaults to it) | warn and skip that branch |
| HUMAnN | `--humann_db` (+ `--metaphlan_db`) | warn and skip that branch |
| CheckM2 | `--checkm2_db` (~3 GB) | error (or `--skip_mag_qc`) |
| GTDB-Tk | `--gtdbtk_db` (R220+, ~110 GB) | error (or `--skip_taxonomy`) |
| DIAMOND | `--diamond_db` (NR `.dmnd`) | error (or `--skip_diamond`) |
| eggNOG-mapper | `--eggnog_db` (~40+ GB) | error (or `--skip_eggnog`) |
| RGI / CARD | `--card_db` (card.json) | error (or `--skip_rgi`) |
| Pathway (optional) | `--pathway_db` (KO→pathway) | column left empty |

dRep, Prodigal and CoverM have no external database dependency. The full
inventory (sizes, acquisition, `--db_dir` layout) is in
[docs/database.md](docs/database.md).

## Testing & reproducibility

### Testing & validation

Two primary regression commands (details and a per-stage matrix in
[docs/workflow.md](docs/workflow.md) and `STATUS.md`):

```bash
# Stub full pipeline — 62/62 tasks, no real tool execution, no database
# dependency; verifies channel topology, parameter guards and output structure
nextflow run main.nf -profile test -stub-run --batch_id demo_stub \
  --kraken2_db test/data/dummy_dbs/kraken2 --bracken_db test/data/dummy_dbs/kraken2 \
  --humann_db test/data/dummy_dbs/humann \
  --metaphlan_db test/data/dummy_dbs/humann/metaphlan \
  --checkm2_db /fake --gtdbtk_db /fake --diamond_db /fake \
  --eggnog_db /fake --card_db /fake

# Real full pipeline — 34/34 tasks, requires the nf-meta environment; read-based
# and database-dependent stages are skipped (databases unavailable locally)
nextflow run main.nf -profile test --batch_id demo_real \
  --skip_read_based --skip_mag_qc --skip_dereplication \
  --skip_taxonomy --skip_annotation
```

- Test data is a ~46 kb synthetic metagenome (2 paired-end samples, reproducible
  via `test/data/make_test_data.py`).
- Parser/orchestration scripts are covered by unit tests
  (`test/test_merge_read_based.py`, `test/test_plot_results.py`, and per-script
  suites referenced in `STATUS.md`).
- `-stub-run` runs every process stub (parsing processes run their real parser),
  so channel wiring and output schemas are exercised without tools or databases.

### Reproducibility & deployment

- **Version pinning:** `environment.yml` pins every tool (Nextflow 26.04.4 —
  the version all verification actually ran on); each process also
  declares a pinned conda version and container tag.
- **`-resume` idempotence:** batch aggregation uses `collectFile(sort)` +
  `toSortedList` for deterministic task hashes (verified cache hits in stub and
  real runs).
- **Fixed seed:** `metabat2_seed=42` keeps MAG IDs stable.
- **`versions.yml`** is emitted per process and aggregated by MultiQC.
- **No hardcoded paths:** all database/index paths come from parameters; the
  gitignored `conf/local.config` is the only place for personal local paths.
- **Profiles:** local (no profile, tools from PATH) / `conda` / `docker` /
  `singularity` / `slurm` / `test` — see [conf/](conf/) and
  [docs/architecture.md](docs/architecture.md).

## Limitations

- **Database-dependent branches not verified with real runs:** read-based
  (Kraken2 standard DB, tens of GB), CheckM2 (~3 GB), dRep, GTDB-Tk (~110 GB),
  and the three annotation branches (NR / eggNOG / CARD) were verified by
  stub-run only; no results were fabricated.
- **Cross-sample merge logic and plotting** (except the CoverM abundance
  heatmap) are verified with stub/synthetic tables, pending real databases.
- **metaSPAdes** has not been run on real data (MEGAHIT has).
- **Container paths not verified on real hardware:** the development machine has
  no Docker/Apptainer/Singularity/sbatch; image-tag existence and profile-stack
  syntax were checked, actual pulls/in-container runs are pending.
- **Test scale:** 46 kb synthetic data, 2 samples, 1 bin — it verifies pipeline
  correctness and numerical cross-consistency, not biological conclusions.
- **Not yet implemented:** `assembly_mode = coassembly` (explicit error;
  channels already reserve `meta.assembly_mode` / `meta.samples`) and binners
  other than MetaBAT2 (`binner` parameter reserved for V2).

## Tools & citations

This workflow orchestrates third-party tools; please cite their original
publications when using results in your own work. No new methods are introduced.

| Tool | Role | Primary citation |
| ---- | ---- | ---------------- |
| [FastQC](https://github.com/s-andrews/FastQC) | read QC | Andrews, 2010 |
| [fastp](https://github.com/OpenGene/fastp) | read trimming / QC | Chen et al., 2018, Bioinformatics |
| [Bowtie2](https://github.com/BenLangmead/bowtie2) | host removal + mapping | Langmead & Salzberg, 2012, Nat. Methods |
| [samtools](https://github.com/samtools/samtools) | BAM handling | Danecek et al., 2021, GigaScience |
| [Kraken2](https://github.com/DerrickWood/kraken2) | taxonomic classification | Wood, Lu & Langmead, 2019, Genome Biol. |
| [Bracken](https://github.com/jenniferlu717/Bracken) | abundance estimation | Lu et al., 2017, PeerJ CS |
| [HUMAnN](https://github.com/biobakery/humann) | functional profiling | Franzosa et al., 2018, Nat. Methods |
| [MEGAHIT](https://github.com/voutcn/megahit) | assembly | Li et al., 2015, Bioinformatics |
| [metaSPAdes](https://github.com/ablab/spades) | assembly | Nurk et al., 2017, Genome Res. |
| [QUAST](https://github.com/ablab/quast) | assembly QC | Gurevich et al., 2013, Bioinformatics |
| [MetaBAT2](https://bitbucket.org/berkeleylab/metabat) | binning + depth | Kang et al., 2019, PeerJ |
| [CheckM2](https://github.com/chklovski/CheckM2) | MAG QC | Chklovski et al., 2023, Nat. Methods |
| [dRep](https://github.com/MrOlm/drep) | dereplication | Olm et al., 2017, ISME J |
| [GTDB-Tk](https://github.com/Ecogenomics/GTDBTk) | MAG taxonomy | Chaumeil et al., 2022, Bioinformatics |
| [Prodigal](https://github.com/hyattpd/Prodigal) | gene prediction | Hyatt et al., 2010, BMC Bioinformatics |
| [DIAMOND](https://github.com/bbuchfink/diamond) | protein alignment | Buchfink et al., 2021, Nat. Methods |
| [eggNOG-mapper](https://github.com/eggnogdb/eggnog-mapper) | functional annotation | Cantalapiedra et al., 2021, Mol. Biol. Evol. |
| [RGI](https://github.com/arpcard/rgi) / CARD | antibiotic resistance | Alcock et al., 2023, Nucleic Acids Res. |
| [CoverM](https://github.com/wwood/CoverM) | MAG abundance | Aroney et al., 2024, Bioinformatics |
| [MultiQC](https://github.com/MultiQC/MultiQC) | report aggregation | Ewels et al., 2016, Bioinformatics |

## Documentation

| Document | Contents |
| --- | --- |
| [docs/architecture.md](docs/architecture.md) | layered architecture, design decisions, workflow DAG |
| [docs/workflow.md](docs/workflow.md) | per-stage flows, mermaid diagram, skip-interaction semantics |
| [docs/channels.md](docs/channels.md) | channel contracts and operator patterns |
| [docs/input.md](docs/input.md) | input format and validation |
| [docs/output.md](docs/output.md) | per-file output inventory |
| [docs/database.md](docs/database.md) | database inventory, sizes, `--db_dir` layout |
| [docs/parameters.md](docs/parameters.md) | full parameter reference (grouped by stage) |

## License

Released under the MIT License for learning and teaching purposes only; see [LICENSE](LICENSE).
