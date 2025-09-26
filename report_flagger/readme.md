# File QC Report Flagger 

Extracts key sequencing QC metrics from your vault of per-file QC artifacts into a **flat CSV** that’s easy to query, pivot, and aggregate (e.g., per dataset). It understands **FASTQ**, **BAM/CRAM**, and **VCF** QC outputs and can **optionally** decrypt and parse encrypted `stats.txt.openssl.gz` files to compute BAM/CRAM GC%.


## Why this exists

Your QC pipeline produces multiple artifacts per EGAF file (gzipped JSONs, FastQC ZIPs, encrypted stats.txt files). When you want to analyze quality trends across **millions of files** or summarize per dataset, you first need a tidy table of metrics. This tool does exactly that:
* walks the standard archive layout by EGAF id
* reads **only** what’s needed from each artifact
* outputs rows like

```
EGAF00001234567,bamcram,mapq,17.42
```

You can then threshold, pivot, and score in downstream steps.

This script has latter been used as the basis for the report_parser script.

## Features

* **FASTQ** (`stdin_fastqc.zip`)
  * %GC from *Basic Statistics* (`gc_content`)
  * duplicate reads % from *Sequence Duplication Levels* (`duplicate_reads`)
  * fraction of reads with mean per-sequence Q≥30 as a proxy quality metric (`quality_reads`)
* **BAM/CRAM** (`*_report.json.gz`)
  * fraction of **alignments** with MAPQ <30 from histogram (`mapq`)
  * duplicate reads % (`duplicate_reads`)
  * % unaligned (from MappedReads) (`unaligned`)
  * optional **GC%** average from encrypted `stats.txt.openssl.gz` (`gc_content`) when `--include-crypt`
* **VCF** (`*_report.json.gz`)
  * Ts/Tv ratio (`tstv_ratio`)
  * average variant QUAL (`avg_qual`)
* **Fast, scalable I/O**
  * process-fork parallelism with configurable workers
  * batch scheduling and periodic GC to keep memory steady
* **No plaintext on disk**
  * encrypted `stats.txt.openssl.gz` is decrypted **in-memory only**


## Archive layout assumptions

The script expects the EGA-like layout rooted at:

```
BASE_PATH = /gpfs/projects/slc00/SL/slc00474/execution-qc/vault/archive
```

Per-file artifacts live under:

```
BASE_PATH/EGAFxxxxxxx/yyy/zzz/execution/
  ├── EGAF..._report.json.gz         # QC JSON (BAM/CRAM or VCF)
  ├── stdin_fastqc.zip               # FASTQ FastQC report
  ├── stats.txt.openssl.gz           # encrypted samtools stats (optional)
  └── key                            # passphrase file for stats decryption (optional)
```

If both `*_report.json.gz` and `stdin_fastqc.zip` exist, the JSON decides the filetype (VCF if it has VCFVersion; otherwise BAM/CRAM). If only the FastQC ZIP exists, it’s treated as FASTQ.

## Installation

* **Python**: 3.8+ (tested on 3.9)
* **Python libs**: standard library plus `cryptography`

```bash
pip install cryptography
```

* Place the script in your repo (e.g., `report_flagger_extract.py`).

No other dependencies or external binaries are required.

## Usage

```bash
python report_flagger_extract.py \
  --file egaf_ids.txt \
  --output qc_metrics.csv \
  --threads 32
```

### CLI arguments

```
--egaf EGAF_ID            analyze a single EGAF id
--file PATH               newline-separated list of EGAF ids (mutually exclusive with --egaf)
--output PATH             CSV output path (will be created/overwritten)
--threads N               number of worker processes (0 = serial, default: 4)
--include-crypt           also decrypt & parse stats.txt.openssl.gz for BAM/CRAM GC%
```

* Use `--threads 0` to run serially (useful for debugging).
* Omit `--include-crypt` unless you specifically need BAM/CRAM GC% from encrypted stats, as it adds CPU cost and can reduce throughput.
* Use `--threads 0` to run serially (useful for debugging).
* Omit `--include-crypt` unless you specifically need BAM/CRAM GC% from encrypted stats, as it adds CPU cost and can reduce throughput.

## Examples

**Single EGAF:**
```bash
python report_flagger_extract.py \
  --egaf EGAF00005001953 \
  --output one.csv \
  --threads 0
```

**Batch via file list:**
```bash
# ids.txt contains one EGAF per line
python report_flagger_extract.py \
  --file ids.txt \
  --output qc_metrics.csv \
  --threads 48
```

**Include encrypted GC% for BAM/CRAM:**
```bash
python report_flagger_extract.py \
  --file ids.txt \
  --output qc_with_gc.csv \
  --threads 32 \
  --include-crypt
```
## Output schema

A single CSV with **one row per metric per file**:

| column   | type   | description                                                              |
|:--------:|:------:|---------------------------------------------------------------------------|
| `egaf_id`  | string | EGAF identifier                                                          |
| `filetype` | enum   | `fastq` \| `bamcram` \| `vcf` \| `__error__`                             |
| `flag_key` | string | metric key (see below)                                                   |
| `value`    | float  | numeric value (units noted below)                                        |

### Keys by filetype

* **FASTQ**
  * `gc_content` — % GC (0–100) from FastQC *Basic Statistics*
  * `duplicate_reads` — % duplicates (0–100)
  * `quality_reads` — % of reads whose **per-sequence mean Q** ≥ 30 (0–100)
* **BAM/CRAM**
  * `mapq` — fraction (%) of **alignments** with MAPQ < 30 (0–100)
  * `duplicate_reads` — % duplicates (0–100)
  * `unaligned` — % reads unaligned = 100 − `MappedReads[0]`×100 (0–100)
  * `gc_content` — **optional** average GC% computed from encrypted stats.txt (0–100)
* **VCF**
  * `tstv_ratio` — transitions/transversions ratio (≈1.5–3.3 for human WGS/WES)
  * `avg_qual` — mean VCF QUAL across records (scale/tool dependent)

### Error sentinel rows

If nothing recognizable is found for an EGAF:
```
filetype="__error__", flag_key="no_qc", value=-1
```
*(no JSON or FastQC ZIP present)*

If the worker crashes unexpectedly:
```
filetype="__error__", flag_key="worker_exception" or future_exception, value=-1
```

Error rows are rare and are included so you can audit coverage (e.g., how many EGAFs had no QC artifacts).

## Interpreting metrics (suggested thresholds)

This extractor **does not** apply thresholds; it only emits values. Common cutoffs I’ve been using:

* **FASTQ**
  * duplicates > **30%**
  * GC outside **35–60%**
  * reads with mean Q < 20 > **30%**  
  *(our `quality_reads` is % with mean Q ≥ 30; derive %<20 in downstream analysis if needed)*
* **BAM/CRAM**
  * fraction of alignments with MAPQ < 30 > **20%**
  * duplicates > **30%**
  * GC outside **35–60%**
* **VCF**
  * Ts/Tv outside **1.5–3.3**
  * average QUAL < **30**

## Performance & scaling

* **Parallelism**: `ProcessPoolExecutor` with `--threads N` (workers = N).
* **Batching**: the driver submits work in mini-batches (`2 * threads`) to avoid flooding the executor queue and to bound memory.
* **GC**: explicit `gc.collect()` after each batch to keep RSS flatter on very large runs.
* **Crypto**: decryption is CPU-bound and single-threaded inside a process. Use `--include-crypt` only when needed.
* **I/O**: many tiny files on shared storage = latency. For maximum throughput:
  * run where the archive is **locally mounted** or on a fast network FS,
  * consider splitting the EGAF list and running **several jobs** in parallel.

## Security & privacy

* `stats.txt.openssl.gz` is **never** written to disk in plaintext.
* Decryption happens entirely in memory using OpenSSL-compatible AES-256-CBC with salted header (`Salted__`) and `EVP_BytesToKey` (MD5).
* The passphrase is read from the `key` file in the same directory.

## Troubleshooting

### “Salted__ header missing”

The encrypted payload didn’t start with the OpenSSL magic. Causes:
* file isn’t base64-encoded OpenSSL output
* wrong key file
* file corrupted

Result: `gc_content` for that file is skipped; other metrics still emit if available.

### UTF-8 decode errors from encrypted stats

The decrypted bytes didn’t parse as text (`stats.txt`); same root causes as above.

### Only `no_qc` rows emitted

That EGAF directory lacked both the JSON report and the FastQC ZIP. Verify archive layout and `BASE_PATH`.

### Slow run times
* Increase `--threads` if CPU-bound.
* Avoid `--include-crypt` unless necessary.
* Ensure storage isn’t the bottleneck; co-locate compute near data if possible.
* Split the workload and run multiple jobs in parallel.

### “Killed” / “bus error”

Likely OOM or a node limit. Reduce `--threads`, split input, or run on a node with more RAM.

## How metrics are computed (implementation notes)

* **FASTQ**
  * `gc_content`: from `>>Basic Statistics` line `%GC\tNN`.
  * `duplicate_reads`: `100 − “Total Deduplicated Percentage”`.
  * `quality_reads`: parses `>>Per sequence quality scores`, computes % of reads with **mean Q ≥ 30** (note: you can convert to “% with mean Q < 20” later if that’s your downstream cutoff).

* **BAM/CRAM**
  * `mapq`: reads `Data.MappingQualityDistribution` (list of `[MAPQ, count]`), sums counts ≤29 as numerator, divides by total counts (denominator). This measures **fraction of alignments** with low MAPQ and avoids >100% artifacts.
  * `unaligned`: uses `Data.MappedReads[0]` (fraction mapped), so `%unaligned = 100 − mapped*100`.
  * `duplicate_reads`: `Data.Duplicates[0] * 100`.
  * `gc_content` (optional): decrypts `stats.txt.openssl.gz` → parses `GCF/GCL` rows → **weighted mean** of GC% by count.

* **VCF**
  * `tstv_ratio`: `Data.TsTvRatio`.
  * `avg_qual`: `AvgQuality` (top-level) or fallback to `Data.AvgQuality` if present.

## Extending

To add a metric:
1. Add logic to the relevant analyzer (`analyze_fastq`, `analyze_bamcram`, `analyze_vcf`) returning `(flag_key, value)`.
2. Keep `flag_key` short and snake_case (e.g., `avg_depth`).
3. The worker already emits whatever the analyzer returns; the CSV schema remains unchanged.

## Example end-to-end

**Input `ids.txt`:**
```
EGAF00005001953
EGAF00005001965
EGAF00005001437
```

**Run:**
```bash
python report_flagger_extract.py --file ids.txt --output qc.csv --threads 16
```

**Sample `qc.csv`:**
```csv
egaf_id,filetype,flag_key,value
EGAF00005001953,bamcram,mapq,27.28
EGAF00005001953,bamcram,duplicate_reads,20.50
EGAF00005001953,bamcram,unaligned,12.84
EGAF00005001965,fastq,quality_reads,45.13
EGAF00005001965,fastq,gc_content,22.00
EGAF00005001965,fastq,duplicate_reads,54.87
EGAF00005001437,fastq,duplicate_reads,38.56
EGAF00005001437,fastq,gc_content,41.00
EGAF00005001437,fastq,quality_reads,92.10
```

## Code structure (for maintainers)

```
report_flagger_extract.py
├─ get_file_paths()               # resolves per-EGAF artifact paths
├─ decrypt_stats_file()           # OpenSSL-compatible in-memory decrypt
├─ analyze_bamcram()              # parses JSON + optional stats text
├─ analyze_fastq()                # parses fastqc_data.txt inside ZIP
├─ analyze_vcf()                  # parses VCF JSON
├─ process_file()                 # worker: detect type, call analyzer, emit rows
├─ read_egaf_file(), chunked()    # input helpers
└─ main()                         # cli, batching, process pool, csv writer
```

* Crypto libs are imported once at module load (cheap in workers created via fork).
* Decryption is **optional** and only attempted when `--include-crypt` and both stats and key exist.
* Errors while reading a metric are swallowed for that metric; other metrics still emit.


## License & acknowledgments

* Internal tool for QC analytics. Adjust `BASE_PATH` to your environment if needed.
* Thanks to the pipeline authors producing the JSON & FastQC artifacts this tool consumes.
