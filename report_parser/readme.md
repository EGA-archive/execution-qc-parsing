# File QC Report Parser

Extracts per-file QC metrics from EGAF-style archive directories into a tidy CSV.

This python script parses key sequencing QC metrics from the execution directory of per-file QC artifacts into a **flat CSV** that’s easy to query, pivot, and aggregate (e.g., per dataset). It understands **FASTQ**, **BAM/CRAM**, and **VCF** QC outputs and can **optionally** decrypt and parse encrypted stats.txt.openssl.gz files to compute BAM/CRAM GC%.

## Why this exists

The File QC Report pipeline produces multiple artifacts per EGAF file (gzipped JSONs, FastQC ZIPs, encrypted stats.txt files). When you want to analyze quality trends across **millions of files** or summarize per dataset, you first need a tidy table of metrics. This tool does exactly that:
* walks the standard archive layout by EGAF id
* reads **only** what’s needed from each artifact
* outputs rows like

```
⠀EGAF00001234567,bamcram,mapq,17.42
```
You can then threshold, pivot, and score in downstream steps.

## Contents

- [Features](#features)
- [What it reads](#what-it-reads)
- [Output schema](#output-schema)
- [Requierements](#requierements)
- [Usage](#usage)
- [Examples](#examples)
- [Details & assumptions](#details--assumptions)
- [Performance notes](#performance-notes)
- [Error handling](#error-handling)
- [Security & privacy](#security--privacy)
- [FAQ](#faq)
- [Acknowledgments](#acknowledgments)

## Features

- **Supports three file types**: FASTQ, BAM/CRAM, and VCF.
- **Outputs a tidy CSV** with one row per metric per file:
  - `egaf_id, filetype, flag_key, value`
- **HPC-friendly**: optional multi-process execution; batch submission ready.
- **Crypto isolated**: optional `--include-crypt` enables GC extraction from encrypted `stats.txt.openssl.gz`; crypto is isolated and skipped otherwise.
- **Robust parsing**: defensive I/O and lenient parsers keep runs going even with partially corrupted QC files.

## What it reads

The script expects the standard EGAF-like directory layout under a fixed base path:

```
BASE_PATH/<egaf[:9]>/<egaf[9:12]>/<egaf[12:15]>/execution/
    ├── {EGAF}_report.json.gz         # BAM/CRAM QC JSON or VCF QC JSON
    ├── stdin_fastqc.zip              # FastQC report
    ├── stats.txt.openssl.gz          # optional, encrypted samtools stats (BAM/CRAM)
    └── key                           # passphrase for the encrypted stats file
```

Base path in this repo’s code:

```
BASE_PATH = Path("/gpfs/projects/slc00/SL/slc00474/execution-qc/vault/archive")
```

> If the archive lives elsewhere, update that constant or overlay a symlink.


If both *_report.json.gz and stdin_fastqc.zip exist, the JSON decides the filetype (VCF if it has VCFVersion; otherwise BAM/CRAM). If only the FastQC ZIP exists, it’s treated as FASTQ.


## Output schema

CSV columns:
- `egaf_id` – EGAF identifier (string)
- `filetype` – one of `fastq`, `bamcram`, `vcf`, or `__error__`
- `flag_key` – metric key (see below)
- `value` – numeric value for the metric (float). For error rows this is `-1`

### Metric keys by file type

**FASTQ** (from `stdin_fastqc/fastqc_data.txt`):
- `duplicate_reads` – percent duplicate reads (100 − *Total Deduplicated %*)
- `gc_content` – `%GC` from “Basic Statistics”
- `quality_reads` – fraction of reads whose **per‑sequence mean Q ≥ 30** (derived from *Per sequence quality scores*)

**BAM/CRAM** (from `{EGAF}_report.json.gz` and optional decrypted `stats.txt`):
- `unaligned` – percent reads unaligned = 100 − (Data.MappedReads[0] × 100)
- `mapq` – **fraction of alignments** with MAPQ \< 30, computed from `Data.MappingQualityDistribution` (sum of counts with q ≤ 29 divided by total counts)
- `duplicate_reads` – percent duplicate reads from `Data.Duplicates[0] × 100`
- `gc_content` – **average GC%** computed from decrypted `stats.txt` histograms (`GCF` and `GCL`)
- `insert_size` – insert/fragment size in base pairs from top-level JSON field `InsertSize` (reported value is taken as-is)

**VCF** (from `{EGAF}_report.json.gz` containing `VCFVersion`):
- `tstv_ratio` – `Data.TsTvRatio`
- `avg_qual` – `AvgQuality` (top-level; some reports also mirror this under `Data`)

**Error sentinel**
- filetype = `__error__`, flag_key = one of:
  - `no_qc` – no recognizable QC file found for the EGAF
  - `worker_exception` – unexpected failure in the worker
  - `future_exception` – failure raised by the executor future


## Requierements

Python 3.8+ is requiered (tested on 3.9)

```bash
module load python/3.9.10
```

**Python libs**: standard library plus cryptography
  
'cryptography' python module is requiered, if missing you may contact support@bsc.es for them to install it on the specific python version you will use

> If you will **not** use encrypted stats (`--include-crypt`), the script still imports `cryptography` but won’t use it unless asked. 


## Usage

```bash
python3.9 report_parser.py \
  --file egaf_ids.txt \
  --output metrics.csv \
  --threads 32 \
  --include-crypt
```

### Arguments

- `--egaf EGAF_ID`  
  process a single egaf id

- `--file PATH`  
  newline-separated text file with many EGAF ids (one per line)

- `--output PATH` **(required)**  
  output CSV path

- `--threads N` *(default: 4)*  
  number of worker processes. use `0` for serial execution.  
  tip: start with `N = min(32, n_cores / 2)` and tune

- `--include-crypt`  
  enable decryption/parsing of `stats.txt.openssl.gz` to compute **BAM/CRAM GC%**.  
  this significantly increases runtime; leave off if GC% is not needed.


## Examples

### Single file (FASTQ)
```bash
python3.9 report_parser.py \
  --egaf EGAF00001234567 \
  --output out.csv \
  --threads 0
```

### Many EGAFs with parallel workers
```bash
python3.9 report_parser.py \
  --file all_egafs.txt \
  --output metrics.csv \
  --threads 48
```

### Enable encrypted GC% for BAM/CRAM
```bash
python3.9 report_parser.py \
  --file all_bams.txt \
  --output bams_with_gc.csv \
  --threads 24 \
  --include-crypt
```

### Sample CSV output
```csv
egaf_id,filetype,flag_key,value
EGAF00005001437,fastq,duplicate_reads,38.56
EGAF00005001965,fastq,quality_reads,54.87
EGAF00005001965,fastq,gc_content,22.00
EGAF00005001953,bamcram,mapq,27.28
EGAF00002250174,vcf,tstv_ratio,2.00
EGAF00002250174,vcf,avg_qual,344.90
```

## Details & assumptions

### FASTQ
- `%GC` is taken directly from the “Basic Statistics” module (`%GC` line).
- `duplicate_reads` is computed as `100 − Total Deduplicated %` from *Sequence Duplication Levels*.
- `quality_reads` evaluates the *Per sequence quality scores* block:
  - the FastQC table has rows of `mean_quality<TAB>count`
  - we sum counts with `mean_quality ≥ 30` and divide by total counts to get a percentage of high‑quality reads.

### BAM/CRAM
- Unaligned reads fraction comes from `Data.MappedReads[0]` (a ratio in [0,1]) → `100 − ratio×100`.
- Low MAPQ fraction uses the **alignment distribution**, not `TotalReads`:
  - total = sum of counts in `MappingQualityDistribution`
  - low = sum of counts with `q ≤ 29`
  - metric = `low / total × 100`  
  This avoids mis-scaling (e.g., when mates/secondary alignments inflate counts).
- `duplicate_reads` uses `Data.Duplicates[0] × 100` when present.
- Average GC% is derived from decrypted `stats.txt` by weighted mean over `GCF` and `GCL` histograms:
  - each line `GCX <percent> <count>` contributes `percent × count`
  - result = sum / total_count
- Insert size:
  - Read from the top-level JSON field InsertSize (not under Data).
	- Reported in base pairs; the script emits the numeric value unchanged as insert_size.
	- 	If the field is missing or non-numeric, the metric is omitted for that EGAF.

### VCF
- Ts/Tv ratio read from `Data.TsTvRatio` when present.
- Average QUAL read from top-level `AvgQuality` (most VCF qc JSONs expose it there).


## Performance notes

- **I/O** is the bottleneck. Keep `BASE_PATH` on a high‑throughput filesystem.
- Use **moderate parallelism**: `--threads 16..48` typically scales well; more than that can cause filesystem contention.
- The script batches futures per `threads × 2` to keep workers fed without overloading memory.
- `--include-crypt` is **computing expensive**. Only use it when you need BAM/CRAM GC%.


## Error handling

- If **no QC file** exists for an EGAF, a single row is written:
  - `egaf_id,__error__,no_qc,-1`
- Unexpected worker failures are recorded as:
  - `egaf_id,__error__,worker_exception,-1` or `future_exception`
- Decryption errors are **not fatal**; the file’s other metrics are still emitted (without GC%).


## Security & privacy

- Encrypted stats files are decrypted **in-memory** only; no plaintext is written to disk. The passphrase is read from the key file in the same directory.
- The decryption routine uses OpenSSL-compatible EVP_BytesToKey (MD5) + AES‑CBC and strips PKCS#7 padding.
- The `--include-crypt` flag is opt‑in to avoid unnecessary key handling in most runs.


## FAQ

**Q: Why are `mapq` percentages sometimes high?**  
A: The metric is the *fraction of alignments* with MAPQ \< 30, not the fraction of reads. When secondary/multi-mapped alignments are present the denominator is larger than read count, by design. This makes it more conservative for QC.

**Q: What if a FastQC module is missing?**  
A: The script skips that metric for that EGAF; no row is written for that key.

**Q: How do I change the base path?**  
A: Edit the `BASE_PATH` constant at the top of the script or create a symlink matching the expected layout.

**Q: Can I merge this CSV with dataset membership?**  
A: Yes—join on `egaf_id` with your dataset mapping to build dataset-level quality scores.

**Q: “Salted__ header missing”**  
A: The encrypted payload didn’t start with the OpenSSL magic. Causes:
* file isn’t base64-encoded OpenSSL output
* wrong key file
* file corrupted
Result: gc_content for that file is skipped; other metrics still emit if available.

**Q: UTF-8 decode errors from encrypted stats**  
A: The decrypted bytes didn’t parse as text (stats.txt); same root causes as above.

**Q: Only no_qc rows emitted**  
A: That EGAF directory lacked both the JSON report and the FastQC ZIP. Verify archive layout and BASE_PATH.

**Q: Slow run times** 
A: * Increase --threads if CPU-bound.
* Avoid --include-crypt unless necessary.
* Split the workload and run multiple jobs in parallel.

**Q: “Killed” / “bus error”** 
A: Likely OOM or a node limit. Reduce --threads, split input, or run on a node with more RAM.

## How metrics are computed (implementation notes)
* **FASTQ**
  * gc_content: from >>Basic Statistics line %GC\tNN.
  * duplicate_reads: 100 − “Total Deduplicated Percentage”.
  * quality_reads: parses >>Per sequence quality scores, computes % of reads with **mean Q ≥ 30** (note: you can convert to “% with mean Q < 20” later if that’s your downstream cutoff).

* **BAM/CRAM**
  * mapq: reads Data.MappingQualityDistribution (list of [MAPQ, count]), sums counts ≤29 as numerator, divides by total counts (denominator). This measures **fraction of alignments** with low MAPQ and avoids >100% artifacts.
  * unaligned: uses Data.MappedReads[0] (fraction mapped), so %unaligned = 100 − mapped*100.
  * duplicate_reads: Data.Duplicates[0] * 100.
  * gc_content (optional): decrypts stats.txt.openssl.gz → parses GCF/GCL rows → **weighted mean** of GC% by count.
  * insert_size: ins = data.get("InsertSize"); if ins is int/float, emit (flag_key="insert_size", value=float(ins)).

* **VCF**
  * tstv_ratio: Data.TsTvRatio.
  * avg_qual: AvgQuality (top-level) or fallback to Data.AvgQuality if present.

## Extending

To add a metric:
1. Add logic to the relevant analyzer (analyze_fastq, analyze_bamcram, analyze_vcf) returning (flag_key, value).
2. Keep flag_key short and snake_case (e.g., avg_depth).
3. The worker already emits whatever the analyzer returns; the CSV schema remains unchanged.


**Code structure (for maintainers)**
```
report_parser.py
├─ get_file_paths()               # resolves per-EGAF artifact paths
├─ decrypt_stats_file()           # OpenSSL-compatible in-memory decrypt
├─ analyze_bamcram()              # parses JSON + optional stats text
├─ analyze_fastq()                # parses fastqc_data.txt inside ZIP
├─ analyze_vcf()                  # parses VCF JSON
├─ process_file()                 # worker: detect type, call analyzer, emit rows
├─ read_egaf_file(), chunked()    # input helpers
└─ main()                         # CLI, batching, process pool, CSV writer
```

* Crypto libs are imported once at module load (cheap in workers created via fork).
* Decryption is **optional** and only attempted when --include-crypt and both stats and key exist.
* Errors while reading a metric are swallowed for that metric; other metrics still emit.

## Acknowledgments
* Internal tool for QC analytics. Adjust BASE_PATH to your environment if needed.
* Thanks to the pipeline authors producing the JSON & FastQC artifacts this tool consumes.
