# QC Error parser, for extracting error execution stats

An ultra-fast, I/O-friendly Python utility to scan **millions** of EGAF entries in the N4 execution directory, detect **missing QC reports**, and infer likely causes from `pipeErrors` logs. Designed based on the error_flagger script, for HPC clusters (SLURM), with bounded memory, robust parallelism, and extensible error pattern matching.


## Features

- **Scales to Millions**: Processes very large EGAF lists with stable memory usage.  
- **HPC-Optimized**: Thread-pooled I/O, batched submissions, a dedicated writer thread.  
- **Two Outputs**:
  - `--output`: Per-EGAF detailed error excerpts, ready for triage.
  - `--missing-error-output`: EGAF IDs that have **missing QC** and **no** error logs, one per line.
- **Broad Error Coverage**:
  - `samtools`/HTSlib assertions, truncation, “unsorted”, “failed to open”, format errors.
  - `FastQC` memory range violations and per-tile anomalies.
  - System/runtime issues (OOM, I/O errors, permission problems, disk space, etc.).
- **Extensible**: Add new regex patterns safely and quickly.
- **Clear Summary**: Terminal output shows high-level counts and error-file breakdowns.
- **Zero Dependencies**: Standard library only.

## Requirements

- Python ≥ 3.8  
- Access to the archive tree (update if needed in `error_parser.py`):
  ```
  /slgpfs/projects/slc00/slc00474/execution-qc/vault/archive
  ```
## What it reads

The script expects the standard EGAF-like directory layout under a fixed base path:

```
BASE_PATH/<egaf[:9]>/<egaf[9:12]>/<egaf[12:15]>/execution/
    ├── {EGAF}_report.json.gz         # BAM/CRAM QC JSON or VCF QC JSON
    ├── stdin_fastqc.zip              # FastQC report
    ├── stats.txt.openssl.gz          # optional, encrypted samtools stats (BAM/CRAM)
    └── ...
```

Base path in this repo’s code:

```
BASE_PATH = Path("/gpfs/projects/slc00/SL/slc00474/execution-qc/vault/archive")
```

> If the archive lives elsewhere, update that constant or overlay a symlink.



## Usage

```bash
./error_parser.py \
  --file egaf_list.txt \
  --output results.txt \
  --missing-error-output missing_error.txt \
  --threads 32 \
  --max-threads 32 \
  --batch-size 1000
```

### Arguments

- `--file` *(string)*: Path to a text/CSV/TSV file with one EGAF ID per line.  
  *(Alternatively, use `--egaf EGAF0000...` to analyze a single ID.)*
- `--output` *(string, required)*: Per-EGAF detailed error excerpts file.  
- `--missing-error-output` *(string, required)*: EGAF IDs with **missing QC** and **no** error logs; one ID per line.  
- `--threads` *(int, default: 16)*: Worker threads for I/O-bound processing.  
- `--max-threads` *(int, default: 16)*: Safety cap for worker threads (min of the two is used).  
- `--batch-size` *(int, default: 1000)*: EGAF submissions per scheduling batch.

### Input Format

`egaf_list.txt`:
```
EGAF00000000001
EGAF00000000002
EGAF00000000003
...
```

## Expected Outputs

### 1) `--output` (Detailed Excerpts)

Only EGAF entries with detected issues are printed, prefixed with `>> EGAF...`:

```
>> EGAF00008XXXXXX
Error in index.e:
[W::hts_close] EOF marker is absent. The input is probably truncated.

Error in unzip.e:
unzip:  cannot find or open /out/xxxxx_fastqc.zip, /out/xxxxx_fastqc.zip.zip or /out/xxxxx_fastqc.zip.ZIP.
```

### 2) `--missing-error-output` (IDs Only)

EGAF IDs with **missing QC reports** and **no error lines** in any `pipeErrors/*.{e,o}`:

```
EGAF00003456789
EGAF00002345678
...
```

### Terminal Summary (Example)

```
Finished checking 1,387,016 files.
1,247,594 files (89.9%) have a missing QC report.
10,502 files have error logs indicating reasons for missing QC reports.
Counts of files per error log file:
- validateSummary.e: 1,478 files
- validateProperPairSummary.e: 971 files
- uncompress.e: 720 files
- fastqc.e: 2,204 files
- index.e: 6,945 files
- stats.e: 11 files
- fastqc.o: 25 files
- unzip.e: 4 files
- finalize.e: 9 files
Total errors detected across all files: 19,386
X files have missing QC report and no error logs. EGAF IDs are written to missing_error.txt

Time elapsed: 1234.56 seconds.
```

## Directory Conventions

For an EGAF ID `EGAF00...`, the script infers:

```
BASE_PATH / EGAF00000 / 123 / 456 / execution /
  ├─ EGAF00000abcdef_report.json.gz      # bam/cram/vcf QC artifact (if present)
  ├─ stdin_fastqc.zip                    # fastq QC artifact (if present)
  └─ pipeErrors/                         # error logs
       *.e, *.o
```

**File Type Inference**

- Presence of `*_report.json.gz` → BAM/CRAM/VCF QC exists.  
- Presence of `stdin_fastqc.zip` → FASTQ QC exists.  
- If neither exists → **missing QC report**; then `pipeErrors` is scanned.

## Error Detection

The function `is_error_line(line: str) -> bool` tests each line against a **union** of regex patterns.

### Key Families Covered

**SAMtools/HTSlib**
- `bgzf_read\(\) on reference file: Success`
- `cram/cram_io\.c: Assertion .* failed`
- `samtools index: "-" is in a format that cannot be usefully indexed`
- `[W::hts_close] EOF marker is absent. The input is probably truncated`
- `[E::hts_hopen] fail to open file`
- `[E::hts_open_format] fail to open file`
- `samtools index: failed to open`
- `unsorted positions`, `corrupted or unsorted`, `unexpected end of file`

**FastQC**
- `Too many tiles \(>1000\)`
- `Memory value .* was outside the allowed range`
- `Terminating due to java.lang.OutOfMemoryError`

**System/Runtime**
- `Exception|ERROR|Error|failed`
- `cannot allocate memory`, `No space left on device`
- `Segmentation fault|Bus error|Aborted|core dumped`
- `Input/output error|permission denied|Broken pipe`
- `No such file or directory|Too many open files|Disk quota exceeded`
- And many others.

> Regular FastQC progress messages (e.g., “Approx 5% complete …”) are **ignored** as non-errors. The script specifically flags “Too many tiles (>1000) …” and out-of-range memory, among others.

### Extending Patterns

1. Open `script.py`, find `def is_error_line`.
2. Append new raw-string regexes to `error_patterns`:
   ```python
   r'my new error phrase here'
   ```
3. Prefer specific tool-scoped prefixes (e.g., `samtools:`) to prevent false positives.


## Performance Tuning

### Knobs

- `--threads`: set up to available cores on the node.
- `--batch-size`: number of EGAFs per scheduling batch (default 1000).
  - Lower if you see memory pressure or metadata storms.
  - Higher if node memory and storage throughput allow.

### Tips for 3M+ IDs

- Prefer **`sbatch`** over interactive runs for reliability and scheduling efficiency.
- If the filesystem is hot, reduce `--threads` and `--batch-size`.
- Split the list and run multiple jobs over different input files in parallel.
- Keep outputs on a fast/high-throughput filesystem.


## Troubleshooting

- **“Bus Error” on Very Large Runs**
  - Decrease `--threads` and/or `--batch-size`.
  - Use `sbatch` with higher `--mem`.
  - The script already avoids accumulating futures; keep batches moderate.

- **“No Files Were Checked”**
  - Verify input list formatting (one EGAF per line).
  - Confirm permissions to read the archive path.
  - Ensure `BASE_PATH` matches your environment.

- **False Positives in FastQC**
  - Only specific FastQC anomalies are treated as errors.
  - Tighten or add patterns as needed.

- **Duplicate EGAF Headers**
  - The script prints `>> EGAF…` once per EGAF block. If you added prints, ensure they don’t re-emit IDs.


## FAQ

**How are file types determined?**  
By presence of QC artifacts in `execution/`:  
- `*_report.json.gz` → BAM/CRAM/VCF QC present  
- `stdin_fastqc.zip` → FASTQ QC present  
If neither is present → **missing QC**; the script then scans `pipeErrors`.

**Does the tool parse QC contents?**  
No. The goal is to explain **missing QC** or concurrent errors during QC generation, not to score QC quality.

**Can I add new error patterns?**  
Yes. Append precise regexes to `error_patterns` in `is_error_line`. Prefer tool-scoped, specific phrases.


