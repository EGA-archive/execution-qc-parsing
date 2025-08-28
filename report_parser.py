from __future__ import annotations

import argparse
import base64
import csv
import gzip
import json
import time
import zipfile
from collections import defaultdict
from io import BytesIO
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Tuple
import concurrent.futures
import gc

# heavy crypto libs
from cryptography.hazmat.backends import default_backend
from cryptography.hazmat.primitives import hashes
from cryptography.hazmat.primitives.ciphers import Cipher, algorithms, modes

# ----------------------------------------------------------------------------
# Configuration
# ----------------------------------------------------------------------------
BASE_PATH = Path("/slgpfs/projects/slc00/slc00474/execution-qc/vault/archive")

# ----------------------------------------------------------------------------
# Helper functions
# ----------------------------------------------------------------------------

def _if_exists(p: Path) -> Optional[Path]:
    return p if p.is_file() else None

def get_file_paths(egaf: str) -> Dict[str, Optional[Path]]:
    pre, mid, suf = egaf[:9], egaf[9:12], egaf[12:15]
    base = BASE_PATH / pre / mid / suf / "execution"
    return {
        "json": _if_exists(base / f"{egaf}_report.json.gz"),
        "zip": _if_exists(base / "stdin_fastqc.zip"),
        "stats": _if_exists(base / "stats.txt.openssl.gz"),
        "key": _if_exists(base / "key"),
    }

# ----------------------------------------------------------------------------
# Decryption helpers (OpenSSL EVP_BytesToKey + AES‑CBC)
# ----------------------------------------------------------------------------

def _evp_bytes_to_key(pwd: bytes, salt: bytes, klen: int, ivlen: int) -> Tuple[bytes, bytes]:
    d = b""; dtot = b""
    while len(dtot) < klen + ivlen:
        md5 = hashes.Hash(hashes.MD5(), backend=default_backend())
        md5.update(d + pwd + salt)
        d = md5.finalize()
        dtot += d
    return dtot[:klen], dtot[klen:klen + ivlen]

def decrypt_stats_file(enc: bytes, pwd: bytes) -> str:
    with gzip.GzipFile(fileobj=BytesIO(enc)) as gz:
        b64_payload = gz.read()
    raw = base64.b64decode(b64_payload)
    if not raw.startswith(b"Salted__"):
        raise ValueError("Salted__ header missing")
    salt, payload = raw[8:16], raw[16:]
    key, iv = _evp_bytes_to_key(pwd, salt, 32, 16)
    cipher = Cipher(algorithms.AES(key), modes.CBC(iv), backend=default_backend()).decryptor()
    padded = cipher.update(payload) + cipher.finalize()
    return padded[:-padded[-1]].decode()

# ----------------------------------------------------------------------------
# QC analyzers – each returns list[(flag_key, value)]
# ----------------------------------------------------------------------------

def _mean_gc_from_stats(text: str) -> Optional[float]:
    total = count = 0
    for parts in (l.split() for l in text.splitlines() if l.startswith("GC")):
        if len(parts) == 3:
            gc_pct, n = float(parts[1]), int(parts[2])
            total += gc_pct * n
            count += n
    return total / count if count else None

def analyze_bamcram(jpath: Path, dec_stats: Optional[str]) -> List[Tuple[str, float]]:
    out: List[Tuple[str, float]] = []
    try:
        with gzip.open(jpath, "rt") as fh:
            data = json.load(fh)

        ins = data.get("InsertSize")
        if isinstance(ins, (int, float)):
            out.append(("insert_size", float(ins)))

        d = data.get("Data", {})

        if "MappedReads" in d:
            out.append(("unaligned", 100 - d["MappedReads"][0] * 100))

        if "MappingQualityDistribution" in d:
            dist = d["MappingQualityDistribution"]
            tot = sum(c for _, c in dist)
            low = sum(c for q, c in dist if q <= 29)
            if tot:
                out.append(("mapq", low / tot * 100))

        if "Duplicates" in d:
            out.append(("duplicate_reads", d["Duplicates"][0] * 100))

        if dec_stats:
            gc_val = _mean_gc_from_stats(dec_stats)
            if gc_val is not None:
                out.append(("gc_content", gc_val))

    except Exception:
        pass

    return out

def analyze_fastq(zpath: Path) -> List[Tuple[str, float]]:
    out: List[Tuple[str, float]] = []
    try:
        with zipfile.ZipFile(zpath) as z:
            with z.open("stdin_fastqc/fastqc_data.txt") as fh:
                lines = (l.decode().rstrip() for l in fh)
                it = iter(lines)
                dup = gc_pct = q20 = None
                for ln in it:
                    if ln.startswith(">>Basic Statistics"):
                        # … unchanged …
                        for ln in it:
                            if ln.startswith("%GC"):
                                gc_pct = float(ln.split("\t")[1])
                                break
                            if ln.startswith(">>END_MODULE"):
                                break

                    elif ln.startswith(">>Sequence Duplication Levels"):
                        # … unchanged …
                        for ln in it:
                            if ln.startswith("#Total Deduplicated Percentage"):
                                dup = 100 - float(ln.split("\t")[1])
                                break
                            if ln.startswith(">>END_MODULE"):
                                break

                    elif ln.startswith(">>Per sequence quality scores"):
                        # skip header
                        for ln in it:
                            if ln.startswith("#Quality"):
                                break

                        tot = below30 = 0.0
                        for ln in it:
                            if ln.startswith(">>END_MODULE"):
                                break
                            q, c = map(float, ln.split("\t"))
                            tot += c
                            # count reads with mean Q ≥ 30
                            if q >= 30:
                                below30 += c

                        if tot:
                            # overwrite q20 with % reads mean-Q ≥ 30
                            q20 = below30 / tot * 100

                if dup is not None:
                    out.append(("duplicate_reads", dup))
                if gc_pct is not None:
                    out.append(("gc_content", gc_pct))
                if q20 is not None:
                    out.append(("quality_reads", q20))
    except Exception:
        pass

    return out

def analyze_vcf(jpath: Path) -> List[Tuple[str, float]]:
    out: List[Tuple[str, float]] = []
    try:
        with gzip.open(jpath, "rt") as fh:
            data = json.load(fh)
        d = data.get("Data", {})
        if (ratio := d.get("TsTvRatio")) is not None:
            out.append(("tstv_ratio", ratio))
        if (qual := data.get("AvgQuality", d.get("AvgQuality"))) is not None:
            out.append(("avg_qual", qual))
    except Exception:
        pass
    return out

# ----------------------------------------------------------------------------
# Worker – always returns a list (may include __error__ sentinel rows)
# ----------------------------------------------------------------------------

def process_file(args: Tuple[str, bool]) -> List[Tuple[str, str, str, float]]:
    egaf, include_crypt = args
    rows: List[Tuple[str, str, str, float]] = []
    try:
        p = get_file_paths(egaf)
        dec_stats = None
        if include_crypt and p["stats"] and p["key"]:
            try:
                dec_stats = decrypt_stats_file(p["stats"].read_bytes(), p["key"].read_bytes().strip())
            except Exception:
                pass  # decryption failed – just skip GC metric
        if p["json"]:
            with gzip.open(p["json"], "rt") as fh:
                j = json.load(fh)
            if "VCFVersion" in j:
                for k, v in analyze_vcf(p["json"]):
                    rows.append((egaf, "vcf", k, v))
            else:
                for k, v in analyze_bamcram(p["json"], dec_stats):
                    rows.append((egaf, "bamcram", k, v))
        elif p["zip"]:
            for k, v in analyze_fastq(p["zip"]):
                rows.append((egaf, "fastq", k, v))
        else:
            rows.append((egaf, "__error__", "no_qc", -1))
    except Exception:
        rows.append((egaf, "__error__", "worker_exception", -1))
    return rows

# ----------------------------------------------------------------------------
# Utility: reading and batching EGAF IDs
# ----------------------------------------------------------------------------

def read_egaf_file(path: Path) -> List[str]:
    return [ln.strip() for ln in path.read_text().splitlines() if ln.strip()]

def chunked(seq: List[str], n: int):
    for i in range(0, len(seq), n):
        yield seq[i : i + n]

# ----------------------------------------------------------------------------
# Main
# ----------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description="Extract QC metrics to CSV")
    grp = ap.add_mutually_exclusive_group(required=True)
    grp.add_argument("--egaf")
    grp.add_argument("--file")
    ap.add_argument("--output", required=True)
    ap.add_argument("--threads", type=int, default=4, help="0 = serial")
    ap.add_argument("--include-crypt", action="store_true")
    args = ap.parse_args()

    # Build list of EGAF IDs
    ids = [args.egaf] if args.egaf else read_egaf_file(Path(args.file))

    rows: List[Tuple[str, str, str, float]] = []
    t0 = time.time()

    if args.threads == 0 or len(ids) == 1:
        # Serial path (debug‑friendly)
        for eg in ids:
            rows.extend(process_file((eg, args.include_crypt)))
    else:
        batch_size = args.threads * 2
        with concurrent.futures.ProcessPoolExecutor(max_workers=args.threads) as ex:
            for batch in chunked(ids, batch_size):
                futs = {ex.submit(process_file, (eg, args.include_crypt)): eg for eg in batch}
                for fut in concurrent.futures.as_completed(futs):
                    try:
                        rows.extend(fut.result())
                    except Exception:
                        # record failure for that EGAF
                        rows.append((futs[fut], "__error__", "future_exception", -1))
                gc.collect()

    # Write CSV
    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    with out.open("w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["egaf_id", "filetype", "flag_key", "value"])
        w.writerows(rows)

    print(f"Wrote {len(rows)} rows from {len(ids)} files to {out} in {time.time() - t0:.1f}s")

if __name__ == "__main__":
    main()
