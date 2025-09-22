import argparse
import base64
import concurrent.futures
import gzip
import json
import logging
import time
import zipfile
from collections import defaultdict
from io import BytesIO
from pathlib import Path
from typing import List, Tuple, Dict, Iterator

# Cryptography library is used for decryption
from cryptography.hazmat.primitives.ciphers import Cipher, algorithms, modes
from cryptography.hazmat.primitives import hashes
from cryptography.hazmat.backends import default_backend

# Base path to the archive directory
BASE_PATH = Path("/slgpfs/projects/slc00/slc00474/execution-qc/vault/archive")

def get_file_paths(file_name: str) -> Tuple[Path, Path, Path, Path]:
    """
    Construct paths for JSON, TXT, ZIP, and stats files based on the file name.
    """
    # Extract parts of the file name to construct the directory path
    prefix = file_name[:9]  # e.g., 'egaf00005'
    middle = file_name[9:12]  # e.g., '432'
    suffix = file_name[12:15]  # e.g., '457'
    base_dir = BASE_PATH / prefix / middle / suffix / "execution"
    # Paths to the required files
    json_path = base_dir / f"{file_name}_report.json.gz"
    txt_path = base_dir / "input_screen_v2.txt"
    zip_path = base_dir / "stdin_fastqc.zip"
    stats_encrypted_path = base_dir / "stats.txt.openssl.gz"
    return json_path, txt_path, zip_path, stats_encrypted_path

def evp_bytes_to_key(password: bytes, salt: bytes, key_len: int, iv_len: int) -> Tuple[bytes, bytes]:
    """
    Derive key and IV using OpenSSL's EVP_BytesToKey algorithm with MD5 hash.
    """
    dtot = b''
    while len(dtot) < key_len + iv_len:
        d = hashes.Hash(hashes.MD5(), backend=default_backend())
        d.update(dtot + password + salt)
        dtot += d.finalize()
    key = dtot[:key_len]
    iv = dtot[key_len:key_len+iv_len]
    return key, iv

def decrypt_stats_file(encrypted_data: bytes, password: bytes) -> str:
    """
    Decrypt and decompress the stats file data in memory.
    """
    # Decompress the gzip data to get the base64-encoded encrypted data
    with gzip.GzipFile(fileobj=BytesIO(encrypted_data)) as gz:
        base64_encrypted_data = gz.read()
    # Decode the base64 data to get the raw encrypted data
    encrypted_data = base64.b64decode(base64_encrypted_data)
    # Check if the data starts with 'Salted__', indicating a salt is used
    if encrypted_data[:8] != b'Salted__':
        raise ValueError("Salted magic header not found in encrypted data.")
    salt = encrypted_data[8:16]
    encrypted_payload = encrypted_data[16:]
    # Derive key and IV from the password and salt
    key, iv = evp_bytes_to_key(password, salt, key_len=32, iv_len=16)
    # Decrypt the data using AES-256-CBC
    cipher = Cipher(
        algorithms.AES(key),
        modes.CBC(iv),
        backend=default_backend()
    )
    decryptor = cipher.decryptor()
    decrypted_padded = decryptor.update(encrypted_payload) + decryptor.finalize()
    # Remove PKCS7 padding
    padding_length = decrypted_padded[-1]
    decrypted_data = decrypted_padded[:-padding_length]
    return decrypted_data.decode('utf-8')

def check_species(txt_path: Path, warn_missing_file: bool
                  ) -> Tuple[str, bool, str]:
    """
    Check species information from the species file.
    """
    if not txt_path.exists():
        if warn_missing_file:
            return (
                f"Warning: species file not found at {txt_path}\n",
                True,
                "species file missing",
            )
        return "", False, ""
    try:
        with txt_path.open('r') as file:
            lines = file.readlines()
        if len(lines) >= 3:
            header_line = lines[1].strip()
            headers = header_line.split("\t")
            try:
                species_idx = headers.index('Genome')
                one_hit_idx = headers.index('#One_hit_one_genome')
                percent_one_hit_idx = headers.index('%One_hit_one_genome')
            except ValueError as e:
                return (
                    f"Error: required columns not found: {e}\n",
                    True,
                    "species error",
                )
            max_value = -1
            max_species = ""
            max_percent = 0.0
            for line in lines[2:]:
                cells = line.strip().split("\t")
                if len(cells) <= max(one_hit_idx, percent_one_hit_idx):
                    continue
                species = cells[species_idx]
                try:
                    one_hit = int(cells[one_hit_idx])
                    percent_one_hit = float(cells[percent_one_hit_idx])
                except ValueError:
                    continue
                if one_hit > max_value:
                    max_value = one_hit
                    max_species = species
                    max_percent = percent_one_hit
            if max_species != "Human":
                return (
                    f"Warning: most reads mapped to {max_species} "
                    f"({max_value} reads).\n",
                    True,
                    "sp not human",
                )
            elif max_percent < 5.0:
                return (
                    f"Warning: species is unknown "
                    f"({max_percent:.2f}% reads mapped to Human).\n",
                    True,
                    "species unknown",
                )
            else:
                return "", False, ""
    except Exception as e:
        return (
            f"Error processing species file: {str(e)}\n",
            True,
            "species error",
        )
    return "", False, ""

def analyze_bam_cram(json_path: Path, stats_encrypted_path: Path, per_file_counts: Dict
                     ) -> Tuple[str, bool]:
    """
    Analyze BAM/CRAM QC reports and check GC content using the encrypted stats file.
    """
    output = ""
    warnings_found = False
    if not json_path.exists():
        per_file_counts["bam_cram_warnings"]["QC report missing"] += 1
        output += f"Error: BAM/CRAM QC report not found at {json_path}\n"
        warnings_found = True
        return output, warnings_found
    try:
        with gzip.open(json_path, 'rt') as file:
            data = json.load(file)
        total_reads = data.get("TotalReads", 0)
        data_section = data.get("Data", {})
        # Check 'MappedReads'
        if 'MappedReads' in data_section:
            mapped_ratio = data_section["MappedReads"][0]
            reads_unaligned = 100 - (mapped_ratio * 100)
            if reads_unaligned > 40.0:
                output += (
                    f"Warning: reads unaligned ({reads_unaligned:.2f}%) exceeds 40%.\n"
                    "High percentage of unaligned reads indicates poor sequencing quality, contamination, or issues with the reference genome alignment.\n"
                )
                per_file_counts["bam_cram_warnings"][
                    "% reads unaligned >40"
                ] += 1
                warnings_found = True
        else:
            output += "Warning: 'MappedReads' not found in QC report.\n"
            warnings_found = True
        # Check 'MappingQualityDistribution'
        if 'MappingQualityDistribution' in data_section:
            map_qual_dist = data_section["MappingQualityDistribution"]
            low_mapq_sum = sum(
                count for quality, count in map_qual_dist if quality <= 29
            )
            low_mapq_percentage = (
                (low_mapq_sum / total_reads) * 100 if total_reads > 0 else 0
            )
            if low_mapq_percentage > 20.0:
                output += (
                    f"Warning: map quality <30 ({low_mapq_percentage:.2f}%) exceeds 20%.\n"
                    "Low average MAPQ values suggest poor alignment confidence, which can lead to false variant calls due to reads mapping to multiple locations.\n"
                )
                per_file_counts["bam_cram_warnings"][
                    "% reads map qual <30 >20"
                ] += 1
                warnings_found = True
        else:
            output += (
                "Warning: 'MappingQualityDistribution' not found in QC report.\n"
            )
            warnings_found = True
        # Check 'Duplicates'
        if 'Duplicates' in data_section:
            duplicates_percentage = data_section["Duplicates"][0] * 100
            if duplicates_percentage > 20:
                output += (
                    f"Warning: duplicate reads ({duplicates_percentage:.2f}%) exceed 20%.\n"
                    "High levels of duplicate reads suggest PCR over-amplification, which can skew allele frequency estimates and reduce data complexity.\n"
                )
                per_file_counts["bam_cram_warnings"][
                    "% duplicate >20"
                ] += 1
                warnings_found = True
        else:
            output += "Warning: 'Duplicates' not found in QC report.\n"
            warnings_found = True
    except Exception as e:
        output += f"Error processing BAM/CRAM QC report: {str(e)}\n"
        warnings_found = True
    # Process encrypted stats.txt.openssl.gz to check average GC% content
    if stats_encrypted_path.exists():
        try:
            # Read the encryption key from the 'key' file in the same directory
            key_path = stats_encrypted_path.parent / "key"
            if not key_path.exists():
                output += f"Error: Encryption key file not found at {key_path}\n"
                per_file_counts["bam_cram_warnings"]["key file missing"] += 1
                warnings_found = True
                return output, warnings_found
            with key_path.open('rb') as key_file:
                password = key_file.read().strip()
            # Read the encrypted and compressed stats file into memory
            with stats_encrypted_path.open('rb') as enc_file:
                encrypted_data = enc_file.read()
            # Decrypt and decompress the stats file in memory
            decrypted_stats = decrypt_stats_file(encrypted_data, password)
            # Parse the decrypted stats content
            gcf_gc_content = []
            gcf_counts = []
            gcl_gc_content = []
            gcl_counts = []
            lines = decrypted_stats.splitlines()
            for line in lines:
                line = line.strip()
                if line.startswith('GCF'):
                    parts = line.split('\t')
                    if len(parts) == 2:
                        gc_percent = float(parts[0][3:])
                        count = int(parts[1])
                        gcf_gc_content.append(gc_percent)
                        gcf_counts.append(count)
                elif line.startswith('GCL'):
                    parts = line.split('\t')
                    if len(parts) == 2:
                        gc_percent = float(parts[0][3:])
                        count = int(parts[1])
                        gcl_gc_content.append(gc_percent)
                        gcl_counts.append(count)
            # Compute total counts and GC content
            total_gcf_counts = sum(gcf_counts)
            total_gcl_counts = sum(gcl_counts)
            total_counts = total_gcf_counts + total_gcl_counts
            total_gc_content = sum(gc_percent * count for gc_percent, count in zip(gcf_gc_content, gcf_counts)) + \
                               sum(gc_percent * count for gc_percent, count in zip(gcl_gc_content, gcl_counts))
            if total_counts > 0:
                average_gc_percent = total_gc_content / total_counts
                if not (35 <= average_gc_percent <= 55):
                    output += (
                        f"Warning: average GC content ({average_gc_percent:.2f}%) is outside the acceptable range (35%-55%).\n"
                        "Values outside the expected range of 40-60% in humans may indicate coverage bias or inefficient sequencing in GC-rich or GC-poor regions, or contamination with other species.\n"
                    )
                    per_file_counts["bam_cram_warnings"][
                        "% GC outside 35-55"
                    ] += 1
                    warnings_found = True
            else:
                output += "Warning: total base count is zero, unable to calculate average GC content.\n"
                warnings_found = True
        except Exception as e:
            output += f"Error processing encrypted stats file: {str(e)}\n"
            warnings_found = True
            per_file_counts["bam_cram_warnings"]["stats.txt error"] += 1
    else:
        output += f"Warning: stats.txt.openssl.gz not found at {stats_encrypted_path}\n"
        warnings_found = True
        per_file_counts["bam_cram_warnings"]["stats.txt missing"] += 1
    return output, warnings_found

def analyze_fastq(zip_path: Path, txt_path: Path, warn_missing_file: bool,
                  per_file_counts: Dict) -> Tuple[str, bool]:
    """
    Analyze FASTQ QC reports and check GC content, duplicate reads, and quality scores.
    """
    output = ""
    warnings_found = False
    if not zip_path.exists():
        per_file_counts["fastq_warnings"]["QC report missing"] += 1
        output += f"Error: FASTQ QC report not found at {zip_path}\n"
        warnings_found = True
        return output, warnings_found
    try:
        with zipfile.ZipFile(zip_path, 'r') as zip_ref:
            with zip_ref.open('stdin_fastqc/fastqc_data.txt') as fastqc_file:
                lines = [line.decode('utf-8') for line in fastqc_file.readlines()]
            # Initialize variables
            lines_iter = iter(lines)
            for line in lines_iter:
                line = line.strip()
                if line.startswith(">>Per sequence GC content"):
                    for line in lines_iter:
                        if line.startswith("#GC Content"):
                            break
                    gc_values = []
                    gc_counts = []
                    for line in lines_iter:
                        if line.startswith(">>END_MODULE"):
                            break
                        parts = line.strip().split("\t")
                        if len(parts) == 2:
                            gc_content = float(parts[0])
                            count = float(parts[1])
                            gc_values.append(gc_content)
                            gc_counts.append(count)
                    total_gc_content = sum(gc_counts)
                    if total_gc_content > 0:
                        gc_in_range = sum(
                            count for gc, count in
                            zip(gc_values, gc_counts)
                            if 35 <= gc <= 55
                        )
                        gc_content_percentage = (
                            (gc_in_range / total_gc_content) * 100
                        )
                        if not (35 <= gc_content_percentage <= 55):
                            output += (
                                f"Warning: GC content ({gc_content_percentage:.2f}%) is out of acceptable range (35%-55%).\n"
                                "Values outside the expected range of 40-60% in humans may indicate coverage bias or inefficient sequencing in GC-rich or GC-poor regions, or contamination with other species.\n"
                            )
                            per_file_counts["fastq_warnings"][
                                "% GC outside 35-55"
                            ] += 1
                            warnings_found = True
                    else:
                        output += (
                            "Warning: total GC content is zero, unable to calculate GC content percentage.\n"
                        )
                        warnings_found = True
                elif line.startswith(">>Sequence Duplication Levels"):
                    for line in lines_iter:
                        if line.startswith(
                            "#Total Deduplicated Percentage"
                        ):
                            parts = line.strip().split("\t")
                            if len(parts) == 2:
                                dedup_percentage = float(parts[1])
                                duplicate_reads_percentage = (
                                    100 - dedup_percentage
                                )
                                if duplicate_reads_percentage > 20:
                                    output += (
                                        f"Warning: duplicate reads ({duplicate_reads_percentage:.2f}%) exceed the acceptable threshold (20%).\n"
                                        "High levels of duplicate reads suggest PCR over-amplification, which can skew allele frequency estimates and reduce data complexity.\n"
                                    )
                                    per_file_counts["fastq_warnings"][
                                        "% duplicate >20"
                                    ] += 1
                                    warnings_found = True
                            break
                elif line.startswith(">>Per sequence quality scores"):
                    for line in lines_iter:
                        if line.startswith("#Quality"):
                            break
                    low_quality_count = 0.0
                    total_quality_count = 0.0
                    for line in lines_iter:
                        if line.startswith(">>END_MODULE"):
                            break
                        parts = line.strip().split("\t")
                        if len(parts) == 2:
                            quality = float(parts[0])
                            count = float(parts[1])
                            total_quality_count += count
                            if quality < 20:
                                low_quality_count += count
                    if total_quality_count > 0:
                        low_quality_percentage = (
                            (low_quality_count / total_quality_count)
                            * 100
                        )
                        if low_quality_percentage > 20:
                            output += (
                                f"Warning: reads with quality score <20 ({low_quality_percentage:.2f}%) exceed the acceptable threshold (20%).\n"
                                "Low average base quality suggests sequencing errors or poor-quality reads.\n"
                            )
                            per_file_counts["fastq_warnings"][
                                "% reads qual <20 >20"
                            ] += 1
                            warnings_found = True
                    else:
                        output += (
                            "Warning: total quality score count is zero, unable to calculate percentage.\n"
                        )
                        warnings_found = True
            # Check species information
            species_warning, _, species_issue = check_species(
                txt_path, warn_missing_file
            )
            if species_warning:
                output += species_warning
                warning_key = species_issue
                if warning_key not in per_file_counts["fastq_warnings"]:
                    per_file_counts["fastq_warnings"][warning_key] = 0
                per_file_counts["fastq_warnings"][warning_key] += 1
                warnings_found = True
    except Exception as e:
        output += f"Error processing FASTQ QC report: {str(e)}\n"
        warnings_found = True
    return output, warnings_found

def analyze_vcf(json_path: Path, per_file_counts: Dict) -> Tuple[str, bool]:
    """
    Analyze VCF QC reports and check Ts/Tv ratio and average quality.
    """
    output = ""
    warnings_found = False
    if not json_path.exists():
        per_file_counts["vcf_warnings"]["QC report missing"] += 1
        output += f"Error: VCF QC report not found at {json_path}\n"
        warnings_found = True
        return output, warnings_found
    try:
        with gzip.open(json_path, 'rt') as file:
            data = json.load(file)
        # Check Ts/Tv ratio
        ts_tv_ratio = data.get("TsTvRatio", None)
        if ts_tv_ratio is not None:
            if not (1.9 <= ts_tv_ratio <= 3.3):
                output += (
                    f"Warning: Ts/Tv ratio ({ts_tv_ratio:.2f}) is outside the acceptable range (1.9 - 3.3).\n"
                    "A Ts/Tv ratio far from the expected human range (~2.0-3.3) may indicate poor variant calling or sequencing errors, with an excess of transversions suggesting potential noise.\n"
                )
                per_file_counts["vcf_warnings"]["Ts/Tv ratio outside 1.9-3.3"] += 1
                warnings_found = True
        else:
            output += "Warning: Ts/Tv ratio not found in QC report.\n"
            per_file_counts["vcf_warnings"]["Ts/Tv ratio missing"] += 1
            warnings_found = True
        # Check average variant quality
        avg_quality = data.get("AvgQuality", None)
        if avg_quality is not None:
            if avg_quality < 30:
                output += (
                    f"Warning: average variant quality ({avg_quality:.2f}) is below the acceptable threshold (30).\n"
                    "Low variant QUAL scores indicate low confidence in variant calls, potentially reflecting sequencing errors or insufficient read support for the variant.\n"
                )
                per_file_counts["vcf_warnings"]["Avg QUAL <30"] += 1
                warnings_found = True
        else:
            output += "Warning: average variant quality not found in QC report.\n"
            per_file_counts["vcf_warnings"]["Avg QUAL missing"] += 1
            warnings_found = True
    except Exception as e:
        output += f"Error processing VCF QC report: {str(e)}\n"
        per_file_counts["vcf_warnings"]["QC report error"] += 1
        warnings_found = True
    return output, warnings_found

def process_file(args: Tuple[str, bool]) -> Tuple[str, Dict, str]:
    """
    Process a single file by determining its type and analyzing accordingly.
    """
    file_name, warn_missing_file = args
    output = ""
    file_type = ""
    warnings_found = False
    # Initialize per-file counts
    per_file_counts = {
        "warnings_found": False,
        "files_with_warnings": 0,
        "files_no_qc_report": 0,
        "fastq_warnings": defaultdict(int),
        "bam_cram_warnings": defaultdict(int),
        "vcf_warnings": defaultdict(int),
    }
    try:
        json_path, txt_path, zip_path, stats_encrypted_path = get_file_paths(file_name)
        if json_path.exists():
            try:
                with gzip.open(json_path, 'rt') as file:
                    data = json.load(file)
            except json.JSONDecodeError as e:
                output += (
                    f"Error: failed to parse JSON QC report at {json_path}: "
                    f"{e}\n"
                )
                per_file_counts["files_no_qc_report"] += 1
                warnings_found = True
                data = {}
            if 'VCFVersion' in data:
                file_type = "VCF"
                result, vcf_warnings_found = analyze_vcf(
                    json_path, per_file_counts
                )
                output += result
                warnings_found = warnings_found or vcf_warnings_found
            elif data:
                file_type = "BAM/CRAM"
                result, bam_warnings_found = analyze_bam_cram(
                    json_path, stats_encrypted_path, per_file_counts
                )
                output += result
                warnings_found = warnings_found or bam_warnings_found
            else:
                output += (
                    f"Error: invalid or empty QC report for {file_name}.\n"
                )
                warnings_found = True
        elif zip_path.exists():
            file_type = "FASTQ"
            result, fastq_warnings_found = analyze_fastq(
                zip_path, txt_path, warn_missing_file, per_file_counts
            )
            output += result
            warnings_found = warnings_found or fastq_warnings_found
        else:
            output += "Error: no valid QC report found.\n"
            per_file_counts["files_no_qc_report"] += 1
            warnings_found = True
        if output.strip():
            if file_type:
                output = f">> {file_name} | {file_type}\n" + output
            else:
                output = f">> {file_name} | unknown\n" + output
        per_file_counts["warnings_found"] = warnings_found
        if warnings_found:
            per_file_counts["files_with_warnings"] += 1
        return output, per_file_counts, file_type
    except Exception as e:
        logging.error(f"Exception in processing file {file_name}: {e}")
        output += f"Exception in processing file {file_name}: {e}\n"
        warnings_found = True
        per_file_counts["warnings_found"] = warnings_found
        per_file_counts["files_with_warnings"] += 1
        return output, per_file_counts, file_type

def read_egaf_from_file(file_path: Path) -> List[str]:
    """
    Read EGAF IDs from a file.
    """
    with file_path.open('r') as f:
        return [line.strip() for line in f if line.strip()]

def chunks(iterable: List, size: int) -> Iterator[List]:
    """
    Yield successive chunks from the iterable.
    """
    for i in range(0, len(iterable), size):
        yield iterable[i:i + size]

def main():
    """
    Main function to parse arguments and initiate processing of files.
    """
    parser = argparse.ArgumentParser(
        description=(
            "Analyze BAM/CRAM, FASTQ, or VCF files for quality control issues."
        ),
        epilog=(
            "Example usage:\n  python script.py --file egaf_list.txt "
            "--output results.txt"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument('--egaf', type=str, help="Single EGAF ID to analyze")
    parser.add_argument(
        '--file',
        type=str,
        help="File (txt/csv/tsv) containing multiple EGAF IDs",
    )
    parser.add_argument(
        '--output',
        type=str,
        required=True,
        help="Output file to append results",
    )
    parser.add_argument(
        '--no-species-warning',
        action='store_true',
        help="Do not warn if the species file is missing",
    )
    parser.add_argument(
        '--threads',
        type=int,
        default=4,
        help="Number of threads for concurrent processing",
    )
    args = parser.parse_args()
    warn_missing_file = not args.no_species_warning
    # Configure logging
    logging.basicConfig(level=logging.INFO, format='%(message)s')
    start_time = time.time()
    output_file = Path(args.output)
    if args.egaf:
        egaf_list = [args.egaf]
    elif args.file:
        egaf_list = read_egaf_from_file(Path(args.file))
    else:
        parser.error("Provide either --egaf or --file argument.")
    # Ensure output file exists and is empty
    output_file.touch()
    output_file.write_text('')
    # Initialize counters
    total_files_checked = 0
    files_with_warnings = 0
    files_no_qc_report = 0
    files_with_qc_report = 0
    fastq_warnings = defaultdict(int)
    bam_cram_warnings = defaultdict(int)
    vcf_warnings = defaultdict(int)
    total_fastq_files = 0
    total_bam_cram_files = 0
    total_vcf_files = 0
    files_with_warnings_and_qc = 0
    batch_size = 10000  # Adjust based on system resources
    args_list = [(egaf, warn_missing_file) for egaf in egaf_list]
    with concurrent.futures.ProcessPoolExecutor(
        max_workers=args.threads
    ) as executor, output_file.open('a') as f:
        # Process files in batches
        for batch in chunks(args_list, batch_size):
            futures = [executor.submit(process_file, arg) for arg in batch]
            for future in concurrent.futures.as_completed(futures):
                try:
                    output, per_file_counts, file_type = future.result()
                    total_files_checked += 1
                    if per_file_counts["warnings_found"]:
                        files_with_warnings += 1
                    files_no_qc_report += per_file_counts["files_no_qc_report"]
                    if per_file_counts["files_no_qc_report"] == 0:
                        files_with_qc_report += 1
                        if per_file_counts["warnings_found"]:
                            files_with_warnings_and_qc += 1
                    # Count file types
                    if file_type == "FASTQ":
                        total_fastq_files += 1
                    elif file_type == "BAM/CRAM":
                        total_bam_cram_files += 1
                    elif file_type == "VCF":
                        total_vcf_files += 1
                    # Accumulate warnings
                    for key in per_file_counts["fastq_warnings"]:
                        fastq_warnings[key] += per_file_counts[
                            "fastq_warnings"
                        ][key]
                    for key in per_file_counts["bam_cram_warnings"]:
                        bam_cram_warnings[key] += per_file_counts[
                            "bam_cram_warnings"
                        ][key]
                    for key in per_file_counts["vcf_warnings"]:
                        vcf_warnings[key] += per_file_counts[
                            "vcf_warnings"
                        ][key]
                    # Write output to the file
                    f.write(output)
                except Exception as e:
                    logging.error(f"Exception in processing: {e}")
    elapsed_time = time.time() - start_time
    # Print final summary
    logging.info(f"\nFinished checking {total_files_checked} files.")
    if total_files_checked > 0:
        missing_qc_percentage = (
            (files_no_qc_report / total_files_checked) * 100
        )
        logging.info(
            f"{files_no_qc_report} files "
            f"({missing_qc_percentage:.1f}%) have a missing QC report."
        )
        if files_with_qc_report > 0:
            warnings_percentage = (
                (files_with_warnings_and_qc / files_with_qc_report) * 100
            )
            logging.info(
                f"{files_with_qc_report} files with report identified, "
                f"{files_with_warnings_and_qc} have shown some warning "
                f"({warnings_percentage:.1f}%)."
            )
        else:
            logging.info("No files with QC report identified.")
    else:
        logging.info("No files were checked.")
    logging.info("\nSummary of warnings:")
    if total_fastq_files > 0:
        logging.info(f"- FASTQ ({total_fastq_files} total files):")
        for warning_key, count in fastq_warnings.items():
            if count > 0:
                logging.info(f" {count} files with warning: {warning_key}")
    if total_bam_cram_files > 0:
        logging.info(f"- BAM/CRAM ({total_bam_cram_files} total files):")
        for warning_key, count in bam_cram_warnings.items():
            if count > 0:
                logging.info(f" {count} files with warning: {warning_key}")
    if total_vcf_files > 0:
        logging.info(f"- VCF ({total_vcf_files} total files):")
        for warning_key, count in vcf_warnings.items():
            if count > 0:
                logging.info(f" {count} files with warning: {warning_key}")
    logging.info(f"\nTime elapsed: {elapsed_time:.2f} seconds.")

if __name__ == "__main__":
    main()
