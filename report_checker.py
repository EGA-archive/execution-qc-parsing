import argparse
import concurrent.futures
import gzip
import json
import logging
import tempfile
import time
import zipfile
from pathlib import Path
from typing import List, Tuple, Dict

# base path
BASE_PATH = Path("/slgpfs/projects/slc00/slc00474/execution-qc/vault/archive")

def get_file_paths(file_name: str) -> Tuple[Path, Path, Path]:
    prefix = file_name[:9]  # egaf00005
    middle = file_name[9:12]  # 432
    suffix = file_name[12:15]  # 457
    base_dir = BASE_PATH / prefix / middle / suffix / "execution"
    json_path = base_dir / f"{file_name}_report.json.gz"
    txt_path = base_dir / "input_screen.txt"
    zip_path = base_dir / "stdin_fastqc.zip"
    return json_path, txt_path, zip_path

def check_species(
    txt_path: Path, warn_missing_file: bool
) -> Tuple[str, bool]:
    if not txt_path.exists():
        if warn_missing_file:
            return (
                f"Warning: Species file input_screen.txt not found at "
                f"{txt_path}\n",
                True,
            )
        return "", False  # do not warn if file is missing and flag is set

    try:
        with txt_path.open('r') as file:
            lines = file.readlines()

        # process third line (first matching species)
        if len(lines) >= 3:
            third_line = lines[2].strip()
            cells = third_line.split("\t")

            # extract species and percentage mapped
            if len(cells) >= 4:
                percent_unmapped = float(cells[3])
                percent_mapped = 100 - percent_unmapped

                # check if the first cell contains "Human"
                first_cell = cells[0]
                if first_cell == "Human":
                    return "", False  # no warning
                else:
                    return (
                        f"Warning: Most reads mapped to {first_cell} "
                        f"({percent_mapped:.2f}%).\n",
                        True,
                    )
    except Exception as e:
        return (
            f"Error processing input_screen.txt species file: {str(e)}\n",
            True,
        )

    return "", False

def analyze_bam_cram(
    json_path: Path, per_file_counts: Dict
) -> Tuple[str, bool]:
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
        mapped_ratio = data["Data"]["MappedReads"][0]
        reads_unaligned = 100 - (mapped_ratio * 100)

        # check aligned reads percentage
        if reads_unaligned > 40.0:
            output += (
                f"Warning: Reads unaligned ({reads_unaligned:.2f}%) "
                f"exceeds 40%.\n"
            )
            per_file_counts["bam_cram_warnings"]["% reads unaligned >40"] += 1
            warnings_found = True

        # sum map quality values from 0 to 29
        map_qual_dist = data["Data"]["MappingQualityDistribution"]
        low_mapq_sum = sum(
            count for quality, count in map_qual_dist if quality <= 29
        )

        low_mapq_percentage = (
            (low_mapq_sum / total_reads) * 100 if total_reads > 0 else 0
        )

        # check map quality percentage
        if low_mapq_percentage > 20.0:
            output += (
                f"Warning: Map quality <30 ({low_mapq_percentage:.2f}%) "
                f"exceeds 20%.\n"
            )
            per_file_counts["bam_cram_warnings"][
                "% reads map qual <30 >20"
            ] += 1
            warnings_found = True

        # check duplicate reads percentage
        duplicates_percentage = data["Data"]["Duplicates"][0] * 100
        if duplicates_percentage > 20:
            output += (
                f"Warning: Duplicate reads ({duplicates_percentage:.2f}%) "
                f"exceed 20%.\n"
            )
            per_file_counts["bam_cram_warnings"]["% duplicate >20"] += 1
            warnings_found = True

    except (json.JSONDecodeError, KeyError, OSError) as e:
        output += f"Error processing BAM/CRAM QC report: {str(e)}\n"
        warnings_found = True

    return output, warnings_found

def analyze_fastq(
    zip_path: Path,
    txt_path: Path,
    warn_missing_file: bool,
    per_file_counts: Dict,
) -> Tuple[str, bool]:
    output = ""
    warnings_found = False

    if not zip_path.exists():
        per_file_counts["fastq_warnings"]["QC report missing"] += 1
        output += f"Error: FASTQ QC report not found at {zip_path}\n"
        warnings_found = True
        return output, warnings_found

    try:
        with zipfile.ZipFile(zip_path, 'r') as zip_ref:
            with tempfile.TemporaryDirectory() as tmpdirname:
                zip_ref.extractall(tmpdirname)
                fastqc_data_path = (
                    Path(tmpdirname) / 'stdin_fastqc' / 'fastqc_data.txt'
                )

                if not fastqc_data_path.exists():
                    output += (
                        f"Error: fastqc_data.txt not found in {zip_path}\n"
                    )
                    warnings_found = True
                    return output, warnings_found

                with fastqc_data_path.open('r') as file:
                    lines = file.readlines()

                # initialize variables
                gc_content_percentage = None
                duplicate_reads_percentage = None

                # parse the fastqc_data.txt file
                lines_iter = iter(lines)
                for line in lines_iter:
                    line = line.strip()
                    if line.startswith(">>Per sequence GC content"):
                        # skip until "#GC Content" line
                        for line in lines_iter:
                            if line.startswith("#GC Content"):
                                break
                        gc_values = []
                        gc_counts = []
                        # read GC content data
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
                                count
                                for gc, count in zip(gc_values, gc_counts)
                                if 35 <= gc <= 55
                            )
                            gc_content_percentage = (
                                (gc_in_range / total_gc_content) * 100
                            )
                            if not (35 <= gc_content_percentage <= 55):
                                output += (
                                    f"Warning: GC content "
                                    f"({gc_content_percentage:.2f}%) is out "
                                    f"of acceptable range (35%-55%).\n"
                                )
                                per_file_counts["fastq_warnings"][
                                    "% GC outside 35-55"
                                ] += 1
                                warnings_found = True
                        else:
                            output += (
                                "Warning: Total GC content is zero, unable "
                                "to calculate GC content percentage.\n"
                            )
                            warnings_found = True

                    elif line.startswith(">>Sequence Duplication Levels"):
                        # skip until "#Total Deduplicated Percentage" line
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
                                            f"Warning: Duplicate reads "
                                            f"({duplicate_reads_percentage:.2f}%) "
                                            "exceed the acceptable "
                                            "threshold (20%).\n"
                                        )
                                        per_file_counts["fastq_warnings"][
                                            "% duplicate >20"
                                        ] += 1
                                        warnings_found = True
                                break

                # check species information
                species_warning, species_flag = check_species(
                    txt_path, warn_missing_file
                )
                if species_flag:
                    output += species_warning
                    per_file_counts["fastq_warnings"]["sp not-human"] += 1
                    warnings_found = True

    except zipfile.BadZipFile as e:
        output += f"Error processing FASTQ QC report: {str(e)}\n"
        warnings_found = True

    return output, warnings_found

def process_file(
    file_name: str, warn_missing_file: bool
) -> Tuple[str, Dict]:
    output = ""
    file_type = ""
    warnings_found = False

    per_file_counts = {
        "warnings_found": False,
        "files_with_warnings": 0,
        "files_no_qc_report": 0,
        "fastq_warnings": {
            "QC report missing": 0,
            "sp not-human": 0,
            "% duplicate >20": 0,
            "% GC outside 35-55": 0,
        },
        "bam_cram_warnings": {
            "QC report missing": 0,
            "% reads unaligned >40": 0,
            "% reads map qual <30 >20": 0,
            "% duplicate >20": 0,
        },
    }

    json_path, txt_path, zip_path = get_file_paths(file_name)

    if json_path.exists():
        file_type = "BAM/CRAM"
        result, warnings_found = analyze_bam_cram(
            json_path, per_file_counts
        )
        output += result
    elif zip_path.exists():
        file_type = "FASTQ"
        result, warnings_found = analyze_fastq(
            zip_path, txt_path, warn_missing_file, per_file_counts
        )
        output += result
    else:
        output += "Error: No valid QC report found.\n"
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

    return output, per_file_counts

def read_egaf_from_file(file_path: Path) -> List[str]:
    with file_path.open('r') as f:
        return [line.strip() for line in f if line.strip()]

def main():
    parser = argparse.ArgumentParser(
        description=(
            "analyze BAM/CRAM or FASTQ files for quality control issues."
        ),
        epilog=(
            "example usage:\n  python script.py --file egaf_list.txt "
            "--output results.txt"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument('--egaf', type=str, help="single EGAF ID to analyze")
    parser.add_argument(
        '--file',
        type=str,
        help="file (txt/csv/tsv) containing multiple EGAF IDs",
    )
    parser.add_argument(
        '--output',
        type=str,
        required=True,
        help="output file to append results",
    )
    parser.add_argument(
        '--no-species-warning',
        action='store_true',
        help="do not warn if the species file is missing",
    )
    parser.add_argument(
        '--threads',
        type=int,
        default=4,
        help="number of threads for concurrent processing",
    )
    args = parser.parse_args()

    warn_missing_file = not args.no_species_warning

    # configure logging
    logging.basicConfig(level=logging.INFO, format='%(message)s')

    start_time = time.time()

    output_file = Path(args.output)

    if args.egaf:
        egaf_list = [args.egaf]
    elif args.file:
        egaf_list = read_egaf_from_file(Path(args.file))
    else:
        parser.error("Provide either --egaf or --file argument.")

    # ensure output file exists and is empty
    output_file.touch()
    output_file.write_text('')

    total_files_checked = 0
    files_with_warnings = 0
    files_no_qc_report = 0

    fastq_warnings = {
        "QC report missing": 0,
        "sp not-human": 0,
        "% duplicate >20": 0,
        "% GC outside 35-55": 0,
    }
    bam_cram_warnings = {
        "QC report missing": 0,
        "% reads unaligned >40": 0,
        "% reads map qual <30 >20": 0,
        "% duplicate >20": 0,
    }

    # use ProcessPoolExecutor for parallel processing
    with concurrent.futures.ProcessPoolExecutor(
        max_workers=args.threads
    ) as executor:
        futures = {
            executor.submit(process_file, egaf, warn_missing_file): egaf
            for egaf in egaf_list
        }

        with output_file.open('a') as f:
            for future in concurrent.futures.as_completed(futures):
                egaf = futures[future]
                try:
                    output, per_file_counts = future.result()
                    total_files_checked += 1
                    if per_file_counts["warnings_found"]:
                        files_with_warnings += 1
                    files_no_qc_report += per_file_counts[
                        "files_no_qc_report"
                    ]

                    for key in fastq_warnings:
                        fastq_warnings[key] += per_file_counts[
                            "fastq_warnings"
                        ][key]

                    for key in bam_cram_warnings:
                        bam_cram_warnings[key] += per_file_counts[
                            "bam_cram_warnings"
                        ][key]

                    # write output to the file
                    f.write(output)
                except Exception as e:
                    logging.error(f"Error processing file {egaf}: {e}")

    elapsed_time = time.time() - start_time

    # print final summary
    logging.info(f"\nFinished checking for {total_files_checked} files.")
    if total_files_checked > 0:
        warning_percentage = (
            (files_with_warnings / total_files_checked) * 100
        )
        logging.info(
            f"{files_with_warnings} files "
            f"({warning_percentage:.1f}%) have shown some warning."
        )
        if files_no_qc_report > 0:
            missing_qc_percentage = (
                (files_no_qc_report / total_files_checked) * 100
            )
            logging.info(
                f"{files_no_qc_report} files "
                f"({missing_qc_percentage:.1f}%) have a missing QC report."
            )
    else:
        logging.info("No files were checked.")

    logging.info("\nSummary of warnings:")

    if any(count > 0 for count in fastq_warnings.values()):
        logging.info("- FASTQ:")
        for warning, count in fastq_warnings.items():
            if count > 0:
                logging.info(f"  {count} {warning}")

    if any(count > 0 for count in bam_cram_warnings.values()):
        logging.info("- BAM/CRAM:")
        for warning, count in bam_cram_warnings.items():
            if count > 0:
                logging.info(f"  {count} {warning}")

    logging.info(f"\nTime elapsed: {elapsed_time:.2f} seconds.")

if __name__ == "__main__":
    main()