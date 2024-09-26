import json
import os
import zipfile
import argparse
import gzip

# get paths for BAM/CRAM or FASTQ files
def get_file_paths(file_name):
    base_path = "/slgpfs/projects/slc00/slc00474/execution-qc/vault/archive"
    prefix = file_name[:9]  # EGAF00005
    middle = file_name[9:12]  # 432
    suffix = file_name[12:15]  # 457
    json_path = os.path.join(base_path, prefix, middle, suffix, "execution", "{}_report.json.gz".format(file_name))
    txt_path = os.path.join(base_path, prefix, middle, suffix, "execution", "input_screen.txt")
    zip_path = os.path.join(base_path, prefix, middle, suffix, "execution", "stdin_fastqc.zip")
    return json_path, txt_path, zip_path


# check species information (common for both BAM/CRAM and FASTQ)
def check_species(txt_path, warn_missing_file=True):
    if not os.path.exists(txt_path):
        if warn_missing_file:
            return "Warning: Species file input_screen.txt not found at {}\n".format(txt_path)
        return ""  # Do not warn if file is missing and flag is set

    try:
        with open(txt_path, 'r') as file:
            lines = file.readlines()

            # process third line (first matching sp)
            if len(lines) >= 3:
                third_line = lines[2].strip()
                cells = third_line.split("\t")

                # extract species and percentage mapped
                if len(cells) >= 4:
                    percent_unmapped = float(cells[3])
                    percent_mapped = 100 - percent_unmapped

                    # check if the first cell in the third line contains "Human"
                    first_cell = cells[0]
                    if first_cell == "Human":
                        return ""  # No warning, pass the OK message
                    else:
                        return "Warning: Most reads mapped to {} ({:.2f}%).\n".format(first_cell, percent_mapped)
    except Exception as e:
        return "Error processing input_screen.txt species file: {}\n".format(str(e))

    return ""


# analyze BAM/CRAM file (JSON.GZ)
def analyze_bam_cram(json_path, txt_path, warn_missing_file):
    output = ""
    if not os.path.exists(json_path):
        return "Error: BAM/CRAM QC report not found at {}\n".format(json_path)

    try:
        # Open and read the .json.gz file
        with gzip.open(json_path, 'rt') as file:
            data = json.load(file)

        total_reads = data["TotalReads"]
        mapped_reads_ratio = data["Data"]["MappedReads"][0]
        reads_unaligned = 100 - (mapped_reads_ratio * 100)

        # check aligned reads percentage
        if reads_unaligned > 40.0:
            output += "Warning: Reads unaligned ({:.2f}%) exceeds 40%.\n".format(reads_unaligned)

        # sum map quality values from 0 to 29
        map_quality_distribution = data["Data"]["MappingQualityDistribution"]
        low_map_quality_sum = sum(count for quality, count in map_quality_distribution if quality <= 29)

        low_map_quality_percentage = (low_map_quality_sum / total_reads) * 100

        # check map quality percentage
        if low_map_quality_percentage > 5.0:
            output += "Warning: Map quality <30 ({:.2f}%) exceeds 5%.\n".format(low_map_quality_percentage)

        # check duplicate reads percentage
        duplicates_percentage = data["Data"]["Duplicates"][0] * 100
        if duplicates_percentage > 20:
            output += "Warning: Duplicate reads ({:.2f}%) exceed 20%.\n".format(duplicates_percentage)

        # run species function
        species_warning = check_species(txt_path, warn_missing_file)
        if species_warning:
            output += species_warning

    except (json.JSONDecodeError, KeyError, OSError) as e:
        output += "Error processing BAM/CRAM QC report: {}\n".format(str(e))

    return output


# analyze FASTQ file (ZIP and TXT extracting)
def analyze_fastq(zip_path, txt_path, warn_missing_file):
    output = ""

    # extract qc data from zip
    if not os.path.exists(zip_path):
        output += "Error: FASTQ QC report not found at {}\n".format(zip_path)
        return output

    try:
        with zipfile.ZipFile(zip_path, 'r') as zip_ref:
            zip_ref.extractall("/tmp/raul")  # Extract to a temporary directory
            fastqc_data_path = "/tmp/raul/stdin_fastqc/fastqc_data.txt"

        if not os.path.exists(fastqc_data_path):
            output += "Error: fastqc_data.txt not found in {}\n".format(zip_path)
            return output

        # analyze fastqc report
        with open(fastqc_data_path, 'r') as file:
            lines = file.readlines()

            # initialize variables
            gc_content_found = False
            duplicate_found = False
            gc_content_percentage = None
            duplicate_reads_percentage = None

            # looping line by line
            for i, line in enumerate(lines):
                # look for GC content
                if ">>Per sequence GC content" in line:
                    gc_content_found = True
                if gc_content_found and line.startswith("#GC Content"):
                    gc_values = []
                    for gc_line in lines[i + 1:]:
                        if ">>END_MODULE" in gc_line:
                            break
                        _, count = gc_line.split("\t")
                        gc_values.append(float(count))

                    # calculate GC content percentage
                    total_gc_content = sum(gc_values)
                    if total_gc_content > 0:
                        gc_content_percentage = (sum([gc_values[i] for i in range(35, 45)]) / total_gc_content) * 100
                    else:
                        gc_content_percentage = None  # Skip GC content percentage calculation if it's zero
                        output += "Warning: Total GC content is zero, unable to calculate GC content percentage.\n"
                    gc_content_found = False

                # look for duplicate reads
                if ">>Sequence Duplication Levels" in line:
                    duplicate_found = True
                if duplicate_found and line.startswith("#Total Deduplicated Percentage"):
                    duplicate_reads_percentage = 100 - float(line.split("\t")[1])
                    duplicate_found = False

        # output GC content result
        if gc_content_percentage is not None:
            if not (35 <= gc_content_percentage <= 55):
                output += "Warning: GC content ({:.2f}%) is out of acceptable range (35%-55%).\n".format(gc_content_percentage)
        else:
            output += "Error: Unable to find GC content percentage.\n"

        # output duplicate reads result
        if duplicate_reads_percentage is not None:
            if duplicate_reads_percentage > 20:
                output += "Warning: Duplicate reads ({:.2f}%) exceed the acceptable threshold (20%).\n".format(duplicate_reads_percentage)
        else:
            output += "Error: Unable to find duplicate percentage.\n"

        # check species function and append
        species_warning = check_species(txt_path, warn_missing_file)
        if species_warning:
            output += species_warning

        return output

    except zipfile.BadZipFile as e:
        output += "Error processing FASTQ QC report: {}\n".format(str(e))
        return output

# main function to infer if EGAF is BAM/CRAM or FASTQ
def analyze_file(file_name, output_file, warn_missing_file):
    output = ""

    json_path, txt_path, zip_path = get_file_paths(file_name)

    # check if is BAM/CRAM (JSON) or FASTQ (ZIP)
    if os.path.exists(json_path):
        output += analyze_bam_cram(json_path, txt_path, warn_missing_file)
    elif os.path.exists(zip_path):
        output += analyze_fastq(zip_path, txt_path, warn_missing_file)
    else:
        output += "Error: No valid QC report found.\n"

    # write output only if there are warnings
    if output:
        output = ">> {}\n".format(file_name) + output
        with open(output_file, 'a') as f:
            f.write(output)


# parse EGAFs to look from file (txt/csv/tsv)
def read_egaf_from_file(file_path):
    with open(file_path, 'r') as f:
        return [line.strip() for line in f.readlines() if line.strip()]


# argument parser
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Analyze BAM/CRAM or FASTQ files.")
    parser.add_argument('--egaf', type=str, help="Single EGAF ID to analyze")
    parser.add_argument('--file', type=str, help="File (txt/csv/tsv) containing multiple EGAF IDs")
    parser.add_argument('--output', type=str, required=True, help="Output file to append results")
    parser.add_argument('--no-species-warning', action='store_true', help="Do not warn if the species file is missing")

    args = parser.parse_args()

    warn_missing_file = not args.no_species_warning

    if args.egaf:
        # Analyze a single EGAF
        analyze_file(args.egaf, args.output, warn_missing_file)
    elif args.file:
        # Read EGAFs from the provided file and analyze each
        egaf_list = read_egaf_from_file(args.file)
        for egaf in egaf_list:
            analyze_file(egaf, args.output, warn_missing_file)
    else:
        print("Provide either --egaf or --file argument.")