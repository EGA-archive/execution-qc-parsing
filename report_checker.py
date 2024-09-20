import json
import os
import zipfile
import argparse

# get paths for BAM/CRAM or FASTQ files
def get_file_paths(file_name):
    base_path = "/slgpfs/projects/slc00/slc00474/execution-qc/vault/archive/vault/archive"
    prefix = file_name[:9]  # EGAF00005
    middle = file_name[9:12]  # 432
    suffix = file_name[12:15]  # 457
    json_path = os.path.join(base_path, prefix, middle, suffix, "execution", f"{file_name}_report.json")
    txt_path = os.path.join(base_path, prefix, middle, suffix, "execution", "input_screen.txt")
    zip_path = os.path.join(base_path, prefix, middle, suffix, "execution", "stdin_fastqc.zip")
    return json_path, txt_path, zip_path


# check species information (common for both BAM/CRAM and FASTQ)
def check_species(txt_path):
    if not os.path.exists(txt_path):
        return f" Error: Species file input_screen.txt not found at {txt_path}\n"

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
                        return f" Warning: Most reads mapped to {first_cell} ({percent_mapped:.2f}%).\n"
    except Exception as e:
        return f"Error processing input_screen.txt species file: {e}\n"

    return ""


# analyze BAM/CRAM file (JSON)
def analyze_bam_cram(json_path, txt_path):
    output = ""
    if not os.path.exists(json_path):
        return f"Error: BAM/CRAM QC report not found at {json_path}\n"

    try:
        with open(json_path, 'r') as file:
            data = json.load(file)

        total_reads = data["TotalReads"]
        mapped_reads_ratio = data["Data"]["MappedReads"][0]
        reads_unaligned = 100 - (mapped_reads_ratio * 100)

        # check aligned reads percentage
        if reads_unaligned > 40.0:
            output += f"Warning: Reads unaligned ({reads_unaligned:.2f}%) exceeds 40%.\n"

        # sum map quality values from 0 to 29
        map_quality_distribution = data["Data"]["MappingQualityDistribution"]
        low_map_quality_sum = sum(count for quality, count in map_quality_distribution if quality <= 29)

        low_map_quality_percentage = (low_map_quality_sum / total_reads) * 100

        # check map quality percentage
        if low_map_quality_percentage > 5.0:
            output += f"Warning: Map quality <30 ({low_map_quality_percentage:.2f}%) exceeds 5%.\n"

        # run species function
        species_warning = check_species(txt_path)
        if species_warning:
            output += species_warning

    except (json.JSONDecodeError, KeyError) as e:
        output += f"Error processing BAM/CRAM QC report: {e}\n"

    return output


# analyze FASTQ file (ZIP and TXT extracting)
def analyze_fastq(zip_path, txt_path):
    output = ""

    # extract qc data from zip
    if not os.path.exists(zip_path):
        output += f"Error: FASTQ QC report not found at {zip_path}\n"
        return output

    try:
        with zipfile.ZipFile(zip_path, 'r') as zip_ref:
            zip_ref.extractall("/tmp")  # Extract to a temporary directory
            fastqc_data_path = "/tmp/stdin_fastqc/fastqc_data.txt"

        if not os.path.exists(fastqc_data_path):
            output += f"Error: fastqc_data.txt not found in {zip_path}\n"
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
                    if gc_values:
                        total_gc_content = sum(gc_values)
                        gc_content_percentage = (sum([gc_values[i] for i in range(35, 45)]) / total_gc_content) * 100
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
                output += f"Warning: GC content ({gc_content_percentage:.2f}%) is out of acceptable range (35%-55%).\n"
        else:
            output += "Error: Unable to find GC content percentage.\n"

        # output duplicate reads result
        if duplicate_reads_percentage is not None:
            if duplicate_reads_percentage > 20:
                output += f"Warning: Duplicate reads ({duplicate_reads_percentage:.2f}%) exceed the acceptable threshold (20%).\n"
        else:
            output += "Error: Unable to find duplicate percentage.\n"

        # check species function and append
        species_warning = check_species(txt_path)
        if species_warning:
            output += species_warning

        return output

    except zipfile.BadZipFile as e:
        output += f"Error processing FASTQ QC report: {e}\n"
        return output


# main function to infer if EGAF is BAM/CRAM or FASTQ
def analyze_file(file_name, output_file):
    output = ""

    json_path, txt_path, zip_path = get_file_paths(file_name)

    # check if is BAM/CRAM (JSON) or FASTQ (ZIP)
    if os.path.exists(json_path):
        output += analyze_bam_cram(json_path, txt_path)
    elif os.path.exists(zip_path):
        output += analyze_fastq(zip_path, txt_path)
    else:
        output += "Error: No valid QC report found.\n"

    # write output only if there are warnings
    if output:
        output = f">> {file_name}\n" + output
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

    args = parser.parse_args()

    if args.egaf:
        # Analyze a single EGAF
        analyze_file(args.egaf, args.output)
    elif args.file:
        # Read EGAFs from the provided file and analyze each
        egaf_list = read_egaf_from_file(args.file)
        for egaf in egaf_list:
            analyze_file(egaf, args.output)
    else:
        print("Provide either --egaf or --file argument.")