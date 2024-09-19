import json
import os
import zipfile
import re

# get paths for BAM/CRAM or FASTQ files
def get_file_paths(file_name):
#   base_path = "/slgpfs/projects/slc00/slc00474/execution-qc/vault/archive/vault/archive"
    base_path = "/Users/raul/Library/CloudStorage/OneDrive-CRG-CentredeRegulacioGenomica/ega.nosync/bioteam/slc00474/execution-qc/vault/archive"
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
        print(f"⚠️ Warning: Species file input_screen.txt not found at {txt_path}")
        return

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
                        print(f"✅ Most reads mapped to Human ({percent_mapped:.2f}%).")
                    else:
                        print(f"🚨 Caution: Most reads mapped to {first_cell} ({percent_mapped:.2f}%).")

    except Exception as e:
        print(f"Error processing input_screen.txt species file: {e}")


# analyze BAM/CRAM file (JSON)
def analyze_bam_cram(json_path, txt_path):
    if not os.path.exists(json_path):
        print(f"⚠️ BAM/CRAM QC report not found at {json_path}")
        return

    print("✅ BAM/CRAM QC report loaded.")

    try:
        with open(json_path, 'r') as file:
            data = json.load(file)

        total_reads = data["TotalReads"]
        mapped_reads_ratio = data["Data"]["MappedReads"][0]
        reads_unaligned = 100 - (mapped_reads_ratio * 100)

        # check aligned reads percentage
        if reads_unaligned > 40.0:
            print(f"⚠️ Reads unaligned ({reads_unaligned:.2f}%) exceeds 40%.")
        else:
            print(f"✅ Reads unaligned ({reads_unaligned:.2f}%) is acceptable (<40%).")

        # sum map quality values from 0 to 29
        map_quality_distribution = data["Data"]["MappingQualityDistribution"]
        low_map_quality_sum = sum(count for quality, count in map_quality_distribution if quality <= 29)

        low_map_quality_percentage = (low_map_quality_sum / total_reads) * 100

        # check map quality percentage
        if low_map_quality_percentage > 5.0:
            print(f"⚠️ Map quality <30 ({low_map_quality_percentage:.2f}%) exceeds 5%.")
        else:
            print(f"✅ Map quality >30 ({low_map_quality_percentage:.2f}%) is acceptable (≤5%).")

        # check species function
        check_species(txt_path)

    except (json.JSONDecodeError, KeyError) as e:
        print(f"Error processing BAM/CRAM QC report: {e}")


# analyze FASTQ file (ZIP and TXT extracting)
def analyze_fastq(zip_path, txt_path):
    if not os.path.exists(zip_path):
        print(f"⚠️ FASTQ QC report not found at {zip_path}")
        return

    try:
        with zipfile.ZipFile(zip_path, 'r') as zip_ref:
            zip_ref.extractall("/tmp")  # extract to tmp directory
            fastqc_data_path = "/tmp/stdin_fastqc/fastqc_data.txt"

        if not os.path.exists(fastqc_data_path):
            print(f"⚠️ fastqc_data.txt not found in {zip_path}")
            return

        print("✅ FASTQ QC report loaded.")

        with open(fastqc_data_path, 'r') as file:
            lines = file.readlines()

        # parse duplicate percentage
        try:
            dup_module_start = next(i for i, line in enumerate(lines) if ">>Sequence Duplication Levels" in line)
            dup_data = lines[dup_module_start:dup_module_start + 10]  # Adjust as necessary
            dup_percentage_line = next(line for line in dup_data if "Total Deduplicated Percentage" in line)
            dup_match = re.search(r'([\d.]+)', dup_percentage_line)

            if dup_match:
                dup_percentage = 100 - float(dup_match.group(1))  # Inverse because it's deduplicated
                if dup_percentage > 20.0:
                    print(f"⚠️ Duplicate reads ({dup_percentage:.2f}%) exceed 20%.")
                else:
                    print(f"✅ Duplicate reads ({dup_percentage:.2f}%) are acceptable (≤20%).")
            else:
                print(f"⚠️ Unable to find duplicate percentage.")
        except StopIteration:
            print(f"⚠️ Sequence Duplication Levels section not found in fastqc_data.txt.")

        # parse GC content
        try:
            gc_module_start = next(i for i, line in enumerate(lines) if ">>Per sequence GC content" in line)
            gc_data = lines[gc_module_start:gc_module_start + 100]
            gc_content_values = [float(line.split()[0]) for line in gc_data if line.strip() and line[0].isdigit()]
            avg_gc_content = sum(gc_content_values) / len(gc_content_values) if gc_content_values else None

            if avg_gc_content:
                if avg_gc_content < 35.0 or avg_gc_content > 55.0:
                    print(f"⚠️ GC content ({avg_gc_content:.2f}%) is outside the 35-55% range.")
                else:
                    print(f"✅ GC content ({avg_gc_content:.2f}%) is within the acceptable range (35-55%).")
            else:
                print(f"⚠️ Unable to find GC content percentage.")
        except StopIteration:
            print(f"⚠️ Per sequence GC content section not found in fastqc_data.txt.")

        # check species function
        check_species(txt_path)

    except zipfile.BadZipFile as e:
        print(f"Error processing FASTQ QC report: {e}")

# main function to infer if EGAF is BAM/CRAM or FASTQ
def analyze_file(file_name):
    # print EGAF identifier
    print(f">> {file_name}")

    json_path, txt_path, zip_path = get_file_paths(file_name)

    # check if it's BAM/CRAM (JSON) or FASTQ (ZIP)
    if os.path.exists(json_path):
        analyze_bam_cram(json_path, txt_path)
    elif os.path.exists(zip_path):
        analyze_fastq(zip_path, txt_path)
    else:
        print("⚠️ No valid QC report found.")


# ID usage. # example usage, modify the file to check the flags for that file, this will change to an option parser
#file_name = "EGAF00005432457" #bam
file_name = "EGAF00002237576" #fastq
analyze_file(file_name)