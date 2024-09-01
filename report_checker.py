import json
import os

def get_file_paths(file_name):
    base_path = "/slgpfs/projects/slc00/slc00474/execution-qc/vault/archive/vault/archive"
    prefix = file_name[:9]  # EGAF00005
    middle = file_name[9:12]  # 432
    suffix = file_name[12:15]  # 457
    json_path = os.path.join(base_path, prefix, middle, suffix, "execution", f"{file_name}_report.json")
    txt_path = os.path.join(base_path, prefix, middle, suffix, "execution", "input_screen.txt")
    return json_path, txt_path

def analyze_file(file_name):
    json_path, txt_path = get_file_paths(file_name)

    # check if the JSON report file exists
    if not os.path.exists(json_path):
        print(f"⚠️ Warning: QC report not found at {json_path}")
        return

    print("📊 QC report loaded.")

    try:
        with open(json_path, 'r') as file:
            data = json.load(file)

        # extract  mapped reads ratio
        total_reads = data["TotalReads"]
        mapped_reads_ratio = data["Data"]["MappedReads"][0]
        reads_unaligned = 100 - (mapped_reads_ratio * 100)

        # check unaligned reads percentage
        if reads_unaligned > 40.0:
            print(f"⚠️ Warning: Reads unaligned out of total reads ({reads_unaligned:.2f}%) exceeds 40%")
        else:
            print(f"✅ Reads unaligned out of total reads ({reads_unaligned:.2f}%) is acceptable (<40%).")

        # sum the map quality values from 0 to 29
        map_quality_distribution = data["Data"]["MappingQualityDistribution"]
        low_map_quality_sum = sum(count for quality, count in map_quality_distribution if quality <= 29)

        # Calculate percentage
        low_map_quality_percentage = (low_map_quality_sum / total_reads) * 100

        # Check if mapq percentage is above 5%
        if low_map_quality_percentage > 5.0:
            print(
                f"⚠️ Warning: Percentage of reads with map quality <30 ({low_map_quality_percentage:.2f}%) exceeds 5%")
        else:
            print(
                f"✅ Percentage of reads with map quality higher than 30 ({low_map_quality_percentage:.2f}%) is acceptable (≤5%).")

    except json.JSONDecodeError:
        print(f"Error: Unable to parse JSON file at {json_path}")
    except KeyError as e:
        print(f"Error: Missing key in JSON data: {e}")

    # check if species input_screen.txt file exists
    if not os.path.exists(txt_path):
        print(f"⚠️ Warning: Fasq_screen specie's file input_screen.txt not found at {txt_path}")
        return

    try:
        with open(txt_path, 'r') as file:
            lines = file.readlines()

            # process the third (first) line
            if len(lines) >= 3:
                third_line = lines[2].strip()
                cells = third_line.split("\t")

                # calculate percentage from the third cell
                if len(cells) >= 4:
                    percent_unmapped = float(cells[3])
                    percent_mapped = 100 - percent_unmapped

                    # check if the first cell in the third line contains "Human"
                    first_cell = cells[0]
                    if first_cell == "Human":
                        print(f"✅ Most reads mapped to Human ({percent_mapped:.2f}%).")
                    else:
                        print(f"️🚨 Caution: Most reads mapped to {first_cell} ({percent_mapped:.2f}%).")


    except Exception as e:
        print(f"Error processing input_screen.txt species file: {e}")

# example usage, modify the file to check the flags for that file, this will change to an option parser
file_name = "EGAF00005432457"
analyze_file(file_name)
