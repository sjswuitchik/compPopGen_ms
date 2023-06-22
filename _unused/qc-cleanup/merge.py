import sys
import csv

# Column headers to keep
headers = ["Sample","Total_Reads","Percent_mapped","Num_duplicates","Percent_properly_paired","Fraction_reads_pass_filter","NumReadsPassingFilters","File"]

# Initialize empty list to store rows of data
rows = []

# Read list of filenames from file
filenames_file = sys.argv[1]
with open(filenames_file, "r") as f:
    filenames = f.read().splitlines()

# Read each input file and store data in rows list
for filename in filenames:
    with open(filename, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            # Keep only the columns with the specified headers
            row_data = {}
            for header in headers:
                if header in row:
                    row_data[header] = row[header]
                elif header == "Total_Reads" and "Total_alignments" in row:
                    row_data[header] = row["Total_alignments"]
                elif header == "Num_duplicates" and "Percent_duplicates" in row:
                    row_data[header] = float(row["Percent_duplicates"]) * float(row["Total_alignments"])
                elif header == "NumReadsPassingFilters" and "Num_filtered_reads" in row:
                     row_data[header] = row["Num_filtered_reads"]                   
                else:
                    row_data[header] = "-"
            row_data["File"] = filename
            rows.append(row_data)

# Write merged data to stdout
writer = csv.DictWriter(sys.stdout, fieldnames=headers)
writer.writeheader()
writer.writerows(rows)

