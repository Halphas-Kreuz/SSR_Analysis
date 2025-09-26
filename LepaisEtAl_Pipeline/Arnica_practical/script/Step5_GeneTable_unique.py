import os

input_file = "../GenotypicTable_nSSR_FullLength_ParameterSet1_sa50_sb10_m15_n20.txt"
output_file = "../new_output/GenotypicTable_fused.txt"

with open(input_file, "r") as fin, open(output_file, "w") as fout:
    lines = fin.readlines()
    headers = lines[0].strip().split("\t")
    fout.write("\t".join(headers) + "\n")
    for line in lines[1:]:
        cells = line.strip().split("\t")
        # Pad cells to match headers length if needed
        if len(cells) < len(headers):
            cells += ["NA"] * (len(headers) - len(cells))
        new_cells = [cells[0]]  # keep the first column (sample/individual name)
        for i in range(1, len(headers)):
            # Find the header for the current pair
            pair_header_idx = i if i % 2 == 1 else i - 1
            pair_header = headers[pair_header_idx]
            value = cells[i]
            if value != "NA" and value != "":
                new_cells.append(f"{pair_header}_{value}")
            else:
                new_cells.append("NA")
        fout.write("\t".join(new_cells) + "\n")
print(f"Transformed table written to {output_file}")