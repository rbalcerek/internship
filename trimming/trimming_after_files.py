import os
import subprocess
import re

input_dir = "."
output_dir = "cleaned_files"
os.makedirs(output_dir, exist_ok=True)

file_pairs = {}
pattern = re.compile(r"^([\w-]+)_([12])\.fq$")
for file_name in os.listdir(input_dir):
    match = pattern.match(file_name)
    if match:
        sample_id, pair_id = match.groups()
        if sample_id not in file_pairs:
            file_pairs[sample_id] = {}
        file_pairs[sample_id][pair_id] = os.path.join(input_dir, file_name)

for sample_id, files in file_pairs.items():
    if "1" in files and "2" in files:
        input_r1 = files["1"]
        input_r2 = files["2"]
        final_output_r1 = os.path.join(output_dir, f"trimmed_{sample_id}_1.fastq")
        final_output_r2 = os.path.join(output_dir, f"trimmed_{sample_id}_2.fastq")

        command = [
            "bbduk.sh",
            f"in1={input_r1}",
            f"in2={input_r2}",
            f"out1={final_output_r1}",
            f"out2={final_output_r2}",
            "ref=adapters.fa",
            "ktrim=f",
            "k=23",
            "hdist=1",
            "minlen=50",
            "overwrite=true",
            "threads=50"
        ]
        subprocess.run(command)
        print(f"trimming done for {sample_id}")
    else:
        print(f"Skipping {sample_id}, incomplete pair.")
