import os
import shutil
import gzip

raw_data_dir = "before"
output_dir = "before_data"

os.makedirs(output_dir, exist_ok=True)

for sample_name in os.listdir(raw_data_dir):
    sample_path = os.path.join(raw_data_dir, sample_name)
    
    if os.path.isdir(sample_path):  
        for file in os.listdir(sample_path):
            if file.endswith("_1.fq.gz") or file.endswith("_2.fq.gz"):
                
                original_file = os.path.join(sample_path, file)
                new_file_name = f"{sample_name}_{file}"
                new_file_path = os.path.join(output_dir, new_file_name)
                
                shutil.copy(original_file, new_file_path)
                print(f"Copied: {original_file} → {new_file_path}")
                
                decompressed_file = new_file_path[:-3]  
                with gzip.open(new_file_path, 'rb') as f_in:
                    with open(decompressed_file, 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
                os.remove(new_file_path)  
                print(f"Unzipped: {decompressed_file}")

print("Done.")
