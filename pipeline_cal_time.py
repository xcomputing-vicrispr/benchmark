from pydantic import BaseModel
import os
import glob
import pandas as pd


from collections import defaultdict
from experiment_functions_500 import *

import os, json, warnings
warnings.filterwarnings("ignore")


df = pd.read_csv("./data/cfdtable.csv")
df[["rna", "dna"]] = df["pair"].str.extract(r"r([ATGC]):d([ATGC])")
df["pos"] = df["pos"].astype(int)


DATA_DIR = "./data"
OUTPUT_DIR = "./data"

class IndexComputingSession(BaseModel):
    idfile: str

class GeneInfo:
    def __init__(self, id: str, strand: str, start: int, end: int, seq: str):
        self.id = id
        self.strand = strand
        self.start = start
        self.end = end
class Data(BaseModel):
    gene_name: str
    species: str

IUPAC_MAP = {
    "A": ["A"],
    "C": ["C"],
    "G": ["G"],
    "T": ["T"],
    "W": ["A", "T"],
    "R": ["A", "G"],
    "M": ["A", "C"],
    "Y": ["C", "T"],
    "S": ["G", "C"],
    "K": ["G", "T"],
    "B": ["C", "G", "T"],
    "D": ["A", "G", "T"],
    "H": ["A", "C", "T"],
    "V": ["A", "C", "G"],
    "N": ["A", "C", "G", "T"]
}


primerDefault = "AAAAAAAAAAAAAAAAAA"
mlseqDefault = "CCACCAGGTGGTTGGTGATTTTGGCGGGGG"
lindelDefault = "CAATCATCGCCACCAGGTGGTTGGTGATTTTGGCGGGGGCAGAGAGGACGGTGGCCACCT"
scaffold_default = "GUUUUAGAGCUAGAAAUAGCAAGUUAAAAUAAGGCUAGUCCGUUAUCAACUUGAAAAAGUGGCACCGAGUCGGUGCUUUU"
CONST_GAP = 500





from collections import defaultdict

def process_exon_list(file_path):
    gene_groups = defaultdict(list)
    
    try:
        with open(file_path, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                
                parts = line.split()
                if len(parts) >= 5:
                    gene_name = parts[0]
                    chrom = parts[1]
                    start = int(parts[3]) + 10000000
                    end = int(parts[4]) + 10000000
                    
                    formatted_coord = f"{chrom}:{start}-{end}"
                    gene_groups[gene_name].append(formatted_coord)
                    
        return gene_groups

    except FileNotFoundError:
        print(f"File not found at {file_path}")
        return {}

file_path = "./exon_list_500k.txt"
#file_path = "./exon_list_1m.txt"
#file_path = "./exon_list_5m.txt"
# file_path = "./exon_list_full.txt"
grouped_data = process_exon_list(file_path)



def convert_json_to_csv(json_data, output_filename):
    try:
        with open(json_data, 'r', encoding='utf-8') as f:
            data = json.load(f)
    except FileNotFoundError:
        return
    
    metadata = data[0]
    sequences = data[1:]
    rows = []
    combined_list = []
    for seq in sequences:
        if seq.get("GC Content", 0) < 25:
            continue

        if seq.get("GC Content", 0) > 75:
            continue

        if seq.get("mm0", 0) > 0:
            continue
        
        if seq.get("mm1", 0) > 0:
            continue

        if seq.get("mfe2", 0) <= -24.5:
            continue
        
        if seq.get("mlScore", 0) <= 0.5:
            continue
        
        if seq.get("rs3", 0) <= 0:
            continue

        seq_str = seq.get("sequence")


        location_str = seq.get("location", "")
        start_val = int(location_str.split(":")[1]) - 10000000
        end_val = int(start_val) + 22
        print(location_str)
        selected_data = {
            "name": "vicrispr",
            "sequence": seq.get("sequence")[:-3] + "...",
            "start": start_val,
            "end": end_val,
            "strand": seq.get("strand")
        }
        rows.append(selected_data)

    if rows:
        df = pd.DataFrame(rows)
        df = df[["name", "sequence", "start", "end", "strand"]]
        df.to_csv(output_filename, index=False, encoding='utf-8-sig')
        print(f"Done, saved at {output_filename}")

def merge_csv_files(folder_path, output_file):

    all_files = glob.glob(os.path.join(folder_path, "data*.csv"))
    
    if not all_files:
        print("Not found file")
        return

    print(f"Found {len(all_files)} file. Merging...")

    list_df = []
    for filename in all_files:
        df = pd.read_csv(filename, index_col=None, header=0, encoding='utf-8-sig')
        list_df.append(df)

    merged_df = pd.concat(list_df, axis=0, ignore_index=True)
    merged_df.to_csv(output_file, index=False, encoding='utf-8-sig')
    print(f"Done! Saved at: {output_file}")


k = 0 
for gene, coords in grouped_data.items():
    print(f"Processing Gene: {gene} with {len(coords)} exons")
    CoordinateComputing(gene, coords)
    convert_json_to_csv(f"./data/vcp{gene}.json", f"./data/data{gene}.csv")
    target_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "./data"))

    pattern = os.path.join(target_dir, f"{gene}*")
    files_to_delete = glob.glob(pattern)

    print(f" Deleting gene: {gene}...")

    for file_path in files_to_delete:
        if os.path.isfile(file_path):
            try:
                os.remove(file_path)
                print(f"--- Deleted: {file_path}")
            except Exception as e:
                print(f"--- Cannot delete file {file_path}: {e}")
    break

folder_path = './data'
output_name = './data/data500k_filter.csv'

merge_csv_files(folder_path, output_name)

def convert_to_normalise(input_csv, output_file):
    df = pd.read_csv(input_csv)

    if 'sequence' in df.columns:
        df['sequence'] = df['sequence'].astype(str).str.upper()
    
    output_df = df[['name', 'sequence', 'start', 'end', 'strand']]

    output_df.to_csv(output_file, header=False, index=False, encoding='utf-8')
    
    print(f"convert to .normalised file done, ready for Bradford J et al. script (mm10db-rejects-accepted-by-other-tools.py): {output_file}")

convert_to_normalise("./data/data500k_filter.csv", "./normalised/exon-only/500k/vicrispr.normalised")

