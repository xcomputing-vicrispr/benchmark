from pydantic import BaseModel
import subprocess
import sys
from pathlib import Path
from itertools import product
import shutil
import os, math
import glob, math
import pandas as pd

from collections import defaultdict

import re, os, json, random, string, subprocess, warnings
warnings.filterwarnings("ignore")
from rs3.seq import predict_seq
from Bio.Seq import Seq
import re, os, json, random, string




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
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "W": "[AT]",
    "R": "[AG]",
    "M": "[AC]",
    "Y": "[CT]",
    "S": "[GC]",
    "K": "[GT]",
    "B": "[CGT]",
    "D": "[AGT]",
    "H": "[ACT]",
    "V": "[ACG]",
    "N": "[ACGT]"
}

IUPAC_MAP_WBT = {
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




def micro(seq: str) -> float:

    res = 0
    l = len(seq)
    for i in range(l):
        if seq[i] in ['A', 'T']:
            res += 1
        else:
            res += 2
    return res



def get_fasta_from_twobit(twobit_file: str, chromosome: str, chrstart: int, chrstop: int) -> str:

    PARENT_DIR = os.path.dirname((os.path.abspath(__file__)))
    DATA_DIR = os.path.join(PARENT_DIR, "data")

    twoBitToFa_path = os.path.join(DATA_DIR, "twoBitToFa")

    if(int(chrstart) > int(chrstop)):
        tmp = chrstart
        chrstart = chrstop
        chrstop = tmp

    command_list = [
        twoBitToFa_path,
        f"{twobit_file}:{chromosome}:{chrstart}-{chrstop}",
        "/dev/stdout"
    ]

    result = subprocess.run(
        command_list,
        shell=False,
        capture_output=True,
        text=True,
        check=True,
        cwd=DATA_DIR,
    )
    parts = result.stdout.split('\n', 1)
    return parts[1].replace('\n', '')

def random_string():
    chars = string.ascii_letters + string.digits
    return ''.join(random.choices(chars, k=10))


def save_sgRNA_list(idd: str, data: dict, gene_name: str, spec: str, pam: str, sgRNA_len: int,
                        status: str = "processing", log = "", stage = 1, gene_strand="+"):
    if stage == 0:
        filename = "vcp" + idd + ".json"
        path = os.path.join(DATA_DIR, filename)

        meta_obj = {"name": gene_name, "spec": spec, "pam": pam, "sgRNA_len": sgRNA_len, "gene_strand": gene_strand, 
                     "status": status, "log": log}
        
        with open(path, 'r', encoding='utf-8') as f:
            data_in_file = json.load(f)
            data_in_file[0] = meta_obj
        with open(path, 'w', encoding='utf-8') as f:
            json.dump(data_in_file, f, indent=4, ensure_ascii=False)
        return

    if idd == "unk":
        x = random_string()
        filename = "vcp" + x + ".json"
        path = os.path.join(DATA_DIR, filename)

        meta_obj = {"name": gene_name, "spec": spec, "pam": pam, "sgRNA_len": sgRNA_len, "gene_strand": gene_strand,
                     "status": status, "log": log}
        data_with_meta = [meta_obj] + data
    elif stage == 1:
        x = idd
        filename = "vcp" + x + ".json"
        path = os.path.join(DATA_DIR, filename)

        meta_obj = {"name": gene_name, "spec": spec, "pam": pam, "sgRNA_len": sgRNA_len, "gene_strand": gene_strand,
                     "status": status, "log": log}
        data_with_meta = [meta_obj] + data
    
    with open(path, "w") as f:
        json.dump(data_with_meta, f, indent=2)
    return x


def consensus(regions):
    if not regions:
        return []
    if len(regions) == 1:
        regions = list(regions)
        return [(regions[0][0], regions[0][1], regions[0][2])]    

    regions = sorted(list(regions), key=lambda x: (x[0], x[1], x[2], x[3]))

    merged = []
    curr_chr, curr_start, curr_end = regions[0][0], regions[0][1], regions[0][2]

    for c, s, e, _ in regions[1:]:
        if c == curr_chr and s <= curr_end:  
            curr_end = max(curr_end, e)
        else:
            merged.append((curr_chr, curr_start, curr_end))
            curr_chr, curr_start, curr_end = c, s, e
    merged.append((curr_chr, curr_start, curr_end))
    return merged


def count_permu_IUPAC(s: str) -> int:
    iupac_counts = {
        "A": 1, "C": 1, "G": 1, "T": 1,
        "R": 2, "Y": 2, "S": 2, "W": 2,
        "K": 2, "M": 2,
        "B": 3, "D": 3, "H": 3, "V": 3,
        "N": 4
    }
    res = 1
    for ch in s.upper():
        res *= iupac_counts.get(ch, 1)
    return res

def xuly(line: str, datafile, pam_name: str, ll: int, off_target: int = 0, num_of_mismatches: int = 0):

    res = []
    parts = line.strip().split('\t')

    permu_num = count_permu_IUPAC(pam_name)

    x = int(int(parts[0]) / permu_num)
    if len(parts) < 8:
        datafile[x]["mm0"] += 1
        return -1, -1

    mismatch_info = parts[7]
    mismatch_info_parse = mismatch_info.split(",")
    positions = []
    
    for item in mismatch_info_parse[0:]:
        pos_str = item.split(":")[0]  
        positions.append(int(pos_str))
    if any(pos >= ll for pos in positions):
        return -1, -1
    
    if off_target == 1:
        mismatch_count_in_seed = sum(1 for pos in positions if pos > 9)
        if mismatch_count_in_seed > num_of_mismatches:
            return -1, -1

    if not mismatch_info:
        mis_match_num = 0
    else:
        mis_match_num = mismatch_info.count(":")

    res = f"{parts[2]}:{parts[3]}", parts[4], mismatch_info
    print(f"{x},{parts[1]},{parts[2]},{parts[3]},{parts[4]},{mismatch_info}")

    mutated_seq = ""

    res = f"{parts[2]}:{parts[3]}", mutated_seq, mismatch_info
    if mis_match_num == 1:
        datafile[x]["mm1"] += 1
    elif mis_match_num == 2:
        datafile[x]["mm2"] += 1
    elif mis_match_num == 3:
        datafile[x]["mm3"] += 1

    return x, res


def getMMDetails(bed_file_path: str, json_final_file: str):
    mm_details = {}

    with open(bed_file_path, 'r', encoding='utf-8') as bed_file:
        for line in bed_file:
            fields = line.strip().split('\t')
            if len(fields) < 5:
                continue
            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])
            region_type = fields[4]

            key = f"{chrom}:{start}"
            mm_details[key] = region_type

    with open(json_final_file, 'r', encoding='utf-8') as f:
        data = json.load(f)
    header = data[0]  
    data = data[1:]
    for entry in data:
        bowtie_details = entry.get('bowtie_details', '')
        
        if not bowtie_details or bowtie_details.strip() == "":
            entry['mismatch_region'] = ""
            continue
        
        regions = []
        print(bowtie_details)
        offtargets = bowtie_details.strip().rstrip(';').split(';')
        for offtarget in offtargets:
            offtarget = offtarget.strip()
            if not offtarget:
                continue
            
            chrom_pos = offtarget.split(',,')[0].strip()
            region = mm_details.get(chrom_pos, "unknown")
            regions.append(region)
        
        # Format: "exon, intron, intergenic"
        print(regions)
        entry['mismatch_region'] = ";".join(regions)


    final_data = [header] + data
    with open(json_final_file, 'w', encoding='utf-8') as f:
        json.dump(final_data, f, indent=4)
    

    return


def getMMRegion(name: str, idfile: str):

    raw_bed_file = os.path.join(DATA_DIR, f"{idfile}_raw.bed")
    bed_file = f"{idfile}_sorted.bed"    


    possible_ext = [".gff3", ".gff", ".gtf"]
    found_ext = None

    for ext in possible_ext:
        file_path = os.path.join(DATA_DIR, f"{name}{ext}")
        if os.path.exists(file_path):
            found_ext = ext
            break

    if not found_ext:
        print(f"File annotation not found: {name}.gff3 / .gff / .gtf")
        return None

    exon_file = f"{name}_exons.sorted{found_ext}"
    gene_file = f"{name}_genes.sorted{found_ext}"
    result_file = f"{idfile}_mm_annotation.bed"

    exon_file = os.path.join(DATA_DIR, exon_file)
    gene_file = os.path.join(DATA_DIR, gene_file)
    bed_file = os.path.join(DATA_DIR, bed_file)
    result_file = os.path.join(DATA_DIR, result_file)



    cmds = [

        f"sort -k1,1V -k2,2n {raw_bed_file} > {bed_file}",


        f"bedtools intersect -a {bed_file} -b {exon_file} -wa | "
        f"awk 'BEGIN{{OFS=\"\\t\"}} {{print $1, $2, $3, $4, \"exon\"}}' > annotated.part1.exonic",

        f"bedtools intersect -a {bed_file} -b {exon_file} -v > remaining_regions.1.bed",

        f"bedtools intersect -a remaining_regions.1.bed -b {gene_file} -wa | "
        f"awk 'BEGIN{{OFS=\"\\t\"}} {{print $1, $2, $3, $4, \"intron\"}}' > annotated.part2.intronic",

        f"bedtools intersect -a remaining_regions.1.bed -b {gene_file} -v > remaining_regions.2.bed",

        f"awk 'BEGIN{{OFS=\"\\t\"}} {{print $1, $2, $3, $4, \"intergenic\"}}' remaining_regions.2.bed > annotated.part3.intergenic",

        f"cat annotated.part1.exonic annotated.part2.intronic annotated.part3.intergenic | "
        f"sort -k4,4V -k2,2n -u > {result_file} && "
        f"rm remaining_regions.*.bed annotated.part* && "
        f"rm {bed_file}"
    ]

    for cmd in cmds:
        try:
            subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")
        except subprocess.CalledProcessError as e:
            print(f" Error: {e}")
            return None

    json_file_path = os.path.join(DATA_DIR, f"vcp{idfile}.json")
    getMMDetails(result_file, json_file_path)

    return


def gc_content(seq: str, len_without_pam: int) -> float:
    seq = seq[0:len_without_pam]
    seq = seq.upper()
    gc_count = seq.count('G') + seq.count('C')
    return (gc_count / len(seq)) * 100

def fold_rna(seq: str):
    result = subprocess.run(
        ["RNAfold"],
        input=seq.encode(),
        stdout=subprocess.PIPE,
        stderr=subprocess.DEVNULL
    )
    output = result.stdout.decode().splitlines()
    if len(output) < 2:
        return None, None
    structure = output[1].split(" ")[0]
    mfe = float(output[1].split("(")[-1].replace(")", ""))
    return structure, mfe

def calMicroScore(seq1: str, seq2: str) -> float:

    out_frame_score = 0.0
    frame_score = 0.0
    l = len(seq1)
    sad = set()
    num_out_frame_score = 0
    total_posiblity = 0

    for i in range(l):
        for j in range(3, l - i):
            tmp = seq1[i:i + j]
            for k in range(l - j):
                if tmp == seq2[k:k + j]:

                    ext = j
                    while (i + ext < l and k + ext < l and seq1[i + ext] == seq2[k + ext]):
                        ext += 1
                    tmp = seq2[k:k+ext]

                    key = (i, k, tmp)
                    if key in sad:
                        continue
                    sad.add(key)

                    delta = k - i + l
                    mh_score = math.exp(-delta / 20) * micro(tmp) * 100
                    frame_score += mh_score
                    if (delta) % 3 != 0:
                        out_frame_score += mh_score
                        num_out_frame_score = num_out_frame_score + 1
                    total_posiblity = total_posiblity + 1
                    print(i, k, delta,micro(tmp), tmp, mh_score)
    return num_out_frame_score / total_posiblity
def gop(a: str, b:str) -> str:
    return f"{a}:{b}"


def get_ml_score(seqlist):
    seq_str = json.dumps(seqlist)
    
    conda_env = os.getenv("ML_CONDA_ENV", "env2_conda")    
    
    CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
    TEST_FILE = os.path.abspath(os.path.join(CURRENT_DIR, "get_rs2.py"))
    TEST_DIR = os.path.dirname(TEST_FILE)

    result = subprocess.Popen(
        ["conda", "run", "-n", conda_env, "python", TEST_FILE, seq_str],
        cwd=TEST_DIR,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE
    )
    # ===========================

    stdout, stderr = result.communicate()
    
    print("=== ML Script STDOUT ===")
    print(stdout.decode())
    print("=== ML Script STDERR ===")
    print(stderr.decode())
    
    try:
       # return json.loads(stdout.decode().strip())
        stdout_lines = stdout.decode().splitlines()
        for line in reversed(stdout_lines):
            line = line.strip()
            if line.startswith("[") and line.endswith("]"):
                return json.loads(line)
    except json.JSONDecodeError as e:
        print(f"JSON decode error: {e}")
        print("Raw output that caused error:")
        print(stdout.decode())
        raise



def get_ml_score_azi3(seqlist):
    res = predict_seq(seqlist, sequence_tracr='Hsu2013')
    return res.tolist()


def get_percent_active(pos, rna, dna):
    row = df[(df["pos"] == pos) & (df["rna"] == rna) & (df["dna"] == dna)]
    if not row.empty:
        return float(row.iloc[0]["pa"])
    else:
        return 1

def get_rna(seq):
    ax = {
        'A': 'U',
        'C': 'G',
        'G': 'C',
        'T': 'A'
    }
    return ax.get(seq, '?')


def pam_to_regex(pam):
    return "".join(IUPAC_MAP.get(base, base) for base in pam.upper())
    
def find_pam_positions(seq, pam):
    regex_mo = pam_to_regex(pam)
    pattern = re.compile(f"(?={regex_mo})")
    return [m.start() for m in pattern.finditer(seq.upper())]


def get_cfd_score(line, pam_name: str, ll: int):
    parts = line.strip().split('\t')

    permu_num = count_permu_IUPAC(pam_name)

    x = int(int(parts[0]) / permu_num)
    if len(parts) < 8:
        return x, 0
    
    mismatch_info = parts[7]
    mismatch_info_parse = mismatch_info.split(",")
    positions = []
    for item in mismatch_info_parse[1:]:
        pos_str = item.split(":")[0]  
        positions.append(int(pos_str))
    if any(pos >= ll for pos in positions):
        return -1, 0

    items = mismatch_info.split(',')

    result = []
    posl = []
    for item in items:
        pos_part, change = item.split(':')
        dna, rna = change.split('>')
        result.append((int(pos_part), dna, get_rna(rna)))
        posl.append(int(pos_part))

    scr = 1
    l = len(posl)
    for item in result:
        k = get_percent_active(item[0], item[2], item[1])
        if k == 0:
            k = 0.0025
        scr = scr * k

    scr = scr / l

    return x, scr


def write_sgrna_to_fasta_with_IUPAC(sgrnas, pam, idfile):
    index = 0 
    bases_per_position = [IUPAC_MAP_WBT[base] for base in pam]
    all_variants = list(product(*bases_per_position))
    base_dir="./data"
    output_filename = f"{base_dir}/{idfile}_sgrna_output.fa"
    with open(output_filename, "w") as f:
        for seq in sgrnas:
            prefix = seq[:-len(pam)]
            for all_posi in all_variants:
                pam_sequence = "".join(all_posi)
                full_sequence = prefix + pam_sequence
                f.write(f">{index}\n{full_sequence}\n")
                index += 1  


def indexComputing(idfile: str, off_target: bool = 0, num_of_mismatches: int = 3):
    try:
        seq_list = []
        seq_list_ml = []
        seq_list_lindel = []

        filename = "vcp" + idfile + ".json"
        file_path = os.path.join(OUTPUT_DIR, filename)
        with open(file_path, 'r', encoding='utf-8') as f:
            datafile = json.load(f)

        pam_name = datafile[0]["pam"]
        bowtie_index_file = datafile[0]["spec"] + "_index"
        spec_name = datafile[0]["spec"]
        sgRNA_len = datafile[0]["sgRNA_len"]

        tmp = datafile[0]
        datafile = datafile[1:]

        for i, row in enumerate(datafile):
            seq_list.append(row["sequence"])
            seq_list_ml.append(row["mlseq"])
            seq_list_lindel.append(row["lindel"])

        write_sgrna_to_fasta_with_IUPAC(seq_list, pam_name, idfile)

        if (pam_name == "NGG"):
            ml_score = get_ml_score(seq_list_ml)
            rs3_score = get_ml_score_azi3(seq_list_ml)
        else:
            ml_score = ["N/S"] * len(datafile)
            rs3_score = ["N/S"] * len(datafile)


        pol = 0
        limit_num = 1000
        ss_map = {19: 5000, 20: 4000, 21: 3000, 22: 2000, 23: 50000}
        ss = ss_map.get(len(pam_name) + sgRNA_len, 700)
        pp = count_permu_IUPAC(pam_name)
        limit_num = max(20, ss/pp)
        sg_file = f"{idfile}_sgrna_output.fa"

        para_mm = 3
        if off_target == 0:
            para_mm = num_of_mismatches
        command = [
        "bowtie",
        "-v", str(1),
        "-a",  
        #"-k", str(limit_num),
        "-f",
        "-x", bowtie_index_file,
        sg_file,
        "/dev/stdout"
    ]

        process = subprocess.Popen(
            command,
            cwd=DATA_DIR,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            bufsize=1
        )

        pos_list = []
        count_dict = {}
        max_hits_per_id = 5

        danh_dau = set()
        vcl = []
        m = defaultdict(int)
        ll = len(datafile[0]['sequence']) - len(pam_name)
        for line in process.stdout:
            pol = pol + 1
            line = line.strip()
            idseq, ttmm = xuly(line, datafile, pam_name, ll, off_target, num_of_mismatches)
            id_in_fa = int(int(line.strip().split('\t')[0]))
            m[id_in_fa] += 1
            if (m[id_in_fa] >= limit_num):
                danh_dau.add(id_in_fa/pp)
                vcl.append(id_in_fa)
            if idseq == -1:
                continue
            datafile[idseq]["bowtie_details"] += ",".join(map(str, ttmm)) + "; "



            chrom = line.strip().split('\t')[2]
            start = int(line.strip().split('\t')[3])
            end = start + sgRNA_len - 1
            pos_list.append((chrom, start, end, idseq))

            x, scr = get_cfd_score(line, pam_name, ll)
            if x == -1:
                continue
            datafile[x]["cfdScore"] += scr

            print(str(x) + " " + str(datafile[x]["cfdScore"]))


        process.stdout.close()
        process.wait()

        for i, row in enumerate(datafile):
            if i in danh_dau:
                if str(row.get("mm3", "0")) != "0":
                    row["mm3"] = f">={row['mm3']}"
                elif str(row.get("mm2", "0")) != "0":
                    row["mm2"] = f">={row['mm2']}"
                elif str(row.get("mm2", "0")) != "0":
                    row["mm1"] = f">={row['mm1']}"
            print(
                f"{i}: {row.get('sequence', '')}, "
                f"location={row.get('location', '')}, "
                f"mm0={row.get('mm0', 0)}, "
                f"mm1={row.get('mm1', 0)}, "
                f"mm2={row.get('mm2', 0)}, "
                f"mm3={row.get('mm3', 0)}"
            )

        for i in range(len(datafile)):
            datafile[i]["mlScore"] = ml_score[i]
            datafile[i]["rs3"] = rs3_score[i]

            datafile[i]["cfdScore"] = round(100 / (100 + datafile[i]["cfdScore"]), 2)

        
        
        
        real_datafile = [tmp] + datafile
        with open(file_path, 'w', encoding='utf-8') as f:
            json.dump(real_datafile, f, ensure_ascii=False, indent=2)

        grouped = defaultdict(list)
        for chrom, start, end, idseq in pos_list:
            grouped[idseq].append((chrom, start, end))

        raw_bed_dir = os.path.join(DATA_DIR, f"{idfile}_raw.bed")
        with open(raw_bed_dir, "w") as f:
            for idseq, regions in grouped.items():
                for chrom, start, end in regions:
                    f.write(f"{chrom}\t{start}\t{end}\t{idseq}\n")

        getMMRegion(spec_name, idfile)

        print(limit_num)
        print(pol)
        print(ll)

        return
    except Exception as e:
        print(f"Error in indexComputing: {str(e)}")
        raise


def CoordinateComputing(idd, bo_string):
    print("CoordinateComputing checkpoint")
    GAP = CONST_GAP
    try:
        results = []
        auke = []
        cutting_sites = []
        query = idd
        spec = "nmd_0_mm10"

        PAM = "NGG"
        len_without_pam = 20
        REV_PAM = Seq(PAM)
        REV_PAM = str(REV_PAM.reverse_complement())
        scaffold = scaffold_default
        
        idd = save_sgRNA_list(idd, results, query, spec, PAM, len_without_pam, status="Finding", log="Finding sgRNA candidates")


        twobit_file = spec + ".2bit"
        final_st_end = bo_string
        pam_size = len(PAM)

        final_st_end = []
        for idx, item in enumerate(bo_string):
            chrom, coords = item.split(':')
            start_str, end_str = coords.split('-')
            final_st_end.append((chrom, int(start_str), int(end_str) + 1, idx + 1))

        aval = 0
        for s in final_st_end:
            seq = get_fasta_from_twobit(twobit_file, s[0], (int(s[1]) - 1), int(s[2]))                
            l = len(seq)        
            template = None
            try:
                check_id = int(s[1]) - 1 - 500
                if check_id < 0:
                    check_id = 0
                    aval = 1
                template = get_fasta_from_twobit(twobit_file, s[0], str(check_id), str(int(s[2]) - 1 + 500))
            except Exception as e:
                template = seq

            x = find_pam_positions(seq, REV_PAM)
            for id in x:
                if id + int(s[1]) not in cutting_sites:         
                    cutting_sites.append(id + int(s[1]))    
                else:
                    continue

                pam_seq = seq[id: min(l - 1, id + pam_size + len_without_pam)]        
                if len(pam_seq) != (pam_size + len_without_pam):                  
                    continue     

                pam_seq = Seq(pam_seq)
                rev_comp = str(pam_seq.reverse_complement())


                gcc = gc_content(rev_comp, len_without_pam)
                if gcc < 10 or gcc > 90:
                    continue

                (ss, mfe) = fold_rna(rev_comp)
                (ss2, mfe2) = fold_rna(str(rev_comp + scaffold))

                microScore = primer = mlseq = lindel = None
                clea1 = clea2 = ""
                try:
                    idl = id + GAP + pam_size + 3 - len_without_pam - 10
                    clea1 = template[idl : id + GAP + pam_size + 3]
                    clea2 = template[id + GAP + pam_size + 3 : id + GAP + pam_size + 3 + len_without_pam + 10]
                    microScore = calMicroScore(clea1, clea2)

                    spe_seq = template[id + GAP - 3 : id + GAP + 27]
                    spe_seq = Seq(spe_seq)
                    spe_seq = str(spe_seq.reverse_complement())

                    lindel = Seq(template[id + GAP + pam_size + 3 - 30 : id + GAP + pam_size + 3 + 30])
                    lindel = str(lindel.reverse_complement())

                    primer = template[id + GAP - 260 : id + GAP + 281]

                except Exception as e:
                    microScore = -999999
                    spe_seq = mlseqDefault
                    lindel = lindelDefault
                    primer = primerDefault
                
                if len(spe_seq) != 30:
                    microScore = -999999
                    spe_seq = mlseqDefault
                    lindel = lindelDefault
                    primer = primerDefault
                if id + GAP + pam_size + 3 - 30 < 0 or id + GAP + pam_size + 3 + 30 > len(template):
                    microScore = -999999
                    lindel = lindelDefault
                    primer = primerDefault
                    spe_seq = mlseqDefault

                if primer == "":
                    microScore = -999999
                    lindel = lindelDefault
                    primer = primerDefault
                    spe_seq = mlseqDefault
                

                results.append({
                    "sequence": rev_comp,
                    "location": gop(s[0], str(id + int(s[1]))),
                    "strand": "-",
                    "GC Content": gcc,
                    "Self-complementary": round(mfe, 2),
                    "Primer": primer,
                    "mlseq": spe_seq.upper(),
                    "mm0": -1,
                    "mm1": 0,
                    "mm2": 0,
                    "mm3": 0,
                    "cfdScore": 0,
                    "mlScore": 0,
                    "microScore": microScore,
                    "mmejpre": clea1 + "," + clea2,
                    "Secondary structure with scaffold": f"{ss2}, ({round(mfe2, 2)} kcal/mol)",
                    "bowtie_details": "",
                    "mismatch_region":"",
                    "lindel": lindel,
                    "rs3": "",
                    "mfe2": round(mfe2, 2)
                })
                auke.append(rev_comp)

            x = find_pam_positions(seq, PAM)
            for id in x:
                if id - len_without_pam + int(s[1]) not in cutting_sites:    
                    cutting_sites.append(id - len_without_pam + int(s[1]))    
                else:
                    continue

                pam_seq = seq[max(0, id - len_without_pam):id + pam_size]         
                if len(pam_seq) != (pam_size + len_without_pam):                  
                    continue     

                gcc = gc_content(pam_seq, len_without_pam)
                if gcc < 10 or gcc > 90:
                    continue

                (ss, mfe) = fold_rna(pam_seq)
                (ss2, mfe2) = fold_rna(pam_seq + scaffold)
                microScore = primer = mlseq = lindel = None
                clea1 = clea2 = ""
                try:
                    idl = id + GAP - len_without_pam - 3 - 10
                    clea1 = template[idl: id + GAP - 3]
                    clea2 = template[id + GAP - 3: id + GAP - 3 + len_without_pam]
                    microScore = calMicroScore(clea1, clea2)

                    primer = template[id + GAP - 280: id + GAP + 261]
                    mlseq = template[id + GAP - 20 - 4: id + GAP - 19 + 25]

                    lindel = template[id + GAP - 3 - 30: id + GAP + pam_size + 30]

                except Exception as e:
                    print(e)
                    microScore = -999999
                    mlseq = mlseqDefault
                    lindel = lindelDefault
                    primer = primerDefault
                if len(mlseq) != 30 or aval == 1:
                    microScore = -999999
                    mlseq = mlseqDefault
                    lindel = lindelDefault
                    primer = primerDefault
                if id + GAP - 3 - 30 < 0 or id + GAP + pam_size + 30 > len(template):
                    microScore = -999999
                    lindel = lindelDefault
                    primer = primerDefault
                    spe_seq = mlseqDefault

                if primer == "":
                    microScore = -999999
                    lindel = lindelDefault
                    primer = primerDefault
                    spe_seq = mlseqDefault

                results.append({
                        "sequence": pam_seq,
                        "location": gop(s[0], str(id - len_without_pam + int(s[1]))),
                        "strand": "+",
                        "GC Content": gcc,
                        "Self-complementary": round(mfe,2),
                        "Primer": primer,
                        "mlseq": mlseq.upper(),
                        "mm0": -1,
                        "mm1": 0,
                        "mm2": 0,
                        "mm3": 0,
                        "cfdScore": 0,
                        "mlScore": 0,
                        "microScore": microScore,
                        "mmejpre": clea1 + "," + clea2,
                        "Secondary structure with scaffold": f"{ss2}, ({round(mfe2, 2)} kcal/mol)",
                        "bowtie_details": "",
                        "mismatch_region":"",
                        "lindel": lindel,
                        "rs3": "",
                        "mfe2": round(mfe2, 2)  
                    })
                auke.append(pam_seq)

        idd = save_sgRNA_list(idd, results, query, spec, PAM, len_without_pam, status="calculating-index_processing", log="indexing")

        if len(results) == 0:
            idd = save_sgRNA_list(idd, results, query, spec, PAM, len_without_pam, status="no_result", log="No result available, check your gene name or region", stage=0)
            return

        indexComputing(idd)
        idd = save_sgRNA_list(idd, results, query, spec, PAM, len_without_pam, status="success", log="done", stage=0)
    except Exception as e:
        print(f"Error in GeneNameComputing: {str(e)}")
        raise 
    return