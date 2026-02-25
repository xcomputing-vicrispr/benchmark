import numpy as np
from pydantic import BaseModel
from Bio.Seq import Seq
import re, os, faiss, pickle
from Bio.Seq import Seq
from Bio import SeqIO
from collections import defaultdict
import pandas as pd
import re, os, faiss, pickle

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


DATA_DIR = "./data"

def get_paths(genome_name: str):
    data_dir = DATA_DIR
    
    paths = {
        "data_dir": data_dir,
        "fasta_path": os.path.join(data_dir, f"{genome_name}.fa"),
        "anno_path": os.path.join(data_dir, f"{genome_name}.gff3"),
        "filtered_anno_path": os.path.join(data_dir, f"{genome_name}_only_genes.gff3"),
        "pkl_path": os.path.join(data_dir, f"{genome_name}.pkl"),
        "ori_pkl_path": os.path.join(data_dir, f"{genome_name}ori.pkl"),
        "name_file": os.path.join(data_dir, f"gw_{genome_name}.csv"),
    }
    return paths

def find_sgRNAs_with_PAM(seq: str, chrom: str, pam: str, seed_len=9, seq_len=20):
    res = []
    pam_pattern_xuoi = pam_to_regex(pam)
    pattern_xuoi = re.compile(f"(?={pam_pattern_xuoi})")
    seq_upper = seq.upper()
    for m in pattern_xuoi.finditer(seq_upper):
        pam_start = m.start() + 1
        kx = seq_upper[m.start() - seq_len : m.start() + len(pam)]
        if kx.startswith("CC") and kx.endswith("GG"):
           continue
        res.append({
            "chrom": f'{chrom}:{pam_start-seq_len}-{pam_start-1}',
            "strand": "+",
            "seq": seq_upper[m.start() - seq_len : m.start() + len(pam)],
            "seq_no_pam": seq_upper[m.start() - seq_len : m.start()],
            "seed": seq_upper[m.start() - seed_len : m.start()]
        })
    print(len(res))
    pam_pattern_nguoc = pam_to_regex(str(Seq(pam).reverse_complement()))
    pattern_nguoc = re.compile(f"(?={pam_pattern_nguoc})")
    
    for m in pattern_nguoc.finditer(seq_upper):
        pam_start = m.start() + 1
        seed = seq_upper[m.start() + len(pam): m.start() + len(pam) + seed_len]
        seed = str(Seq(seed).reverse_complement())

        seq = str(Seq(seq_upper[m.start() : m.start() + len(pam) + seq_len]).reverse_complement())
        seq_no_pam = str(Seq(seq_upper[m.start() + len(pam): m.start() + len(pam) + seq_len]).reverse_complement())
        

        res.append({
            "chrom": f'{chrom}:{pam_start + len(pam)}-{pam_start + len(pam) + seq_len - 1}',
            "strand": "-",
            "seq": seq,
            "seq_no_pam": seq_no_pam,
            "seed": seed
        })
    print(len(res))
    return res

NUC2ONEHOT = {
    'A': 0b1000,
    'C': 0b0100,
    'G': 0b0010,
    'T': 0b0001
}

def seq_to_bits(seq):
    bits = 0
    for base in seq:
        bits = (bits << 4) | NUC2ONEHOT[base]
    
    # Convert to uint8 array
    bit_length = len(seq) * 4
    byte_length = (bit_length + 7) // 8
    return np.frombuffer(bits.to_bytes(byte_length, byteorder='big'), dtype=np.uint8)

mapping = {
    'A': np.array([1,0,0,0], dtype=np.float32),
    'C': np.array([0,1,0,0], dtype=np.float32),
    'G': np.array([0,0,1,0], dtype=np.float32),
    'T': np.array([0,0,0,1], dtype=np.float32)
}

def hamming(s1: str, s2: str):
    x = 0
    return sum(char1 != char2 for char1, char2 in zip(s1, s2))

def pam_to_regex(pam):
    return "".join(IUPAC_MAP.get(base, base) for base in pam.upper())

def one_hot_encode(seq):
    return np.concatenate([mapping[base] for base in seq.upper()])

def find_sgRNAs(seq: str, chrom: str, pam="NGG", guide_len=20):
    res = []
    pam_pattern_xuoi = pam_to_regex(pam)
    pattern_xuoi = re.compile(f"(?={pam_pattern_xuoi})")
    seq_upper = seq.upper()
    for m in pattern_xuoi.finditer(seq_upper):
        pam_start = m.start() + 1
        res.append({
            "chrom": f'{chrom}:{pam_start-guide_len}-{pam_start-1}',
            "strand": "+",
            "seq_no_pam": seq_upper[m.start() - guide_len : m.start()],
        })

    print(len(res))

    pam_pattern_nguoc = pam_to_regex(str(Seq(pam).reverse_complement()))
    pattern_nguoc = re.compile(f"(?={pam_pattern_nguoc})")
    
    for m in pattern_nguoc.finditer(seq_upper):
        pam_start = m.start() + 1
        seq_no_pam = str(Seq(seq_upper[m.start() + len(pam): m.start() + len(pam) + guide_len]).reverse_complement())
        res.append({
            "chrom": f'{chrom}:{pam_start + len(pam)}-{pam_start + len(pam) + guide_len - 1}',
            "strand": "-",
            "seq_no_pam": seq_no_pam,
        })
    print(len(res))
    return res

def load_all_metadata_from_pkl(file_path):
    all_metadata = []
    with open(file_path, 'rb') as f:
        while True:
            try:
                batch_metadata = pickle.load(f)
                all_metadata.extend(batch_metadata)
            except EOFError:
                break
    return all_metadata

def load_filtered_genes(gff_path):
    genes_list = []
    with open(gff_path, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split('\t')
            if len(parts) < 9:
                continue
            
            genes_list.append({
                'chrom': parts[0],
                'start': int(parts[3]),
                'end': int(parts[4]),
                'strand': parts[6],
                'attributes': parts[8]
            })
    return genes_list

def buildFaissIndex(genome_name: str, PAM: str, sgRNA_length: int):

    paths = get_paths(genome_name)
    fasta_path = paths["fasta_path"]
    pkl_path = paths["pkl_path"]
    ori_pkl_path = paths["ori_pkl_path"]
    anno_path = paths["anno_path"]
    filtered_anno_path = paths["filtered_anno_path"]

    with open(anno_path, 'r') as f_in, open(filtered_anno_path, 'w') as f_out:
        for line in f_in:
            if line.startswith('#'):
                continue
            
            parts = line.strip().split('\t')
            if len(parts) > 2 and parts[2] == 'gene':
                f_out.write(line)

    name = f"{genome_name}"

    faiss_path = os.path.join(DATA_DIR, f"{name}.faiss")

    sgRNA_len = sgRNA_length
    dim = sgRNA_len * 4

    # crated basic index + bọc IndexIDMap
    dim_bits = sgRNA_length * 4  # 2 bit/nu (A,C,G,T)
    index_flat = faiss.IndexBinaryFlat(dim_bits)
    index = faiss.IndexBinaryIDMap(index_flat)

    all_sgRNAs_metadata = []
    current_id = 0

    # add theo batch
    print("Faiss index building started")
    with open(ori_pkl_path, 'ab') as f_meta:
        for record in SeqIO.parse(fasta_path, "fasta"):
            chrom = record.id
            seq = str(record.seq)

            sgRNAs_for_chrom = find_sgRNAs(seq, chrom, PAM, sgRNA_length)
            valid_sgRNAs = [g for g in sgRNAs_for_chrom if not any(nu not in 'ACGT' for nu in g["seq_no_pam"]) and len(g["seq_no_pam"]) == sgRNA_length]
            
            if not valid_sgRNAs:
                print(f"Not found any sgRNAs at {chrom}")
                continue

            # Save this batch metadata 
            pickle.dump(valid_sgRNAs, f_meta) 
            
            # Encode
            vectors_for_chrom = np.vstack([seq_to_bits(g["seq_no_pam"]) for g in valid_sgRNAs]).astype(np.uint8)

            # Unique ID
            ids_for_chrom = np.arange(current_id, current_id + len(valid_sgRNAs)).astype(np.int64)
            
            # Add vector and ID batch vào
            index.add_with_ids(vectors_for_chrom, ids_for_chrom)
            
            # Update ID for next batch
            current_id += len(valid_sgRNAs)

            print(f"Done {chrom}. {len(valid_sgRNAs)}")
        
    # --- 3. Save Index and Metadata ---
    print(f"{index.ntotal}")
    faiss.write_index_binary(index, faiss_path)

    print(f"faiss index {faiss_path}")

    index = None

    #Start merge metadata
    all_sgRNAs_merged = load_all_metadata_from_pkl(ori_pkl_path)
    print(f"Merge {len(all_sgRNAs_merged)} bản")

    #New pkl
    with open(pkl_path, 'wb') as f:
        pickle.dump(all_sgRNAs_merged, f)
    print(f"New file {pkl_path}")

    os.remove(ori_pkl_path)
    print(f"Delete file {ori_pkl_path}")

    return

def queryFaissIndex(genome_name: str, PAM: str, sgrna_length: int, seed_length: int, hamming_distance: int, flank_up: int, flank_down: int):

    paths = get_paths(genome_name)
    fasta_path = paths["fasta_path"]
    pkl_path = paths["pkl_path"]
    output_file_path = paths["name_file"]
    filtered_anno_path = paths["filtered_anno_path"]

    name = f"{genome_name}"
    faiss_path = os.path.join(DATA_DIR, f"{name}.faiss")
    index_loaded = faiss.read_index_binary(faiss_path)


    with open(pkl_path, 'rb') as f:
        sgRNAs_loaded = pickle.load(f)

    PAM = PAM
    seed_len = seed_length
    all_sgRNAs = []
    for record in SeqIO.parse(fasta_path, "fasta"):
            chrom = record.id
            seq = str(record.seq)

            sgRNAs_for_chrom = find_sgRNAs_with_PAM(seq, chrom, PAM, seed_len)
            all_sgRNAs.extend(sgRNAs_for_chrom)
    print("Tat ca la co:", len(all_sgRNAs))
   # sgRNAs_for_chrom = find_sgRNAs_with_PAM_v2(seq, chrom, PAM, seed_len)
    seed_strand_to_indices = defaultdict(list)

    for i, g in enumerate(all_sgRNAs):
        key = (g["seed"])
        seed_strand_to_indices[key].append(i)

    filtered = []
    for seed, idxs in seed_strand_to_indices.items():
        if len(idxs) == 1:
            only_index = idxs[0]
            filtered.append(all_sgRNAs[only_index])

    print("After seed filter", len(filtered))

    seq_list = [
        {
            "seq_no_pam": g["seq_no_pam"],
            "chrom": g["chrom"],
            "strand": g["strand"],
            "seq": g["seq"],
            "kc": ""
        }
        for g in filtered
    ]
    print(len(seq_list))
    count_cc_gg = sum(
        1 for g in seq_list 
        if g["seq"].startswith("CC") and g["seq"].endswith("GG")
    )
    print(len(seq_list))
    i = 0
    encoded_queries = []
    total = len(seq_list)
    for query in seq_list:
        sequence = query["seq_no_pam"]
        if any(nuc not in "ACGT" for nuc in sequence):
            continue
        if len(sequence) != sgrna_length:
            continue
        vector = seq_to_bits(sequence)
        encoded_queries.append(vector)
        i = i + 1
        if i % 100 == 0:
            print(f"Progress embedding: {i}/{total} ({i/total:.2%})", end='\r')
    query_vectors = np.vstack(encoded_queries).astype(np.uint8)

    print("encoded xong")
    sgRNAs_final = []
    dl = []
    vt = set()
    D, I = index_loaded.search(query_vectors, k=3)
    hammingDistance = hamming_distance

    print("Done")
    newdt = []
    i = 0
    for q, idx_row, dist_row in zip(seq_list, I, D):
        idx = idx_row[1]

        s1 = sgRNAs_loaded[idx_row[0]]  # xâu gần nhất
        s2 = sgRNAs_loaded[idx_row[1]]  # xâu thứ 2
        s3 = sgRNAs_loaded[idx_row[2]]  # xâu thứ 3
        smss1 = s1["seq_no_pam"]
        smss2 = s2["seq_no_pam"]
        smss3 = s3["seq_no_pam"]
        dist0 = hamming(q["seq_no_pam"], smss1)
        dist1 = hamming(q["seq_no_pam"], smss2)
        dist2 = hamming(q["seq_no_pam"], smss3)

        if (dist1 < hammingDistance):
            continue
        sgRNAs_final.append(
            {
                "seq_no_pam": q["seq_no_pam"],
                "chrom": q["chrom"],
                "strand": q["strand"],
                "seq": q["seq"],
                "kc": str(dist0) + "/" + str(dist1) + "/" + str(dist2),
                "details": str(s1) + "," + str(s2) + ", " + str(s3)
            }
        )
        i = i + 1
        if i % 100 == 0:
            print(f"Progress filter distance: {i}/{total} ({i/total:.2%})", end='\r')

    print("sgRNAs after filter:", len(sgRNAs_final))
    output_file = os.path.join(DATA_DIR, 'testing.csv')
    print("Done")

    genes_data = load_filtered_genes(filtered_anno_path)
    results = []
    for j, sgRNA in enumerate(sgRNAs_final, 1):

        print(j)
        chrom_info = sgRNA["chrom"].split(':')
        chrom = chrom_info[0]
        position = chrom_info[1]
        sg_start = int(position.split('-')[0])
        sg_end = int(position.split('-')[1])
        
        for gene in genes_data:
            if gene['chrom'] != chrom:
                continue
                
            if gene['strand'] == '+':
                tss_start = gene['start'] - flank_up
                tss_end = gene['start'] + flank_down
            else:
                tss_start = gene['end'] - flank_down
                tss_end = gene['end'] + flank_up

            if not (sg_end < tss_start or sg_start > tss_end):
                print(f"Match: {sgRNA['seq']} | Pos: {gene['start']}-{gene['end']} | {j}/{len(sgRNAs_final)}")
                
                results.append({
                    'sgRNA_seq': sgRNA["seq"],
                    'sgRNA_loc': sgRNA["chrom"],
                    'Strand': sgRNA["strand"],
                    'gene_start': gene['start'],
                    'gene_end': gene['end'],
                    'gene_strand': gene['strand'],
                    'kc': sgRNA["kc"],
                    'details': sgRNA["details"],
                    'gene_data': gene['attributes']
                })                
                break
    df = pd.DataFrame(results)
    df.to_csv(output_file_path, index=False)
    print("Saved results")
    print(count_cc_gg)

def cleanFaissIndex(genome_name: str):
    paths = get_paths(genome_name)
    pkl_path = paths["pkl_path"]

    name = f"{genome_name}"
    faiss_path = os.path.join(DATA_DIR, f"{name}.faiss")

    if os.path.exists(faiss_path):
        os.remove(faiss_path)
        print(f"Deleted Faiss index file: {faiss_path}")

    if os.path.exists(pkl_path):
        os.remove(pkl_path)
        print(f"Deleted pickle file: {pkl_path}")

genome_name = "ecolik12"

#genome name, PAM, sgRNA length, length seed region, mismatches, TSS up flanking allowed gene , TSS down flanking allowed gene
buildFaissIndex(genome_name, "NGG", 20)
queryFaissIndex(genome_name, "NGG", 20, 11, 3, 500, 500)
cleanFaissIndex(genome_name)