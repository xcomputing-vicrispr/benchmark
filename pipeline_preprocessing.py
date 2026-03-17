'''
This script read FASTA file and ANNOTATION file and generates optimized intermediate files 
required for querying within the web application.

Notes:
    - Make sure you put a fasta file and annotation file in ./data/preprocessing, we put there ecoli-o127 such an example
    - You can change your input genome at line 34 and 35
    - You can use SBS (github.com/jakeb1996/SBS) to measure the time this script executed

Run:
    - python3 pipeline_preprocessing.py

Input:
    - Genome sequence files (.fna or .fa).
    - Annotation files (.gff, .gff3, or .gtf).
    
Output:
    - A structured set of indexed files and databases ready for webapp integration
'''



import subprocess, time, gffutils, os
import signal
from experiment_functions_500 import *


DATA_DIR = "./data/preprocess"
OUTPUT_DIR = "./data/preprocess"

# faname = "mm10.fna"
# annoname = "mm10.gff"

faname = "ecolio127.fna"
annoname = "ecolio127.gff"

# faname = "osativa70.fa"
# annoname = "osativa70.gff3"


def createMMRegion(anno_file):
    # get .anno. (.gff, .gtf, .gff3)
    ext = os.path.splitext(anno_file)[1]

    pt = anno_file.split(".")[0]

    # output filename for  exon and gene
    exon_out = f"{pt}_exons.sorted{ext}"
    gene_out = f"{pt}_genes.sorted{ext}"


    exon_out = os.path.join(DATA_DIR, exon_out)
    gene_out = os.path.join(DATA_DIR, gene_out)

    anno_file = os.path.join(DATA_DIR, anno_file)

    #  exon
    cmd_exon = f"awk 'BEGIN{{OFS=\"\\t\"}} $3 == \"exon\"' {anno_file} | sort -k1,1V -k4,4n > {exon_out}"
    #  gene
    cmd_gene = f"awk 'BEGIN{{OFS=\"\\t\"}} $3 == \"gene\"' {anno_file} | sort -k1,1V -k4,4n > {gene_out}"

    try:
        subprocess.run(cmd_exon, shell=True, check=True, executable='/bin/bash')
        subprocess.run(cmd_gene, shell=True, check=True, executable='/bin/bash')
        print(f"Successfully:\n  {exon_out}\n  {gene_out}")
    except subprocess.CalledProcessError as e:
        print(" Error:", e)

def nonModel_processing(fa_name, anno_name):
    processes = []
    cleanup_done = False
    task_cancelled = False
    
    def cleanup_all_processes(signum=None, frame=None):
        nonlocal cleanup_done, task_cancelled
        if cleanup_done:
            return
        cleanup_done = True
        task_cancelled = True
        
        print(f"Task cancelled - cleaning up {len(processes)} processes...")
        for proc in processes:
            try:
                if proc.poll() is None:
                    print(f"Terminating subprocess PID {proc.pid}")
                    try:
                        pgid = os.getpgid(proc.pid)
                        os.killpg(pgid, signal.SIGTERM)
                    except ProcessLookupError:
                        pass
                    
                    try:
                        proc.wait(timeout=3)
                    except subprocess.TimeoutExpired:
                        try:
                            os.killpg(pgid, signal.SIGKILL)
                            proc.wait(timeout=1)
                        except:
                            pass
            except Exception as e:
                print(f"Error killing subprocess: {e}")
    
    
    try:
        start = time.time()
        print("sap ghi xong fasta")

        name_fa = fa_name.split(".")[0]
        name_anno = anno_name.split(".")[0]

        nonmdFA = f"nmd_0_{name_fa}.fa"
        nonmdAN = f"nmd_0_{name_anno}.gff3"

        print(name_fa, name_anno, nonmdFA, nonmdAN)
        
        
        temp_fasta = None
        temp_anno = None
        found = False
        for ext in [".fa", ".fna"]:
            temp_fasta = os.path.join(DATA_DIR, f"{name_fa}{ext}")
            if os.path.exists(temp_fasta):
                found = True
                break
        if not found:
            raise FileNotFoundError("Not tmp found")
        
        found = False
        for ext in [".gtf", ".gff", ".gff3"]:
            temp_anno  = os.path.join(DATA_DIR, f"{name_anno}{ext}")
            if os.path.exists(temp_anno):
                found = True
                break
        if not found:
                raise FileNotFoundError("Not found tmp")


        fasta_path = os.path.join(DATA_DIR, nonmdFA)
        anno_path = os.path.join(DATA_DIR, nonmdAN)

        print(fasta_path)
        os.rename(temp_fasta, fasta_path)


        cmd = [
            "conda", "run", "-n", "agat",
            "agat_convert_sp_gxf2gxf.pl",
            "--gff", temp_anno,
            "-o", anno_path,
            "-v", "2"
        ]
        print("Start converting annotaion file")
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, preexec_fn=os.setpgrp)
        processes.append(proc)
        stdout, stderr = proc.communicate()
        
        
        print("AGAT stdout:", stdout)
        print("AGAT stderr:", stderr)
        
        if proc.returncode != 0:
            raise Exception(f"AGAT conversion failed with code {proc.returncode}: {stderr}")
        
        if not os.path.exists(anno_path):
            raise FileNotFoundError(f"AGAT output file not created: {anno_path}")
        
        if os.path.getsize(anno_path) == 0:
            raise Exception(f"AGAT created empty file: {anno_path}")
        
        print(f"AGAT conversion successful. Output size: {os.path.getsize(anno_path)} bytes")

        for temp_file in [temp_fasta, temp_anno]:
            if os.path.exists(temp_file):
                os.remove(temp_file)
                
        cmd_2bit = "../faToTwoBit" + " " + nonmdFA + " " + nonmdFA.split(".")[0] + ".2bit"
        cmd_bowtie_build = "bowtie-build" + " " + nonmdFA + " " + nonmdFA.split(".")[0] + "_index"

        proc = subprocess.Popen(cmd_2bit, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, cwd=DATA_DIR, preexec_fn=os.setpgrp)
        processes.append(proc)
        stdout, stderr = proc.communicate()
        proc.wait()
                
        if proc.returncode != 0:
            raise Exception(f"2bit conversion failed: {stderr}")

        end_2bit = time.time()
        print("Converted to .2bit", end_2bit - start)

        proc = subprocess.Popen(cmd_bowtie_build, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, cwd=DATA_DIR, preexec_fn=os.setpgrp)
        processes.append(proc)
        stdout, stderr = proc.communicate()
        proc.wait()
        
        
        if proc.returncode != 0:
            raise Exception(f"Bowtie build failed: {stderr}")

        end_bowtie = time.time()
        print("Bowtie build successfully", end_bowtie - start)
        
        
        db_path = os.path.join(DATA_DIR, f"nmd_0_{name_anno}.db")
        gff3_file = os.path.join(DATA_DIR, nonmdAN)
        
        print(f"Creating gffutils database from {gff3_file}")
        if not os.path.exists(gff3_file):
            raise FileNotFoundError(f"GFF3 file not found: {gff3_file}")
        
        file_size = os.path.getsize(gff3_file)
        print(f"GFF3 file size: {file_size} bytes")
        if file_size == 0:
            raise Exception(f"GFF3 file is empty: {gff3_file}")
        
        db = gffutils.create_db(
            gff3_file,
            dbfn=db_path,
            merge_strategy="create_unique",
            keep_order=True,
            disable_infer_transcripts=True,
            disable_infer_genes=True,
        )
        print(".db created successfully")
        end_gff = time.time()
        createMMRegion(nonmdAN)

        print("da tao gff3", end_gff - start)
    
    except Exception as e:
        print(f"Task failed: {e}")
        cleanup_all_processes()
        raise



nonModel_processing(faname, annoname)