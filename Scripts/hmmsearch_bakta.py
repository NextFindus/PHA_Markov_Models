#USAGE: python3 hmmsearch_bakta.py hmm_dir bakta_dir output_dir
#bakta_dir is unedited bakta output: */*.faa

from pathlib import Path
import subprocess
import sys
import os

# get hmms from specified directory
directory = Path(sys.argv[1])
list_hmms = [f for f in directory.iterdir() if f.is_file()]
print(f"Access to:\n{list_hmms}")

# get genomes from specified directory, only */*.faa files (due to output of bakta)
metadirectory = Path(sys.argv[2])
print(metadirectory)
list_meta=[]
list_meta = [f for f in metadirectory.glob('**/*') if f.is_file() and f.name.endswith('.faa')]
print(f"Access to:\n{list_meta}")

# determine outdirectory
outdir = Path(sys.argv[3])
print(f"Output to:\n{outdir}")

def hmmsearch(hmm, meta, outdir):
    # getting filenames
    hmm_path = hmm.resolve()
    meta_path = meta.resolve()
    # pretty output filenames
    filename = f"{hmm.stem}_{meta.stem}"
    tblfilename = f"{hmm.stem}_{meta.stem}.tbl"
    outfile = f"{outdir}/{filename}"
    tbloutfile = f"{outdir}/{tblfilename}"
    # use hmmsearch in terminal
    cmd = f"hmmsearch --cpu 20 -o {outfile} --tblout {tbloutfile} {hmm_path} {meta_path}"
    print(cmd)
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    return result

# iterate through hmms and genome files, scanning hmms against all of them
for meta in list_meta:
    for hmm in list_hmms:
        print(hmm)
        result = hmmsearch(hmm, meta, outdir)
        print(result)
