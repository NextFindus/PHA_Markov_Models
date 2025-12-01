#Usage: python3 hmmsearch_metafile.py hmm_directory metagenome_directory output_directory

from pathlib import Path
import subprocess
import sys
import os

# get hmms from specified directory
directory = Path(sys.argv[1])
list_hmms = [f for f in directory.iterdir() if f.is_file()]
print(f"Access to:\n{list_hmms}")

# get genomes from specified directory
meta = Path(sys.argv[2])
print(f"Access to:\n{meta}")

# determine outdirectory
outdir = Path(sys.argv[3])
print(f"Output to:\n{outdir}")

def hmmsearch(hmm, metafile, outdir):
    # get filenames
    hmm_path = hmm.resolve()
    meta_path = meta.resolve()
    # pretty outputfilenames
    filename = f"{hmm.stem}_{meta.stem}"
    tblfilename = f"{hmm.stem}_{meta.stem}.tbl"
    outfile = f"{outdir}/{filename}"
    tbloutfile = f"{outdir}/{tblfilename}"
    # use hmmsearch in terminal
    cmd = f"hmmsearch --cpu 20 -o {outfile} --tblout {tbloutfile} {hmm_path} {meta_path}"
    print(cmd)
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    return result

# run all hmms over metagenomefile
for hmm in list_hmms:
    print(hmm)
    result = hmmsearch(hmm, meta, outdir)
    print(result)
