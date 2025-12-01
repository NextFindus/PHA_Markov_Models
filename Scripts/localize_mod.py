#USAGE: python3 localize_mod.py organism_name
# INFO: score files as .tbl in /proteins, output goes in /locations

import sys
import pandas as pd
from pathlib import Path
import glob
import os



# get name of the organism to identify files
name = sys.argv[1]
print(name)
score_files_pattern = f"./proteins/*{name}.tbl"
output_file = f'./locations/{name}_locations.csv'

# Read all score files and concatenate them into a single DataFrame
score_files = glob.glob(score_files_pattern)
score_df_list = []
for files in score_files:
    with open(files, 'r') as f:
        lines = f.readlines()
    cleaned_lines = []
    for line in lines:
        if not line.startswith('#'):
            cleaned_line = line.replace('#', '').strip()
            cleaned_lines.append(cleaned_line)
    df = pd.DataFrame([x.split() for x in cleaned_lines], columns=[
        "contig", "accession", "query_name", "query_accession",
        "E-value_full", "score_full", "bias_full",
        "E-value_best", "score_best", "bias_best",
        "exp", "reg", "clu", "ov", "env", "dom", "rep", "inc", 'start', 'stop',
        'strain', 'description'
    ])
    # clean up column names
    df['HMM'] = os.path.basename(files).split("_")[1]
    df['contig'] = pd.to_numeric(df['contig'].apply(lambda x: x.split('_')[1]))
    df['ID'] = df['description'].apply(lambda x: str(x).split(';')[0])
    # sort out low values
    df = df[pd.to_numeric(df['score_full']) >= 90]
    score_df_list.append(df)

score_combined_df = pd.concat(score_df_list, ignore_index=True)

# Rename Score column
final_df = pd.DataFrame(score_combined_df)
final_df = final_df.rename(columns={'score_full': 'score'})

# Sort by 'Sequence Id' and 'Start'
final_df['start'] = pd.to_numeric(final_df['start'])
final_df_sorted = final_df.sort_values(by=['contig', 'start'])

# calculate distance between genes in genome
final_df_sorted['jump'] = pd.to_numeric(final_df_sorted['start'].shift(-1)) - pd.to_numeric(final_df_sorted['stop'])

# rearrange columns
final_df_sorted = final_df_sorted[['contig', 'ID', 'HMM', 'score', 'start', 'stop', 'jump']]

# Now we have a table with all hits!
final_df_sorted.to_csv(f"{output_file}", index=False)
print(f"Data successfully written to {output_file}")
