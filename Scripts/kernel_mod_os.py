#USAGE: python3 kernel_mod_os.py name general_cutoff trusted_cutoff 
# INFO: HMMer files in /tbl_files, name=*_vs_swiss/trembl.tbl, output in /output

import pandas as pd
import sys
import seaborn as sns
from matplotlib import pyplot as plt
import numpy as np
import copy



# gamma transformation
def gamma_transf(gamma):
    """Return a pair of functions computing the gamma-power and its inverse"""
    def forward_gamma(x):
        return np.power(x, gamma)
    def inverse_gamma(x):
        return np.power(x, 1/gamma)
    return forward_gamma, inverse_gamma

# plotting the graph
def plot_hit_kde(df, ax, label, stack_os=False):
    if stack_os:
        col_name = "os_" + label
        # get all organisms from table
        print('extracting...')
        df[col_name] = df.description_of_target.apply(extract_genus_from_desc)
        # keep the 10 most abundant
        print('condensing...')
        df[col_name] = condense_organisms(df, col_name, max_n=18, keep_best_n=10)
        os_order = list(df[col_name].groupby(df[col_name]).count().sort_values(ascending=False).index)
        # actual plotting part
        print('plotting...')
        sns.kdeplot(x="score", data=df, ax=ax, legend=True, common_norm=True, common_grid=True,
                    bw_adjust=0.7**(np.log10(len(df))),
                    hue=col_name, palette="tab20", hue_order=os_order, multiple="stack", lw=0.3)
        lg = ax.get_legend()
        patches = lg.get_patches()
        labels = [t.get_text() for t in lg.get_texts()]
    else:
        sns.kdeplot(x="score", data=df, ax=ax, label=label, legend=True, common_norm=True,
                    common_grid=True, lw=0.9, bw_adjust=0.6**(np.log10(len(df))))
        patches, labels = copy.copy(ax.get_legend_handles_labels())
    return patches, labels

def extract_genus_from_desc(s):
    if isinstance(s, str) and "OS=" in s:
        i = 3 + s.find("OS=")
        return s[i:].split()[0]
    return "other"

def condense_organisms(df, os_label, max_n=15, keep_best_n=10):
    max_score = df.score.max()
    os_col = df[os_label]
    i = 0
    while os_col.nunique() > max_n:
        i += 1
        cnt = df.groupby(os_col)["score"].apply(lambda x: ((x / max_score) ** 4).sum())
        for s in cnt.nsmallest(len(cnt) - keep_best_n).index:
            if s.startswith("uncultured"):
                new = "uncultured"
            elif s == "other":
                continue
            else:
                new = "other"
            os_col = os_col.replace(s, new)
    df[os_label] = os_col.apply(italic_os)
    return df[os_label]

# put all names in italics
def italic_os(s):
    if not s.startswith("uncultured") and s != "other":
        s = f"$\\it{{{s}}}$"
    return s

# Read command line arguments
name = sys.argv[1]
general_cutoff = float(sys.argv[2])
trusted_cutoff = float(sys.argv[3])

# Read data from trembl file
data = []
with open(f"tbl_files/{name}_vs_trembl.tbl") as f:
    for line in f:
        if not line.startswith('#'):
            row = line.strip().split()
            row = row[:18] + [" ".join(row[18:])]  # Consolidate remaining columns into description
            data.append(row)

# Convert list to DataFrame
columns = [
    'target_name', 'target_accession', 'query_name', 'query_accession', 'evalue',
    'score', 'bias', 'dom_evalue', 'dom_score', 'dom_bias', 'exp', 'reg', 'clu',
    'ov', 'env', 'dom', 'rep', 'inc', 'description_of_target'
]
df = pd.DataFrame(data, columns=columns)
df['score'] = pd.to_numeric(df['score'])
max_score = df['score'].max()

# read data from swiss file
with open(f'tbl_files/{name}_vs_swiss.tbl', 'r') as i:
    lines2 = i.readlines()
del lines2[-10:]
del lines2[:3]

scores_swiss = []
for x in lines2:
    scores_swiss.append(x.split()[5])
scores_swiss = pd.to_numeric(scores_swiss)

# Plotting
sns.set_style('whitegrid')
plt.figure(figsize=(10, 6))
ax = plt.gca()
patches, labels = plot_hit_kde(df, ax, label=name, stack_os=True)

# adding the other scores
swiss = sns.kdeplot(scores_swiss, ax=ax, label='Swiss', lw=0.9, bw_adjust=0.6**(np.log10(len(scores_swiss))))
patches.append(swiss.get_lines()[0])
labels.append('Swiss')
print('finalizing...')

# adding cutoffs
ax.set_xlim(max_score, 0)
gen = ax.axvline(x=general_cutoff, color='r', linestyle='--', label='general cutoff')
tru = ax.axvline(x=trusted_cutoff, color='b', linestyle='--', label='trusted cutoff')
patches.extend([gen, tru])
labels.extend(['general cutoff', 'trusted cutoff'])

# final modifications
forward_gamma, inverse_gamma = gamma_transf(0.1)
plt.yscale('function', functions=(forward_gamma, inverse_gamma))
plt.title(f'Kernel Density Plot of {name} in UniProtKB')
plt.xlabel('Score')
plt.ylabel('Density')

# Adjusting the legend location
plt.legend(patches, labels, bbox_to_anchor=(1.05, 1), loc='upper left', borderaxespad=0.)
plt.savefig(f"output/{name}_kdeplot.png", format='png', bbox_inches='tight')
print(f"saved to output/{name}_kdeplot.png")
