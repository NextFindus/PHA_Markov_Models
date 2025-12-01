#Usage: python3 repeat_finder.py sequences.fasta
#INFO: Uses fasta files of the phaR+100 bp upstream. ATG at 101 is marked as bold and all inverted repeats between 5 and 20 bp within -30 bp and their respective inverted sequences are randomly colored. If a repeat is found within -7 bp, it is colored red.

from Bio import SeqIO
from Bio.Seq import Seq
import random
import sys
from collections import defaultdict

#random color generator 
def generate_color():
    return f"#{random.randint(0, 0xFFFFFF):06x}"


def find_inverted_repeats(sequences, min_length=5, max_length=20, min_gap=0, max_gap=100):
    inverted_repeats = defaultdict(list)
    
    for record in sequences:
        seq_str = str(record.seq).upper()
        seq_len = len(seq_str)
        #creating list of kmers 
        kmer_dict = defaultdict(list)
        for k in range(min_length, max_length + 1):
            for i in range(seq_len - k + 1):
                kmer = seq_str[i:i+k]
                kmer_dict[kmer].append(i)
        #compare kmers to sequence and save metadata about possible repeats
        for k in range(min_length, max_length + 1):
            for i in range(seq_len - k + 1):
                left = seq_str[i:i+k]
                rev_comp = str(Seq(left).reverse_complement())
                
                for j in kmer_dict[rev_comp]:
                    if i + k + min_gap <= j <= i + k + max_gap:
                        inverted_repeats[record.id].append({
                            'left_start': i,
                            'left_end': i + k - 1,
                            'right_start': j,
                            'right_end': j + k - 1,
                            'repeat': left,
                            'color': generate_color()
                        })
    
    return inverted_repeats

#go through all repeats of every sequence to find common ones
def find_common_inverted_repeats(inverted_repeats, num_sequences):
    repeat_counts = defaultdict(set)
    for seq_id, repeats in inverted_repeats.items():
        for repeat in repeats:
            repeat_counts[repeat['repeat']].add(seq_id)
    
    return {repeat for repeat, sequences in repeat_counts.items() if len(sequences) == num_sequences}

#confirm ATG at postion 101
def has_atg_at_101(sequence):
    return len(sequence) >= 103 and sequence[100:103] == 'ATG'

#only keep inverted repeats near ATG
def filter_repeats_near_atg(inverted_repeats, atg_position=101, upstream_distance=30):
    filtered_repeats = defaultdict(list)
    for seq_id, repeats in inverted_repeats.items():
        for repeat in repeats:
            if (atg_position - upstream_distance <= repeat['left_end'] < atg_position or
                atg_position - upstream_distance <= repeat['right_end'] < atg_position):
                filtered_repeats[seq_id].append(repeat)
    return filtered_repeats

#define right color for every base
def color_sequence(seq_str, repeats, common_repeats):
    repeat_regions = []
    for repeat in repeats:
        left_start, left_end = repeat['left_start'], repeat['left_end']
        right_start, right_end = repeat['right_start'], repeat['right_end']
        color = repeat['color']
        is_common = repeat['repeat'] in common_repeats
        
        # Calculate distance to ATG for both left and right regions
        distance_to_atg = min(
            min(abs(100 - i) for i in range(left_start, left_end + 1)),
            min(abs(100 - i) for i in range(right_start, right_end + 1))
        )
        if distance_to_atg <= 7:
            color = '#ff0000'
        repeat_regions.append((left_start, left_end, right_start, right_end, color, is_common, distance_to_atg))
    #Sort by distance to ATG, nearer=priority during coloring
    repeat_regions.sort(key=lambda x: x[6])  

    colored_positions = [None] * len(seq_str)

    for left_start, left_end, right_start, right_end, color, is_common, _ in repeat_regions:
        for i in range(left_start, left_end + 1):
            if colored_positions[i] is None:
                colored_positions[i] = (color, is_common)
        for i in range(right_start, right_end + 1):
            if colored_positions[i] is None:
                colored_positions[i] = (color, is_common)

    return colored_positions

#make svg with fasta accession and sequence, mark ATG and inverted repeats
def generate_svg_output(sequences, inverted_repeats, common_repeats, output_file):
    char_width = 10
    char_height = 20
    line_spacing = 10
    max_width = max(len(str(record.seq)) for record in sequences)
    svg_width = max(max_width * char_width, 800)
    svg_height = (len(sequences) * (char_height + line_spacing)) + 50

    svg = [f'<svg width="{svg_width}" height="{svg_height}" xmlns="http://www.w3.org/2000/svg">']
    svg.append('<style>')
    svg.append('text { font-family: monospace; font-size: 16px; }')
    svg.append('.common { font-weight: bold; text-decoration: underline; }')
    svg.append('.big { font-size: 24px; font-weight: bold; }')  # Style for big letters
    svg.append('</style>')

    y_offset = 20
    for record in sequences:
        svg.append(f'<text x="10" y="{y_offset}">{record.id}:</text>')
        y_offset += char_height + 5

        seq_str = str(record.seq).upper()
        has_atg = has_atg_at_101(seq_str)

        colored_positions = color_sequence(seq_str, inverted_repeats[record.id], common_repeats)

        for i, base in enumerate(seq_str):
            x = i * char_width + 10
            text_color = 'black'
            bg_color = colored_positions[i][0] if colored_positions[i] else None

            if has_atg and 100 <= i < 103:
                bg_color = 'yellow'

            if bg_color:
                    svg.append(f'<rect x="{x}" y="{y_offset - char_height}" width="{char_width}" height="{char_height}" fill="{bg_color}" />')

            if colored_positions[i] and colored_positions[i][1]:  # is_common
                text_color = 'red'

                svg.append(f'<text x="{x}" y="{y_offset}" fill="{text_color}">{base}</text>')

        y_offset += line_spacing

    svg.append('</svg>')

    with open(output_file, 'w') as f:
        f.write('\n'.join(svg))

#same as html
def generate_html_output(sequences, inverted_repeats, common_repeats, output_file):
    with open(output_file, 'w') as f:
        f.write("<html><head><style>")
        f.write("body { font-family: monospace; white-space: pre; }")
        f.write(".common { font-weight: bold; text-decoration: underline; }")
        f.write(".start-codon { background-color: yellow; font-weight: bold; }")
        f.write(".big { font-size: 24px; font-weight: bold; }")  # Style for big letters
        f.write("</style></head><body>")
        
        for record in sequences:
            f.write(f"<h2>{record.id}</h2>")
            seq_str = str(record.seq).upper()
            has_atg = has_atg_at_101(seq_str)
            
            colored_positions = color_sequence(seq_str, inverted_repeats[record.id], common_repeats)
            
            colored_seq = []
            for i, base in enumerate(seq_str):
                if has_atg and 100 <= i < 103:
                    colored_seq.append(f"<span class='start-codon'>{base}</span>")
                elif colored_positions[i]:
                    color, is_common = colored_positions[i]
                    class_name = " common" if is_common else ""
                    colored_seq.append(f"<span class='{class_name}' style='background-color: {color};'>{base}</span>")
                else:
                    colored_seq.append(base)
            
            f.write("".join(colored_seq))
            f.write("<br><br>")
        
        f.write("</body></html>")

# Main
if __name__ == "__main__":
    #get sequences from file
    input_file = sys.argv[1]
    sequences = list(SeqIO.parse(input_file, "fasta"))
    #find all inverted repeats
    all_inverted_repeats = find_inverted_repeats(sequences)
    #filter them to extract ones near ATG
    filtered_repeats = filter_repeats_near_atg(all_inverted_repeats)
    #identify common repeats in all sequences
    common_repeats = find_common_inverted_repeats(filtered_repeats, len(sequences))
    #make html and svg of them
    generate_html_output(sequences, filtered_repeats, common_repeats, f"{input_file}.html")
    generate_svg_output(sequences, filtered_repeats, common_repeats, f"{input_file}.svg")

    total_repeats = sum(len(repeats) for repeats in filtered_repeats.values())
    print(f"\nTotal number of sequences processed: {len(sequences)}")
    print(f"Total number of inverted repeats found within 30 bases upstream of ATG: {total_repeats}")
    print(f"Number of common inverted repeats: {len(common_repeats)}")

