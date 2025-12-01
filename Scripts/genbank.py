#USAGE: python3 genbank.py 

import argparse
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

def extract_regions_by_locus_tag(input_file, output_file, locus_tags):
    extracted_records = []
    # Parse the input GenBank file
    records = SeqIO.parse(input_file, "genbank")

    for record in records:
        all_features = []
        loci_positions = []
        for feature in record.features:
            if 'locus_tag' in feature.qualifiers and feature.qualifiers['locus_tag'][0] in locus_tags:
                loci_positions.append((feature.location.start, feature.location.end))
        if loci_positions:
            # Determine the minimum start and maximum end positions from the locus tags
            start_pos = min(start for start, end in loci_positions)
            end_pos = max(end for start, end in loci_positions)
            # Extracting all features within the start and end positions
            for feature in record.features:
                if feature.location.start >= start_pos and feature.location.end <= end_pos:
                    all_features.append(feature)
            # Extracting the sequence within the start and end positions
            sub_seq = record.seq[start_pos:end_pos]
            # Adjusting the feature locations
            adjusted_features = []
            for feature in all_features:
                new_start = feature.location.start - start_pos
                new_end = feature.location.end - start_pos
                new_location = FeatureLocation(new_start, new_end, strand=feature.location.strand)
                new_feature = SeqFeature(new_location, type=feature.type, qualifiers=feature.qualifiers)
                adjusted_features.append(new_feature)
            # Ensure annotations are copied
            annotations = record.annotations.copy()
            annotations["molecule_type"] = annotations.get("molecule_type", "DNA")
            sub_record = SeqRecord(
                sub_seq, id=record.id, name=record.name, description=record.description,
                features=adjusted_features, annotations=annotations
            )
            extracted_records.append(sub_record)
    if extracted_records:
        # Extracted regions to output file
        SeqIO.write(extracted_records, output_file, "genbank")
    else:
        print("No matching regions found.")

def main():
    parser = argparse.ArgumentParser(description="Extract specified regions from a GenBank file by locus tag.")
    parser.add_argument("input_file", help="Path to the input GenBank file.")
    parser.add_argument("output_file", help="Path to the output GenBank file.")
    parser.add_argument("locus_tags", nargs=2, help="Two locus tags, region to extract")

    args = parser.parse_args()

    extract_regions_by_locus_tag(args.input_file, args.output_file, args.locus_tags)

if __name__ == "__main__":
    main()
