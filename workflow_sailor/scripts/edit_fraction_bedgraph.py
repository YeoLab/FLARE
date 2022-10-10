import argparse
import os
import pandas as pd
import sys

def split_coverage_col(edits_coverage):
    splitted = edits_coverage.split(',')
    edits = int(splitted[0])
    coverage = int(splitted[1])
    edit_fraction = float(edits)/float(coverage)
    
    return edits, coverage, edit_fraction

def get_new_columns(ranked_sites_bed_path):
    ranked_sites_df = pd.read_csv(ranked_sites_bed_path, sep='\t', names=['chrom', 'start', 'end', 'score', 'edits,coverage', 'strand'])   
    if len(ranked_sites_df) > 0:
        ranked_sites_df['edits'], ranked_sites_df['coverage'], ranked_sites_df['edit_fraction'] = zip(*ranked_sites_df['edits,coverage'].apply(split_coverage_col))

        return ranked_sites_df[['chrom', 'start', 'end', 'edit_fraction']]
    else:
        return pd.DataFrame(columns=['chrom', 'start', 'end', 'edit_fraction'])

def main():
    parser = argparse.ArgumentParser(description='Generate bedgraph for site-level edit fractions')
    parser.add_argument('ranked_sites_bed', type=str)
    parser.add_argument('output_bedgraph', type=str)
    
    args = parser.parse_args()
    ranked_sites_bed_path = args.ranked_sites_bed
    output_bedgraph_path = args.output_bedgraph
    
    new_columns = get_new_columns(ranked_sites_bed_path)
    new_columns.to_csv(output_bedgraph_path, sep='\t', header=False, index=False)

if __name__ == '__main__':
    main()

