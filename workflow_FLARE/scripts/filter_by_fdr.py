import argparse
import os
import pandas as pd
import sys

parser = argparse.ArgumentParser(description='Split up regions file')
parser.add_argument('peaks_with_fdr_input_file', type=str)
parser.add_argument('fdr_threshold', type=float)
parser.add_argument('peaks_filtered_by_fdr_output_file', type=str)

args = parser.parse_args()
peaks_with_fdr_input_file = args.peaks_with_fdr_input_file
fdr_threshold = args.fdr_threshold
peaks_filtered_by_fdr_output_file = args.peaks_filtered_by_fdr_output_file


df = pd.read_csv(peaks_with_fdr_input_file, sep='\t')

filtered_df = df[df.p_adj < fdr_threshold]
filtered_df[['chrom', 'start', 'end', 'mean_depth', 'p_adj', 'strand', 'overall_region_id']].to_csv(peaks_filtered_by_fdr_output_file, sep='\t', header=False, index=False)
    
