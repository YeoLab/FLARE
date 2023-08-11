import argparse
import os
import pandas as pd
import sys

# chr1:879944-879974(-)|NOC2L:exon_0

parser = argparse.ArgumentParser(description='Add label index to merged peaks')
parser.add_argument('peaks_filtered_by_fdr_input_file', type=str)
parser.add_argument('peaks_filtered_by_fdr_reindexed_output_file', type=str)

args = parser.parse_args()
peaks_filtered_by_fdr_input_file = args.peaks_filtered_by_fdr_input_file
peaks_filtered_by_fdr_reindexed_output_file = args.peaks_filtered_by_fdr_reindexed_output_file

df = pd.read_csv(peaks_filtered_by_fdr_input_file, sep='\t', 
            names=['chrom', 'start', 'end', 'mean_depth', 'region_adj_p', 'strand', 'overall_region_id'])

            
def reindex(r):
    return '{}:{}-{}({})|{}'.format(r.chrom, r.start, r.end, r.strand, r.overall_region_id)

        
df['region_id'] = df.apply(reindex, axis=1)
df['subregion'] = [i.split('|')[1].split(',')[0].split(':')[1] for i in df['region_id']]
df['region'] = [i.split('|')[1].split(',')[0].split(':')[1].split('_')[0] for i in df['region_id']]
df['gene'] = [i.split('|')[1].split(',')[0].split(':')[0] for i in df['region_id']]

df[['region_id', 'chrom', 'start', 'end', 'strand', 'subregion', 'region', 'gene']].to_csv(peaks_filtered_by_fdr_reindexed_output_file, sep='\t', index=False, header=True)