import argparse
import os
import pandas as pd
import sys    
from scipy.stats import poisson
import numpy as np

parser = argparse.ArgumentParser(description='Get window significance using control background')
parser.add_argument('edit_fraction_target', type=str)
parser.add_argument('edit_fraction_background', type=str)

args = parser.parse_args()

edit_fraction = {}
edit_fractions['target'] = args.edit_fraction_target
edit_fractions['background'] = args.edit_fraction_background

print("Target file specified: {}".format(edit_fraction_target))
print("Background file specified: {}".format(edit_fraction_background))




def apply_input_normalization_to_site(r):   
    fraction_in = r.fraction_in
    subregion_fraction_in = r.subregion_fraction_in
    region_fraction_in = r.region_fraction_in
    gene_fraction_in = r.gene_fraction_in
        
    site_coverage = int(r.coverage_ip)
    site_edits = int(r.conversions_ip)   
    
    region_fraction_ip = r.region_fraction_ip
    subregion_fraction_ip = r.subregion_fraction_ip
    gene_fraction_ip = r.gene_fraction_ip
    
    
    if gene_fraction_in > 0 and gene_fraction_in <= 1 and fraction_in > 0 and fraction_in <= 1 and gene_fraction_ip > 0 and gene_fraction_ip <= 1:
        editing_delta = gene_fraction_ip / gene_fraction_in
        window_background = editing_delta * np.mean([fraction_in, subregion_fraction_in, region_fraction_in])
    else:
        window_background = np.mean([subregion_fraction_ip, region_fraction_ip])
    
    baseline_edits = window_background * site_coverage
        
    # Calculate a Poisson statistic based on extrapolated context editing as lambda
    return 1 - poisson.cdf(site_edits, baseline_edits, loc=0)




target_and_background = {}

for label, filename in edit_fractions.items():
    edits = pd.read_csv(filename, sep='\t', index_col=0)
    for prefix in ['', 'subregion_', 'region_', 'gene_']:
        if prefix != '':
            edits['{}fraction'.format(prefix)] = (edits['{}conversions'.format(prefix)]-edits['conversions'])/(edits['{}coverage'.format(prefix)] - edits['coverage'])
        else:
            edits['{}fraction'.format(prefix)] = (edits['{}conversions'.format(prefix)])/(edits['{}coverage'.format(prefix)])
    
    edits = edits.fillna(0)
    target_and_background[label] = edits
    

ip_edits_with_input_ratios = target_and_background.get('target').join(target_and_background.get('background'), lsuffix='_ip', rsuffix='_in')
    

print("Applying pvalue to windows based on input distribution...")
ip_edits_with_input_ratios['pvalue'] = edits.apply(apply_input_normalization_to_site, axis=1)



        
        
df.to_csv(peaks_with_scores, sep='\t', index=False, header=True)