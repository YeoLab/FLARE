from multiprocessing import Pool
from scipy.stats import poisson
import pandas as pd
import numpy as np
import json
import sys
from collections import defaultdict
from statsmodels.stats.multitest import fdrcorrection
from argparse import ArgumentParser
import os
import math


def write(text):
    sys.stdout.write('{}\n'.format(text))

    
def get_window_from_depth(depth, depth_windows, mean_fraction_by_depth_windows):
    index = math.floor(depth/10)
    if index > len(depth_windows) - 2:
        index = len(depth_windows) - 2
    return mean_fraction_by_depth_windows.get(index)

def zero_truncated_poisson_cdf(site_edits, baseline_num_edits):
    numerator = poisson.cdf(site_edits, baseline_num_edits, loc=0) - poisson.cdf(0, baseline_num_edits, loc=0)
    denominator = 1 - poisson.cdf(0, baseline_num_edits, loc=0)
        
    return numerator/denominator


def apply_exon_edit_c_pvalue_to_site(r):   
    site_coverage = int(r.target_bases)
    window_background = float(r.background_rate)
    site_edits = int(r.edited_bases)
            
    # Calculate a Poisson statistic
    baseline_num_edits = (window_background * site_coverage)
    return 1 - zero_truncated_poisson_cdf(site_edits, baseline_num_edits)

    
def poisson_filter(edit_c_df):
    region_sites = edit_c_df.copy()

    prefix = 'subregion_'
    region_sites['{}fraction'.format(prefix)] = (region_sites['{}conversions'.format(prefix)]-region_sites['edited_bases'])/(region_sites['{}coverage'.format(prefix)] - region_sites['target_bases'])
    region_sites.replace([np.inf, -np.inf], np.nan, inplace=True)
    region_sites = region_sites.fillna(0)
    
    depth_windows = [0, 10, 20, 30, 40, 50, 1000000]
    mean_fraction_by_depth_windows = {}
    for i, depth_window in enumerate(depth_windows):
        if i < len(depth_windows) -1:
            print(i,  depth_windows[i],  depth_windows[i+1])
            mean_fraction_by_depth_windows[i] = region_sites[(region_sites.mean_depth < depth_windows[i+1])].subregion_fraction.mean()
    print('Mean fraction by depth:\n', '\t', mean_fraction_by_depth_windows)
    
    region_sites['background_rate'] = region_sites.mean_depth.apply(get_window_from_depth, args=(depth_windows, mean_fraction_by_depth_windows,))
    
    # If there were any edits at all in this region, then
    # calculate the p-value for this site
    region_sites['poisson_p'] = region_sites.apply(
        apply_exon_edit_c_pvalue_to_site, axis=1)

    region_sites['poisson_p'] = region_sites.poisson_p.fillna(1)
    
    region_sites['p_adj'] = fdrcorrection(list(region_sites.poisson_p), alpha=0.1)[1]

    return region_sites


def load_json_info(input_json_filepath):
    # Load parameters/filepaths
    with open(input_json) as f:
        data = json.load(f)

    output_folder = data.get('output_folder')
    label = data.get('label')

    write("Label: {}".format(label))
    return output_folder,label

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# SCRIPT STARTS RUNNING HERE ... 
parser = ArgumentParser(description='Multipoisson peak calling')
parser.add_argument('input_json', type=str)
parser.add_argument('--edit_c_regions', type=str, default=None)

args = parser.parse_args()
input_json = args.input_json
edit_c_regions = args.edit_c_regions


write('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\ninput is {}'.format(input_json))

output_folder,label= load_json_info(input_json)

    

output_file = '{}/peak_calling/{}/{}_all_regions_with_fdr.tsv'.format(output_folder, label, label)
write("Output folder is {}\n".format(output_folder))
write("Outputfile is {}\n".format(output_file))


write("Window Edit Fractions File is {}\n".format(edit_c_regions))
edit_c_regions = pd.read_csv(edit_c_regions, sep='\t', index_col=0)
#print(edit_c_regions.columns)
    
# Add poisson statistic and FDR
write("\n\tPoisson statistic calculation...")
def output_blank():
    blank_fdr_output = pd.DataFrame(columns=[['chrom', 'start', 'end', 'overall_region_id', 'coverage', 'strand', 'conversions',
           'edit_c', 'chrom_region', 'start_region', 'end_region',
           'edit_c_region', 'strand_region', 'coverage_region',
           'conversions_region', 'poisson_p', 'p_adj']])
    blank_fdr_output.to_csv(output_file, sep='\t', header=True, index=False)

if len(edit_c_regions) == 0:
    write('\nNo values, outputting empty file\n')
    output_blank()
    sys.exit(0)
    
with_poisson = poisson_filter(edit_c_regions)
if len(with_poisson) > 0:
    # Benjamini-Hochberg FDR correction
    with_poisson['region_id'] = with_poisson.index
    overall_region_id = [':'.join(i.split(':')[:2]) for i in with_poisson.region_id]
    with_poisson['overall_region_id'] = overall_region_id
    with_poisson.to_csv(output_file, sep='\t', header=True, index=False)


else:
    write('\nNo values, outputting empty file\n')
    output_blank()
