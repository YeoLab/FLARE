import argparse
import os
import pandas as pd
import sys    
import numpy as np
from scipy.stats import nbinom

parser = argparse.ArgumentParser(description='Add statistical score to peaks')
parser.add_argument('peaks_with_edit_fraction', type=str)
parser.add_argument('peaks_with_scores', type=str)
parser.add_argument('all_windows_with_fdr')

args = parser.parse_args()
peaks_with_edit_fraction = args.peaks_with_edit_fraction
peaks_with_scores = args.peaks_with_scores
all_windows_with_fdr = args.all_windows_with_fdr 

all_windows = pd.read_csv(all_windows_with_fdr, sep='\t')
all_windows.replace([np.inf, -np.inf], np.nan, inplace=True)
all_windows = all_windows.fillna(0)
background_rate = all_windows.subregion_fraction.mean()

print("Input file specified: {}".format(peaks_with_edit_fraction))
print("Input file specified: {}".format(peaks_with_scores))
print("Background rate for negative binomial peak scoring calculated as: {}".format(background_rate))


def get_confidence(conversions, coverage, frac=background_rate):
    # # Determine nbinom CDF
    # (nbinom.cdf(nbr_failures,nbr_successes,prob_success))
    confidence = 1 - nbinom.cdf(coverage-conversions, conversions, frac)
    return confidence


def conf_for_peak(r):
    coverage = r.target_bases
    conversions = r.edited_bases
    
    return get_confidence(conversions, coverage)


df = pd.read_csv(peaks_with_edit_fraction, sep='\t')
df['score'] = df.apply(conf_for_peak, axis=1)

if 'overall_region_id' in df.columns:
    df = df.drop('overall_region_id', axis=1)
    
df.to_csv(peaks_with_scores, sep='\t', index=False, header=True)