#!/usr/bin/env python
from pyfaidx import Fasta
import argparse
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
import concurrent.futures
from argparse import ArgumentParser
from argparse import RawDescriptionHelpFormatter
import pyBigWig
import numpy as np
import pandas as pd
import os
import glob
import pybedtools
import gffutils
from Bio import SeqIO
from tqdm import trange
from collections import defaultdict, OrderedDict
pd.options.mode.chained_assignment = None  # default='warn'
import sys
import pandas as pd
import re
import json
from multiprocessing import Pool
import time
import datetime
import pysam

class ReadDensity:
    """
    ReadDensity class
    Attributes:
        self.pos(positive *.bw file)
        self.neg(negative *.bw file)
    """

    def __init__(self, pos, neg, name=None):
        try:
            self.pos = pyBigWig.open(pos)
            self.neg = pyBigWig.open(neg)
            self.name = name if name is not None else pos.replace(
                'fwd', '*'
            ).replace(
                'rev', '*'
            )

        except Exception as e:
            print("couldn't open the bigwig files!")
            print(e)

    def values(self, chrom, start, end, strand, reverse=False):
        """
        Parameters
        ----------
        chrom : basestring
            (eg. chr1)
        start : int
            0-based start (first position in chromosome is 0)
        end : int
            1-based end (last position is not included)
        strand : str
            either '+' or '-'
        Returns
        -------
        densites : list
            values corresponding to density over specified positions.
        """
        chrom = str(chrom)
        if reverse:
            if strand == '+':
                strand = '-'
            elif strand == '-':
                strand = '+'
        
        try:
            if strand == "+":
                return list(pd.Series(self.pos.values(chrom, start, end)).fillna(0))
            elif strand == "-":
                return list(pd.Series(self.neg.values(chrom, start, end)).fillna(0))
            else:
                print("Strand neither + or -")
                return 1
        except RuntimeError as e:
            print(e)
            # usually occurs when no chromosome exists in the bigwig file
            return list(pd.Series([np.NaN] * abs(start - end)).fillna(0))
        
def find(s, ch):
    return [i for i, ltr in enumerate(s) if ltr == ch]

def get_c_positions_and_coverage_in_window(chrom, start, end, strand, rdd, edit_type, reverse=False):
    """
    Given region specified by chrom start end and strand, 
    return the positions (0-based) of target bases (or complement of target bases on neg strand).
    """ 
    d = {}
    
    chrom = str(chrom)

    FA = Fasta(fasta, rebuild=False)
        
    sequence = ''
    try:
        sequence = FA['chr{}'.format(chrom)][start:end].seq
    except Exception as e:
        sequence = FA['{}'.format(chrom)][start:end].seq
    
    to_find = edit_type[0]
    to_find_rc = COMPLEMENT.get(to_find)

    if strand == '+':
        looking_for = to_find
    elif strand == '-':
        looking_for = to_find_rc
    else:
        print("Strand error: {}".format(strand))
        sys.exit(1)

    # TO REMOVE
    # print('REMOVE: looking for {}, sequence: {}'.format(looking_for, sequence))
    #print(sequence)
    if len(sequence) > 0:
        matches = re.finditer(looking_for, sequence.upper())
    else:
        print('Failed on {} {} {} {}'.format(chrom, start, end, strand))
        return d

    # print("REMOVE matches: {}".format(matches))
    
    relpos = [match.start() for match in matches]
    abspos = ["{}:{}".format(chrom, start + p) for p in relpos]
    coverage = rdd.values(chrom=chrom, start=start, end=end, strand=strand, reverse=reverse)
    coverage = [np.abs(c) for c in coverage]  

    # TO REMOVE
    #print("REMOVE: coverage: {}".format(coverage))
    c_coverage = [coverage[p] for p in relpos]
    for p, c in zip(abspos, c_coverage):
        d[p] = c
    return d

def sum_all_c_coverage_in_window(chrom, start, stop, strand, rdd, edit_type, reverse=False):
    all_coverage = 0
    c_positions = get_c_positions_and_coverage_in_window(chrom, start, stop, strand, rdd, edit_type, reverse=reverse)

    # REMOVE
    #print("REMOVE c_positions: {}".format(c_positions))
    for pos, cov in c_positions.items():
        all_coverage += cov
    return int(all_coverage)


COMPLEMENT = {
    'C': 'G',
    'A': 'T',
    'G': 'C',
    'T': 'A'
}

def get_origin_and_destination_bases(edit_type, strand='+'):
    edited_origin = edit_type[0]
    edited_destination = edit_type[1]

    if strand == '-':
        edited_origin = COMPLEMENT.get(edited_origin)
        edited_destination = COMPLEMENT.get(edited_destination)
    elif strand == '+':
        edited_origin = edited_origin
        edited_destination = edited_destination

    return edited_origin, edited_destination
    

def get_num_edited_reads_and_coverage(chrom, start, end, edited_origin, edited_destination, bam_path, reference_path):
    # potentially start should be +1?
    
    chrom = str(chrom)
    
    reference = pybedtools.BedTool.seq([chrom,start-1,end-1], reference_path)
    pybedtools.cleanup()

    substrate_bases = [i.start()+start for i in re.finditer(edited_origin, reference.upper())]

    set_of_edit_sites = set(stamp_sites[stamp_sites.chrom_stamp.astype(str) == chrom].start_stamp)
    target_bases = set(substrate_bases).intersection(set_of_edit_sites)

    # REMOVE
    #print("REMOVE start: {}".format(start))
    #print('REMOVE edited origin', edited_origin)
    #print('REMOVE reference:', reference.upper())
    #print('REMOVE substrate bases:', substrate_bases)
    #print('REMOVE set of edit sites:', [s for s in set_of_edit_sites if s >= start and s <= end])
    #print('REMOVE target bases:', target_bases)
    
    read_edits = defaultdict(lambda: defaultdict(lambda:[]))
    #print(target_bases, edited_origin, edited_destination)
    
    #print("bam path is ", bam_path)
    samfile = pysam.AlignmentFile(bam_path, "rb" )

    total_target_bases = 0
    total_edited_bases = 0
    coverages = []

    for pileupcolumn in samfile.pileup(chrom, start-1, end-1, stepper='nofilter', min_base_quality=0):
        
        if pileupcolumn.pos >= start-1:
            coverage = pileupcolumn.n
            coverages.append(coverage)

            if pileupcolumn.pos > end-1:
                break 

            # Subtraction takes care of 1-indexing
            adjusted_pileup_pos = pileupcolumn.pos+1
            if adjusted_pileup_pos in target_bases:  
                total_target_bases += coverage
                for pileupread in pileupcolumn.pileups:
                    if not pileupread.is_del and not pileupread.is_refskip:
                        # query position is None if is_del or is_refskip is set.
                        read_name = pileupread.alignment.query_name,
                        read_base = pileupread.alignment.query_sequence[pileupread.query_position].upper()
                        
                        ref_start = pileupread.alignment.reference_start
                        ref_end = pileupread.alignment.reference_end
                        
                        #print("REMOVE", "position: {}, read base: {}".format(adjusted_pileup_pos, read_base))
                        #print("REMOVE", "edited_destination: {}".format(edited_destination))

                        if read_base == edited_destination:
                            read_edits[read_name]['bases'].append(read_base)
                            read_edits[read_name]['positions'].append(adjusted_pileup_pos)
                            
                            total_edited_bases += 1
        
    total_reads_in_region = 0
    for read in samfile.fetch(chrom, start, end):
        total_reads_in_region += 1        
        
    num_edited_reads = len(list(read_edits.keys()))
    
    #print('>>>>>', chrom, start, end, target_bases, pos_edits)
    return num_edited_reads, np.mean(coverages), total_reads_in_region, len(substrate_bases)
    

def total_conversions_in_window(stamp_sites, chrom, start, end, strand, edit_type):
    stamp_sites_in_window = stamp_sites[(stamp_sites.start_stamp >= start) &
                 (stamp_sites.end_stamp <= end) & 
                 (stamp_sites.chrom_stamp.astype(str) == chrom) &
                 (stamp_sites.strand_stamp == strand)
                ]
    #print(chrom, start, end, strand, edit_type)
    #print('\t', stamp_sites_in_window.num_edited_stamp.to_list())
    return stamp_sites_in_window.num_edited_stamp.sum()


def get_edit_c_info(t):
    region_label, chrom, start, end, strand, edit_type, bam_path, reference_path, include_read_level_info, reverse = t
    edit_fraction_info = {}
    
    d = ReadDensity(
        pos=forward_bw,
        neg=reverse_bw
    )
        
    chrom = str(chrom)
    total_conversions = total_conversions_in_window(stamp_sites, chrom, start, end, strand, edit_type)
    total_target_bases = sum_all_c_coverage_in_window(chrom, start, end, strand, d, edit_type, reverse=reverse)

    # TO REMOVE
    #print("REMOVE: total_target_bases: {}".format(total_target_bases))
    if total_target_bases == 0:
        fraction_bases_edited = 0
    else:
        fraction_bases_edited = total_conversions/total_target_bases
       
            
    if not include_read_level_info:
        edit_fraction_info[region_label] = {
                         'start': start,
                         'end': end,
                         'chrom': chrom,
                         'strand': strand,
                         'target_bases': total_target_bases,
                         'edited_bases': total_conversions,
                         'edit_fraction': fraction_bases_edited
                        }
            
    else:
        edited_origin, edited_destination = get_origin_and_destination_bases(edit_type, strand=strand)
        num_edited_reads, mean_depth, total_reads_in_region, num_substrate_bases = get_num_edited_reads_and_coverage(chrom, start, end, edited_origin, edited_destination, bam_path, reference_path)

        try:
            fraction_reads_edited = num_edited_reads/total_reads_in_region
        except Exception as e:
            fraction_reads_edited = None

        edit_fraction_info[region_label] = {
                                 'start': start,
                                 'end': end,
                                 'chrom': chrom,
                                 'strand': strand,
                                 'target_bases': total_target_bases,
                                 'edited_bases': total_conversions,
                                 'edit_fraction': fraction_bases_edited,
                                 'num_edited_reads': num_edited_reads,
                                 'total_reads_in_region': total_reads_in_region,
                                 'fraction_reads_edited': fraction_reads_edited,
                                 'mean_depth': mean_depth,
                                 'num_substrate_bases': num_substrate_bases
                                }
            
    region_level_edit_fraction_info = pd.DataFrame(edit_fraction_info).T
    
    return region_level_edit_fraction_info


def load_stamp_sites(stamp_sites_filepath):
    stamp_sites = pd.read_csv(stamp_sites_filepath,
                                             sep='\t', names=['chrom_stamp', 'start_stamp', 'end_stamp', 'conf_stamp', 'coverage_stamp', 'strand_stamp', 
                                                              'geneid', 'genename', 'region', 'annot'])
    stamp_sites['num_edited_stamp'] = [int(float(i.split(',')[0])) for i in stamp_sites.coverage_stamp]
    stamp_sites['total_coverage_stamp'] = [int(float(i.split(',')[1])) for i in stamp_sites.coverage_stamp]
    return stamp_sites

def clean_edit_type_param(edit_type):
    edit_type = edit_type.upper()
    edit_type = edit_type.replace('I', 'G')
    return edit_type

def load_json_info(input_json_filepath, regions_override=None):
    # Load parameters/filepaths for stamp sites, bigwigs, and regions to probe editC in
    with open(input_json) as f:
        data = json.load(f)

    bam_path = data.get('bam')
    output_folder = data.get('output_folder')

    label = data.get('label')
    stamp_sites_file = data.get('stamp_sites_file')
    forward_bw = data.get('forward_bw')
    reverse_bw = data.get('reverse_bw')
    fasta = data.get('fasta')
    regions_for_edit_c = data.get('regions_for_edit_c')
    edit_type = data.get('edit_type')
    reverse = data.get('reverse', False)
    
    edit_type = clean_edit_type_param(edit_type)
    
    assert(edit_type in ['AG', 'AC', 'AT', 'GT', 'GC', 'GA', 'CT', 'CG', 'CA', 'TC', 'TG', 'TA'])
    
    if regions_override:
        regions_for_edit_c = regions_override

    print("Label: {}".format(label))
    print("STAMP sites: {}".format(stamp_sites_file))
    print(".bam path: {}".format(bam_path))
    print("Forward bigwig: {}".format(forward_bw))
    print("Reverse bigwig: {}".format(reverse_bw))
    print("Reference being used: {}".format(fasta))
    print("Regions being annotated: {}\n".format(regions_for_edit_c))
    print("Edit type?: {}".format(edit_type))
    print("Reverse coverage calculation?: {}".format(reverse))
    
    # Load STAMP Sites and regions to probe
    print('\n\nLoading STAMP sites...')
    stamp_sites = load_stamp_sites(stamp_sites_file)
    # Regions should have columns: 
    # region_id, chrom, start, end, strand
    
    print('\n\nLoading regions to annotate with editC information -- should have following columns: region_id, chrom, start, end, strand')
    regions = pd.read_csv(regions_for_edit_c, sep='\t')
    
    return output_folder,label, stamp_sites, forward_bw, reverse_bw, fasta, regions, edit_type, bam_path, reverse


def keep_regions_with_edits(stamp_sites, regions_for_edit_c):
    original_region_count = len(regions_for_edit_c)
    print('Original region count: {}'.format(original_region_count))
    
    stamp_sites_bed = pybedtools.BedTool.from_dataframe(stamp_sites[['chrom_stamp', 'start_stamp', 'end_stamp', 'strand_stamp', 'strand_stamp', 'strand_stamp']])
    regions_for_edit_c_bed =  pybedtools.BedTool.from_dataframe(regions_for_edit_c[['chrom', 'start', 'end', 'region_id', 'subregion', 'strand', 'region', 'gene']])
    
    regions_to_keep = regions_for_edit_c_bed.intersect(stamp_sites_bed, wa=True, wb=True, s=True).to_dataframe(names=['chrom', 'start', 'end', 'region_id', 'subregion', 'strand', 'region', 'gene', 'chrom_stamp', 'start_stamp', 'end_stamp', 'strand_stamp', 'strand_stamp1', 'strand_stamp2'])
    
    pybedtools.cleanup()

    new_region_count = len(regions_to_keep)
    print('New region count: {}'.format(new_region_count))

    # TO REMOVE
    #print("REMOVE stamp: {}".format(regions_to_keep[['chrom_stamp', 'start_stamp', 'end_stamp', 'strand_stamp', 'strand_stamp1', 'strand_stamp2']]))
    regions_to_keep = regions_to_keep[['region_id', 'chrom', 'start', 'end', 'strand', 'subregion', 'region', 'gene']].drop_duplicates()

    return regions_to_keep


def print_time():
    now = datetime.datetime.now()
    print ("Current date and time : ")
    print (now.strftime("%Y-%m-%d %H:%M:%S"))

    
def edit_fraction(output_folder, label, forward_bw, reverse_bw, fasta, filtered_regions, bam_path, edit_type, include_read_level_info=False, outfile_override=False, cores=32, reverse=False):    
    region_ids = list(filtered_regions.region_id)
    chroms = list(filtered_regions.chrom)
    starts = list(filtered_regions.start)
    ends = list(filtered_regions.end)
    strands = list(filtered_regions.strand)
    edit_type_list = [edit_type for i in range(len(filtered_regions))]  
    bam_path_list = [bam_path for i in range(len(filtered_regions))]
    reference_path_list = [fasta for i in range(len(filtered_regions))]
    include_read_level_info_list = [include_read_level_info for i in range(len(filtered_regions))]
    reverse_list = [reverse for i in range(len(filtered_regions))]
        
    print_time()
     # Set up a multi-processing pool to process each region individually and add editC information 
    print("\n\nCalculating edit fraction information for each region using multi-processing pool... this step could take a while (fastest when ppn = 8)")
    if include_read_level_info:
        print("\nIncluding read level information")
    else:
        print("Not including read level information")
    num_processes = cores

    p = Pool(num_processes)
    region_level_edit_c_dfs = p.map(get_edit_c_info, zip(region_ids, chroms, starts, ends, strands, edit_type_list, bam_path_list, reference_path_list, include_read_level_info_list, reverse_list))
    p.close()
    p.join()

    edit_c_for_all_regions = pd.concat(region_level_edit_c_dfs)
    
    col_order = ['chrom', 'start', 'end', 'edit_fraction', 'strand', 'target_bases', 'edited_bases']
    if include_read_level_info:
        col_order = col_order + ['num_edited_reads', 'total_reads_in_region', 'fraction_reads_edited', 'mean_depth', 'num_substrate_bases']

    edit_c_for_all_regions = edit_c_for_all_regions.fillna(int(0))[col_order]

    print_time()
    return edit_c_for_all_regions


def get_collapsed_genes(gene_and_region_span):
    gene_collapse = gene_and_region_span.groupby('gene').gene.min()
    chrom_collapse = gene_and_region_span.groupby('gene').chrom.min()
    start_collapse = gene_and_region_span.groupby('gene').start.min()
    end_collapse = gene_and_region_span.groupby('gene').end.max()
    strand_collapse = gene_and_region_span.groupby('gene').strand.min()

    collapsed_df = pd.concat([gene_collapse, chrom_collapse, start_collapse, end_collapse, strand_collapse], axis=1)
    collapsed_df.columns = ['region_id', 'chrom', 'start', 'end', 'strand']
    return collapsed_df

def get_collapsed_subregions(gene_and_region_span):
    gene_and_region_collapse = gene_and_region_span.groupby('gene_and_subregion').gene_and_subregion.min()
    chrom_collapse = gene_and_region_span.groupby('gene_and_subregion').chrom.min()
    start_collapse = gene_and_region_span.groupby('gene_and_subregion').start.min()
    end_collapse = gene_and_region_span.groupby('gene_and_subregion').end.max()
    strand_collapse = gene_and_region_span.groupby('gene_and_subregion').strand.min()

    collapsed_df = pd.concat([gene_and_region_collapse, chrom_collapse, start_collapse, end_collapse, strand_collapse], axis=1)
    collapsed_df.columns = ['region_id', 'chrom', 'start', 'end', 'strand']
    return collapsed_df

def get_dict_for_gene_data(gene_edit_fractions):
    gene_counts = defaultdict(lambda:{})
    for r in gene_edit_fractions.iterrows():
        gene = r[0]

        r = r[1]
        #print(gene)

        coverage = r.target_bases
        conversions = r.edited_bases

        gene_counts[gene]['coverage'] = coverage
        gene_counts[gene]['conversions'] = conversions
    return gene_counts

def get_dicts_for_region_and_subregion_data(subregion_edit_fractions):
    subregion_counts = defaultdict(lambda:{})
    region_counts = defaultdict(lambda:defaultdict(lambda:0))
    for r in subregion_edit_fractions.iterrows():
        subregion = r[0]
        region = subregion.split('_')[0]

        r = r[1]
        #print(subregion, region)

        coverage = r.target_bases
        conversions = r.edited_bases

        subregion_counts[subregion]['coverage'] = coverage
        subregion_counts[subregion]['conversions'] = conversions

        region_counts[region]['coverage'] += coverage
        region_counts[region]['conversions'] += conversions
    return region_counts, subregion_counts


def add_region_and_subregion_conversions(r):
    region_id = r.region_id
    subregion = ':'.join(region_id.split(':')[0:2])
    region = subregion.split('_')[0]
    gene = region.split(':')[0]
    
    try:
        subregion_coverage, subregion_conversions = subregion_counts.get(subregion)['coverage'], subregion_counts.get(subregion)['conversions']
    except Exception as e:
        print(e)
        print(region_id, subregion, region, gene)
        subregion_coverage = 0
        subregion_conversions = 0
    
    return subregion_coverage, subregion_conversions


parser = argparse.ArgumentParser(description='Calculate edit fraction in given regions')
parser.add_argument('input_json', type=str)
parser.add_argument('--outfile_override', type=str, default=None)
parser.add_argument('--regions_override', type=str, default=None)
parser.add_argument('--keep_all', type=bool, default=False)
parser.add_argument('--final_peaks', type=bool, default=False)
parser.add_argument('--cores', type=int, default=32)

args = parser.parse_args()
input_json = args.input_json
regions_override = args.regions_override
outfile_override = args.outfile_override
keep_all = args.keep_all
final_peaks = args.final_peaks
cores = args.cores

now = datetime.datetime.now()
print ("Current date and time : ")
print (now.strftime("%Y-%m-%d %H:%M:%S"))
    
print('input is {}'.format(input_json))
print('cores: {}'.format(cores))

suffix = ''
if regions_override:
    print('Regions override provided: {}'.format(regions_override))
    suffix = regions_override.split('_')[-1]
    
if keep_all:
    print('Keeping all peaks...')
    final_peaks = True
    
if final_peaks:
    print("Asessing final merged peaks... will not report background editing values")
    
    
print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
print("Calculating Edit-C Region Baseline")
print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~")
output_folder,label, stamp_sites, forward_bw, reverse_bw, fasta, regions, edit_type, bam_path, reverse = load_json_info(input_json, regions_override)

print('Looking at {}>{} edits'.format(edit_type[0], edit_type[1]))

    
if outfile_override:
    print("Outfile override provided: {}".format(outfile_override))
    out_file = outfile_override
else:
    output_folder = '{}/editc_outputs'.format(output_folder)
    print('loaded, calculating editC')
    out_dir = '{}/{}'.format(output_folder, label)
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)
            
    out_file = '{}/{}_edit_c_for_all_regions_{}_all_info.tsv'.format(out_dir, label, suffix)

print('Outputting to {}'.format(out_file)) 

# Keeping only windows with edits
if not keep_all:
    filtered_regions = keep_regions_with_edits(stamp_sites, regions)
    filtered_regions.to_csv('{}.filtered_regions.tsv'.format(out_file.split('.tsv')[0]), sep='\t', header=False, index=False)
elif keep_all:
    # load pre-filtered regions if we are just looking at editing in the control sample but within regions decided from the test sample
    filtered_regions = pd.read_csv(regions_override, sep='\t')
    print("{} regions provided to analyze...".format(len(filtered_regions)))

if len(filtered_regions) == 0:
    print("No edits to analyze")
    with open(out_file, 'w') as f:
            f.write('region_id\tchrom\tstart\tend\tedit_c\tstrand\tcoverage\tconversions\toverall_region_id\tchrom_region\tstart_region\tend_region\tedit_c_region\tstrand_region\tcoverage_region\n')
    sys.exit(0)
    

print("Analyzing edit fraction for {} filtered windows...".format(len(filtered_regions)))

#include_read_level_info = False# True <-- read level info is failing for MARINE outputs, removing for now
include_read_level_info = True
if keep_all:
  include_read_level_info = False  
    
window_edit_fractions = edit_fraction(output_folder,label, forward_bw, reverse_bw, fasta, filtered_regions, bam_path, edit_type=edit_type, outfile_override=outfile_override, include_read_level_info=include_read_level_info, cores=cores, reverse=reverse)
window_edit_fractions['region_id'] = window_edit_fractions.index

if not final_peaks:
    print("Adding contextual gene, region, and subregion editing information to dataframe...")

    regions['gene_and_subregion'] = regions.gene + ':' + regions.subregion
    regions['gene_and_region'] = regions.gene + ':' + regions.region

    filtered_regions['gene_and_region'] = filtered_regions.gene + ':' + filtered_regions.region
    filtered_regions['gene_and_subregion'] = filtered_regions.gene + ':' + filtered_regions.subregion

    unique_filtered_gene_and_subregion = list(filtered_regions['gene_and_subregion'].unique())
    gene_and_region_span = regions[regions.gene_and_subregion.isin(unique_filtered_gene_and_subregion)] 
    gene_and_region_span['chrom'] = gene_and_region_span.chrom.astype(str)

    collapsed_subregions = get_collapsed_subregions(gene_and_region_span)
    print("Analyzing edit fraction for {} collapsed subregions (ie introns and exons)...".format(len(collapsed_subregions)))
    subregion_edit_fractions = edit_fraction(output_folder,label, forward_bw, reverse_bw, fasta, collapsed_subregions, bam_path, edit_type=edit_type, outfile_override=outfile_override, cores=cores, reverse=reverse)
    region_counts, subregion_counts = get_dicts_for_region_and_subregion_data(subregion_edit_fractions)
    
    window_edit_fractions['subregion_coverage'],\
    window_edit_fractions['subregion_conversions'] = zip(*window_edit_fractions.apply(add_region_and_subregion_conversions, axis=1))

window_edit_fractions = window_edit_fractions.drop('region_id', axis=1)   

    
window_edit_fractions.to_csv(out_file, sep='\t', index=True)

sys.exit(0)