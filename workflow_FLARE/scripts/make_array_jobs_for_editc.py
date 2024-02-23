import argparse
from glob import glob
import os
import math
from collections import defaultdict
import subprocess
    

# Chunk regions
def get_split_regions(regions_directory):
    # Region file chunk should end with _{index} ....
    regions = sorted(glob('{}/*_*'.format(regions_directory)))
    #regions = regions[0:10]
    regions = [r.split('{}/'.format(regions_directory))[1] for r in regions]
    len(regions)

    # Chunk out the regions so we don't surpass the 500 task limit per job
    split_regions = defaultdict(lambda:[])

    for i, r in enumerate(regions):
        bucket = math.floor(i / 500)
        split_regions[bucket].append(r)
    
    return split_regions
    

def make_bash_files(sample_id, split_regions, main_directory, scripts_directory, input_json, outer_window, keep_all=False):
    # Make the bash array-task files to qsub
    # make output directory for job err and out files
    if not os.path.exists('{}/outs'.format(main_directory)):
        os.makedirs('{}/outs'.format(main_directory))
    if not os.path.exists('{}/bash_scripts/{}'.format(main_directory, sample_id)):
        os.makedirs('{}/bash_scripts/{}'.format(main_directory, sample_id))

    if keep_all:
        bash_command =  'python {}/calculate_edit_c_for_regions.py {} --regions_override {}/{} --keep_all True\n'
    else:
        bash_command =  'python {}/calculate_edit_c_for_regions.py {} --regions_override {}/{}\n'
    
    for chunk_id, region_chunk in split_regions.items():
        print('\t\tJob:', chunk_id, 'Tasks:', len(region_chunk))
        for region in region_chunk:
            bash_script = '{}/bash_scripts/{}/{}_bash_job_edit_c_{}.sh'.format(main_directory, sample_id, sample_id, region.split('_')[-1])
            print('\t\tMaking bash scripts:', bash_script)

            with open(bash_script, 'w') as f:
                f.write('#!/bin/bash\n')
                c = bash_command.format(scripts_directory, input_json, regions_directory, region)
                f.write(c)


parser = argparse.ArgumentParser(description='Calculate edit-c.')
parser.add_argument('--sample_id', type=str)
parser.add_argument('--input_json', type=str)
parser.add_argument('--regions_directory', type=str)
parser.add_argument('--main_directory', type=str)
parser.add_argument('--scripts_directory', type=str)
parser.add_argument('--outer_window', type=int, default=0)
parser.add_argument('--keep_all', type=str, default="false")


args = parser.parse_args()
input_json = args.input_json
regions_directory = args.regions_directory
main_directory = args.main_directory
scripts_directory = args.scripts_directory
sample_id = args.sample_id
outer_window = args.outer_window
keep_all = args.keep_all.lower() == "true"

print("\tPreparing array jobs to determine edit c across regions efficiently...")
split_regions = get_split_regions(regions_directory)
make_bash_files(sample_id, split_regions, main_directory, scripts_directory, input_json, outer_window, keep_all=keep_all)
