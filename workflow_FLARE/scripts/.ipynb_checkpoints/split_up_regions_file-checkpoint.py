import argparse
import os
import pandas as pd
import sys

parser = argparse.ArgumentParser(description='Split up regions file')
parser.add_argument('regions_filepath', type=str)
parser.add_argument('split_regions_folder', type=str)
parser.add_argument('--chunk_size', type=int, default=5000)

args = parser.parse_args()
regions_filepath = args.regions_filepath
split_regions_folder = args.split_regions_folder
chunk_size = args.chunk_size


df = pd.read_csv(regions_filepath, sep='\t')
list_df = [df[i:i+chunk_size] for i in range(0,df.shape[0],chunk_size)]

sys.stdout.write('{} being split into {} chunks of size {} in {}\n'.format(regions_filepath, len(list_df), chunk_size, split_regions_folder))

basename = os.path.basename(regions_filepath).split('.')[0]

if not os.path.isdir(split_regions_folder):
    os.mkdir(split_regions_folder)
    
for i,chunk in enumerate(list_df):
    filename_for_chunk = '{}/{}_{}'.format(split_regions_folder, basename, i)
    print("Chunk #{} with destination of {}".format(i, filename_for_chunk))

    chunk.to_csv(filename_for_chunk, sep='\t', index=False)

    
    