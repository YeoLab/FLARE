import argparse
import os
from glob import glob
import pandas as pd
import sys    
import pandas as pd
import json




# Validate paths provided
def validate_paths(paths, verbose=False):
    for to_validate in paths:
        if os.path.exists(to_validate):
            if verbose:
                sys.stdout.write('Confirming {} exists...\n'.format(to_validate))
            else:
                pass
        else:
            sys.stdout.write('Error: {} does not exist!\n'.format(to_validate))
            sys.exit(1)
    

def create_jsons(sailor_output_folder, peakcalling_input_folder, peakcalling_output_folder, fasta, edit_type, fdr_filter,
                 cluster_identification_regions_folder=None,
                 edit_fraction_regions_file=None
                ):        
    validate_paths([sailor_output_folder, peakcalling_input_folder, peakcalling_output_folder, peakcalling_regions_folder, fasta], verbose=True)

    if not (cluster_identification_regions_folder or edit_fraction_regions_file):
        sys.stdout.write("Either a regions folder (Cluster Identification Mode) or regions file (Edit Fraction Mode) must be provided!")
        sys.exit(1)
    if (cluster_identification_regions_folder and edit_fraction_regions_file):
        sys.stdout.write("Found both regions folder (Cluster Identification Mode) and regions file (Edit Fraction Mode)! Only provide one.")
        sys.exit(1)
        
    if edit_fraction_regions_file:
        sys.stdout.write("Edit Fraction Mode using {}".format(edit_fraction_regions_file))
        keep_all = True
    else:
        sys.stdout.write("Cluster Identification Moe using {}".format(cluster_identification_regions_folder))
        keep_all = False
    
    stamp_site_files = glob('{}/*bed'.format(sailor_output_folder))
    sys.stdout.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
    sys.stdout.write("Found {} finished SAILOR outputs in {}...\n".format(len(stamp_site_files), sailor_output_folder))
    for s in stamp_site_files:
        sys.stdout.write('\t{}\n'.format(s))
    sys.stdout.write('~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n')
    sys.stdout.write('Writing JSON inputs for FLARE for each SAILOR output...\n')
    input_json_dicts = []
    for stamp_site_file in stamp_site_files:
        sample_id = stamp_site_file.split('.combined')[0].split('/')[-1]

        forward_bw = '{}/8_bw_and_bam/{}.fwd.sorted.bw'.format(sailor_output_folder, sample_id)
        reverse_bw = '{}/8_bw_and_bam/{}.rev.sorted.bw'.format(sailor_output_folder, sample_id)
        bam = '{}/8_bw_and_bam/{}_filtered_merged.sorted.bam'.format(sailor_output_folder, sample_id)

        # Validate presence of all files necessary to run peakcaller 
        validate_paths([forward_bw, reverse_bw, bam])

        input_json_dict = {
            "label": sample_id,
            "output_folder": peakcalling_output_folder,
            "stamp_sites_file": stamp_site_file,
            "forward_bw": forward_bw,
            "reverse_bw": reverse_bw,
            "fasta": fasta,
            "edit_type": edit_type,
            "bam": bam,
            "fdr_threshold": fdr_filter,
            "keep_all": keep_all
        }

        if cluster_identification_regions_folder:
            input_json_dict['regions'] = cluster_identification_regions_folder
            json_filename = '{}/{}_cluster_identification_input.json'.format(peakcalling_input_folder, sample_id)

        elif edit_fraction_regions_file:
            input_json_dict['regions_file'] = edit_fraction_regions_file
            json_filename = '{}/{}_edit_fraction_input.json'.format(peakcalling_input_folder, sample_id)
        
        sys.stdout.write("\tWriting {}...\n".format(json_filename))
        with open(json_filename, 'w') as f:
            f.write(json.dumps(input_json_dict))

    sys.stdout.write('Done!\n')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Given finished SAILOR output folders, generate appropriate peakcalling input configuration jsons')
    parser.add_argument('sailor_output_folder', type=str)
    parser.add_argument('peakcalling_input_folder', type=str)
    parser.add_argument('peakcalling_output_folder', type=str)
    parser.add_argument('fasta', type=str)
    parser.add_argument('edit_type', type=str)
    parser.add_argument('fdr_filter', type=float, default=0.1)
    
    parser.add_argument('--cluster_identification_regions_folder', type=str, default=None)
    parser.add_argument('--edit_fraction_regions_file', type=str, default=None)


    args = parser.parse_args()
    sailor_output_folder = args.sailor_output_folder
    peakcalling_input_folder = args.peakcalling_input_folder
    peakcalling_output_folder = args.peakcalling_output_folder
    fasta = args.fasta
    cluster_identification_regions_folder = args.peakcalling_regions_folder
    edit_fraction_regions_file = args.peakcalling_regions_folder

    edit_type = args.edit_type
    fdr_filter = args.fdr_filter 
    
    create_jsons(sailor_output_folder, peakcalling_input_folder, peakcalling_output_folder, fasta, edit_type, fdr_filter,
                 cluster_identification_regions_folder=cluster_identification_regions_folder,
                 edit_fraction_regions_file=edit_fraction_regions_file
                )