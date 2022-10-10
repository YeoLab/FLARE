import pandas as pd
import pybedtools
import math
import os
from argparse import ArgumentParser


# num rows per chunk -- you can increase if you want to reduce the number of jobs you kick off, but increase the runtime per job,
# or decrease it if you want to increase the number of jobs you kick off but decrease the runtime per job.
CHUNK_SIZE = 7500000


def get_gene(annot):
    return annot.split('gene_name \"')[1].split('\"')[0]

def to_tuple(gene):
    return gene.split(',')

def sort_and_explode_bed(df):
    df_bed = pybedtools.BedTool.from_dataframe(df[['chrom', 'start', 'end', 'annot', 'gene_name', 'strand']])
    df_bed_merged_sorted = df_bed.sort().merge(s=True, c=[4,5,6], o=['first', 'distinct', 'distinct'])
    df_bed_merged_sorted_df = df_bed_merged_sorted.to_dataframe()
    df_bed_merged_sorted_df.columns = ['chrom', 'start', 'end', 'region', 'gene', 'strand']
    df_bed_merged_sorted_df['gene'] = df_bed_merged_sorted_df.gene.apply(to_tuple)
    df_bed_merged_sorted_df_exploded = df_bed_merged_sorted_df.explode('gene')
    #Filter teeny exons
    df_bed_merged_sorted_df_exploded['length'] = df_bed_merged_sorted_df_exploded.end - df_bed_merged_sorted_df_exploded.start
    print('\t{} rows before length filter...'.format(len(df_bed_merged_sorted_df_exploded)))
    df_bed_merged_sorted_df_exploded = df_bed_merged_sorted_df_exploded[df_bed_merged_sorted_df_exploded.length >= 15]
    print('\t{} rows after length filter...'.format(len(df_bed_merged_sorted_df_exploded)))
    df_bed_merged_sorted_df_exploded_bed = pybedtools.BedTool.from_dataframe(df_bed_merged_sorted_df_exploded)
    return df_bed_merged_sorted_df_exploded_bed

def get_genes_and_exons(gtf_path):
    gtf = pd.read_csv(gtf_path, sep='\t', comment='#',
                      names=['chrom', 'source', 'feature', 'start', 'end', '.', 'strand', ',', 'annot'])

    genes = gtf[gtf.feature == 'gene']
    exons = gtf[gtf.feature == 'exon']

    genes['gene_name'] = genes.annot.apply(get_gene)
    exons['gene_name'] = exons.annot.apply(get_gene)
    
    return genes, exons


def get_intron_and_exons_df(genes_bed, exons_bed):
    # Take the difference between genes and exons regions to get intron regions. Then combine intron and exon regions to subdivide the genome by just these categories within gene regions.
    introns_bed = genes_bed.subtract(exons_bed)
    introns_df = introns_bed.to_dataframe()
    introns_df.columns = ['chrom', 'start', 'end', 'region', 'gene', 'strand', 'length']
    introns_df = introns_df.drop('length', axis=1)
    introns_df.region = 'intron'

    exon_df = exons_bed.to_dataframe()
    exon_df.columns = ['chrom', 'start', 'end', 'region', 'gene', 'strand', 'length']
    exon_df = exon_df.drop('length', axis=1)
    exon_df.region = 'exon'


    introns_and_exons_df = (pybedtools.BedTool.from_dataframe(introns_df).cat(pybedtools.BedTool.from_dataframe(exon_df), 
                                                                              postmerge=False)).to_dataframe()
    introns_and_exons_sorted_df = pybedtools.BedTool.from_dataframe(introns_and_exons_df).sort().to_dataframe()
    introns_and_exons_sorted_df.columns = ['chrom', 'start', 'end', 'region', 'gene', 'strand']
    return introns_and_exons_sorted_df


def output_region_files(introns_and_exons_sorted_df, output_dir, window_size):
    all_windows = []

    window_overlap = int(window_size/2)
    
    print('window size:', window_size)
    print('window_overlap:', window_overlap)

    to_iter = introns_and_exons_sorted_df
    print('Overall number of introns and exons: ', len(to_iter))

    overall_chunks = 0

    gene_list = to_iter.gene.unique()
    total_windows_processed = 0

    for i, gene in enumerate(gene_list):
        region_num = 0
        for j,r in enumerate(to_iter[to_iter.gene == gene].iterrows()):
            r = r[1]
            chrom = r.chrom
            strand = r.strand
            start = r.start
            end = r.end
            region = r.region
            gene = r.gene
            subregion = '{}_{}'.format(region, region_num)

            while start + window_size <= end:
                window_id = '{}:{}_{}:{}-{}'.format(gene, region, region_num, start, start+window_size)
                window = [window_id, chrom, start, start+window_size, strand, subregion, region, gene]
                all_windows.append(window)
                start += window_size - window_overlap

            if start < end:
                window_id = '{}:{}_{}:{}-{}'.format(gene, region, region_num, start, end)
                window = [window_id, chrom, start, end, strand, subregion, region, gene]
                all_windows.append(window)

            region_num += 1

        if len(all_windows) >= CHUNK_SIZE or i == len(gene_list) - 1:
            chunk_id = overall_chunks
            overall_chunks += 1
            print('\tchunk id is ', chunk_id)
            print('\toutputting ', len(all_windows), 'rows')
            pd.DataFrame(all_windows, columns=['region_id', 'chrom', 'start', 'end', 'strand', 'subregion', 'region', 'gene']).to_csv('{}/{}_{}'.format(output_dir, output_dir, chunk_id), sep='\t', index=False, header=True)


            total_windows_processed += len(all_windows)
            # Print progress
            print('{}/{} genes processed...'.format(i, len(gene_list)))

            all_windows = []
            

def main(gtf_path, output_dir, window_size):
    print('Extract genes and exons from {}...'.format(gtf_path))
    genes, exons = get_genes_and_exons(gtf_path)
    print('Process exons...')
    exons_bed_merged_sorted_df_exploded_bed = sort_and_explode_bed(exons)
    print('Process genes...')
    genes_bed_merged_sorted_df_exploded_bed = sort_and_explode_bed(genes)
    print("Making a dataframe including all introns and exons, ordered, for each gene...")
    introns_and_exons_sorted_df = get_intron_and_exons_df(genes_bed_merged_sorted_df_exploded_bed, exons_bed_merged_sorted_df_exploded_bed)
    print("Outputing regions files... this will take a while.")
    output_region_files(introns_and_exons_sorted_df, output_dir, window_size)
    print("Done!")


        
if __name__ == "__main__":
    parser = ArgumentParser(description='Calculate edit fraction in given regions')
    parser.add_argument('gtf_path', type=str)
    parser.add_argument('output_dir', type=str)
    parser.add_argument('--window_size', type=int, default=30)

    args = parser.parse_args()
    gtf_path = args.gtf_path
    output_dir = args.output_dir
    window_size = args.window_size

    #gtf_path = '/projects/ps-yeolab3/bay001/annotations/hg19/gencode_v19/gencode.v19.annotation.gtf'
    #output_dir = 'hg19_script'

    print('gtf path:', gtf_path)
    print('output dir:', output_dir)

    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
        
    main(gtf_path, output_dir, window_size)
    