# Generate stranded bigwig (.bw) files from a bam file.
#     Note: make sure you are using the correct chrom sizes file for your genome build

#module load samtools
#module load bedtools
#module load makebigwigfiles;


#Sometimes have to link for new conda environments, if getting error:
##     bedGraphToBigWig: error while loading shared libraries: libcrypto.so.1.0.0: cannot open shared object file: No such file or directory
# Solution:
#ln -s /tscc/nfs/home/ekofman/miniconda3/envs/workhorse/lib/libssl.so.3 /tscc/nfs/home/ekofman/miniconda3/envs/workhorse/lib/libssl.so.1.0.0
#ln -s /tscc/nfs/home/ekofman/miniconda3/envs/workhorse/lib/libcrypto.so.3 /tscc/nfs/home/ekofman/miniconda3/envs/workhorse/lib/libcrypto.so.1.0.0

# INPUTS:

# This should be the full path to a bam file from which you wish to make stranded bigwigs files.
input_bam=$1
samplename=$2

# Reference chrom sizes data is generated from index of fasta used for alignment.
# e.g. /projects/ps-yeolab3/bay001/annotations/hg19/hg19.chrom.sizes
# Input here is reference fasta used for aligment, used to create the chrom.sizes file
fasta=$3


# The directory where all outputs from this script should be placed -- a folder will be generated within this folder with the name of each sample,
# and all files for that particular sample will be place within this sample-specific subfolder.
output_dir=$4




echo "Make bigwigs for $input_bam..."
echo "Sample name is $samplename... will output in $output_dir..."

mkdir -p $output_dir

echo "Splitting by strand into two begraphs using genomecov..."
bedtools genomecov -split -strand - -bg -ibam $input_bam > $output_dir/$samplename.fwd.bg;
echo "Done with fwd..."
bedtools genomecov -split -strand + -bg -ibam $input_bam > $output_dir/$samplename.rev.bg;
echo "Done with rev..."


echo "Sorting each split file..."
bedtools sort -i $output_dir/$samplename.fwd.bg > $output_dir/$samplename.fwd.sorted.bg;
echo "Done with fwd..."
bedtools sort -i $output_dir/$samplename.rev.bg > $output_dir/$samplename.rev.sorted.bg;
echo "Done with rev..."


echo "Removing intermediate files..."
rm $output_dir/$samplename.fwd.bg
rm $output_dir/$samplename.rev.bg

echo "Making chrom.sizes file from $fasta_index ($fasta.fai)"
echo "cut -f 1,2 $fasta.fai > $output_dir/chrom.sizes"

cut -f 1,2 $fasta.fai > $output_dir/chrom.sizes

echo "Bedgraph to bigwig"
bedGraphToBigWig $output_dir/$samplename.fwd.sorted.bg $output_dir/chrom.sizes $output_dir/$samplename.fwd.sorted.bw
echo "Done with fwd..."
bedGraphToBigWig $output_dir/$samplename.rev.sorted.bg $output_dir/chrom.sizes $output_dir/$samplename.rev.sorted.bw
echo "Done with rev..."

echo "Done with $samplename!"

