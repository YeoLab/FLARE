# FLARE

### FLagging Areas of RNA-editing Enrichment (FLARE) 

We present FLagging Areas of RNA-editing Enrichment (FLARE),  a Snakemake-based pipeline that uses a statistical approach to determine regions of enriched RNA editing, using SAILOR-derived editing sites as a starting point. FLARE is configurable for use with any type of base pair change – we include with this release of FLARE an update of SAILOR to enable detection of all edit types.

# Requirements

* Your system must be at least Linux Centos7

* Make sure that the environment your STAMP processing pipelines will be running on have snakemake installed (https://snakemake.readthedocs.io/en/v5.6.0/getting_started/installation.html). 

* You also will need to have Singularity installed for several steps of the pipeline to work. We have created singularity images containing all necessary python packages that will automatically be loaded for you in the course of running the pipeline, as long as your system has singularity installed.

# Running the SAILOR snakemake pipeline

## Before you start:

* To get the known snps bedfile for all chromosome, download the appropriate individual chromosome bedfiles, for example from here: https://ftp.ncbi.nih.gov/snp/organisms/, then combine them all (remembering to remove any header lines, and retain only chrom start end columns), ie:

    ```
    for b in $(ls *.bed); do echo $b; tail -n+2 $b | cut -f1,2,3  | sort >> mm10_dbsnp_combined.bed3; done
    ```

    Make sure that your chromosome nomenclature is the same as in the fasta file you are using ("chr1" vs "1")!
    This file should then contain lots of lines like this, tab-separated, without headers:
    
    1   1334235  1334236
    ...

## Parameters
All SAILOR configuration information must be saved in a .json file with the following contents:
```
{
  "samples_path": "/path/to/aligned/bams/",
  "samples": [
    "sample1.sortedByCoord.out.bam",
    "sample2.sortedByCoord.out.bam",
    ....,
    etc.
  ],
  "reverse_stranded": true,
  "reference_fasta": "/path/to/fasta/used/to/align/genome.fa",
  "known_snps": "/path/to/known/common/snps/file/for/organism/in/chrom/start/end/format/b151_GRCh38p7_common.bed3",
  "edit_type": "CT", (or "AG", "TC", etc.)
  "output_dir": "/path/to/output/directory"
}
```
Create your .json config file and call it something sensible based on your experiment, for example 'sailor.json'. Multiple .bam files contained in one directory can be processed with one run of the SAILOR pipeline.

In order to run a snakemake pipeline, there a few parameters that snakemake needs to know about. The first is which Snakefile to use -- the Snakefile contains the instructions for running each step of the pipeline, and for the SAILOR pipeline will be found at your local version of /STAMP/workflow_sailor/Snakefile. The second is which config file to use -- this is the config file you just created, which contains the parameters particular to your run of the pipeline. Always absolute paths. You will also need to tell snakemake to use singularity, and specify singularity arugments allowing the virtual environments to have access to your local filesystem -- the "bind" parameters should reflect locations that the pipeline should have access to, for example folders containing relevant input bams, fastas, gtfs or dbsnp files. So, an example snakemake run could like like the following:

`
snakemake --snakefile /full/path/to/STAMP/workflow_sailor/Snakefile --configfile /full/path/to/your/config/file/sailor.json --verbose --use-singularity     --singularity-args '--bind /home --bind /projects'  -j1
`

However, this will launch the snakemake pipeline on your head node, and all subsequent snakemake jobs (which can number into the hundreds depending on how many samples you are processing) will also be run there. 

## Instructions on HPCC (cluster with submission):
If you have access to a high performance compute cluster, you probably want jobs to be automatically submitted to this system instead so that you can take full advantage of the parallelization built into this pipeline, especially if you are analyzing many samples. Cluster submission information can be placed into a "profile" file. In this case, you can model your profile file on the file at /full/path/to/STAMP/profiles/tscc_sailor/config.yaml, which by default has the following contents:

```
cluster: "qsub -N {rule}.{wildcards} -l nodes=1:ppn={params.threads},walltime={params.run_time} -A yeo-group -q home -V -t 0"
verbose: true
notemp: false
latency: 300
printshellcmds: true
directory: .
snakefile: /path/to/STAMP/workflow_sailor/Snakefile
use-singularity: true
singularity-args: '--bind /oasis --bind /projects --bind /home'
jobs: 8
skip-script-cleanup: true
singularity-prefix: /projects/ps-yeolab4/software/stamp/0.99.0/bin/.singularity
nolock: true
```

* Change "cluster" to be relevant to your cluster system, specifically changing the names of the parameters to match your system's requirements.
    * When in doubt, just hardcode the nodes (in this example 'ppn') to be equal to 1 and run time (in this example 'walltime') values at 5 hours. You can increase the walltime value if your job is running out of time. 
* Change "directory" to be the full path to the directory where you want log files for each step to be deposited during the snakemake run
* Change "snakefile" to the absolute path for your version of the workflow_sailor Snakefile, as mentioned earlier
* Change "singularity-args" to include the correct directory binding relevant to your system so snakemake can find all necessary files
* Change "singularity-prefix" to reflect the absolute path to where you want the singularity images used for the run to be stored (should have a lot of space)
* Change "jobs" to reflect the maximum number of jobs you want submitted to your cluster simulataneously

Put your profile configuration .yaml file in a new folder you can call /full/path/to/STAMP/profiles/my_profile/, for example.
With more run information tucked away into the profile file, the snakemake launch command becomes simpler as it can reference the parameters from this profile file (note that it is actually the folder containing the profile.yaml file that is specified, not the file itself):

`
snakemake --profile /full/path/to/STAMP/profiles/my_profile/ --configfile /full/path/to/your/config/file/sailor.json
`

Running that command should launch your SAILOR run.
An example set of completed outputs from a successful SAILOR run, using the small .bam file and config inputs found in the "examples" folder, looks like this. Note that these folders may not necessarily appear chronologically in an order matching their number. The final outputs are the .ranked.bed files.

```
1_split_strands
3_index_reads
4_filter_reads
5_pileup_reads
6_vcfs
7_scored_outputs
8_bw_and_bam
9_edit_fraction_bedgraphs
subsampled.bam.combined.readfiltered.formatted.varfiltered.snpfiltered.ranked.bed
```

# Running the peakcalling snakemake pipeline

## Before you start:

To run the peakcalling pipeline, you will first need a file specifying the regions in which peakcalling will occur. To generate this file, use the script in workflow_peakcaller/scripts called generate_peakcalling_regions.py

Copy this script to wherever you'd like to generate the helper files, which will be take up about 8-10 GB of space. Then run the script by typing   `generate_peakcalling_regions.py <full/path/to/your/genome/gtf/file> <genome_name>_peakcalling_regions`

The .gtf file you use should include gene and exon level information, i.e. the third column should at least contain the descriptors "gene" and "exon."

If using the following command, for example:
`generate_peakcalling_regions.py <full/path/to/your/genome/gtf/mm10.gtf > mm10_peakcalling_regions`

Once the script finishes running, you will see a new folder called mm10_peakcalling_regions, and within that folder a slew of files with increasing indices, i.e. mm10_peakcalling_regions_0, mm10_peakcalling_regions_1... 

## Parameters
All peakcaller configuration information must be saved in a .json file with the following contents:
```
{
    "label": "label_for_this_sample",
    "output_folder": "/absolute/path/to/peakcalling/folder/outputs_where_all_samples_will_be_placed/peakcalling_outputs/",
    "stamp_sites_file": "/absolute/path/to/sailor/output/this_sample.bam.combined.readfiltered.formatted.varfiltered.snpfiltered.ranked.bed",
    "forward_bw": "/absolute/path/to/sailor/output/8_bw_and_bam/this_sample.sortedByCoord.out.bam.fwd.sorted.bw",
    "reverse_bw": "/absolute/path/to/sailor/output/8_bw_and_bam/this_sample.sortedByCoord.out.bam.rev.sorted.bw",
    "fasta": "/path/to/fasta/used/to/align/genome.fa",
    "regions": "/absolute/path/to/peakcalling_regions/folder",
    "edit_type": "CT", (or "AG", "TC", etc.-- should be same as what was used in SAILOR run)
    "bam": "/absolute/path/to/sailor/output/8_bw_and_bam/this_sample.sortedByCoord.out.bam_filtered_merged.sorted.bam"
}
```

Note that this specifies parameters for peakcalling of only *one* sample. Four of the inputs to the peakcalling pipeline are SAILOR outputs for this sample,
and of those, three are specifically contained in the 8_bw_and_bam folder from that SAILOR run. Other than that, you must specify the edit type (should be the same as whatever was specified for the SAILOR run), a label for your sample, the overall output folder where you want the peakcalling folder for this sample to be generated, and the fasta that was used to align the genome (same sample as was used for the SAILOR run). 

### Use the make_peakcalling_jsons.py to set up configuration .json files for multiple samples
If you are calling peaks on many samples, it will be too time consuming to manually create .json files for each sample. Instead, use the provided make_peakcalling_jsons.py script, which can be found in workflow_peakcaller/scripts, to generate all .json files simultaneously.

* First create a new folder where these .jsons will be stored. 
* Also create a folder where all peakcalling outputs will be stored.
* Then, call the script using the following parameters:

    `
    python make_peakcalling_jsons.py <path_to_sailor_output_folder> <path_to_new_folder_where_jsons_will_be_placed> <path_to_new_folder_where_peakcalling_outputs_will_be_placed> <path_to_peakcalling_regions_folder> <path_to_fasta_used_for_alignment> <edit_type (i.e. CT or AG)> <fdr_threshold (default 0.1)>
    `

Once this completes, you will have one .json file for each sample successfully processed by SAILOR.

### Launch the snakemake job

Using the same config.yaml file you made previously to launch your sailor run, you can then launch the peak-calling pipeline similarly for each sample:
`
snakemake --profile /full/path/to/STAMP/profiles/my_profile/ --configfile /full/path/to/your/config/file/peakcalling_info_for_this_sample.json
`

This will create a new folder within the folder you specified at the "output_folder" parameter. Within that folder, a series of subfolders will be created as the job runs. The ultimate output you should expect will end up looking like this:

```
/absolute/path/to/peakcalling/folder/outputs_where_all_samples_will_be_placed/peakcalling_outputs/
    >bash_scripts
        > label_for_this_sample
            > label_for_this_sample_bash_job_edit_c_0.sh 
            > ...
            > label_for_this_sample_bash_job_edit_c_number_of_regions.sh 
    >bedgraphs
        > label_for_this_sample_editc.bedgraph
    >editc_outputs
        > label_for_this_sample
            > label_for_this_sample_edit_c_for_all_regions_0_all_info.filtered_regions.tsv
            > label_for_this_sample_edit_c_for_all_regions_0_all_info.tsv
            > ...
            > label_for_this_sample_edit_c_for_all_regions_number_of_regions_all_info.filtered_regions.tsv
            > label_for_this_sample_edit_c_for_all_regions_number_of_regions_0_all_info.tsv
    >peak_calling
        > label_for_this_sample
            > label_for_this_sample.fdr_0.1.d_15.edit_fraction.filtered_regions.tsv
            > label_for_this_sample.fdr_0.1.d_15.edit_fraction.tsv
            > label_for_this_sample.fdr_0.1.d_15.merged.bed
            > label_for_this_sample.fdr_0.1.d_15.regions
        label_for_this_sample.fdr_0.1.d_15.scored.tsv
```

The final scored peaks are found within the peak_calling folder, in this example *label_for_this_sample.fdr_0.1.d_15.scored.tsv*.

