# MalariaGEN Parasite SNP & Indel Calling Pipeline
A NextFlow implementation of the MalariaGEN *Plasmodium falciparum* release SNP & indel calling pipeline. This pipeline was used to create Pf8's SNP & indel calls, and will be used to do the same for Pf9. 

This pipeline was modified once Pf8 was finished post hoc. For finding the version of the codebase that was used for generating Pf8, please refer to the commit history. 

## Contents
- [MalariaGEN Parasite SNP \& Indel Calling Pipeline](#malariagen-parasite-snp--indel-calling-pipeline)
  - [Contents](#contents)
  - [Pipeline Overview](#pipeline-overview)
      - [1. Setup](#1-setup)
      - [2. Generating BAMs and GVCFs workflow](#2-generating-bams-and-gvcfs-workflow)
      - [3. Joint genotyping workflow](#3-joint-genotyping-workflow)
- [Running the Pipeline](#running-the-pipeline)
    - [Step 1: Setting up the environment](#step-1-setting-up-the-environment)
    - [Step 2: Generate BAM and GVCF files](#step-2-generate-bam-and-gvcf-files)
    - [Step 3: Joint genotyping](#step-3-joint-genotyping)
  - [Output](#output)
  - [Options](#options)
    - [Inputs](#inputs)
    - [Outputs](#outputs)
    - [Additional](#additional)
    - [Clean up](#clean-up)
    - [Run modes](#run-modes)
    - [Notes](#notes)
  - [Creating a manifest](#creating-a-manifest)
  - [Dependencies](#dependencies)
- [To Do List to Finish Repo](#to-do-list-to-finish-repo)

<a name="summary"></a>
## Pipeline Overview

The `parasite-snp-indel-calling` pipeline is composed of two workflows that must be run sequentially: generating BAM & GVCF files, and joint genotyping

#### 1. Setup
1. Index reference FASTA
2. Create picard dictionary of reference FASTA

#### 2. Generating BAMs and GVCFs workflow
3. Pull CRAM files from iRODs
4. Remap CRAM files to reference
5. Mark duplicates
6. BQSR
7. Merge BAMs
8. Index BAMs
9. Get BAM stats
10. Generate GVCFs 

#### 3. Joint genotyping workflow
10. Haplotype calling
11. Import GVCFs into GenomicDB
12. Genotyping
13. Collect interval targets of all samples
14. Variant recalibration of SNPs
15. Variant recalibration of INDELs
16. Apply recalibration models to each sample
17. Merge recalibrated and genotyped VCFs for each sample
18. Variant function annotation for each sample
19. Variant filtration for each sample

---

<a name="users"></a>
# Running the Pipeline 

<a name="setup"></a>
### Step 1: Setting up the environment

Instructions in this section are specific to the Wellcome Sanger Institute's compute cluster, 'farm22'. 

**1. Login to the farm as user malpipe**
```
ssh malpipe@farm22
```
You may need to be given access to this user account. 

**2. Get iRODs authentication**
See 'Getting Connected' in [this guide](https://ssg-confluence.internal.sanger.ac.uk/display/FARM/iRODS#iRODS-GettingConnected) for further information (available only to Wellcome Sanger Institute employees).

Ensure you specify a long enough time with `iinit` to enable the pipeline to run. In this example we choose 120 hours.
```
module load ISG/IRODS/1.0
iinit --ttl 120
```
Enter the malpipe user iRODs password when prompted (available upon request).

**3. Go into an appropriate lustre directory and clone this repo**

**4. Ensure sample manifests are in place**

Due to the large number of files created in a pipeline run, the first stage (BAM and GVCF generation) must be run in batches to avoid exceeding file quotas. Batch sizes should be approximately 5,000 samples. A manifest  is required for each batch of samples. For guidance on creating a batch manifest see [Creating a manifest](#creating-a-manifest). 

---

<a name="pipeline"></a>
### Step 2: Generate BAM and GVCF files

Per sample BAM and GVCFs with their corresponding index files are generated in the first stage in batches of ~5000 samples. We generate GVCFs per genomic interval, per sample. Totalling 16 GVCFs per sample.

This stage must be run once for each batch of samples. 
Each batch must be run from within a dedicated batch directory.
**Batches cannot be run concurrently**.

After each batch, the `work/` directory can be deleted and another batch started.
See  [Clean up](#clean-up) for guidance on removing the work directory following each batch run.



```
# Navigate to pipeline stage 1 batch 1 run location
cd <PATH_TO_PIPELINE_RUN_DIRECTORY>

# Set up for a batch
mkdir bams_gvcfs
cd bams_gvcfs
mkdir run_<BATCH_NUMBER>

cd <PATH_TO_PIPELINE_RUN_DIRECTORY>/bams_gvcfs/run_<BATCH_NUMBER>
```

Create a shell script file (e.g., `run.sh`) within the batch directory, which will contain instructions for the pipeline run. Use the following code:

```
#!/bin/bash

# Set up proxies/parameters and load relevant modules

export http_proxy=http://wwwcache.sanger.ac.uk:3128
export https_proxy=http://wwwcache.sanger.ac.uk:3128
export NXF_OPTS='-Xms512M -Xmx2G'

module load nextflow-23.10.0
module load ISG/singularity/3.11.4
module load ISG/IRODS/1.0

export PIPELINE_CODEBASE=<PATH_TO_DESIRED_REPO_LOCATION>
export PIPELINE_MANIFESTS=<PATH_TO_DESIRED_REPO_LOCATION>/batched_manifests

BATCH_NUMBER=$1

bsub \
    -J snp-indel_${BATCH_NUMBER} \
    -o snp-indel_${BATCH_NUMBER}.o \
    -e snp-indel_${BATCH_NUMBER}.e \
    -n 1 \ 
    -M 8000 \ # memory allocation
    -R "select[mem>=8000] rusage[mem=8000] span[hosts=1]" \
    -q basement \
    "nextflow run ${PIPELINE_CODEBASE}/main.nf \
    --generate_bams \
    --generate_gvcfs \
    --manifest ${PIPELINE_MANIFESTS}/<CURRENT_BATCH_MANIFEST_FILE> \
    --results_dir 'output' \
    -resume \
    -with-trace -with-report -profile sanger_lsf"
```
- The unpadded intervals file specifies the genomic coordinates in predetermined chunk sizes for performing joint genotyping and is provided with this repository
- The resume flag enables the pipeline to be re-run from last saved point in case of failure to complete

Ensure permissions are added to the run.sh file:
```
chmod +x run.sh
```

Run the batch by supplying the batch number as the only argument:
```
# example for batch 1
./run.sh 1
```

Check on pipeline progress:
```
tail -n 50 <JOB_NAME>_${BATCH_NUMBER}.o
```

When a batch finishes:
1. Ensure BAMS and GVCFs in `--results_dir` are all present and formed as expected 
2. Delete batch work directory (see [Clean up](#clean-up))
3. Launch next batch (repeat steps in this section for batches 2, 3, etc.)

---

### Step 3: Joint genotyping

To enable the next stage of the pipeline, a GVCF map file (`.gvcf_map`) for each chromosome needs to be created. Each file takes the format: `<SAMPLE ID>` `<Path/to/chromosome/gvcf>`

Use the [`generate_gvcf_maps.sh`](https://github.com/malariagen/parasite-snp-indel-calling/blob/master/resources/generate_gvcf_maps.sh) file contained within `resources/` to generate these files.
The newly created `gvcf_maps/` directory can be use in the next stage.

This stage performs variant calling across all samples for each chromosome. The output of this stage is 16 filtered VCFs; one per region containing positions from all samples. 

Use the following to create a shell script:
```
#!/bin/bash

# Set up proxies/parameters and load relevant modules

export http_proxy=http://wwwcache.sanger.ac.uk:3128
export https_proxy=http://wwwcache.sanger.ac.uk:3128
export NXF_OPTS='-Xms512M -Xmx2G'

module load nextflow-23.10.0
module load ISG/singularity/3.11.4
module load ISG/IRODS/1.0

bsub \
    -J <JOB_NAME> \
    -o <JOB_NAME>.o \
    -e <JOB_NAME>.e \
    -n 1 \
    -M 8000 \
    -R "select[mem>=8000] rusage[mem=8000] span[hosts=1]" \
    -q basement \
    "nextflow run main.nf \
    --joint_genotyping \
    --gvcf_map_directory <INSERT_PATH_HERE>/gvcf_maps \
    --results_dir 'output' \
    -resume \
    -with-trace -with-report -profile sanger_lsf"
```
Execute the `.sh` file (no arguments needed) to start joint genotyping.

---

<a name="output"></a>
## Output

For each sample, the following files are produced and copied into the results directory specified using `--results_dir` after the generate BAMs and GVCFs stage:

File | Description
---- | -----------
\<sample\>.bam | Output merged bam file
\<sample\>.bam.bai | Output index of bam file
\<sample\>.bamstats | Output from `samtools stats` on (intermediate) bam file
\<sample\>.\<chromosome\>.vcf.gz| Output variant file for each chromosome
\<sample\>.\<chromosome\>.vcf.gz.tbi| Output index of variant file for each chromosome

For each chromosome, the following files are produced after joint genotyping:

File | Description
---- | -----------
\<chromosome\>.filt.vcf.gz | All sample filtered VCF file
\<chromosome\>.filt.vcf.gz.tbi | Index of all sample filtered VCF file



<a name="options"></a>
## Options
### Inputs
    Workflows:
    --generate_bams                   Run sample BAM generation 
    --generate_gvcfs                  Run sample GVCF generation
    --joint_genotyping                Run joint genotyping

    Inputs:
    --manifest <PATH>                 Path to sample manifest file (Must be used if --generate_bams is specified)
    --bam_files <PATH>                Path to sample bam files (Must be used if --generate_bams is not specified but --generate_gvcfs is specified)
    --gvcf_map_directory <PATH>       Path to gvcf map directory (Must be used if --joint_genotyping is specified)

### Outputs
    --results_dir                     Path to a directory for output files

### Additional
    -N <email address>               For email notifications when the pipeline completes/fails
    -resume                           If specified, pipeline can be restarted from last saved point
    -with-trace                       Enables the generation of a trace.txt file, containing detailed execution information for all processes in the pipeline. Useful for monitoring and debugging workflow performance. 
    -with-report                      Generates an HTML report that provides a high-level summary of the pipeline execution
    --help                            If true, print help information then quit (Default: false)

<a name="cleanup_for_users"></a>
### Clean up

The `work` directory keeps the intermediate files of the pipeline. For running stage 1, the work directory needs to be removed from a batch directory after successful completion of the batch. You can use `rm -rf work` but only when no other pipelines are running (more dangerous). The run names and session ids can be found by using `nextflow log -q`.

### Run modes
The default way to run this pipeline is to generate BAMs and GVCFs in one pipeline run (key flags: `--generate_bams`, `--generate_gvcfs`, `--manifest`), then perform joint genotyping in a second pipeline run (key flags: `--joint_genotyping`, `--gvcf_map_directory`). However, it is also possible to split the first run into two runs: firstly only generating BAMs (key flags: `--generate_bams`, `--manifest`), then generating GVCFs (key flags: `--generate_gvcfs`, `--bam_files`), in which case `bam_files` must be a path to a tab-separated file containing sample ID, path to .bam, and path to .bam.bai as its only three columns, with no header. 

<a name="note_for_users"></a>
### Notes
- In the event of failure, Nextflow provides a detailed description of which process(es) failed, along with a stack trace from the process itself. It is informative to visit the process subdirectory work and viewing the .command.log and .command.sh files. Once the problem has been rectified, restart the pipeline using the commands above, but add the `-resume` flag to the nextflow command.
- By default up to 1000 jobs are submitted at a time, but this can be changed by specifying `--queue_size` e.g., to run 200 jobs at a time, add `--queue_size 200`.
- Other parameters are specified in the nextflow.config. For example, in `bsqr.nf` within `modules`, there is a call to `BaseRecalibrator` that looks like:
	```bash
	${gatk} --java-options ${jvm_args} BaseRecalibrator -R ${reference_file} -I ${bam_file} -O ${gatk_recalibration_report} ${gatk_base_recalibrator_options}
	```
	Any additional option that can be used with GATK4's `BaseRecalibrator` (with the exception of options already used `-R`, `-I`, `-O`) can be specified using `--gatk_base_recalibrator_options`.
	
<a name="manifest"></a>
## Creating a manifest

To run the pipeline on a set of samples, you need to supply a manifest file as input. This manifest file can be generated by running the python script [`create_manifest.py`](https://github.com/malariagen/parasite-snp-indel-calling/blob/master/resources/manifest_creation/create_manifest.py), available in `resources/`. For details on how to successfully run this script please see the associated README.

<a name="dependencies"></a>
## Dependencies
The pipeline uses executables that are installed in a Singularity image specified by `process.container` in the nextflow.config file. Versions can be derived from running `singularity exec <image_path> <software executable> -v` on farm22.