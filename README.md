# Snakemake workflow: rna-seq-star-deseq2

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4737358.svg)](https://doi.org/10.5281/zenodo.4737358)
[![Snakemake](https://img.shields.io/badge/snakemake-≥6.1.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/snakemake-workflows/rna-seq-star-deseq2/workflows/Tests/badge.svg?branch=master)](https://github.com/snakemake-workflows/rna-seq-star-deseq2/actions?query=branch%3Amaster+workflow%3ATests)


<!-- TABLE OF CONTENTS -->
<details open="open">
  <summary><h2 style="display: inline-block">Table of Contents</h2></summary>
  <ol>
    <li>
      <a href="#about">About</a>
    </li>
    <li>
      <a href="#updates">Updates</a>
    </li>
    <li>
      <a href="#setting-up">Setting up</a>
    </li>
    <li>
      <a href="#workflow-dag">Workflow DAG</a>
    </li>
  </ol>
</details>




## About

This workflow performs a ExomeSeq analysis with bwa, Mutect2, Mutect, Strelka, Varscan, Sequenza.
It has been customized to work on the HPC4Health Slurm cluster and includes extra analysis to do genotype matching.

Link to the original snakemake workflow: [RNAseq-Exome-snakemake/Snakemake/exome/](https://github.com/pmcc-modil/pipelines/tree/main/RNAseq-Exome-snakemake/Snakemake/exome)

## Updates

* [TODO]: Add a script to initiate a blank working copy for building envs
* [TODO]: Implemented first working copy of genotype identification

## Setting up

As states above, this workflow is set up to work on the H4H cluster. As such, there are major limitations when working on this system. For instance, only the home directory contains internet access while typically, all analysis needs to be conducted on the project directories which can only be accessed on nodes with no internet access. Additionally, the home directory contains a very limited amount of storage space.

### 0. Install Snakemake as a conda env (Build node)
More detailed information can be found on the [snakemake docs](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html).
```
conda install -n base -c conda-forge mamba
mamba create -c conda-forge -c bioconda -n snakemake snakemake
```

### 1. Clone the workflow (Build node)
We first want to clone the workflow to our workflow directory in our home. This requires internet access, so be sure to log in to the build node first and activate your snakemake conda env:
```
salloc --partition=build -c 1 -t 2:0:0 --mem 2G
conda activate snakemake

cd ~/workflows
git clone git@github.com:elsamah/ExomeSeq.git
```

### 2. Setup a blank analysis directory (Build/Project node)
The basic idea of this workflow is that we set up all the conda-envs using a set output path in our home directory. However, to run the `--create-env-only` command, we need to have a directory set up with all the files needed to execute the entire workflow. In the future, I will add a script to do this automatically, but for now, it must be made manually:

```
wflowdir='~/workflows/ExomeSeq'
mkdir -p ~/workflows/intialize/ExomeSeq
cd ~/workflows/intialize/ExomeSeq

mkdir config data resources ref
ln -s ${wflowdir}/workflow/scripts .
ln -s ${wflowdir}/slurm .

cp ${wflowdir}/config/* config/
```
At this step, you will want to populate your `config/units.tsv` and `config/samples.tsv` with some junk nonsensical data, then you'll want to `touch` all the files into creation in the `data/` directory. For example, if you choose not to modify the `units.tsv` or `samples.tsv` (it should be fine), you can touch the following files in your data dir:

```
cd data
touch  A.1.fastq.gz A.2.fastq.gz A2.1.fastq.gz A2.2.fastq.gz B.1.fastq.gz B.2.fastq.gz
```

While not listed here, many of your reference files (e.g. `genoma.fa`) will be in your Project directories, which you cannot access from the build node.  Be sure to reconfigure the `config/config.yaml` to point to reference files that are in the home directory (`resources/` directory). As nothing will actually be run, you can simply `touch` these files for the sake of satisfying the snakemake file checker.

In `ref/` directory
```
touch 1000G_phase1.snps.high_confidence.hg38.vcf dbsnp_144.hg38.vcf Mills_and_1000G_gold_standard.indels.hg38.vcf BWAgenome.fa genome.fa
```

In `resources/` directory
```
touch genome.fa xenomeidx-both.kmers.high-bits xenomeidx-graft.kmers.low-bits.lwr INS-B-014-SB.processed.bam xenomeidx-both.kmers.low-bits.lwr xenomeidx-graft.kmers.low-bits.upr S04380110_Covered.headless.bed xenomeidx-both.kmers.low-bits.upr xenomeidx-host.header xengsortidx.hash xenomeidx-both.lhs-bits xenomeidx-host.kmers-d0 xengsortidx.info xenomeidx-both.rhs-bits xenomeidx-host.kmers-d1 xengsort.sif xenomeidx-graft.header xenomeidx-host.kmers.header xenomeidx-both.header xenomeidx-graft.kmers-d0 xenomeidx-host.kmers.high-bits xenomeidx-both.kmers-d0 xenomeidx-graft.kmers-d1 xenomeidx-host.kmers.low-bits.lwr xenomeidx-both.kmers-d1 xenomeidx-graft.kmers.header xenomeidx-host.kmers.low-bits.upr xenomeidx-both.kmers.header xenomeidx-graft.kmers.high-bits
```

In `data/` directory
```
touch 1428_PAR_S8_L004_R1_001.fastq.gz B.2.fastq.gz MCF7_PAR_S5_L001_R1_001.fastq.gz 1428_PAR_S8_L004_R2_001.fastq.gz CAMA1_PAR_S3_L001_R1_001.fastq.gz MCF7_PAR_S5_L001_R2_001.fastq.gz 1428_RES_S9_L004_R1_001.fastq.gz CAMA1_PAR_S3_L001_R2_001.fastq.gz MCF7_RES_S1_L003_R1_001.fastq.gz 1428_RES_S9_L004_R2_001.fastq.gz CAMA1_RES_S4_L001_R1_001.fastq.gz MCF7_RES_S1_L003_R2_001.fastq.gz 361_PAR_S1_L001_R1_001.fastq.gz CAMA1_RES_S4_L001_R2_001.fastq.gz T47D_PAR_S4_L003_R1_001.fastq.gz 361_PAR_S1_L001_R2_001.fastq.gz KPL1_PAR_S2_L003_R1_001.fastq.gz T47D_PAR_S4_L003_R2_001.fastq.gz 361_RES_S2_L001_R1_001.fastq.gz KPL1_PAR_S2_L003_R2_001.fastq.gz T47D_RES_S5_L004_R1_001.fastq.gz 361_RES_S2_L001_R2_001.fastq.gz KPL1_RES_S3_L003_R1_001.fastq.gz T47D_RES_S5_L004_R2_001.fastq.gz A.1.fastq.gz KPL1_RES_S3_L003_R2_001.fastq.gz ZR_PAR_S10_L004_R1_001.fastq.gz A2.1.fastq.gz LY2_PAR_S6_L004_R1_001.fastq.gz ZR_PAR_S10_L004_R2_001.fastq.gz A2.2.fastq.gz LY2_PAR_S6_L004_R2_001.fastq.gz ZR_RES_S11_L004_R1_001.fastq.gz A.2.fastq.gz LY2_RES_S7_L004_R1_001.fastq.gz ZR_RES_S11_L004_R2_001.fastq.gz B.1.fastq.gz LY2_RES_S7_L004_R2_001.fastq.gz
```

### 3. Setup your project directory (Project node)
Similary to Step 2, you will want to set up your project directory in your group directory using the same structure. However, this time you will need to add a meaningful `config/units.tsv`, `config/samples.tsv`, `config/config.yaml`, and populate your `data/` directory with symlinks to all the fastqs you will be working with.

### 4. Download the snakemake-wrappers (Build node)
Again, while in the build node, you can only access the home directory. Because of this, you will need to have a local copy of `snakemake-wrappers` in your home directory in order to build your `conda-envs`

This workflow requires: `v0.75.0` and `0.59.2`.

Visit the the git repo for [snakemake-wrappers](https://github.com/snakemake/snakemake-wrappers). Switch to the **tag** that you want to download (e.g. v0.75.0), then download the zip file (e.g. `wget 'https://github.com/snakemake/snakemake-wrappers/archive/refs/tags/v0.75.0.zip'`) to your snakemake-wrappers directory. Once unzipped, be sure to relabel it as `v0.75.0` or whatever tag you had it.

Basically, you will need to add a path that will have the following structure:
```
/path/to/snakemake-wrappers/[TAG]/bio/
```

### 5. Building the conda-envs (Build node)
Now that you have your "intialize" and your "project" directory set up, you can start building your conda environments required for the workflow.

Before you build your conda-envs, you'll need to switch the `workdir:` path in your `workflow/Snakefile`:
```
> workflow/Snakefile
workdir: "/path/to/workflows/intialize/ExomeSeq"
```

Once the Snakefile workdir path has been changed, you can run the entire snakemake workflow using `conda-create-envs-only`:

```
cd ~/workflows/ExomeSeq
condaprefix=$(readlink -f .)"/.snakemake/conda"

snakemake \
--use-conda \
--use-singularity \
--conda-create-envs-only \
--conda-frontend conda \
--conda-prefix ${condaprefix} \
--wrapper-prefix 'file:////cluster/home/selghamr/v0.77.0/' \
--cores 4
```
**Note:** the `--wrapper-prefix` is labelled as `file:///` with three `/`'s.  It WILL throw a warning/error during building/running the workflow, but it MUST be set this way or it will not work.

The idea behind this is that snakemake will install the conda envs to `.snakemake/conda`. It creates a hash to label that environment it builds. However, the hash is generated based on the `conda-prefix` and the hash of the `env.yaml` that it is building. As long as the `conda-prefix` and `env.yaml` remain unchanged, it will use the pre-existing environment found in `.snakemake/conda/[HASH]`

### 6. Run your workflow (Build/Project node)
Once the workflow has been set up and the environments have been created, you can finally run your pre-configured workflow on your project directory. You will need to reconfigure some paths in your `scheduler.sh` script. Be sure to activate your `snakemake` env before queuing the scheduler.

Before you run your workflow on your project directory, you'll need to switch the `workdir:` path in your `workflow/Snakefile`:
```
> workflow/Snakefile
workdir: "/path/to/group_directory/project"
```

The `scheduler.sh` script runs the following command in a 5-day long job, with its main purpose to track the job queue and job submission on the slurm cluster.
```
cd ~/workflows/ExomeSeq
condaprefix='/cluster/home/selghamr/workflows/ExomeSeq/.snakemake/conda'

snakemake \
--jobs 6 \
--profile slurm \
--cluster-config slurm/cluster.json \
--conda-frontend conda \
--use-conda \
--use-singularity \
--conda-prefix ${condaprefix} \
--wrapper-prefix 'file:////cluster/home/selghamr/v0.77.0/' \
--rerun-incomplete
```

You can submit the scheduler as followed:
```
conda activate snakemake
sbatch scheduler.sh
```

## Workflow DAG
