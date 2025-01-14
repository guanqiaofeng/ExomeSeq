
```
mamba create -c conda-forge -c bioconda -n snakemake6153 snakemake=6.15.3

salloc --partition=build -c 1 -t 2:0:0 --mem 2G

mamba activate snakemake6153

cd ~
mkdir workflow
cd workflow

git clone git@github.com:elsamah/ExomeSeq.git

wflowdir='~/workflow/ExomeSeq'
mkdir -p ~/workflow/intialize/ExomeSeq
cd ~/workflow/intialize/ExomeSeq

mkdir config data resources ref
ln -s ${wflowdir}/workflow/scripts .
ln -s ${wflowdir}/slurm .

cp ~/workflow/ExomeSeq/config/* config/

```

```
cd data
touch  A.1.fastq.gz A.2.fastq.gz A2.1.fastq.gz A2.2.fastq.gz B.1.fastq.gz B.2.fastq.gz
touch 1428_PAR_S8_L004_R1_001.fastq.gz B.2.fastq.gz MCF7_PAR_S5_L001_R1_001.fastq.gz 1428_PAR_S8_L004_R2_001.fastq.gz CAMA1_PAR_S3_L001_R1_001.fastq.gz MCF7_PAR_S5_L001_R2_001.fastq.gz 1428_RES_S9_L004_R1_001.fastq.gz CAMA1_PAR_S3_L001_R2_001.fastq.gz MCF7_RES_S1_L003_R1_001.fastq.gz 1428_RES_S9_L004_R2_001.fastq.gz CAMA1_RES_S4_L001_R1_001.fastq.gz MCF7_RES_S1_L003_R2_001.fastq.gz 361_PAR_S1_L001_R1_001.fastq.gz CAMA1_RES_S4_L001_R2_001.fastq.gz T47D_PAR_S4_L003_R1_001.fastq.gz 361_PAR_S1_L001_R2_001.fastq.gz KPL1_PAR_S2_L003_R1_001.fastq.gz T47D_PAR_S4_L003_R2_001.fastq.gz 361_RES_S2_L001_R1_001.fastq.gz KPL1_PAR_S2_L003_R2_001.fastq.gz T47D_RES_S5_L004_R1_001.fastq.gz 361_RES_S2_L001_R2_001.fastq.gz KPL1_RES_S3_L003_R1_001.fastq.gz T47D_RES_S5_L004_R2_001.fastq.gz A.1.fastq.gz KPL1_RES_S3_L003_R2_001.fastq.gz ZR_PAR_S10_L004_R1_001.fastq.gz A2.1.fastq.gz LY2_PAR_S6_L004_R1_001.fastq.gz ZR_PAR_S10_L004_R2_001.fastq.gz A2.2.fastq.gz LY2_PAR_S6_L004_R2_001.fastq.gz ZR_RES_S11_L004_R1_001.fastq.gz A.2.fastq.gz LY2_RES_S7_L004_R1_001.fastq.gz ZR_RES_S11_L004_R2_001.fastq.gz B.1.fastq.gz LY2_RES_S7_L004_R2_001.fastq.gz
```

```
cd ../ref
touch 1000G_phase1.snps.high_confidence.hg38.vcf dbsnp_144.hg38.vcf Mills_and_1000G_gold_standard.indels.hg38.vcf BWAgenome.fa genome.fa
```

```
cd ~/workflows/ExomeSeq/workflow/rules
for file in *; do
    if [ -f "$file" ]; then
        sed -i 's/selghamr/t135250uhn/g' "$file"
    fi
done
```
```
cd ../resources
touch genome.fa xenomeidx-both.kmers.high-bits xenomeidx-graft.kmers.low-bits.lwr INS-B-014-SB.processed.bam xenomeidx-both.kmers.low-bits.lwr xenomeidx-graft.kmers.low-bits.upr S04380110_Covered.headless.bed xenomeidx-both.kmers.low-bits.upr xenomeidx-host.header xengsortidx.hash xenomeidx-both.lhs-bits xenomeidx-host.kmers-d0 xengsortidx.info xenomeidx-both.rhs-bits xenomeidx-host.kmers-d1 xengsort.sif xenomeidx-graft.header xenomeidx-host.kmers.header xenomeidx-both.header xenomeidx-graft.kmers-d0 xenomeidx-host.kmers.high-bits xenomeidx-both.kmers-d0 xenomeidx-graft.kmers-d1 xenomeidx-host.kmers.low-bits.lwr xenomeidx-both.kmers-d1 xenomeidx-graft.kmers.header xenomeidx-host.kmers.low-bits.upr xenomeidx-both.kmers.header xenomeidx-graft.kmers.high-bits
```



Error1:
```
more slurm-14530035.out 
Building DAG of jobs...
Using shell: /usr/bin/bash
Provided cluster nodes: 12
Traceback (most recent call last):
  File "/cluster/home/t135250uhn/miniconda3/envs/snakemake6153/lib/python3.10/site-packages/snakemake/__init__.py", line 701, in snakemake
    success = workflow.execute(
  File "/cluster/home/t135250uhn/miniconda3/envs/snakemake6153/lib/python3.10/site-packages/snakemake/workflow.py", line 1060, in execute
    logger.run_info("\n".join(dag.stats()))
  File "/cluster/home/t135250uhn/miniconda3/envs/snakemake6153/lib/python3.10/site-packages/snakemake/dag.py", line 2191, in stats
    yield tabulate(rows, headers="keys")
  File "/cluster/home/t135250uhn/miniconda3/envs/snakemake6153/lib/python3.10/site-packages/tabulate/__init__.py", line 2048, in tabulate
    list_of_lists, headers = _normalize_tabular_data(
  File "/cluster/home/t135250uhn/miniconda3/envs/snakemake6153/lib/python3.10/site-packages/tabulate/__init__.py", line 1471, in _normaliz
e_tabular_data
    rows = list(map(lambda r: r if _is_separating_line(r) else list(r), rows))
  File "/cluster/home/t135250uhn/miniconda3/envs/snakemake6153/lib/python3.10/site-packages/tabulate/__init__.py", line 1471, in <lambda>
    rows = list(map(lambda r: r if _is_separating_line(r) else list(r), rows))
  File "/cluster/home/t135250uhn/miniconda3/envs/snakemake6153/lib/python3.10/site-packages/tabulate/__init__.py", line 107, in _is_separa
ting_line
    (len(row) >= 1 and row[0] == SEPARATING_LINE)
  File "/cluster/home/t135250uhn/miniconda3/envs/snakemake6153/lib/python3.10/site-packages/snakemake/rules.py", line 1138, in __eq__
    return self.name == other.name and self.output == other.output
AttributeError: 'str' object has no attribute 'name'
```
https://github.com/snakemake/snakemake/issues/1899
Solved by 
```
mamba install 'tabulate=0.8.10'
```
Error2
```
$ more deconvolutexengsort_14530333.err
Building DAG of jobs...
Using shell: /usr/bin/bash
Provided cores: 4
Rules claiming more threads will be scaled down.
Select jobs to execute...

[Mon Jan  6 14:20:54 2025]
rule deconvolutexengsort:
    input: /cluster/projects/cesconlab/workspace/Guanqiao/workflow/WES_2024Dec/data/BPTO93T_5_S5_R1_001.fastq.gz, /cluster/projects/cesconlab/worksp
ace/Guanqiao/workflow/WES_2024Dec/data/BPTO93T_5_S5_R2_001.fastq.gz
    output: results/xengsort/BPTO93T-graft.1.fq.gz, results/xengsort/BPTO93T-graft.2.fq.gz, results/xengsort/BPTO93T-neither.1.fq, results/xengsort/
BPTO93T-neither.2.fq, results/xengsort/BPTO93T-both.1.fq, results/xengsort/BPTO93T-both.2.fq, results/xengsort/BPTO93T-ambiguous.1.fq, results/xengs
ort/BPTO93T-ambiguous.2.fq, results/xengsort/BPTO93T-host.1.fq, results/xengsort/BPTO93T-host.2.fq
    jobid: 0
    wildcards: sample=BPTO93T
    threads: 4
    resources: mem_mb=7124, disk_mb=7124, tmpdir=/tmp

mkdir: cannot create directory ‘tmp’: File exists
[Mon Jan  6 14:20:57 2025]
Error in rule deconvolutexengsort:
    jobid: 0
    output: results/xengsort/BPTO93T-graft.1.fq.gz, results/xengsort/BPTO93T-graft.2.fq.gz, results/xengsort/BPTO93T-neither.1.fq, results/xengsort/
BPTO93T-neither.2.fq, results/xengsort/BPTO93T-both.1.fq, results/xengsort/BPTO93T-both.2.fq, results/xengsort/BPTO93T-ambiguous.1.fq, results/xengs
ort/BPTO93T-ambiguous.2.fq, results/xengsort/BPTO93T-host.1.fq, results/xengsort/BPTO93T-host.2.fq
    shell:
        
    module load apptainer/1.0.2
    module load pigz/2.6 
    
    mkdir -p results/xengsort
    mkdir tmp
    zcat /cluster/projects/cesconlab/workspace/Guanqiao/workflow/WES_2024Dec/data/BPTO93T_5_S5_R1_001.fastq.gz > tmp/BPTO93T_R1.fastq
    zcat /cluster/projects/cesconlab/workspace/Guanqiao/workflow/WES_2024Dec/data/BPTO93T_5_S5_R2_001.fastq.gz > tmp/BPTO93T_R2.fastq
    
    apptainer run /cluster/projects/cesconlab/envs/containers/xengsort/xengsort.sif     xengsort classify     --index /cluster/projects/cesconlab/Re
ferences/xengsort/idx_grcm38_grch38/xengsortidx     --fastq tmp/BPTO93T_R1.fastq     --pairs tmp/BPTO93T_R2.fastq     --prefix results/xengsort/BPTO
93T     --compression none     -T 4     --progress     --filter
    
    rm tmp/BPTO93T_R1.fastq
    rm tmp/BPTO93T_R2.fastq
    pigz -4 results/xengsort/BPTO93T-graft.1.fq.gz
    pigz -4 results/xengsort/BPTO93T-graft.2.fq.gz
    
        (one of the commands exited with non-zero exit code; note that snakemake uses bash strict mode!)

Shutting down, this might take some time.
Exiting because a job execution failed. Look above for error message
```
Fix:
/cluster/home/t135250uhn/workflow/ExomeSeq/workflow/rules/align.smk
```
rule deconvolutexengsort:
  input:
    f1 =  get_r1,
    f2 =  get_r2
  output:
    graftf1 =  "results/xengsort/{sample}-graft.1.fq.gz",
    graftf2 =  "results/xengsort/{sample}-graft.2.fq.gz",
    neitherf1 =  temp("results/xengsort/{sample}-neither.1.fq"),
    neitherf2 =  temp("results/xengsort/{sample}-neither.2.fq"),
    bothf1 =  temp("results/xengsort/{sample}-both.1.fq"),
    bothf2 =  temp("results/xengsort/{sample}-both.2.fq"),
    ambiguousf1 =  temp("results/xengsort/{sample}-ambiguous.1.fq"),
    ambiguousf2 =  temp("results/xengsort/{sample}-ambiguous.2.fq"),
    hostf1 =  temp("results/xengsort/{sample}-host.1.fq"),
    hostf2 =  temp("results/xengsort/{sample}-host.2.fq")
  params:
    xengsortidx=config["ref"]["xengsortidx"],
    xengsortcontainer=config['env']['xengsort'],
    sampleid="{sample}"
  threads: 4
  shell:
    """
    module load apptainer/1.0.2
    module load pigz/2.6

    mkdir -p results/xengsort
    mkdir tmp
    zcat {input.f1} > tmp/{params.sampleid}_R1.fastq
    zcat {input.f2} > tmp/{params.sampleid}_R2.fastq

    apptainer run {params.xengsortcontainer} \
    xengsort classify \
    --index {params.xengsortidx} \
    --fastq tmp/{params.sampleid}_R1.fastq \
    --pairs tmp/{params.sampleid}_R2.fastq \
    --prefix results/xengsort/{params.sampleid} \
    --compression none \
    -T {threads} \
    --progress \
    --filter

    rm tmp/{params.sampleid}_R1.fastq
    rm tmp/{params.sampleid}_R2.fastq
    pigz -{threads} {output.graftf1}
    pigz -{threads} {output.graftf2}
    """
```
`mkdir tmp` --> `mkdir -p tmp`

Error3:
```
$ more deconvolutexengsort_14530381.err 
Building DAG of jobs...
Using shell: /usr/bin/bash
Provided cores: 4
Rules claiming more threads will be scaled down.
Select jobs to execute...

[Mon Jan  6 14:27:25 2025]
rule deconvolutexengsort:
    input: /cluster/projects/cesconlab/workspace/Guanqiao/workflow/WES_2024Dec/data/DCBTO14N_14_S14_R1_001.fastq.gz, /cluster/projects/ces
conlab/workspace/Guanqiao/workflow/WES_2024Dec/data/DCBTO14N_14_S14_R2_001.fastq.gz
    output: results/xengsort/DCBTO14N-graft.1.fq.gz, results/xengsort/DCBTO14N-graft.2.fq.gz, results/xengsort/DCBTO14N-neither.1.fq, resu
lts/xengsort/DCBTO14N-neither.2.fq, results/xengsort/DCBTO14N-both.1.fq, results/xengsort/DCBTO14N-both.2.fq, results/xengsort/DCBTO14N-am
biguous.1.fq, results/xengsort/DCBTO14N-ambiguous.2.fq, results/xengsort/DCBTO14N-host.1.fq, results/xengsort/DCBTO14N-host.2.fq
    jobid: 0
    wildcards: sample=DCBTO14N
    threads: 4
    resources: mem_mb=7550, disk_mb=7550, tmpdir=/tmp

pigz: skipping: results/xengsort/DCBTO14N-graft.1.fq.gz does not exist
[Mon Jan  6 14:53:02 2025]
Error in rule deconvolutexengsort:
    jobid: 0
    output: results/xengsort/DCBTO14N-graft.1.fq.gz, results/xengsort/DCBTO14N-graft.2.fq.gz, results/xengsort/DCBTO14N-neither.1.fq, resu
lts/xengsort/DCBTO14N-neither.2.fq, results/xengsort/DCBTO14N-both.1.fq, results/xengsort/DCBTO14N-both.2.fq, results/xengsort/DCBTO14N-am
biguous.1.fq, results/xengsort/DCBTO14N-ambiguous.2.fq, results/xengsort/DCBTO14N-host.1.fq, results/xengsort/DCBTO14N-host.2.fq
    shell:
        
    module load apptainer/1.0.2
    module load pigz/2.6 
    
    mkdir -p results/xengsort
    mkdir -p tmp
    zcat /cluster/projects/cesconlab/workspace/Guanqiao/workflow/WES_2024Dec/data/DCBTO14N_14_S14_R1_001.fastq.gz > tmp/DCBTO14N_R1.fastq
    zcat /cluster/projects/cesconlab/workspace/Guanqiao/workflow/WES_2024Dec/data/DCBTO14N_14_S14_R2_001.fastq.gz > tmp/DCBTO14N_R2.fastq
    
    apptainer run /cluster/projects/cesconlab/envs/containers/xengsort/xengsort.sif     xengsort classify     --index /cluster/projects/ce
sconlab/References/xengsort/idx_grcm38_grch38/xengsortidx     --fastq tmp/DCBTO14N_R1.fastq     --pairs tmp/DCBTO14N_R2.fastq     --prefix
 results/xengsort/DCBTO14N     --compression none     -T 4     --progress     --filter
    
    rm tmp/DCBTO14N_R1.fastq
    rm tmp/DCBTO14N_R2.fastq
    pigz -4 results/xengsort/DCBTO14N-graft.1.fq.gz
    pigz -4 results/xengsort/DCBTO14N-graft.2.fq.gz
    
        (one of the commands exited with non-zero exit code; note that snakemake uses bash strict mode!)

Shutting down, this might take some time.
Exiting because a job execution failed. Look above for error message
```
To fix in "rule deconvolutexengsort:"
```
pigz -4 results/xengsort/DCBTO14N-graft.1.fq.gz
pigz -4 results/xengsort/DCBTO14N-graft.2.fq.gz
```
change to 
```
pigz -4 results/xengsort/DCBTO14N-graft.1.fq
pigz -4 results/xengsort/DCBTO14N-graft.2.fq
```
Error 4
```
$ more deconvolutexengsort_14531730.err
Building DAG of jobs...
Using shell: /usr/bin/bash
Provided cores: 4
Rules claiming more threads will be scaled down.
Select jobs to execute...

[Tue Jan  7 10:37:36 2025]
rule deconvolutexengsort:
    input: /cluster/projects/cesconlab/workspace/Guanqiao/workflow/WES_2024Dec/data/BPTO51_1_S1_R1_001.fastq.gz, /cluster/projects/cesconlab/
workspace/Guanqiao/workflow/WES_2024Dec/data/BPTO51_1_S1_R2_001.fastq.gz
    output: results/xengsort/BPTO51-graft.1.fq.gz, results/xengsort/BPTO51-graft.2.fq.gz, results/xengsort/BPTO51-neither.1.fq, results/xengs
ort/BPTO51-neither.2.fq, results/xengsort/BPTO51-both.1.fq, results/xengsort/BPTO51-both.2.fq, results/xengsort/BPTO51-ambiguous.1.fq, result
s/xengsort/BPTO51-ambiguous.2.fq, results/xengsort/BPTO51-host.1.fq, results/xengsort/BPTO51-host.2.fq
    jobid: 0
    wildcards: sample=BPTO51
    threads: 4
    resources: mem_mb=7888, disk_mb=7888, tmpdir=/tmp

Waiting at most 60 seconds for missing files.
MissingOutputException in line 1 of /cluster/home/t135250uhn/workflow/ExomeSeq/workflow/rules/align.smk:
Job Missing files after 60 seconds:
results/xengsort/BPTO51-neither.1.fq
results/xengsort/BPTO51-neither.2.fq
results/xengsort/BPTO51-both.1.fq
results/xengsort/BPTO51-both.2.fq
results/xengsort/BPTO51-ambiguous.1.fq
results/xengsort/BPTO51-ambiguous.2.fq
results/xengsort/BPTO51-host.1.fq
results/xengsort/BPTO51-host.2.fq
This might be due to filesystem latency. If that is the case, consider to increase the wait time with --latency-wait.
Job id: 0 completed successfully, but some output files are missing. 0
Removing output files of failed job deconvolutexengsort since they might be corrupted:
results/xengsort/BPTO51-graft.1.fq.gz, results/xengsort/BPTO51-graft.2.fq.gz
Shutting down, this might take some time.
Exiting because a job execution failed. Look above for error message
```
To fix in "rule deconvolutexengsort:"
```
    apptainer run {params.xengsortcontainer} \
    xengsort classify \
    --index {params.xengsortidx} \
    --fastq tmp/{params.sampleid}_R1.fastq \
    --pairs tmp/{params.sampleid}_R2.fastq \
    --prefix results/xengsort/{params.sampleid} \
    --compression none \
    -T {threads} \
    --progress \
    --filter
```
remove `--filter`

Error 5
```
[Tue Jan  7 14:20:41 2025]
rule picardMarkDuplicates:
    input: results/alignment/BPNO51/BPNO51_sorted.bam, results/alignment/BPNO51/BPNO51_sorted.bam.bai
    output: results/alignment/BPNO51/BPNO51_sorted.dedup.bam, results/alignment/BPNO51/BPNO51_picardmetrics.txt
    jobid: 0
    wildcards: sample=BPNO51
    threads: 4
    resources: mem_mb=1822, disk_mb=1822, tmpdir=/tmp

Activating conda environment: /cluster/home/t135250uhn/workflow/ExomeSeq/.snakemake/conda/1b527edd50c39f616d401330675ce543
Error: Unable to access jarfile /cluster/home/t135250uhn/workflows/ExomeSeq/.snakemake/conda/9b770440ff173434e53ee101c7452a0a/share/picard-2.
26.0-0/picard.jar
[Tue Jan  7 14:20:47 2025]
Error in rule picardMarkDuplicates:
    jobid: 0
    output: results/alignment/BPNO51/BPNO51_sorted.dedup.bam, results/alignment/BPNO51/BPNO51_picardmetrics.txt
    conda-env: /cluster/home/t135250uhn/workflow/ExomeSeq/.snakemake/conda/1b527edd50c39f616d401330675ce543
    shell:
        
    java -Xmx12g -jar /cluster/home/t135250uhn/workflows/ExomeSeq/.snakemake/conda/9b770440ff173434e53ee101c7452a0a/share/picard-2.26.0-0/pic
ard.jar MarkDuplicates INPUT=results/alignment/BPNO51/BPNO51_sorted.bam OUTPUT=results/alignment/BPNO51/BPNO51_sorted.dedup.bam METRICS_FILE=
results/alignment/BPNO51/BPNO51_picardmetrics.txt ASSUME_SORTED=true MAX_RECORDS_IN_RAM=100000 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true
 USE_JDK_DEFLATER=true USE_JDK_INFLATER=true
    
        (one of the commands exited with non-zero exit code; note that snakemake uses bash strict mode!)

Shutting down, this might take some time.
Exiting because a job execution failed. Look above for error message
```
To fix in rule picardMarkDuplicates of `rules/align.smk`
```
rule picardMarkDuplicates:
  input:
    bam="results/alignment/{sample}/{sample}_sorted.bam",
    bai="results/alignment/{sample}/{sample}_sorted.bam.bai",
  output:
    dedup="results/alignment/{sample}/{sample}_sorted.dedup.bam",
    metrics="results/alignment/{sample}/{sample}_picardmetrics.txt"
  params:
    picard="/cluster/home/t135250uhn/workflows/ExomeSeq/.snakemake/conda/9b770440ff173434e53ee101c7452a0a/share/picard-2.26.0-0"
  threads: 4
  conda:
    "/cluster/home/t135250uhn/workflows/ExomeSeq/workflow/envs/bwa.yaml",
  shell:
    """
    java -Xmx12g -jar {params.picard}/picard.jar MarkDuplicates INPUT={input.bam} OUTPUT={output.dedup} METRICS_FILE={output.metrics} ASSUME_SORTED=
true MAX_RECORDS_IN_RAM=100000 VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true USE_JDK_DEFLATER=true USE_JDK_INFLATER=true
    """
```
update `picard` by find its correct path using:
```
conda activate /cluster/home/t135250uhn/workflow/ExomeSeq/.snakemake/conda/1b527edd50c39f616d401330675ce543
conda list | grep picard
find /cluster/home/t135250uhn/workflow/ExomeSeq/.snakemake/conda/1b527edd50c39f616d401330675ce543 -name "picard.jar"

/cluster/home/t135250uhn/workflow/ExomeSeq/.snakemake/conda/1b527edd50c39f616d401330675ce543/share/picard-2.26.0-0/picard.jar
```
now update to 
```
picard="/cluster/home/t135250uhn/workflow/ExomeSeq/.snakemake/conda/1b527edd50c39f616d401330675ce543/share/picard-2.26.0-0"
```
Error 6
```
$ more Strelka_14533257.err
Building DAG of jobs...
Using shell: /usr/bin/bash
Provided cores: 4
Rules claiming more threads will be scaled down.
Select jobs to execute...

[Tue Jan  7 17:18:13 2025]
rule Strelka:
    input: results/alignment/BPNO51/BPNO51.realigned.recal.bam, ref/genome.fa, config/strelka_config_bwa.ini, resources/INS-B-014-SB.processe
d.bam
    output: results/Strelka/BPNO51/BPNO51.myAnalysis, results/Strelka/BPNO51/BPNO51.myAnalysis/results/variants/BPNO51_Slk_somatic.indels.vcf
.gz, results/Strelka/BPNO51/BPNO51.myAnalysis/results/variants/BPNO51_Slk_somatic.snvs.vcf.gz, results/Strelka/BPNO51/BPNO51.myAnalysis/resul
ts/variants/BPNO51_Slk_somatic.indels.vcf.gz.tbi, results/Strelka/BPNO51/BPNO51.myAnalysis/results/variants/BPNO51_Slk_somatic.snvs.vcf.gz.tb
i
    jobid: 0
    wildcards: sample=BPNO51
    threads: 4
    resources: mem_mb=29222, disk_mb=29222, tmpdir=/tmp

Activating conda environment: /cluster/home/t135250uhn/workflow/ExomeSeq/.snakemake/conda/ca496e3542c49ad4d5d8124e6f30f37f
/usr/bin/bash: line 2: /cluster/home/t135250uhn/workflows/ExomeSeq/.snakemake/conda/236aa367b1347b9561439ce4facd36c0/share/strelka-2.9.10-1/b
in/configureStrelkaSomaticWorkflow.py: No such file or directory
[Tue Jan  7 17:18:16 2025]
Error in rule Strelka:
    jobid: 0
    output: results/Strelka/BPNO51/BPNO51.myAnalysis, results/Strelka/BPNO51/BPNO51.myAnalysis/results/variants/BPNO51_Slk_somatic.indels.vcf
.gz, results/Strelka/BPNO51/BPNO51.myAnalysis/results/variants/BPNO51_Slk_somatic.snvs.vcf.gz, results/Strelka/BPNO51/BPNO51.myAnalysis/resul
ts/variants/BPNO51_Slk_somatic.indels.vcf.gz.tbi, results/Strelka/BPNO51/BPNO51.myAnalysis/results/variants/BPNO51_Slk_somatic.snvs.vcf.gz.tb
i
    conda-env: /cluster/home/t135250uhn/workflow/ExomeSeq/.snakemake/conda/ca496e3542c49ad4d5d8124e6f30f37f
    shell:
        
    ## configuration
    /cluster/home/t135250uhn/workflows/ExomeSeq/.snakemake/conda/236aa367b1347b9561439ce4facd36c0/share/strelka-2.9.10-1/bin/configureStrelka
SomaticWorkflow.py     --normal=resources/INS-B-014-SB.processed.bam      --tumor=results/alignment/BPNO51/BPNO51.realigned.recal.bam     --r
ef=ref/genome.fa     --config=config/strelka_config_bwa.ini     --runDir=results/Strelka/BPNO51/BPNO51.myAnalysis

    ## running pipeline
    results/Strelka/BPNO51/BPNO51.myAnalysis/runWorkflow.py -m local -j 20

    echo "Current Dir: "$(pdir)
    mv results/Strelka/BPNO51/BPNO51.myAnalysis/results/variants/somatic.indels.vcf.gz results/Strelka/BPNO51/BPNO51.myAnalysis/results/varia
nts/BPNO51_Slk_somatic.indels.vcf.gz
    mv results/Strelka/BPNO51/BPNO51.myAnalysis/results/variants/somatic.snvs.vcf.gz results/Strelka/BPNO51/BPNO51.myAnalysis/results/variant
s/BPNO51_Slk_somatic.snvs.vcf.gz
    mv results/Strelka/BPNO51/BPNO51.myAnalysis/results/variants/somatic.indels.vcf.gz.tbi results/Strelka/BPNO51/BPNO51.myAnalysis/results/v
ariants/BPNO51_Slk_somatic.indels.vcf.gz.tbi
    mv results/Strelka/BPNO51/BPNO51.myAnalysis/results/variants/somatic.snvs.vcf.gz.tbi results/Strelka/BPNO51/BPNO51.myAnalysis/results/var
iants/BPNO51_Slk_somatic.snvs.vcf.gz.tbi

    
        (one of the commands exited with non-zero exit code; note that snakemake uses bash strict mode!)

Removing output files of failed job Strelka since they might be corrupted:
results/Strelka/BPNO51/BPNO51.myAnalysis
Shutting down, this might take some time.
Exiting because a job execution failed. Look above for error message
```
Solve it:
```
(base) [t135250uhn@node123 ExomeSeq]$ conda activate /cluster/home/t135250uhn/workflow/ExomeSeq/.snakemake/conda/ca496e3542c49ad4d5d8124e6f30f37f
(/cluster/home/t135250uhn/workflow/ExomeSeq/.snakemake/conda/ca496e3542c49ad4d5d8124e6f30f37f) [t135250uhn@node123 ExomeSeq]$ conda list | grep strelka
strelka                   2.9.10               hdfd78af_2    bioconda
(/cluster/home/t135250uhn/workflow/ExomeSeq/.snakemake/conda/ca496e3542c49ad4d5d8124e6f30f37f) [t135250uhn@node123 ExomeSeq]$ find /cluster/home/t135250uhn/workflow/ExomeSeq/.snakemake/conda/ca496e3542c49ad4d5d8124e6f30f37f -name "configureStrelkaSomaticWorkflow.py"
/cluster/home/t135250uhn/workflow/ExomeSeq/.snakemake/conda/ca496e3542c49ad4d5d8124e6f30f37f/bin/configureStrelkaSomaticWorkflow.py
/cluster/home/t135250uhn/workflow/ExomeSeq/.snakemake/conda/ca496e3542c49ad4d5d8124e6f30f37f/share/strelka-2.9.10-2/bin/configureStrelkaSomaticWorkflow.py
```
Now update rule Strelka:
```
params:
    strelka="/cluster/home/t135250uhn/workflows/ExomeSeq/.snakemake/conda/236aa367b1347b9561439ce4facd36c0/share/strelka-2.9.10-1/bin",
```
to 
```
params:
    strelka="/cluster/home/t135250uhn/workflow/ExomeSeq/.snakemake/conda/ca496e3542c49ad4d5d8124e6f30f37f/share/strelka-2.9.10-2/bin",
```
Error 7
same error as Error 5/6, it is for varscan.smk. Correct it as same as above
also update conda with the correct path

Error 8
```
[2025-01-08T15:31:29.755024Z] [node77.h4h.uhnresearch.ca] [28080_1] [TaskManager] [ERROR] 
Failed to complete command task: 'CallGenome+callGenomeSegment_chromId_000_chr1_0008' launched from sub-workflow 'CallGenome', error code: 1, 
command: '
/cluster/home/t135250uhn/workflow/ExomeSeq/.snakemake/conda/ca496e3542c49ad4d5d8124e6f30f37f/share/strelka-2.9.10-2/libexec/strelka2 
--region chr1:94840545-106695612 
--ref /cluster/projects/cesconlab/workspace/Guanqiao/workflow/WES_2024Dec/ref/genome.fa 
--max-indel-size 49 
--min-mapping-quality 20 
--somatic-snv-rate 0.0001 
--shared-site-error-rate 0.0000000005
--shared-site-error-strand-bias-fraction 0.0 
--somatic-indel-rate 0.000001 
--shared-indel-error-factor 2.2 
--tier2-min-mapping-quality 0 
--strelka-snv
-max-filtered-basecall-frac 0.4 
--strelka-snv-max-spanning-deletion-frac 0.75 
--strelka-snv-min-qss-ref 15 
--strelka-indel-max-window-filtered-basecall
-frac 0.3 
--strelka-indel-min-qsi-ref 40 
--ssnv-contam-tolerance 0.15 
--indel-contam-tolerance 0.15 
--somatic-snv-scoring-model-file /cluster/home/t135250uhn/workflow/ExomeSeq/.snakemake/conda/ca496e3542c49ad4d5d8124e6f30f37f/share/strelka-2.9.10-2/share/config/somaticSNVScoringModels.json 
--somatic-indel-scoring-model-file /cluster/home/t135250uhn/workflow/ExomeSeq/.snakemake/conda/ca496e3542c49ad4d5d8124e6f30f37f/share/strelka-2.9.10-2/share/config/somaticIndelScoringModels.json 
--normal-align-file /cluster/projects/cesconlab/workspace/Guanqiao/workflow/WES_2024Dec/resources/INS-B-014-SB.processed.bam 
--tumor-align-file /cluster/projects/cesconlab/workspace/Guanqiao/workflow/WES_2024Dec/results/alignment/BPNO51/BPNO51.realigned.recal.bam 
--somatic-snv-file /cluster/projects/cesconlab/workspace/Guanqiao/workflow/WES_2024Dec/results/Strelka/BPNO51/BPNO51.myAnalysis/workspace/genomeSegment.tmpdir/somatic.snvs.unfiltered.chromId_000_chr1_0008.vcf 
--somatic-indel-file /cluster/projects/cesconlab/workspace/Guanqiao/workflow/WES_2024Dec/results/Strelka/BPNO51/BPNO51.myAnalysis/workspace/genomeSegment.tmpdir/somatic.indels.unfiltered.chromId_000_chr1_0008.vcf 
--stats-file /cluster/projects/cesconlab/workspace/Guanqiao/workflow/WES_2024Dec/results/Strelka/BPNO51/BPNO51.myAnalysis/workspace/genomeSegment.tmpdir/runStats.chromId_000_chr1_0008.xml
--strelka-skip-header 
--strelka-chrom-depth-file /cluster/projects/cesconlab/workspace/Guanqiao/workflow/WES_2024Dec/results/Strelka/BPNO51/BPNO51.myAnalysis/workspace/chromDepth.tsv 
--strelka-max-depth-factor 3.0
'
[2025-01-08T15:31:29.776880Z] [node77.h4h.uhnresearch.ca] [28080_1] [TaskManager] [ERROR] [CallGenome+callGenomeSegment_chromId_000_chr1_0008] Error Message:
[2025-01-08T15:31:29.795050Z] [node77.h4h.uhnresearch.ca] [28080_1] [TaskManager] [ERROR] [CallGenome+callGenomeSegment_chromId_000_chr1_0008] Anomalous task wrapper stderr output. Wrapper signal file: '/cluster/projects/cesconlab/workspace/Guanqiao/workflow/WES_2024Dec/results/Strelka/BPNO51/BPNO51.myAnalysis/workspace/pyflow.data/logs/tmp/taskWrapperLogs/000/107/pyflowTaskWrapper.signal.txt'
[2025-01-08T15:31:29.806699Z] [node77.h4h.uhnresearch.ca] [28080_1] [TaskManager] [ERROR] [CallGenome+callGenomeSegment_chromId_000_chr1_0008] Logging 5 line(s) of task wrapper log output below:
[2025-01-08T15:31:29.824197Z] [node77.h4h.uhnresearch.ca] [28080_1] [TaskManager] [ERROR] [CallGenome+callGenomeSegment_chromId_000_chr1_0008] [taskWrapper-stderr] [2025-01-08T15:31:29.396861Z] [node77.h4h.uhnresearch.ca] [28080_1] [pyflowTaskWrapper:CallGenome+callGenomeSegment_chromId_000_chr1_0008] [wrapperSignal] wrapperStart
[2025-01-08T15:31:29.839917Z] [node77.h4h.uhnresearch.ca] [28080_1] [TaskManager] [ERROR] [CallGenome+callGenomeSegment_chromId_000_chr1_0008] [taskWrapper-stderr] [2025-01-08T15:31:29.451702Z] [node77.h4h.uhnresearch.ca] [28080_1] [pyflowTaskWrapper:CallGenome+callGenomeSegment_chromId_000_chr1_0008] [wrapperSignal] taskStart
[2025-01-08T15:31:29.844013Z] [node77.h4h.uhnresearch.ca] [28080_1] [TaskManager] [ERROR] [CallGenome+callGenomeSegment_chromId_000_chr1_0008] [taskWrapper-stderr] [2025-01-08T15:31:29.554423Z] [node77.h4h.uhnresearch.ca] [28080_1] [pyflowTaskWrapper:CallGenome+callGenomeSegment_chromId_000_chr1_0008] [wrapperSignal] taskExitCode -11
[2025-01-08T15:31:29.856201Z] [node77.h4h.uhnresearch.ca] [28080_1] [TaskManager] [ERROR] [CallGenome+callGenomeSegment_chromId_000_chr1_0008] [taskWrapper-stderr] [2025-01-08T15:31:29.565490Z] [node77.h4h.uhnresearch.ca] [28080_1] [pyflowTaskWrapper:CallGenome+callGenomeSegment_chromId_000_chr1_0008] [wrapperSignal] taskStderrTail 1
[2025-01-08T15:31:29.864091Z] [node77.h4h.uhnresearch.ca] [28080_1] [TaskManager] [ERROR] [CallGenome+callGenomeSegment_chromId_000_chr1_0008] [taskWrapper-stderr] Last 0 stderr lines from task (of 0 total lines):
[2025-01-08T15:31:29.879082Z] [node77.h4h.uhnresearch.ca] [28080_1] [TaskManager] [ERROR] Shutting down task submission. Waiting for remaining tasks to complete.
```
Probably mem request is too low. In `slurm/cluster.json` change 
```
Strelka:
  cpus-per-task : 4
  time : 3-00:00:00
```
to
```
Strelka:
  cpus-per-task : 4
  mem : 20G
  time : 3-00:00:00
```
