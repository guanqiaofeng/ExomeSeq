
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
cd ../resources
touch genome.fa xenomeidx-both.kmers.high-bits xenomeidx-graft.kmers.low-bits.lwr INS-B-014-SB.processed.bam xenomeidx-both.kmers.low-bits.lwr xenomeidx-graft.kmers.low-bits.upr S04380110_Covered.headless.bed xenomeidx-both.kmers.low-bits.upr xenomeidx-host.header xengsortidx.hash xenomeidx-both.lhs-bits xenomeidx-host.kmers-d0 xengsortidx.info xenomeidx-both.rhs-bits xenomeidx-host.kmers-d1 xengsort.sif xenomeidx-graft.header xenomeidx-host.kmers.header xenomeidx-both.header xenomeidx-graft.kmers-d0 xenomeidx-host.kmers.high-bits xenomeidx-both.kmers-d0 xenomeidx-graft.kmers-d1 xenomeidx-host.kmers.low-bits.lwr xenomeidx-both.kmers-d1 xenomeidx-graft.kmers.header xenomeidx-host.kmers.low-bits.upr xenomeidx-both.kmers.header xenomeidx-graft.kmers.high-bits
```
