
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

mkdir config data resources


```

