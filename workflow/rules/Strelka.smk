rule Strelka:
  input:
    tumor="results/alignment/{sample}/{sample}.realigned.recal.bam",
    ref = 'ref/genome.fa',
    conf="config/strelka_config_bwa.ini",
    normal=norm,
  params:
    strelka="/cluster/home/selghamr/workflows/ExomeSeq/.snakemake/conda/236aa367b1347b9561439ce4facd36c0/share/strelka-2.9.10-1/bin",
    outdir="results/Strelka/{sample}/{sample}",
  output: directory("results/Strelka/{sample}/{sample}.myAnalysis")
  threads: 4
  conda:
    "/cluster/home/selghamr/workflows/ExomeSeq/workflow/envs/strelka.yaml",
  shell:
    """
    {params.strelka}/configureStrelkaSomaticWorkflow.py --normal={input.normal}  --tumor={input.tumor} --ref={input.ref} --config={input.conf} --runDir={output}
    make -C {params.outdir}.myAnalysis
    """
