rule Strelka:
  input:
    tumor="results/alignment/{sample}/{sample}.realigned.recal.bam",
    ref = 'ref/genome.fa',
    conf="config/strelka_config_bwa.ini",
    normal=norm,
  params:
    outdir="results/Strelka/{sample}/{sample}",
  output: directory("results/Strelka/{sample}/{sample}.myAnalysis")
  threads: 4
  conda:
    "../envs/strelka.yaml",
  shell:
    """
    configureStrelkaWorkflow.pl --normal={input.normal}  --tumor={input.tumor} --ref={input.ref} --config={input.conf} --output-dir={output}
    make -C {params.outdir}.myAnalysis
    """
