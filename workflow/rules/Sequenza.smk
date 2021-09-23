rule Sequenza:
  input:
    snp="results/Varscan/snv/{sample}/{sample}.snp",
    copynum="results/Varscan/cnv/{sample}/{sample}.vscn.copynumber"
  params:
    script="scripts/SequenzaSingleSample_v2.1_hg38.R",
  output: directory("results/Sequenza/{sample}")
  threads: 4
  conda:
    "/cluster/home/selghamr/workflows/ExomeSeq/workflow/envs/sequenza.yaml",
  shell:
     """
     Rscript {params.script} -s {input.snp} -c {input.copynum} -o {output}
     """
