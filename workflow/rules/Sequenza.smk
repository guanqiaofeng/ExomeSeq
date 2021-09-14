rule Sequenza:
  input:
    snp="{output_dir}/Varscan/snv/{sample}/{sample}.snp",
    copynum="{output_dir}/Varscan/cnv/{sample}/{sample}.vscn.copynumber"
  params:
    script="/cluster/home/amammoli/SequenzaSingleSample_v2.1_hg19.R",
  output: directory("{output_dir}/Sequenza/{sample}")
  threads: 4
  conda:
    "../envs/sequenza.yaml",
  shell:
     """
     Rscript {params.script} -s {input.snp} -c {input.copynum} -o {output}
     """
