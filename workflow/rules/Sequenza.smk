rule Sequenza:
  input:
    bam="results/alignment/{sample}/{sample}.realigned.recal.bam",
    #snp="results/Varscan/snv/{sample}/{sample}.snp",
    #copynum="results/Varscan/cnv/{sample}/{sample}.vscn.copynumber"
    normal=norm,
  params:
    #script="scripts/SequenzaSingleSample_v2.1_hg38.R",
    refgc='ref/GCgenome.wig'
    fasta='ref/genome.fa'
  output: "results/Sequenza/{sample}.gz"
  threads: 4
  conda:
    "/cluster/home/selghamr/workflows/ExomeSeq/workflow/envs/sequenza.yaml",
  shell:
    """
    sequenza-utils \
    bam2seqz -t {input.bam} \
    -n {input.norm} \
    -gc {params.refgc} \
    -F {params.fasta} \
    -o {output}
    """

#     """
#     Rscript {params.script} -s {input.snp} -c {input.copynum} -o {output}
#     """
