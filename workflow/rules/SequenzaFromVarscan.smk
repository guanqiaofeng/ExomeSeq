rule Sequenza:
    input:
        snp="results/Varscan/Varscan/snv/{sample}/{sample}.snp",
        copynum="results/Varscan/Varscan/cnv/{sample}/{sample}.vscn.copynumber"
    params:
        script="scripts/SequenzaSingleSample_v2.1_hg38.R",
    output: directory("results/Sequenza/{sample}")
    threads: 4
    shell:
     """
     module load R/3.3.0
     Rscript {params.script} -s {input.snp} -c {input.copynum} -o {output}
     """
