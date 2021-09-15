
rule haploGATK:
  input:
    bam="results/alignment/{sample}/{sample}.realigned.recal.bam",
    ref = '/cluster/tools/data/genomes/human/hg38/iGenomes/Sequence/WholeGenomeFasta/genome.fa',
  output: "results/Haplotype/{sample}/{sample}.raw.snps.indels.vcf"
  threads: 4
  conda:
    "../envs/gatk.yaml"
  shell:
    """
    java -Xmx8g -jar $gatk_dir/GenomeAnalysisTK.jar -T HaplotypeCaller -nct 4 -R {input.ref} -I {input.bam} -o {output}
    """
