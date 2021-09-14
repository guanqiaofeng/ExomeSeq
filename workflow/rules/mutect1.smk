rule MuTect1:
  input:
    bam = "{output_dir}/alignment/{sample}/{sample}.realigned.recal.bam",
    ref = '/cluster/tools/data/genomes/human/hg38/iGenomes/Sequence/WholeGenomeFasta/genome.fa',
    interval = region
  output:
    vcf="{output_dir}/MuTect1/{sample}/{sample}.mut1.vcf",
    stats="{output_dir}/MuTect1/{sample}/{sample}.call_stats",
    coverage="{output_dir}/MuTect1/{sample}/{sample}.wig.txt"
  threads: 2
  conda:
    "../envs/mutect1.yaml",
  shell:
    """
    java -Xmx10g -jar $mutect_dir/muTect-1.1.4.jar -T MuTect \
    -R {input.ref} \
    -L {input.interval} \
    --input_file:tumor {input.bam} \
    --vcf {output.vcf} \
    --out {output.stats} \
    --coverage_file {output.coverage} \
    --downsampling_type NONE \
    --fraction_contamination 0.02
    """
