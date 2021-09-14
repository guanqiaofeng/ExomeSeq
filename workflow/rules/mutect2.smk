intervals = pd.read_table(
    config["bed_file"]
).set_index(
    "interval", drop=False
)
def get_intervals(wildcards):
    inter = wildcards.interval
    bed = "/cluster/home/selghamr/workflows/ExomeSeq/resources/hg38_bed/" + inter + ".bed"
    return bed



rule MuTect2:
  input:
    bam = "{output_dir}/alignment/{sample}/{sample}.realigned.recal.bam",
    ref = '/cluster/tools/data/genomes/human/hg38/iGenomes/Sequence/WholeGenomeFasta/genome.fa',
  params:
    intervals = get_intervals
  output: "{output_dir}/MuTect2/{sample}/{sample}_{interval}.mut2.vcf"
  threads: 2
  conda:
    "../envs/gatk.yaml",
  shell:
    """
    java -Xmx8g -jar $gatk_dir/GenomeAnalysisTK.jar \
    -T MuTect2  \
    -R {input.ref} \
    -L {params.intervals} \
    --input_file:tumor {input.bam} \
    -o {output} \
    """

rule MuTect2Merge:
  input: directory("{output_dir}/MuTect2/{sample}"),
  params:
    script="/cluster/home/amammoli/concatvcfs",
    out="{output_dir}/MuTect2Merge/{sample}",
    samp="{sample}"
  output: "{output_dir}/MuTect2Merge/{sample}/{sample}_merged_mut2.vcf"
  threads: 2
  shell:
    """
    ff=$(find {input} -type f -name '*.mut2.vcf')
    echo 'sh {params.script}' \
    $ff$' > {output}' > {params.out}/'merge_{params.samp}_VCFs.sh'
    sh {params.out}/merge_{params.samp}_VCFs.sh
    """

rule filterMuTect2:
  input:
    vcf = "{output_dir}/MuTect2Merge/{sample}/{sample}_merged_mut2.vcf",
  params:
    outdirsnv="{output_dir}/MuTect2Merge/{sample}/{sample}.snvs",
    outdirindel="{output_dir}/MuTect2Merge/{sample}/{sample}.indels",
  output:
    snv="{output_dir}/MuTect2Merge/{sample}/{sample}.snvs.recode.vcf",
    indel="{output_dir}/MuTect2Merge/{sample}/{sample}.indels.recode.vcf"
  threads: 3
  conda:
    "../envs/gatk.yaml",
  shell:
    """
    vcftools --vcf {input.vcf} --remove-indels --recode --recode-INFO-all --out {params.outdirsnv} --remove-filtered-all
    vcftools --vcf {input.vcf} --keep-only-indels --recode --recode-INFO-all --out {params.outdirindel} --remove-filtered-all
    """
