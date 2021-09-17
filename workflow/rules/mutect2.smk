intervals = pd.read_table(
    config["bed_file"]
).set_index(
    "interval", drop=False
)
def get_intervals(wildcards):
    inter = wildcards.interval
    bed = "/cluster/home/selghamr/workflows/ExomeSeq/resources/hg38_bed/" + inter + ".bed"
    return bed

def get_MuTect2_output(wildcards):
    res = []
    for i in intervals.itertuples():
        res.append(
            "results/MuTect2/{}/{}_{}.mut2.vcf".format(
                wildcards.sample, wildcards.sample, i.interval
            )
        )
    return res

output_dir = os.environ.get("output_dir")


rule MuTect2:
  input:
    bam = "results/alignment/{sample}/{sample}.realigned.recal.bam",
    ref = 'ref/genome.fa',
  params:
    intervals = get_intervals
  output: "results/MuTect2/{sample}/{sample}_{interval}.mut2.vcf"
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
  input: get_MuTect2_output,
  params:
    script="scripts/concatvcfs",
    out="results/MuTect2Merge/{sample}/",
    samp="{sample}"
  output: "results/MuTect2Merge/{sample}/{sample}_merged_mut2.vcf"
  threads: 2
  shell:
    """
    echo {input} > debug.tmp ;
    ff= {input} ;
    echo 'sh {params.script}' \
    $ff$' > {output}' > {params.out}/'merge_{params.samp}_VCFs.sh'
    sh {params.out}/merge_{params.samp}_VCFs.sh
    """

rule filterMuTect2:
  input:
    vcf = "results/MuTect2Merge/{sample}/{sample}_merged_mut2.vcf",
  params:
    outdirsnv="results/MuTect2Merge/{sample}/{sample}.snvs",
    outdirindel="results/MuTect2Merge/{sample}/{sample}.indels",
  output:
    snv="results/MuTect2Merge/{sample}/{sample}.snvs.recode.vcf",
    indel="results/MuTect2Merge/{sample}/{sample}.indels.recode.vcf"
  threads: 3
  conda:
    "../envs/gatk.yaml",
  shell:
    """
    vcftools --vcf {input.vcf} --remove-indels --recode --recode-INFO-all --out {params.outdirsnv} --remove-filtered-all
    vcftools --vcf {input.vcf} --keep-only-indels --recode --recode-INFO-all --out {params.outdirindel} --remove-filtered-all
    """
