norm="/cluster/projects/bhklab/Data/TFRI_Exome/alignment/Pugh_INS_ex_dil_111_AGTACAAG_L005/Pugh_INS_ex_dil_111_AGTACAAG_L005.realigned.recal.bam"
rule varscanCopyNumber:
  input:
    ref = '/cluster/tools/data/genomes/human/hg38/iGenomes/Sequence/WholeGenomeFasta/genome.fa',
    normal = norm,
    tumor = "{output_dir}/alignment/{sample}/{sample}.realigned.recal.bam",
    bed = region
  params:
    outdir="{output_dir}/Varscan/cnv/{sample}/{sample}.vscn",
  output: "{output_dir}/Varscan/cnv/{sample}/{sample}.vscn.copynumber"
  threads: 3
  conda:
    "../envs/varscan.yaml",
  shell:
    """
    samtools mpileup -B -q 1 -d 1000000 -l {input.bed} -f {input.ref} {input.normal} {input.tumor} | java -Xmx12g -jar $varscan_dir/VarScan.jar copynumber - {params.outdir} --mpileup 1
    """
rule varscanSomatic:
  input:
    ref = '/cluster/tools/data/genomes/human/hg38/iGenomes/Sequence/WholeGenomeFasta/genome.fa',
    normal = norm,
    tumor = "{output_dir}/alignment/{sample}/{sample}.realigned.recal.bam",
    bed = region
  output:
    snp="{output_dir}/Varscan/snv/{sample}/{sample}.snp",
    indel="{output_dir}/Varscan/snv/{sample}/{sample}.indel"
  threads: 3
  conda:
    "../envs/varscan.yaml",
  shell:
    """
    samtools mpileup -B -q 1 -d 1000000 -l {input.bed} -f {input.ref} {input.normal} {input.tumor} | java -Xmx12g -jar $varscan_dir/VarScan.jar somatic --mpileup 1 --output-snp {output.snp} --output-indel {output.indel}
    """

rule varscanProcessSomatic:
  input:
    snp="{output_dir}/Varscan/snv/{sample}/{sample}.snp",
    indel="{output_dir}/Varscan/snv/{sample}/{sample}.indel"
  output:
    snp="{output_dir}/Varscan/snv/{sample}/{sample}.snp.Somatic.hc",
    indel="{output_dir}/Varscan/snv/{sample}/{sample}.indel.Somatic.hc"
  threads: 4
  conda:
    "../envs/varscan.yaml",
  shell:
    """
    java -Xmx8g -jar $varscan_dir/VarScan.jar processSomatic {input.snp} {output.snp}
    java -Xmx8g -jar $varscan_dir/VarScan.jar processSomatic {input.indel} {output.indel}
    """

rule varscanToVCF:
  input:
    snp="{output_dir}/Varscan/snv/{sample}/{sample}.snp.Somatic.hc",
    indel="{output_dir}/Varscan/snv/{sample}/{sample}.indel.Somatic.hc"
  output:
    snp="{output_dir}/Varscan/snv/{sample}/{sample}.snp.Somatic.hc.vcf",
    indel="{output_dir}/Varscan/snv/{sample}/{sample}.indel.Somatic.hc.vcf"
  threads: 4
  conda:
    "../envs/varscan.yaml",
  shell:
    """
    python /cluster/home/amammoli/varscan/hg19/varscan2vcf.py {input.snp} > {output.snp}
    python /cluster/home/amammoli/varscan/hg19/varscan2vcf.py {input.indel} > {output.indel}
    """
