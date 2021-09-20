rule varscanCopyNumber:
  input:
    ref = 'ref/genome.fa',
    normal = norm,
    tumor = "results/alignment/{sample}/{sample}.realigned.recal.bam",
    bed = region
  params:
    outdir="results/Varscan/cnv/{sample}/{sample}.vscn",
  output: "results/Varscan/cnv/{sample}/{sample}.vscn.copynumber"
  threads: 3
  conda:
    "/cluster/home/selghamr/workflows/ExomeSeq/workflow/envs/varscan.yaml"
  shell:
    """
    samtools mpileup -B -q 1 -d 1000000 -l {input.bed} -f {input.ref} {input.normal} {input.tumor} | java -Xmx12g -jar $varscan_dir/VarScan.jar copynumber - {params.outdir} --mpileup 1
    """
rule varscanSomatic:
  input:
    ref = 'ref/genome.fa',
    normal = norm,
    tumor = "results/alignment/{sample}/{sample}.realigned.recal.bam",
    bed = region
  output:
    snp="results/Varscan/snv/{sample}/{sample}.snp",
    indel="results/Varscan/snv/{sample}/{sample}.indel"
  threads: 3
  conda:
    "/cluster/home/selghamr/workflows/ExomeSeq/workflow/envs/varscan.yaml"
  shell:
    """
    samtools mpileup -B -q 1 -d 1000000 -l {input.bed} -f {input.ref} {input.normal} {input.tumor} | java -Xmx12g -jar $varscan_dir/VarScan.jar somatic --mpileup 1 --output-snp {output.snp} --output-indel {output.indel}
    """

rule varscanProcessSomatic:
  input:
    snp="results/Varscan/snv/{sample}/{sample}.snp",
    indel="results/Varscan/snv/{sample}/{sample}.indel"
  output:
    snp="results/Varscan/snv/{sample}/{sample}.snp.Somatic.hc",
    indel="results/Varscan/snv/{sample}/{sample}.indel.Somatic.hc"
  threads: 4
  conda:
    "/cluster/home/selghamr/workflows/ExomeSeq/workflow/envs/varscan.yaml"
  shell:
    """
    java -Xmx8g -jar $varscan_dir/VarScan.jar processSomatic {input.snp} {output.snp}
    java -Xmx8g -jar $varscan_dir/VarScan.jar processSomatic {input.indel} {output.indel}
    """

rule varscanToVCF:
  input:
    snp="results/Varscan/snv/{sample}/{sample}.snp.Somatic.hc",
    indel="results/Varscan/snv/{sample}/{sample}.indel.Somatic.hc"
  output:
    snp="results/Varscan/snv/{sample}/{sample}.snp.Somatic.hc.vcf",
    indel="results/Varscan/snv/{sample}/{sample}.indel.Somatic.hc.vcf"
  threads: 4
  conda:
    "/cluster/home/selghamr/workflows/ExomeSeq/workflow/envs/varscan.yaml"
  shell:
    """
    python scripts/varscan2vcf.py {input.snp} > {output.snp}
    python scripts/varscan2vcf.py {input.indel} > {output.indel}
    """
