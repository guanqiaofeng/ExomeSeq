rule vcfIntersectINDEL:
  input:
    var_vcf="results/Varscan/snv/{sample}/{sample}.indel.Somatic.hc.vcf",
    mut2_vcf="results/MuTect2Merge/{sample}/{sample}.indels.recode.vcf",
    strelka_vcf="results/Strelka/{sample}/{sample}.myAnalysis/results/variants/somatic.indels.vcf.gz",
    ref = 'ref/genome.fa',
    sequenza="results/Sequenza/{sample}.gz"
  params:
    outdir="results/vcfIntersect",
    script= "vcfIntersect.sh",
  shell:
    """
    echo "sh {params.script} {params.outdir}/indels {sample} {sample} {input.var_vcf} {input.mut2_vcf} {input.strelka_vcf}" \
    > {params.outdir}/bash_scripts/{sample}_Indel_overlap.sh
    sh {params.outdir}/bash_scripts/{sample}_Indel_overlap.sh
    """

rule vcfIntersectSNV:
  input:
    var_vcf="results/Varscan/snv/{sample}/{sample}.snvs.Somatic.hc.vcf",
    mut2_vcf="results/MuTect2Merge/{sample}/{sample}.snvs.recode.vcf",
    strelka_vcf="results/Strelka/{sample}/{sample}.myAnalysis/results/variants/somatic.snvs.vcf.gz",
    mut1_vcf="results/MuTect1/{sample}/{sample}.mut1.vcf ",
    ref = 'ref/genome.fa',
  params:
    outdir="results/vcfIntersect",
    script= "vcfIntersect.sh",
  shell:
    """
    echo "sh {params.script} {params.outdir}/snvs {sample} {sample} {input.mut1_vcf} {input.var_vcf} {input.mut2_vcf} {input.strelka_vcf}" \
    > {params.outdir}/bash_scripts/{sample}_snvs_overlap.sh
    sh {params.outdir}/bash_scripts/{sample}_snvs_overlap.sh
    """
