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
    samp="{sample}"
  output:
    bash_indel="results/vcfIntersect/bash_scripts/{sample}_Indel_overlap.sh",
  shell:
    """
    echo "sh {params.script} {params.outdir}/indels {params.samp} {params.samp} {input.var_vcf} {input.mut2_vcf} {input.strelka_vcf}" > {output.bash_indel}
    sh {output.bash_snv}
    """

rule vcfIntersectSNV:
  input:
    var_vcf="results/Varscan/snv/{sample}/{sample}.snp.Somatic.hc.vcf",
    mut2_vcf="results/MuTect2Merge/{sample}/{sample}.snvs.recode.vcf",
    strelka_vcf="results/Strelka/{sample}/{sample}.myAnalysis/results/variants/somatic.snvs.vcf.gz",
    mut1_vcf="results/MuTect1/{sample}/{sample}.mut1.vcf",
    ref = 'ref/genome.fa',
  params:
    outdir="results/vcfIntersect",
    script= "vcfIntersect.sh",
    samp="{sample}"
  output:
    #bash_snv="results/vcfIntersect/bash_scripts/{sample}_snvs_overlap.sh",
    vcf_file="results/vcfIntersect/indel/{sample}/{indel}.vcf",
    vcf_fil="results/vcfIntersect/indel/{sample}/fil_{indel}.vcf"
  shell:
    "sh {params.script} "
    "{params.outdir}/snvs "
    "{params.samp} "
    "{params.samp} "
    "{input.mut1_vcf} "
    "{input.var_vcf} "
    "{input.mut2_vcf} "
    "{input.strelka_vcf} "
