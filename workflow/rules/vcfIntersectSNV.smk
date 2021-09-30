rule vcfIntersectINDEL:
  input:
    var_vcf="results/Varscan/snv/{sample}/{sample}.indel.Somatic.hc.vcf",
    mut2_vcf="results/MuTect2Merge/{sample}/{sample}.indels.recode.vcf",
    strelka_vcf="results/Strelka/{sample}/{sample}.myAnalysis",
    ref = 'ref/genome.fa',
    sequenza="results/Sequenza/{sample}.gz"
  params:
    temp_file="results/vcfIntersectTempINDEL/{sample}",
    outdir="results/vcfIntersectINDEL/{sample}",
    var_final="results/vcfIntersectINDEL/{sample}/{sample}_var_indel.gz",
    mut2_final="results/vcfIntersectINDEL/{sample}/{sample}_mut2_indel.gz",
    strelka_final="results/vcfIntersectINDEL/{sample}/{sample}_strelka_indel.gz",
    gatk="/cluster/home/selghamr/workflows/ExomeSeq/.snakemake/conda/9933b5f3a92c804102746a579b8a499c/opt/gatk-3.8"
  output: directory("results/vcfIntersect/{sample}_intersect_indel")
  threads: 4
  conda:
    "/cluster/home/selghamr/workflows/ExomeSeq/workflow/envs/gatk.yaml",
  shell:
    """
    mkdir -p {params.temp_file}
    mkdir -p {params.outdir}

    ls -la {input.sequenza}
    vcf-sort -c  {input.var_vcf} > {params.temp_file}_var_sorted.vcf
    vcf-sort -c  {input.mut2_vcf} > {params.temp_file}_mut2_sorted.vcf
    vcf-sort -c  {input.strelka_vcf}/results/passed.somatic.indels.vcf > {params.temp_file}_strelka_sorted.vcf


    ## Left align indels
    gatk3 -Xmx12g \
    --unsafe \
    -T LeftAlignAndTrimVariants \
    -R {input.ref} \
    --variant {params.temp_file}_var_sorted.vcf \
    -o {params.temp_file}_var_sorted_left.vcf \
    --dontTrimAlleles

    bgzip -c {params.temp_file}_var_sorted_left.vcf > {params.var_final}
    tabix -p vcf {params.var_final}


    ## Left align indels
    gatk3 -Xmx12g \
    --unsafe \
    -T LeftAlignAndTrimVariants \
    -R {input.ref} \
    --variant {params.temp_file}_mut2_sorted.vcf \
    -o {params.temp_file}_mut2_sorted_left.vcf \
    --dontTrimAlleles

    bgzip -c {params.temp_file}_mut2_sorted_left.vcf > {params.mut2_final}
    tabix -p vcf {params.mut2_final}


    ## Left align indels
    gatk3 -Xmx12g \
    --unsafe \
    -T LeftAlignAndTrimVariants \
    -R {input.ref} \
    --variant {params.temp_file}_strelka_sorted.vcf \
    -o {params.temp_file}_strelka_sorted_left.vcf \
    --dontTrimAlleles

    bgzip -c {params.temp_file}_strelka_sorted_left.vcf > {params.strelka_final}
    tabix -p vcf {params.strelka_final}

    /cluster/projects/pughlab/bin/bcftools/bcftools isec --nfiles +2 {params.outdir}/*indel.gz -p {output}
    """

rule vcfIntersectSNV:
  input:
    mut1_vcf="results/MuTect1/{sample}/{sample}.mut1.vcf",
    var_vcf="results/Varscan/snv/{sample}/{sample}.snp.Somatic.hc.vcf",
    mut2_vcf="results/MuTect2Merge/{sample}/{sample}.snvs.recode.vcf",
    strelka_vcf="results/Strelka/{sample}/{sample}.myAnalysis",
    ref = 'ref/genome.fa',
    sequenza="results/Sequenza/{sample}",
    validate = "results/vcfIntersect/{sample}_intersect_indel",

  params:
    temp_file="results/vcfIntersectTemp/{sample}",
    outdir="results/vcfIntersectSNV/{sample}",
    mut1_final="results/vcfIntersectSNV/{sample}/{sample}_mut1_snv.gz",
    var_final="results/vcfIntersectSNV/{sample}/{sample}_var_snv.gz",
    mut2_final="results/vcfIntersectSNV/{sample}/{sample}_mut2_snv.gz",
    strelka_final="results/vcfIntersectSNV/{sample}/{sample}_strelka_snv.gz",
    snv=get_snv_intersects,
    gatk="/cluster/home/selghamr/workflows/ExomeSeq/.snakemake/conda/9933b5f3a92c804102746a579b8a499c/opt/gatk-3.8"
  output:
    dir=directory("results/vcfIntersect/{sample}_intersect_snv"),
    file="results/vcfIntersect/{sample}_intersect_snv/{snv}.vcf"
  threads: 4
  conda:
    "../envs/varscan.yaml",
  shell:
    """
    mkdir -p {params.temp_file}
    mkdir -p {params.outdir}
    ls -la {input.sequenza}
    ls -la {input.validate}

    vcf-sort -c  {input.mut1_vcf} > {params.temp_file}_mut1_sorted.vcf
    vcf-sort -c  {input.var_vcf} > {params.temp_file}_var_sorted.vcf
    vcf-sort -c  {input.mut2_vcf} > {params.temp_file}_mut2_sorted.vcf
    vcf-sort -c  {input.strelka_vcf}/results/passed.somatic.snvs.vcf > {params.temp_file}_strelka_sorted.vcf

    ## Left align indels (MuTect1)
    gatk3 -Xmx12g \
    --unsafe \
    -T LeftAlignAndTrimVariants \
    -R {input.ref} \
    --variant {params.temp_file}_mut1_sorted.vcf \
    -o {params.temp_file}_mut1_sorted_left.vcf \
    --dontTrimAlleles

    bgzip -c {params.temp_file}_mut1_sorted_left.vcf > {params.mut1_final}
    tabix -p vcf {params.mut1_final}

    ## Left align indels (Varscan)
    gatk3 -Xmx12g \
    --unsafe \
    -T LeftAlignAndTrimVariants \
    -R {input.ref} \
    --variant {params.temp_file}_var_sorted.vcf \
    -o {params.temp_file}_var_sorted_left.vcf \
    --dontTrimAlleles

    bgzip -c {params.temp_file}_var_sorted_left.vcf > {params.var_final}
    tabix -p vcf {params.var_final}


    ## Left align indels (MuTect2)
    gatk3 -Xmx12g \
    --unsafe \
    -T LeftAlignAndTrimVariants \
    -R {input.ref} \
    --variant {params.temp_file}_mut2_sorted.vcf \
    -o {params.temp_file}_mut2_sorted_left.vcf \
    --dontTrimAlleles

    bgzip -c {params.temp_file}_mut2_sorted_left.vcf > {params.mut2_final}
    tabix -p vcf {params.mut2_final}


    ## Left align indels (Strelka)
    gatk3 -Xmx12g \
    --unsafe \
    -T LeftAlignAndTrimVariants \
    -R {input.ref} \
    --variant {params.temp_file}_strelka_sorted.vcf \
    -o {params.temp_file}_strelka_sorted_left.vcf \
    --dontTrimAlleles

    bgzip -c {params.temp_file}_strelka_sorted_left.vcf > {params.strelka_final}
    tabix -p vcf {params.strelka_final}

    /cluster/projects/pughlab/bin/bcftools/bcftools isec --nfiles {params.snv} {params.outdir}/*snv.gz -p {output.dir}
    """
