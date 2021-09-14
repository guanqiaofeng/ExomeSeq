rule vcfIntersectINDEL:
  input:
    var_vcf="{output_dir}/Varscan/snv/{sample}/{sample}.indel.Somatic.hc.vcf",
    mut2_vcf="{output_dir}/MuTect2Merge/{sample}/{sample}.indels.recode.vcf",
    strelka_vcf="{output_dir}/Strelka/{sample}/{sample}.myAnalysis",
    ref = '/cluster/tools/data/genomes/human/hg38/iGenomes/Sequence/WholeGenomeFasta/genome.fa',
    sequenza="{output_dir}/Sequenza/{sample}"
  params:
    temp_file="{output_dir}/vcfIntersectTempINDEL/{sample}",
    outdir="{output_dir}/vcfIntersectINDEL/{sample}",
    var_final="{output_dir}/vcfIntersectINDEL/{sample}/{sample}_var_indel.gz",
    mut2_final="{output_dir}/vcfIntersectINDEL/{sample}/{sample}_mut2_indel.gz",
    strelka_final="{output_dir}/vcfIntersectINDEL/{sample}/{sample}_strelka_indel.gz",
  output: directory("{output_dir}/vcfIntersect/{sample}_intersect_indel")
  threads: 4
  conda:
    "../envs/varscan.yaml",
  shell:
    """
    mkdir -p {params.temp_file}
    mkdir -p {params.outdir}

    ls -la {input.sequenza}
    vcf-sort -c  {input.var_vcf} > {params.temp_file}_var_sorted.vcf
    vcf-sort -c  {input.mut2_vcf} > {params.temp_file}_mut2_sorted.vcf
    vcf-sort -c  {input.strelka_vcf}/results/passed.somatic.indels.vcf > {params.temp_file}_strelka_sorted.vcf


    ## Left align indels
    java -jar $gatk_dir/GenomeAnalysisTK.jar \
    --unsafe \
    -T LeftAlignAndTrimVariants \
    -R {input.ref} \
    --variant {params.temp_file}_var_sorted.vcf \
    -o {params.temp_file}_var_sorted_left.vcf \
    --dontTrimAlleles

    bgzip -c {params.temp_file}_var_sorted_left.vcf > {params.var_final}
    tabix -p vcf {params.var_final}


    ## Left align indels
    java -jar $gatk_dir/GenomeAnalysisTK.jar \
    --unsafe \
    -T LeftAlignAndTrimVariants \
    -R {input.ref} \
    --variant {params.temp_file}_mut2_sorted.vcf \
    -o {params.temp_file}_mut2_sorted_left.vcf \
    --dontTrimAlleles

    bgzip -c {params.temp_file}_mut2_sorted_left.vcf > {params.mut2_final}
    tabix -p vcf {params.mut2_final}


    ## Left align indels
    java -jar $gatk_dir/GenomeAnalysisTK.jar \
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
    mut1_vcf="{output_dir}/MuTect1/{sample}/{sample}.mut1.vcf",
    var_vcf="{output_dir}/Varscan/snv/{sample}/{sample}.snp.Somatic.hc.vcf",
    mut2_vcf="{output_dir}/MuTect2Merge/{sample}/{sample}.snvs.recode.vcf",
    strelka_vcf="{output_dir}/Strelka/{sample}/{sample}.myAnalysis",
    ref = '/cluster/tools/data/genomes/human/hg38/iGenomes/Sequence/WholeGenomeFasta/genome.fa',
    sequenza="{output_dir}/Sequenza/{sample}",
    validate = "{output_dir}/vcfIntersect/{sample}_intersect_indel",

  params:
    temp_file="{output_dir}/vcfIntersectTemp/{sample}",
    outdir="{output_dir}/vcfIntersectSNV/{sample}",
    mut1_final="{output_dir}/vcfIntersectSNV/{sample}/{sample}_mut1_snv.gz",
    var_final="{output_dir}/vcfIntersectSNV/{sample}/{sample}_var_snv.gz",
    mut2_final="{output_dir}/vcfIntersectSNV/{sample}/{sample}_mut2_snv.gz",
    strelka_final="{output_dir}/vcfIntersectSNV/{sample}/{sample}_strelka_snv.gz",
  output: directory("{output_dir}/vcfIntersect/{sample}_intersect_snv")
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

    ## Left align indels
    java -jar $gatk_dir/GenomeAnalysisTK.jar \
    --unsafe \
    -T LeftAlignAndTrimVariants \
    -R {input.ref} \
    --variant {params.temp_file}_mut1_sorted.vcf \
    -o {params.temp_file}_mut1_sorted_left.vcf \
    --dontTrimAlleles

    bgzip -c {params.temp_file}_mut1_sorted_left.vcf > {params.mut1_final}
    tabix -p vcf {params.mut1_final}

    ## Left align indels
    java -jar $gatk_dir/GenomeAnalysisTK.jar \
    --unsafe \
    -T LeftAlignAndTrimVariants \
    -R {input.ref} \
    --variant {params.temp_file}_var_sorted.vcf \
    -o {params.temp_file}_var_sorted_left.vcf \
    --dontTrimAlleles

    bgzip -c {params.temp_file}_var_sorted_left.vcf > {params.var_final}
    tabix -p vcf {params.var_final}


    ## Left align indels
    java -jar $gatk_dir/GenomeAnalysisTK.jar \
    --unsafe \
    -T LeftAlignAndTrimVariants \
    -R {input.ref} \
    --variant {params.temp_file}_mut2_sorted.vcf \
    -o {params.temp_file}_mut2_sorted_left.vcf \
    --dontTrimAlleles

    bgzip -c {params.temp_file}_mut2_sorted_left.vcf > {params.mut2_final}
    tabix -p vcf {params.mut2_final}


    ## Left align indels
    java -jar $gatk_dir/GenomeAnalysisTK.jar \
    --unsafe \
    -T LeftAlignAndTrimVariants \
    -R {input.ref} \
    --variant {params.temp_file}_strelka_sorted.vcf \
    -o {params.temp_file}_strelka_sorted_left.vcf \
    --dontTrimAlleles

    bgzip -c {params.temp_file}_strelka_sorted_left.vcf > {params.strelka_final}
    tabix -p vcf {params.strelka_final}

    /cluster/projects/pughlab/bin/bcftools/bcftools isec --nfiles +2 {params.outdir}/*snv.gz -p {output}
    """
