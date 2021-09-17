
rule vcftoMAFsnv
  input:
    ref = 'ref/genome.fa',
    vcf_inter = expand("results/vcfIntersect/{sample}_intersect_snv", sample = samples["sample"])
  params:
    samp="{sample}",
    snv = get_snvs,
  output: get_maf_output(type=snv),
  threads: 4
  conda:
    "../envs/VCFtoMAF.yaml",
  run:
    if wildcards.snv != '0002':
      shell(
      """
      bcftools view -f PASS {input.vcf_inter}/{params.snv} > {input.vcf_inter}/fil_{params.snv}
      perl scripts/vcf2maf.pl \
        --input-vcf {input.vcf_inter}/fil_{params.snv} \
        --output-maf {output} \
        --vep-forks 4 \
        --species homo_sapiens \
        --buffer-size 1000 \
        --ref-fasta={input.ref} \
        --filter-vcf ref/ExAC_nonTCGA.r1.sites.hg19ToHg38.vep.vcf.gz \
        --tumor-id={params.samp} \
        --ncbi-build GRCh38 \
        --vep-path=ref/98 \
        --vep-data=ref/98""")
    else:
      shell(
      """
      bcftools view -f PASS {input.vcf_inter}/{params.snv} > {input.vcf_inter}/fil_{params.snv}
      perl scripts/vcf2maf.pl \
        --input-vcf {input.vcf_inter}/fil_{params.snv} \
        --output-maf {output} \
        --vep-forks 4 \
        --species homo_sapiens \
        --buffer-size 1000 \
        --ref-fasta={input.ref} \
        --filter-vcf ref/ExAC_nonTCGA.r1.sites.hg19ToHg38.vep.vcf.gz \
        --normal-id unmatched \
        --vcf-tumor-id TUMOR \
        --vcf-normal-id NORMAL \
        --tumor-id={params.samp} \
        --ncbi-build GRCh38 \
        --vep-path=ref/98 \
        --vep-data=ref/98""")
