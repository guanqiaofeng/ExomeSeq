rule vcftoMAFindel:
  input:
    ref = 'ref/genome.fa',
    vcf_inter = "results/vcfIntersect/{sample}_intersect_indel/{indel}.vcf"
  params:
    samp="{sample}",
    indel = get_indels,
  output:  "results/MAF_38_f/indel/{sample}/{indel}.maf",
  threads: 4
  conda:
    "../envs/VCFtoMAF.yaml",
  run:
    if wildcards.indel != '0001':
      """
      bcftools view -f PASS {input.vcf_inter}/{params.indel} > {input.vcf_inter}/fil_{params.indel}
      perl scripts/vcf2maf.pl \
        --input-vcf {input.vcf_inter}/fil_{params.indel} \
        --output-maf {output} \
        --vep-forks 4 \
        --species homo_sapiens \
        --buffer-size 1000 \
        --ref-fasta={input.ref} \
        --tumor-id={params.samp} \
        --ncbi-build GRCh38 \
        --filter-vcf ref/ExAC_nonTCGA.r1.sites.hg19ToHg38.vep.vcf.gz \
        --vep-path=ref/98 \
        --vep-data=ref/98""")
    else:
      """
      bcftools view -f PASS {input.vcf_inter}/{params.indel} > {input.vcf_inter}/fil_{params.indel}
      perl scripts/vcf2maf.pl \
        --input-vcf {input.vcf_inter}/fil_{params.indel} \
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
