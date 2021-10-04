
rule vcftoMAFsnv:
  input:
    ref = 'ref/genome.fa',
    vcf_inter = "results/vcfIntersect/snvs/{sample}",
  params:
    samp="{sample}",
    snvs = "{snv}",
#    snv = get_snvs,
#    intersect = get_snv_intersects,
  output:  "results/MAF_38_f/snv/{sample}/{intersect}.maf",
  threads: 4
  conda:
    "../envs/VCFtoMAF.yaml",
  shell:
    """
    if [ {params.snvs} != '0002' ]; then
      bcftools view -f PASS {input.vcf_inter}/{params.snv}.vcf > {input.vcf_inter}/fil_{params.snv}.vcf;
      perl scripts/vcf2maf.pl \
        --input-vcf {input.vcf_inter}/fil_{params.snv}.vcf \
        --output-maf {output} \
        --vep-forks 4 \
        --species homo_sapiens \
        --buffer-size 1000 \
        --ref-fasta={input.ref} \
        --filter-vcf ref/ExAC_nonTCGA.r1.sites.hg19ToHg38.vep.vcf.gz \
        --tumor-id={params.samp} \
        --ncbi-build GRCh38 \
        --vep-path=ref/98 \
        --vep-data=ref/98
    else
      bcftools view -f PASS {input.vcf_inter}/{params.snv}.vcf > {input.vcf_inter}/fil_{params.snv}.vcf;
      perl scripts/vcf2maf.pl \
        --input-vcf {input.vcf_inter}/fil_{params.snv}.vcf \
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
        --vep-data=ref/98
    fi
        """
