
rule vcftoMAFsnv
  input:
    ref = '/cluster/tools/data/genomes/human/hg38/iGenomes/Sequence/WholeGenomeFasta/genome.fa',
    vcf_inter = "{output_dir}/vcfIntersect/{sample}_intersect_snv"
  params:
    samp="{sample}",
    snv = get_snvs,
    partition="himem"
  output: "{output_dir}/MAF_38_final/snv/{sample}/{snv}.maf",
  threads: 4
  run:
    if wildcards.snv != '0002':
      shell("""module load samtools vep/98
      bcftools view -f PASS {input.vcf_inter}/{params.snv} > {input.vcf_inter}/fil_{params.snv}
      perl /cluster/projects/pughlab/bin/vcf2maf-1.6.17/vcf2maf.pl \
        --input-vcf {input.vcf_inter}/fil_{params.snv} \
        --output-maf {output} \
        --vep-forks 4 \
        --species homo_sapiens \
        --buffer-size 1000 \
        --ref-fasta={input.ref} \
        --filter-vcf /cluster/projects/pughlab/references/VEP_cache/ExAC_nonTCGA.r1.sites.hg19ToHg38.vep.vcf.gz \
        --tumor-id={params.samp} \
        --ncbi-build GRCh38 \
        --vep-path=/cluster/tools/software/centos7/vep/98 \
        --vep-data=/cluster/projects/pughlab/references/VEP_cache/98""")
    else:
      shell("""module load samtools vep/98
      bcftools view -f PASS {input.vcf_inter}/{params.snv} > {input.vcf_inter}/fil_{params.snv}
      perl /cluster/projects/pughlab/bin/vcf2maf-1.6.17/vcf2maf.pl \
        --input-vcf {input.vcf_inter}/fil_{params.snv} \
        --output-maf {output} \
        --vep-forks 4 \
        --species homo_sapiens \
        --buffer-size 1000 \
        --ref-fasta={input.ref} \
        --filter-vcf /cluster/projects/pughlab/references/VEP_cache/ExAC_nonTCGA.r1.sites.hg19ToHg38.vep.vcf.gz \
        --normal-id unmatched \
        --vcf-tumor-id TUMOR \
        --vcf-normal-id NORMAL \
        --tumor-id={params.samp} \
        --ncbi-build GRCh38 \
        --vep-path=/cluster/tools/software/centos7/vep/98 \
        --vep-data=/cluster/projects/pughlab/references/VEP_cache/98""")
