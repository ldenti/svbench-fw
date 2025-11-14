BAM = pjoin(WD, "minimap2.bam")
BAM_HT = pjoin(WD, "minimap2-haplotagged.bam")


rule faidx:
    input:
        REF,
    output:
        REF + ".fai",
    conda:
        "../envs/samtools.yml"
    shell:
        """
        samtools faidx {input}
        """


rule align_reads:
    input:
        fa=REF,
        fq=FQ,
    output:
        bam=BAM,
    threads: workflow.cores
    conda:
        "../envs/minimap2.yml"
    shell:
        """
        minimap2 -ax map-hifi --MD --eqx -Y -R '@RG\\tID:{SAMPLE_NAME}\\tSM:{SAMPLE_NAME}' -t {threads} {input.fa} {input.fq} | samtools view -bS | samtools sort -@ $(({threads}-1)) -T {output.bam} > {output.bam}
        samtools index {output.bam}
        """


rule deepvariant:
    input:
        fa=REF,
        bam=BAM,
    output:
        vcf=pjoin(WD, "deepvariant.vcf.gz"),
    threads: workflow.cores
    conda:
        "../envs/deepvariant.yml"
    shell:
        """
        singularity run --bind /:/wd docker://google/deepvariant:1.8.0 /opt/deepvariant/bin/run_deepvariant --model_type PACBIO --ref /wd/{input.fa} --reads /wd/{input.bam} --output_vcf /wd/{output.vcf} --num_shards {threads} --sample_name {SAMPLE_NAME}
        # tabix -p vcf {output.vcf}
        """


rule whathshap:
    input:
        fa=REF,
        bam=BAM,
        vcf=rules.deepvariant.output.vcf,
    output:
        vcf=pjoin(WD, "deepvariant-phased.vcf.gz"),
    threads: workflow.cores
    conda:
        "../envs/whatshap.yml"
    shell:
        """
        whatshap phase -o {output.vcf} --reference {input.fa} {input.vcf} {input.bam}
        tabix -p vcf {output.vcf}
        """


rule wh_haplo:
    input:
        fa=REF,
        bam=BAM,
        vcf=rules.whathshap.output.vcf,
    output:
        bam=BAM_HT,
    threads: workflow.cores
    conda:
        "../envs/whatshap.yml"
    shell:
        """
        whatshap haplotag --reference {input.fa} {input.vcf} {input.bam} | samtools view -bS | samtools sort > {output.bam}
        samtools index {output.bam}
        """
