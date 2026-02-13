rule faidx:
    input:
        pjoin(WD, "input", "refs", "{ref}.fa"),
    output:
        pjoin(WD, "input", "refs", "{ref}.fa.fai"),
    conda:
        "../envs/samtools.yml"
    shell:
        """
        samtools faidx {input}
        """


rule align_reads:
    input:
        fa=pjoin(WD, "input", "refs", "{ref}.fa"),
        fq=FQ,
    output:
        bam=pjoin(WD, "{ref}", "alignments.bam"),
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
        fa=pjoin(WD, "input", "refs", "{ref}.fa"),
        fai=pjoin(WD, "input", "refs", "{ref}.fa.fai"),
        bam=rules.align_reads.output.bam,
    output:
        vcf=pjoin(WD, "{ref}", "deepvariant.vcf.gz"),
    params:
        # XXX: ugly workaround, this might not work everytime
        bind="/" + WD.split("/")[1],
    threads: workflow.cores
    conda:
        "../envs/deepvariant.yml"
    shell:
        """
        singularity run --bind {params.bind}:{params.bind} docker://google/deepvariant:1.9.0 /opt/deepvariant/bin/run_deepvariant --model_type PACBIO --ref {input.fa} --reads {input.bam} --output_vcf {output.vcf} --num_shards {threads} --sample_name {SAMPLE_NAME}
        # tabix -p vcf {output.vcf}
        """


rule whathshap:
    input:
        fa=pjoin(WD, "input", "refs", "{ref}.fa"),
        fai=pjoin(WD, "input", "refs", "{ref}.fa.fai"),
        bam=rules.align_reads.output.bam,
        vcf=rules.deepvariant.output.vcf,
    output:
        vcf=pjoin(WD, "{ref}", "deepvariant-phased.vcf.gz"),
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
        fa=pjoin(WD, "input", "refs", "{ref}.fa"),
        fai=pjoin(WD, "input", "refs", "{ref}.fa.fai"),
        bam=rules.align_reads.output.bam,
        vcf=rules.whathshap.output.vcf,
    output:
        bam=pjoin(WD, "{ref}", "alignments-ht.bam"),
    threads: workflow.cores
    conda:
        "../envs/whatshap.yml"
    shell:
        """
        whatshap haplotag --reference {input.fa} {input.vcf} {input.bam} | samtools view -bS | samtools sort > {output.bam}
        samtools index {output.bam}
        """
