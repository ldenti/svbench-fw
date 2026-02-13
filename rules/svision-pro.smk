# svision-pro v2.4 from conda wasn't working. I haven't tried with v2.5 but still, cloning it works.
# I had some conflicts between modules compiled using NumPy 1.x that cannot be run in NumPy 2.0.2
rule get_svisionpro:
    output:
        exe=pjoin(WD, "software", "svision-pro", "SVision-pro"),
        repo=directory(pjoin(WD, "software", "svision-pro")),
        # env=pjoin(WD, "software", "svision-pro", "environment.fixed.yml"),
        model=pjoin(
            WD,
            "software",
            "svision-pro",
            "src",
            "pre_process",
            "model_liteunet_256_8_16_32_32_32.pth",
        ),
    shell:
        """
        rm -rf {output.repo}
        git clone https://github.com/songbowang125/SVision-pro.git {output.repo}
        cd {output.repo}
        git checkout 93f98e413065ed2cbd9ada1a4ac4e16ea7b35740
        chmod +x {output.exe}
        """


rule svisionpro:
    input:
        exe=rules.get_svisionpro.output.exe,
        fa=pjoin(WD, "input", "refs", "{ref}.fa"),
        bam=pjoin(WD, "{ref}", "alignments-ht.bam"),
        model=rules.get_svisionpro.output.model,
    output:
        vcf=pjoin(
            WD,
            "{ref}",
            "svisionpro-s{s}-q{q}",
            f"{SAMPLE_NAME}.svision_pro_v2.5.s{{s}}.vcf",
        ),
    params:
        odir=pjoin(WD, "{ref}", "svisionpro-s{s}-q{q}"),
    log:
        time=pjoin(WD, "times", "{ref}", "svisionpro-s{s}-q{q}.time"),
    conda:
        "../envs/svisionpro.yml"  # pjoin(WD, "software", "svision-pro", "environment.fixed.yml"),
    threads: workflow.cores
    shell:
        """
        /usr/bin/time -vo {log.time} {input.exe} --min_mapq {wildcards.q} --min_sv_size 50 --target_path {input.bam} --genome_path {input.fa} --model_path {input.model} --out_path {params.odir} --sample_name {SAMPLE_NAME} --detect_mode germline --process_num {threads} --min_supp {wildcards.s} --preset hifi
        """


rule svisionpro_post:
    input:
        fa=pjoin(WD, "input", "refs", "{ref}.fa"),
        bam=pjoin(WD, "{ref}", "alignments-ht.bam"),
        vcf=rules.svisionpro.output.vcf,
    output:
        vcf=pjoin(WD, "{ref}", "callsets", "svisionpro-s{s}-q{q}.vcf.gz"),
    log:
        time=pjoin(WD, "times", "{ref}", "svisionpro-s{s}-q{q}.hiphase.time"),
    conda:
        "../envs/hiphase.yml"
    threads: workflow.cores
    shell:
        """
        echo {SAMPLE_NAME} > {input.vcf}.sample.txt
        bcftools sort {input.vcf} | bcftools view --exclude "INFO/SVTYPE='BND'" | python3 ./scripts/to_upper.py | bcftools reheader --samples {input.vcf}.sample.txt | python3 ./scripts/fix_svision.py | bgzip -c > {input.vcf}.gz
        tabix -p vcf {input.vcf}.gz
        /usr/bin/time -vo {log.time} hiphase --bam {input.bam} --reference {input.fa} --vcf {input.vcf}.gz --output-vcf {output.vcf} --threads {threads} > {output.vcf}
        """
