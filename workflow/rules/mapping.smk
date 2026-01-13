"""
Mapping rules: read mapping, merging, duplicate marking, and filtering
"""


rule map_init:
    """
    Map raw reads to reference genome using bowtie2.
    Auto-detects single-end vs paired-end based on units.tsv.
    """
    input:
        fq=get_fq,
    output:
        bam=OUT_DIR
        + "/{sample}/bam/{sample}.{lib}.{unit}."
        + PREFIX
        + ".{ref}.init.bam",
    threads: 12
    log:
        OUT_DIR
        + "/{sample}/logs/{sample}.{lib}.{unit}."
        + PREFIX
        + ".{ref}.map_init.log",
    params:
        ref=lambda wildcards: config["ref"][wildcards.ref],
        bt_param=config["bt_param"],
    run:
        # Check if single-end or paired-end
        if len(input.fq) == 1:
            # Single-end
            shell(
                f"(bowtie2 -x {params.ref} -p {threads} --rg-id {wildcards.unit} "
                f"--rg LB:{wildcards.lib} --rg SM:{wildcards.sample} {params.bt_param} "
                f"-U {input.fq[0]} | samtools view -bh - > {output.bam}) &> {log}"
            )
        elif len(input.fq) == 2:
            # Paired-end (PE1 and PE2 only)
            shell(
                f"(bowtie2 -x {params.ref} -p {threads} --rg-id {wildcards.unit} "
                f"--rg LB:{wildcards.lib} --rg SM:{wildcards.sample} {params.bt_param} "
                f"-1 {input.fq[0]} -2 {input.fq[1]} | samtools view -bh - > {output.bam}) &> {log}"
            )
        elif len(input.fq) == 3:
            # Paired-end with collapsed reads (PE1, PE2, and collapsed)
            shell(
                f"(bowtie2 -x {params.ref} -p {threads} --rg-id {wildcards.unit} "
                f"--rg LB:{wildcards.lib} --rg SM:{wildcards.sample} {params.bt_param} "
                f"-1 {input.fq[0]} -2 {input.fq[1]} -U {input.fq[2]} | samtools view -bh - > {output.bam}) &> {log}"
            )


rule merge_lib:
    """Merge multiple units within a library, remove unmapped reads, and sort."""
    input:
        get_bam_lib,
    output:
        bam=temp(TMP_DIR + "/{sample}/{sample}.{lib}." + PREFIX + ".{ref}.srt.bam"),
    threads: 2
    shell:
        """
        samtools merge -@{threads} - {input} |\
        samtools view -bh -F4 - |\
        samtools sort -@{threads} - > {output.bam}
        """


rule mark_duplicates:
    """Mark PCR/optical duplicates using Picard MarkDuplicates."""
    input:
        bam=TMP_DIR + "/{sample}/{sample}.{lib}." + PREFIX + ".{ref}.srt.bam",
    output:
        bam=temp(TMP_DIR + "/{sample}/{sample}.{lib}." + PREFIX + ".{ref}.mrkdup.bam"),
        metrics=OUT_DIR
        + "/{sample}/logs/{sample}.{lib}."
        + PREFIX
        + ".{ref}.mrkdup_metrics.txt",
    shell:
        """
        picard MarkDuplicates I={input} O={output.bam} M={output.metrics}
        """


rule merge_sample:
    """Merge all libraries for a sample."""
    input:
        get_bam_sample,
    threads: 8
    output:
        bam=OUT_DIR + "/{sample}/bam/{sample}." + PREFIX + ".{ref}.mrkdup.bam",
        bai=OUT_DIR + "/{sample}/bam/{sample}." + PREFIX + ".{ref}.mrkdup.bam.bai",
    shell:
        """
        samtools merge -@{threads} - {input} > {output.bam}
        samtools index {output.bam}
        """


rule filter_bam:
    """Filter BAM by mapping quality and read length."""
    input:
        bam=OUT_DIR + "/{sample}/bam/{sample}." + PREFIX + ".{ref}.mrkdup.bam",
    output:
        bam=OUT_DIR + "/{sample}/bam/{sample}." + PREFIX + ".{ref}.filter.bam",
        bai=OUT_DIR + "/{sample}/bam/{sample}." + PREFIX + ".{ref}.filter.bam.bai",
    params:
        qlen=config["bam_filter"]["qlen"],
        rlen=config["bam_filter"]["rlen"],
        mq=config["bam_filter"]["mq"],
    shell:
        """
        samtools view -q{params.mq} -F3844 -e 'qlen >= {params.qlen} && rlen >= {params.rlen}' -bh {input.bam} > {output.bam}
        samtools index {output.bam}
        """
