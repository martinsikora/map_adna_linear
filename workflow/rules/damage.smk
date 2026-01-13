## --------------------------------------------------------------------------------
## Module: Damage Estimation
## Estimate ancient DNA damage patterns using metaDMG and mapDamage.
## --------------------------------------------------------------------------------


rule get_mapdamage:
    input:
        bam=OUT_DIR + "/{sample}/bam/{sample}." + PREFIX + ".{ref}.filter.bam",
        bai=OUT_DIR + "/{sample}/bam/{sample}." + PREFIX + ".{ref}.filter.bam.bai",
    output:
        frag_m=OUT_DIR
        + "/{sample}/mapdamage/{sample}."
        + PREFIX
        + ".{ref}.filter/Fragmisincorporation_plot.pdf",
        frag_l=OUT_DIR
        + "/{sample}/mapdamage/{sample}."
        + PREFIX
        + ".{ref}.filter/Length_plot.pdf",
    params:
        ref=lambda wildcards: config["ref"][wildcards.ref] + ".fasta",
        px=OUT_DIR
        + "/{sample}/mapdamage/{sample}."
        + PREFIX
        + ".{ref}.filter",
    shell:
        """
        mkdir -p {params.px}
        touch {params.px}/Fragmisincorporation_plot.pdf
        touch {params.px}/Length_plot.pdf
        mapDamage -i {input.bam} -r {params.ref} -d {params.px} --no-stats
        """


rule sort_bam_sample:
    input:
        bam=OUT_DIR + "/{sample}/bam/{sample}." + PREFIX + ".{ref}.filter.bam",
        bai=OUT_DIR + "/{sample}/bam/{sample}." + PREFIX + ".{ref}.filter.bam.bai",
    output:
        bam=temp(TMP_DIR + "/{sample}/{sample}.{ref}.sample.filter.srt_name.bam"),
    params:
        mq=config["bam_filter"]["mq"],
    shell:
        """
        samtools sort -n {input.bam} > {output.bam}
        """


rule get_damage_global_sample:
    input:
        bam=TMP_DIR + "/{sample}/{sample}.{ref}.sample.filter.srt_name.bam",
    output:
        dam=OUT_DIR
        + "/{sample}/metadamage/sample_level/{ref}."
        + PREFIX
        + ".filter.global.tsv",
        bdam=OUT_DIR
        + "/{sample}/metadamage/sample_level/{ref}."
        + PREFIX
        + ".filter.global.bdamage.gz",
        dfit=OUT_DIR
        + "/{sample}/metadamage/sample_level/{ref}."
        + PREFIX
        + ".filter.global.dfit.gz",
    log:
        OUT_DIR + "/{sample}/logs/{ref}.get_damage_global_sample.log",
    params:
        px=OUT_DIR
        + "/{sample}/metadamage/sample_level/{ref}."
        + PREFIX
        + ".filter.global",
    shell:
        """
        ({METADMG_CPP} getdamage -p 25 -l 30 -r 0 {input.bam} -o {params.px}) 2> {log} 1> /dev/null
        {METADMG_CPP} print {output.bdam} -howmany 25 > {output.dam}
        {METADMG_CPP} dfit {output.bdam} --showfits 2 --out {params.px}
        """


rule get_damage_local_sample:
    input:
        bam=TMP_DIR + "/{sample}/{sample}.{ref}.sample.filter.srt_name.bam",
    output:
        dam=OUT_DIR
        + "/{sample}/metadamage/sample_level/{ref}."
        + PREFIX
        + ".filter.local.tsv",
        bdam=OUT_DIR
        + "/{sample}/metadamage/sample_level/{ref}."
        + PREFIX
        + ".filter.local.bdamage.gz",
        dfit=OUT_DIR
        + "/{sample}/metadamage/sample_level/{ref}."
        + PREFIX
        + ".filter.local.dfit.gz",
    log:
        OUT_DIR + "/{sample}/logs/{ref}.get_damage_local_sample.log",
    params:
        px=OUT_DIR
        + "/{sample}/metadamage/sample_level/{ref}."
        + PREFIX
        + ".filter.local",
    shell:
        """
        ({METADMG_CPP} getdamage -p 25 -l 30 -r 1 {input.bam} -o {params.px}) 2> {log} 1> /dev/null
        {METADMG_CPP} print {output.bdam} -howmany 25 -bam {input.bam} > {output.dam}
        {METADMG_CPP} dfit {output.bdam} --showfits 2 --bam {input.bam} --out {params.px}
        """


rule filter_sort_bam_library:
    input:
        bam=TMP_DIR + "/{sample}/{sample}.{lib}." + PREFIX + ".{ref}.mrkdup.bam",
    output:
        bam=temp(TMP_DIR + "/{sample}/{sample}.{lib}.{ref}.lib.filter.srt_name.bam"),
    params:
        mq=config["bam_filter"]["mq"],
    shell:
        """
        samtools view -bh -q {params.mq} -F3844 {input.bam} | samtools sort -n > {output.bam}
        """


rule get_damage_global_library:
    input:
        bam=TMP_DIR + "/{sample}/{sample}.{lib}.{ref}.lib.filter.srt_name.bam",
    output:
        dam=OUT_DIR
        + "/{sample}/metadamage/library_level/{lib}/{ref}."
        + PREFIX
        + ".filter.global.tsv",
        bdam=OUT_DIR
        + "/{sample}/metadamage/library_level/{lib}/{ref}."
        + PREFIX
        + ".filter.global.bdamage.gz",
        dfit=OUT_DIR
        + "/{sample}/metadamage/library_level/{lib}/{ref}."
        + PREFIX
        + ".filter.global.dfit.gz",
    log:
        OUT_DIR + "/{sample}/logs/{lib}.{ref}.get_damage_global_library.log",
    params:
        px=OUT_DIR
        + "/{sample}/metadamage/library_level/{lib}/{ref}."
        + PREFIX
        + ".filter.global",
        lib_type=lambda wildcards: "ss" if wildcards.lib.endswith("-ss") else "ds",
    shell:
        """
        ({METADMG_CPP} getdamage -p 25 -l 30 -r 0 {input.bam} -o {params.px}) 2> {log} 1> /dev/null
        {METADMG_CPP} print {output.bdam} -howmany 25 > {output.dam}
        {METADMG_CPP} dfit {output.bdam} --showfits 2 --lib {params.lib_type} --out {params.px}
        """


rule get_damage_local_library:
    input:
        bam=TMP_DIR + "/{sample}/{sample}.{lib}.{ref}.lib.filter.srt_name.bam",
    output:
        dam=OUT_DIR
        + "/{sample}/metadamage/library_level/{lib}/{ref}."
        + PREFIX
        + ".filter.local.tsv",
        bdam=OUT_DIR
        + "/{sample}/metadamage/library_level/{lib}/{ref}."
        + PREFIX
        + ".filter.local.bdamage.gz",
        dfit=OUT_DIR
        + "/{sample}/metadamage/library_level/{lib}/{ref}."
        + PREFIX
        + ".filter.local.dfit.gz",
    log:
        OUT_DIR + "/{sample}/logs/{lib}.{ref}.get_damage_local_library.log",
    params:
        px=OUT_DIR
        + "/{sample}/metadamage/library_level/{lib}/{ref}."
        + PREFIX
        + ".filter.local",
        lib_type=lambda wildcards: "ss" if wildcards.lib.endswith("-ss") else "ds",
    shell:
        """
        ({METADMG_CPP} getdamage -p 25 -l 30 -r 1 {input.bam} -o {params.px}) 2> {log} 1> /dev/null
        {METADMG_CPP} print {output.bdam} -howmany 25 -bam {input.bam} > {output.dam}
        {METADMG_CPP} dfit {output.bdam} --showfits 2 --bam {input.bam} --lib {params.lib_type} --out {params.px}
        """
