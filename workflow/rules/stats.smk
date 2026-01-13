"""
Statistics rules: BAM statistics, library statistics, and summary reports
"""


def get_dmg_files(wildcards):
    """Get all damage analysis files (sample-level and library-level) for a given sample and reference."""
    files = []

    ## sample-level outputs
    files.extend(
        expand(
            OUT_DIR
            + "/{sample}/metadamage/sample_level/{ref}."
            + PREFIX
            + ".filter.{ext}",
            sample=wildcards.sample,
            ref=wildcards.ref,
            ext=[
                "global.dfit.gz",
                "local.dfit.gz",
            ],
        )
    )

    ## library-level outputs
    libraries = unit_df.loc[(wildcards.sample), "lib_id"].unique()
    for lib in libraries:
        files.extend(
            expand(
                OUT_DIR
                + "/{sample}/metadamage/library_level/{lib}/{ref}."
                + PREFIX
                + ".filter.{ext}",
                sample=wildcards.sample,
                lib=lib,
                ref=wildcards.ref,
                ext=[
                    "global.dfit.gz",
                    "local.dfit.gz",
                ],
            )
        )
    return sorted(set(files))


rule get_bam_stats:
    """Calculate BAM statistics for each sample."""
    input:
        bam=OUT_DIR + "/{sample}/bam/{sample}." + PREFIX + ".{ref}.filter.bam",
        dmg=get_dmg_files,
    output:
        tsv=OUT_DIR
        + "/{sample}/tables/{sample}."
        + PREFIX
        + ".{ref}.filter.bam_stats.tsv",
    params:
        mq=config["bam_filter"]["mq"],
    shell:
        """
        Rscript {workflow.basedir}/scripts/bam_stats.R -b {input.bam} -o {output.tsv} --mq {params.mq} {input.dmg}
        """


rule get_lib_stats:
    """Calculate library-level statistics."""
    input:
        OUT_DIR + "/{sample}/tables/{sample}." + PREFIX + ".{ref}.filter.bam_stats.tsv",
    output:
        tsv=OUT_DIR
        + "/{sample}/tables/{sample}."
        + PREFIX
        + ".{ref}.filter.lib_stats.tsv",
    shell:
        """
        Rscript {workflow.basedir}/scripts/lib_stats.R -s {wildcards.sample} -r {wildcards.ref} -p {PREFIX} -o {OUT_DIR}
        """


rule plot_sample_stats:
    """Plot sample statistics"""
    input:
        bam=OUT_DIR + "/{sample}/bam/{sample}." + PREFIX + ".{ref}.filter.bam",
        stats=OUT_DIR
        + "/{sample}/tables/{sample}."
        + PREFIX
        + ".{ref}.filter.bam_stats.tsv",
    output:
        pdf=OUT_DIR + "/{sample}/plots/{sample}." + PREFIX + ".{ref}.filter.summary.pdf",
    shell:
        """
        Rscript {workflow.basedir}/scripts/plot_sample.R -b {input.bam} -s {wildcards.sample} -r {wildcards.ref} -p {PREFIX} -o {OUT_DIR}
        """


rule summary_stats:
    """Generate summary statistics and plots across all samples."""
    input:
        bam_stats=expand(
            OUT_DIR
            + "/{sample}/tables/{sample}."
            + PREFIX
            + ".{{ref}}.filter.bam_stats.tsv",
            sample=SAMPLES,
        ),
        lib_stats=expand(
            OUT_DIR
            + "/{sample}/tables/{sample}."
            + PREFIX
            + ".{{ref}}.filter.lib_stats.tsv",
            sample=SAMPLES,
        ),
    output:
        tsv_bam=OUT_DIR + "/summary/tables/" + PREFIX + ".{ref}.filter.bam_stats.tsv",
        tsv_lib=OUT_DIR + "/summary/tables/" + PREFIX + ".{ref}.filter.lib_stats.tsv",
        pdf_bam=OUT_DIR
        + "/summary/plots/"
        + PREFIX
        + ".{ref}.filter.bam_stats_coverage.pdf",
        pdf_lib=expand(
            OUT_DIR + "/summary/plots/" + PREFIX + ".{{ref}}.filter.{ext}",
            ext=[
                "lib_stats_p_dup.pdf",
                "lib_stats_n_reads_p_dup.pdf",
                "lib_stats_n_map_uniq_p_dup.pdf",
                "lib_stats_n_reads_n_map_uniq.pdf",
                "lib_stats_n_reads_n_map_uniq.pdf",
            ],
        ),
    params:
        px_tsv=OUT_DIR + "/summary/tables/" + PREFIX + ".{ref}.filter",
        px_pdf=OUT_DIR + "/summary/plots/" + PREFIX + ".{ref}.filter",
    shell:
        """
        Rscript {workflow.basedir}/scripts/summarise_stats.R --prefix_out_plots {params.px_pdf} --prefix_out_tables {params.px_tsv} {input.bam_stats}
        """
