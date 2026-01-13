# map_adna_linear

Snakemake workflow for mapping ancient DNA reads to reference genomes, with
damage estimation and summary statistics. It supports single-end, paired-end,
and collapsed reads.

## Repository layout

- `workflow/Snakefile` - main Snakemake workflow
- `workflow/rules/mapping.smk` - mapping and BAM processing rules
- `workflow/rules/damage.smk` - metaDMG and mapDamage rules
- `workflow/rules/stats.smk` - summary statistics and plots
- `dataset_example/config.yml` - example configuration
- `dataset_example/units.tsv` - example sample sheet
- `run_dataset.sh` - run locally
- `run_dataset_slurm.sh` - submit to SLURM

## Requirements

- `snakemake` (9+ recommended)
- `bowtie2`, `samtools`, `picard`
- `mapDamage`, `metaDMG-cpp`
- R (for plotting/stats scripts)
- R packages: `argparse`, `tidyverse`, `ggrepel`, `viridis`, `scales`, `Rsamtools`, `cowplot`, `cigarillo`, `vroom`
- Python with `pandas`

## Configuration

Create a dataset directory with a `config/` folder that contains:

- `config/config.yml`
- `config/units.tsv`

An example dataset lives at `dataset_example/`.

Example `config/config.yml` (pairs with the `units.tsv` example below):

```yaml
prefix: demo
out_dir: results
units: config/units.tsv
ref:
  ref_genome_a: /path/to/bowtie2_prefix
bt_param: "-D 20 -R 3 -N 1 -L 20 -i S,1,0.50 --end-to-end --no-unal"
bam_filter:
  mq: 20
  qlen: 25
  rlen: 25
tmp_dir: /path/to/tmp/dir
metaDMG-cpp: /path/to/metaDMG-cpp
```

`config/units.tsv` columns:

- `sample_id`: sample identifier used in output directory names
- `lib_id`: library identifier (used in library-level outputs)
- `unit_id`: lane/run identifier
- `fq_se`: single-end or collapsed FASTQ (optional if paired-end only)
- `fq_pe1`: paired-end read 1 FASTQ (optional)
- `fq_pe2`: paired-end read 2 FASTQ (optional)

Example:

| sample_id | lib_id | unit_id | fq_se | fq_pe1 | fq_pe2 |
| --- | --- | --- | --- | --- | --- |
| SAMPLE01 | SAMPLE01_LIB01 | SAMPLE01_LIB01_L001 | /data/fastq/SAMPLE01_LIB01_L001.fq.gz | /data/fastq/SAMPLE01_LIB01_L001_R1.fq.gz | /data/fastq/SAMPLE01_LIB01_L001_R2.fq.gz |
| SAMPLE01 | SAMPLE01_LIB01 | SAMPLE01_LIB01_L002 | /data/fastq/SAMPLE01_LIB01_L002.fq.gz | /data/fastq/SAMPLE01_LIB01_L002_R1.fq.gz | /data/fastq/SAMPLE01_LIB01_L002_R2.fq.gz |
| SAMPLE02 | SAMPLE02_LIB01 | SAMPLE02_LIB01_L001 | /data/fastq/SAMPLE02_LIB01_L001.fq.gz | /data/fastq/SAMPLE02_LIB01_L001_R1.fq.gz | /data/fastq/SAMPLE02_LIB01_L001_R2.fq.gz |

## Run locally

From the repo root:

```bash
bash run_dataset.sh /path/to/dataset
bash run_dataset.sh dataset_example
```

You can pass extra Snakemake args after the dataset path:

```bash
bash run_dataset.sh /path/to/dataset -- --cores 8 --rerun-incomplete
```

## Run on SLURM

```bash
bash run_dataset_slurm.sh /path/to/dataset --jobs 50 --partition general
bash run_dataset_slurm.sh dataset_example --dry-run
```

See full options:

```bash
bash run_dataset_slurm.sh --help
```

## Outputs

For each `sample_id`, outputs are written under `${out_dir}/{sample_id}/`:

- `bam/` - merged and duplicate-marked BAMs
- `metadamage/` - sample- and library-level damage estimates
- `mapdamage/` - mapDamage plots
- `tables/` - per-sample statistics tables
- `plots/` - per-sample summary plots (reference similarity, read length distribution and aDNA damage)
- `logs/` - mapping and damage logs

Summary outputs across all samples are written under `${out_dir}/summary/`.

Example output structure:

```text
${out_dir}/
├── SAMPLE01/
│   ├── bam/
│   │   ├── SAMPLE01.SAMPLE01_LIB01.SAMPLE01_LIB01_L001.demo.ref_genome_a.init.bam
│   │   ├── SAMPLE01.SAMPLE01_LIB01.SAMPLE01_LIB01_L002.demo.ref_genome_a.init.bam
│   │   ├── SAMPLE01.demo.ref_genome_a.mrkdup.bam
│   │   ├── SAMPLE01.demo.ref_genome_a.mrkdup.bam.bai
│   │   ├── SAMPLE01.demo.ref_genome_a.filter.bam
│   │   └── SAMPLE01.demo.ref_genome_a.filter.bam.bai
│   ├── logs/
│   │   ├── SAMPLE01.SAMPLE01_LIB01.SAMPLE01_LIB01_L001.demo.ref_genome_a.map_init.log
│   │   ├── SAMPLE01.SAMPLE01_LIB01.SAMPLE01_LIB01_L002.demo.ref_genome_a.map_init.log
│   │   ├── SAMPLE01.SAMPLE01_LIB01.demo.ref_genome_a.mrkdup_metrics.txt
│   │   ├── ref_genome_a.get_damage_global_sample.log
│   │   ├── ref_genome_a.get_damage_local_sample.log
│   │   ├── SAMPLE01_LIB01.ref_genome_a.get_damage_global_library.log
│   │   └── SAMPLE01_LIB01.ref_genome_a.get_damage_local_library.log
│   ├── metadamage/
│   │   ├── sample_level/
│   │   │   ├── ref_genome_a.demo.filter.global.tsv
│   │   │   ├── ref_genome_a.demo.filter.global.bdamage.gz
│   │   │   ├── ref_genome_a.demo.filter.global.dfit.gz
│   │   │   ├── ref_genome_a.demo.filter.local.tsv
│   │   │   ├── ref_genome_a.demo.filter.local.bdamage.gz
│   │   │   └── ref_genome_a.demo.filter.local.dfit.gz
│   │   └── library_level/
│   │       └── SAMPLE01_LIB01/
│   │           ├── ref_genome_a.demo.filter.global.tsv
│   │           ├── ref_genome_a.demo.filter.global.bdamage.gz
│   │           ├── ref_genome_a.demo.filter.global.dfit.gz
│   │           ├── ref_genome_a.demo.filter.local.tsv
│   │           ├── ref_genome_a.demo.filter.local.bdamage.gz
│   │           └── ref_genome_a.demo.filter.local.dfit.gz
│   ├── mapdamage/
│   │   └── SAMPLE01.demo.ref_genome_a.filter/
│   │       ├── Fragmisincorporation_plot.pdf
│   │       └── Length_plot.pdf
│   ├── tables/
│   │   ├── SAMPLE01.demo.ref_genome_a.filter.bam_stats.tsv
│   │   └── SAMPLE01.demo.ref_genome_a.filter.lib_stats.tsv
│   └── plots/
│       └── SAMPLE01.demo.ref_genome_a.filter.summary.pdf
└── summary/
    ├── tables/
    │   ├── demo.ref_genome_a.filter.bam_stats.tsv
    │   └── demo.ref_genome_a.filter.lib_stats.tsv
    └── plots/
        ├── demo.ref_genome_a.filter.bam_stats_coverage.pdf
        ├── demo.ref_genome_a.filter.lib_stats_n_map_uniq_p_dup.pdf
        ├── demo.ref_genome_a.filter.lib_stats_n_reads_n_map_uniq.pdf
        ├── demo.ref_genome_a.filter.lib_stats_n_reads_p_dup.pdf
        └── demo.ref_genome_a.filter.lib_stats_p_dup.pdf
```

## Output tables

`bam_stats.tsv` (`${out_dir}/{sample}/tables/{sample}.{PREFIX}.{ref}.filter.bam_stats.tsv`) columns:

- `sample_id`: sample identifier
- `library_id`: library identifier (`all` indicates sample-level aggregate)
- `contig_id`: reference contig ID (`genome` indicates genome-wide aggregate)
- `contig_l`: contig length
- `n_reads`: number of reads mapped to the contig
- `dp_avg`: average depth of coverage
- `dp_sd`: depth standard deviation
- `coverage_bp`: covered bases (depth > 0)
- `coverage_p`: fraction of covered bases
- `coverage_p_exp`: expected covered fraction (Poisson)
- `coverage_p_ratio`: observed/expected coverage fraction
- `dp_cv`: depth coefficient of variation
- `coverage_evenness_score`: evenness metric for coverage
- `cov_pos_rel_entropy_100`: relative coverage position entropy (100 bp bins)
- `cov_pos_rel_entropy_1000`: relative coverage position entropy (1000 bp bins)
- `read_l_avg`: average read length
- `mq_avg`: average mapping quality
- `n_soft_clip_avg`: average soft-clipped bases per read
- `ani`: average nucleotide identity proxy (1 - mismatches/read length)
- `edit_dist_mode`: modal edit distance (NM)
- `edit_dist_avg`: mean edit distance
- `edit_dist_avg_decay`: mean decay of edit distance histogram
- `edit_dist_decay_end`: edit distance bin where counts begin to increase
- `dmg_f_5p`: 5' damage frequency (metaDMG)
- `dmg_f_3p`: 3' damage frequency (metaDMG)
- `dmg_amp`: damage amplitude (metaDMG A)
- `dmg_decay`: damage decay (metaDMG q)
- `dmg_noise`: damage noise (metaDMG c)
- `dmg_Zfit`: metaDMG Z-fit statistic
- `dmg_ntot_5p`: total 5' counts used for damage fit
- `dmg_ntot_3p`: total 3' counts used for damage fit

`lib_stats.tsv` (`${out_dir}/{sample}/tables/{sample}.{PREFIX}.{ref}.filter.lib_stats.tsv`) columns:

- `sample_id`: sample identifier
- `library_id`: library identifier
- `num_reads_total`: total reads from fastq logs (paired reads counted as two)
- `num_reads_mapped`: reads examined after mapping
- `prop_mapped`: fraction mapped (`num_reads_mapped/num_reads_total`)
- `num_reads_dup`: duplicate reads
- `prop_dup`: fraction duplicates (`num_reads_dup/num_reads_mapped`)
- `num_reads_final`: final mapped reads after filtering
- `prop_final`: fraction retained (`num_reads_final/num_reads_total`)
