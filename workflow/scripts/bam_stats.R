#!/usr/bin/env Rscript
# Optimized bam statistics script with chunked BAM parsing and cached entropy
# generation. Produces the same outputs as bam_stats.R but minimizes repeated
# work and reduces peak memory usage.

suppressPackageStartupMessages({
  library(argparse)
  library(tidyverse)
  library(viridis)
  library(Rsamtools)
  library(cigarillo)
  library(vroom)
})

options(dplyr.summarise.inform = FALSE)

parser <- ArgumentParser()
parser$add_argument("-b", "--bam_file",
  action = "store",
  dest = "bam_file",
  help = "Input BAM file"
)
parser$add_argument("-o", "--out_file",
  action = "store",
  dest = "out_file",
  help = "Output filename"
)

parser$add_argument("--mq",
  action = "store",
  dest = "mq",
  type = "double",
  default = 20,
  help = "MQ cutoff for filtered stats [default %(default)s]"
)
parser$add_argument("--chunk-size",
  action = "store",
  dest = "chunk_size",
  type = "integer",
  default = 500000,
  help = "Number of alignments to process per scanBam chunk [default %(default)s]"
)
parser$add_argument("metadmg_files",
  nargs = "+",
  help = "Files with metadamage dfit data"
)
args <- parser$parse_args()

hdr <- scanBamHeader(args$bam_file)
genome <- tibble(
  contig_id = names(hdr[[1]]$targets),
  contig_l = hdr[[1]]$targets
)
genome_l <- sum(genome$contig_l)

rg_idx <- which(names(hdr[[1]]$text) == "@RG")
rg_info <- map_dfr(hdr[[1]]$text[rg_idx], ~ {
  tibble(
    rg_id = gsub("ID:", "", .x[grepl("^ID", .x)]),
    library_id = gsub("LB:", "", .x[grepl("^LB", .x)]),
    sample_id = gsub("SM:", "", .x[grepl("^SM", .x)])
  )
})

if (!nrow(rg_info)) {
  stop("No read groups detected in BAM header")
}

entropy_max_even <- function(total_reads, n_bins) {
  if (is.na(total_reads) || is.na(n_bins) || total_reads <= 0 || n_bins <= 0) {
    return(NA_real_)
  }
  q <- total_reads %/% n_bins
  r <- total_reads %% n_bins
  p1 <- if (q == 0) 0 else q / total_reads
  p2 <- (q + 1) / total_reads
  term1 <- if (p1 > 0) (n_bins - r) * p1 * log(p1, base = 2) else 0
  term2 <- if (r > 0 && p2 > 0) r * p2 * log(p2, base = 2) else 0
  -(term1 + term2)
}

extract_soft_clip <- function(cigar_vec) {
  if (!length(cigar_vec)) {
    return(integer(0))
  }
  tab <- cigarillo::tabulate_cigar_ops(cigar_vec)
  if (!"S" %in% colnames(tab)) {
    return(integer(nrow(tab)))
  }
  as.integer(tab[, "S", drop = TRUE])
}

coverage_stats <- function(cov_tbl, sample_id, library_id) {
  if (!nrow(cov_tbl)) {
    return(tibble())
  }
  cov_tbl |>
    group_by(contig_id) |>
    summarise(
      contig_l = contig_l[1],
      dp_avg = sum(dp * n) / contig_l,
      dp_sd = sqrt(sum((dp - dp_avg)^2 * n) / sum(n)),
      coverage_bp = sum(n[dp > 0]),
      coverage_p = coverage_bp / contig_l,
      coverage_p_exp = 1 - exp(-dp_avg),
      coverage_p_ratio = coverage_p / coverage_p_exp,
      dp_cv = dp_sd / dp_avg,
      coverage_evenness_score = 1 - sum((ceiling(dp_avg) - dp[dp <= ceiling(dp_avg)]) * n[dp <= ceiling(dp_avg)] / (ceiling(dp_avg) * contig_l)),
      .groups = "drop"
    ) |>
    mutate(sample_id = sample_id, library_id = library_id)
}

run_genomecov <- function(rg_ids, sample_id, library_id) {
  if (!length(rg_ids)) {
    return(tibble())
  }
  rg_flag <- paste(paste("-r", shQuote(rg_ids)), collapse = " ")
  cmd <- sprintf(
    "samtools view -bh -F 1024 %s -q %s %s | bedtools genomecov -ibam stdin",
    rg_flag, args$mq, shQuote(args$bam_file)
  )
  r1 <- read_tsv(pipe(cmd),
    col_types = cols(
      contig_id = col_character(),
      dp = col_double(),
      n = col_double(),
      contig_l = col_double(),
      p = col_double()
    ),
    col_names = c("contig_id", "dp", "n", "contig_l", "p"),
    show_col_types = FALSE
  )
  coverage_stats(r1, sample_id, library_id)
}

summarise_nm_table <- function(df) {
  if (!nrow(df)) {
    return(tibble(
      contig_id = character(),
      library_id = character(),
      sample_id = character(),
      n_reads = double(),
      edit_dist_mode = double(),
      edit_dist_avg = double(),
      edit_dist_avg_decay = double(),
      edit_dist_decay_end = double()
    ))
  }
  df |>
    arrange(contig_id, library_id, sample_id, nm) |>
    group_by(contig_id, library_id, sample_id) |>
    summarise(
      n_reads = sum(n),
      edit_dist_mode = nm[which.max(n)],
      edit_dist_avg = sum(nm * n) / n_reads,
      edit_dist_avg_decay = {
        diffs <- diff(n)
        if (length(diffs)) mean(diffs) / n_reads else 0
      },
      edit_dist_decay_end = {
        diffs <- diff(n)
        idx <- which(diffs > 0)[1]
        ifelse(is.na(idx), max(nm), idx)
      },
      .groups = "drop"
    )
}

finalize_read_stats <- function(df) {
  if (!nrow(df)) {
    return(tibble(
      contig_id = character(),
      library_id = character(),
      sample_id = character(),
      read_l_avg = double(),
      mq_avg = double(),
      n_soft_clip_avg = double(),
      ani = double()
    ))
  }
  df |>
    mutate(
      read_l_avg = if_else(n_reads > 0, sum_qwidth / n_reads, NA_real_),
      mq_avg = if_else(n_reads > 0, sum_mapq / n_reads, NA_real_),
      n_soft_clip_avg = if_else(n_reads > 0, sum_soft_clip / n_reads, NA_real_),
      ani = if_else(sum_qwidth > 0, 1 - (sum_nm / sum_qwidth), NA_real_)
    ) |>
    select(contig_id, library_id, sample_id, read_l_avg, mq_avg, n_soft_clip_avg, ani)
}

compute_entropy <- function(counts_tbl, bin_size, genome_bins_map, genome_bins_total) {
  if (!nrow(counts_tbl)) {
    return(tibble(
      contig_id = character(),
      library_id = character(),
      sample_id = character(),
      bin_size = integer(),
      cov_pos_rel_entropy = double()
    ))
  }
  counts_tbl |>
    group_by(contig_id, library_id, sample_id) |>
    summarise(
      entropy_obs = {
        probs <- n / sum(n)
        -sum(probs * log(probs, base = 2))
      },
      total_reads = sum(n),
      .groups = "drop"
    ) |>
    mutate(
      n_bins = if_else(
        contig_id == "genome",
        genome_bins_total,
        genome_bins_map[contig_id]
      ),
      entropy_max = map2_dbl(total_reads, n_bins, entropy_max_even),
      cov_pos_rel_entropy = if_else(
        is.na(entropy_max) | entropy_max == 0,
        NA_real_,
        abs(entropy_obs / entropy_max)
      ),
      bin_size = bin_size
    ) |>
    select(contig_id, library_id, sample_id, bin_size, cov_pos_rel_entropy)
}

safe_bind_rows <- function(lst) {
  if (!length(lst)) {
    return(tibble())
  }
  bind_rows(lst)
}

cov_scopes_lib <- rg_info |>
  group_by(sample_id, library_id) |>
  summarise(rg_ids = list(rg_id), .groups = "drop")

cov_scopes_sample <- rg_info |>
  group_by(sample_id) |>
  summarise(rg_ids = list(rg_id), .groups = "drop") |>
  mutate(library_id = "all")

cov_scopes <- bind_rows(cov_scopes_lib, cov_scopes_sample)

r_cov <- cov_scopes |>
  mutate(stats = pmap(list(rg_ids, sample_id, library_id), run_genomecov)) |>
  select(stats) |>
  unnest(stats)

bin_sizes <- c(100L, 1000L)
entropy_hist <- setNames(
  vector("list", length(bin_sizes)),
  as.character(bin_sizes)
)
nm_histories <- list()
stats_histories <- list()
has_reads <- FALSE

param <- ScanBamParam(
  what = c("rname", "pos", "mapq", "qwidth", "cigar"),
  tag = c("NM", "RG"),
  flag = scanBamFlag(
    isDuplicate = FALSE,
    isSupplementaryAlignment = FALSE,
    isSecondaryAlignment = FALSE
  ),
  mapqFilter = args$mq
)

bf <- BamFile(args$bam_file, yieldSize = args$chunk_size)
open(bf)
repeat {
  chunk <- scanBam(bf, param = param)[[1]]
  if (length(chunk[["pos"]]) == 0) {
    break
  }
  chunk_tbl <- tibble(
    contig_id = as.character(chunk[["rname"]]),
    pos = as.integer(chunk[["pos"]]),
    mapq = as.numeric(chunk[["mapq"]]),
    qwidth = as.integer(chunk[["qwidth"]]),
    cigar = chunk[["cigar"]],
    nm = as.integer(chunk$tag$NM),
    rg_id = chunk$tag$RG
  ) |>
    left_join(rg_info, by = "rg_id") |>
    filter(!is.na(library_id))

  if (!nrow(chunk_tbl)) {
    next
  }
  has_reads <- TRUE
  chunk_tbl <- chunk_tbl |>
    mutate(n_soft_clip = extract_soft_clip(cigar))

  nm_histories[[length(nm_histories) + 1]] <- chunk_tbl |>
    group_by(contig_id, library_id, sample_id, nm) |>
    summarise(n = n(), .groups = "drop")

  stats_histories[[length(stats_histories) + 1]] <- chunk_tbl |>
    group_by(contig_id, library_id, sample_id) |>
    summarise(
      n_reads = n(),
      sum_qwidth = sum(qwidth),
      sum_mapq = sum(mapq),
      sum_soft_clip = sum(n_soft_clip),
      sum_nm = sum(nm),
      .groups = "drop"
    )

  for (bs in bin_sizes) {
    bin_key <- as.character(bs)
    entropy_hist[[bin_key]] <- bind_rows(
      entropy_hist[[bin_key]],
      chunk_tbl |>
        mutate(bp_bin = as.character(pos %/% bs)) |>
        group_by(contig_id, library_id, sample_id, bp_bin) |>
        summarise(n = n(), .groups = "drop")
    )
  }
}
close(bf)

if (!has_reads) {
  nm_counts <- tibble(
    contig_id = character(),
    library_id = character(),
    sample_id = character(),
    nm = integer(),
    n = integer()
  )
  stats_counts <- tibble(
    contig_id = character(),
    library_id = character(),
    sample_id = character(),
    n_reads = double(),
    sum_qwidth = double(),
    sum_mapq = double(),
    sum_soft_clip = double(),
    sum_nm = double()
  )
  entropy_hist <- setNames(map(bin_sizes, ~ tibble(
    contig_id = character(),
    library_id = character(),
    sample_id = character(),
    bp_bin = character(),
    n = double()
  )), as.character(bin_sizes))
} else {
  nm_counts <- safe_bind_rows(nm_histories) |>
    group_by(contig_id, library_id, sample_id, nm) |>
    summarise(n = sum(n), .groups = "drop")

  stats_counts <- safe_bind_rows(stats_histories) |>
    group_by(contig_id, library_id, sample_id) |>
    summarise(
      n_reads = sum(n_reads),
      sum_qwidth = sum(sum_qwidth),
      sum_mapq = sum(sum_mapq),
      sum_soft_clip = sum(sum_soft_clip),
      sum_nm = sum(sum_nm),
      .groups = "drop"
    )
}

nm_counts_sample <- nm_counts |>
  group_by(contig_id, sample_id, nm) |>
  summarise(n = sum(n), .groups = "drop") |>
  mutate(library_id = "all")

nm_counts_genome_lib <- nm_counts |>
  group_by(library_id, sample_id, nm) |>
  summarise(n = sum(n), .groups = "drop") |>
  mutate(contig_id = "genome")

nm_counts_genome_sample <- nm_counts_sample |>
  group_by(sample_id, nm) |>
  summarise(n = sum(n), .groups = "drop") |>
  mutate(contig_id = "genome", library_id = "all")

nm_all <- bind_rows(
  nm_counts,
  nm_counts_sample,
  nm_counts_genome_lib,
  nm_counts_genome_sample
)

r_nm <- summarise_nm_table(nm_all)

stats_counts_sample <- stats_counts |>
  group_by(contig_id, sample_id) |>
  summarise(
    n_reads = sum(n_reads),
    sum_qwidth = sum(sum_qwidth),
    sum_mapq = sum(sum_mapq),
    sum_soft_clip = sum(sum_soft_clip),
    sum_nm = sum(sum_nm),
    .groups = "drop"
  ) |>
  mutate(library_id = "all")

stats_counts_genome_lib <- stats_counts |>
  group_by(library_id, sample_id) |>
  summarise(
    n_reads = sum(n_reads),
    sum_qwidth = sum(sum_qwidth),
    sum_mapq = sum(sum_mapq),
    sum_soft_clip = sum(sum_soft_clip),
    sum_nm = sum(sum_nm),
    .groups = "drop"
  ) |>
  mutate(contig_id = "genome")

stats_counts_genome_sample <- stats_counts_genome_lib |>
  group_by(sample_id) |>
  summarise(
    n_reads = sum(n_reads),
    sum_qwidth = sum(sum_qwidth),
    sum_mapq = sum(sum_mapq),
    sum_soft_clip = sum(sum_soft_clip),
    sum_nm = sum(sum_nm),
    .groups = "drop"
  ) |>
  mutate(contig_id = "genome", library_id = "all")

r_stats <- bind_rows(
  finalize_read_stats(stats_counts),
  finalize_read_stats(stats_counts_sample),
  finalize_read_stats(stats_counts_genome_lib),
  finalize_read_stats(stats_counts_genome_sample)
)

entropy_tables <- list()
for (bs in bin_sizes) {
  bin_key <- as.character(bs)
  counts_tbl <- entropy_hist[[bin_key]] |>
    group_by(contig_id, library_id, sample_id, bp_bin) |>
    summarise(n = sum(n), .groups = "drop")
  counts_sample <- counts_tbl |>
    group_by(contig_id, sample_id, bp_bin) |>
    summarise(n = sum(n), .groups = "drop") |>
    mutate(library_id = "all")
  counts_genome_lib <- counts_tbl |>
    group_by(library_id, sample_id, bp_bin) |>
    summarise(n = sum(n), .groups = "drop") |>
    mutate(contig_id = "genome")
  counts_genome_sample <- counts_sample |>
    group_by(sample_id, bp_bin) |>
    summarise(n = sum(n), .groups = "drop") |>
    mutate(contig_id = "genome", library_id = "all")
  combined_counts <- bind_rows(
    counts_tbl,
    counts_sample,
    counts_genome_lib,
    counts_genome_sample
  )
  genome_bins_map <- setNames(
    pmax(1, genome$contig_l %/% bs + 1),
    genome$contig_id
  )
  genome_bins_total <- sum(genome_bins_map)
  entropy_tables[[bin_key]] <- compute_entropy(
    combined_counts, bs, genome_bins_map, genome_bins_total
  )
}

r_entropy <- bind_rows(entropy_tables) |>
  pivot_wider(
    names_from = bin_size,
    values_from = cov_pos_rel_entropy,
    names_prefix = "cov_pos_rel_entropy_"
  )
entropy_cols <- paste0("cov_pos_rel_entropy_", bin_sizes)
if (!nrow(r_entropy)) {
  r_entropy <- tibble(
    contig_id = character(),
    library_id = character(),
    sample_id = character()
  )
}
for (col in entropy_cols) {
  if (!col %in% names(r_entropy)) {
    r_entropy[[col]] <- NA_real_
  }
}

read_metadmg <- function(path) {
  dat <- vroom::vroom(path,
    delim = "\t",
    show_col_types = FALSE,
    progress = FALSE
  )
  if (!nrow(dat)) {
    return(tibble())
  }
  dat |>
    mutate(
      contig_id = case_when(
        taxid == 0 ~ "genome",
        TRUE ~ as.character(taxid)
      )
    )
}

metadmg_tbl <- map_dfr(args$metadmg_files, ~ {
  r1 <- read_metadmg(.x)
  f_str <- strsplit(.x, "/")[[1]]
  lib_idx <- match("library_level", f_str)
  library_id <- ifelse(
    !is.na(lib_idx) && length(f_str) >= lib_idx + 1,
    f_str[lib_idx + 1],
    "all"
  )
  if (!nrow(r1)) {
    return(tibble(
      library_id = character(),
      contig_id = character(),
      dmg_f_5p = double(),
      dmg_f_3p = double(),
      dmg_amp = double(),
      dmg_decay = double(),
      dmg_noise = double(),
      dmg_Zfit = double(),
      dmg_ntot_5p = double(),
      dmg_ntot_3p = double()
    ))
  }
  dmg_totals <- r1 |>
    select(
      contig_id, starts_with("fwN"), starts_with("bwN"),
      starts_with("fwK"), starts_with("bwK")
    ) |>
    pivot_longer(
      cols = -contig_id,
      names_to = c("strand", ".value"),
      names_pattern = "(fw|bw)([A-Z]+)"
    ) |>
    group_by(contig_id, strand) |>
    summarise(dmg_ntot = sum(N) + sum(K), .groups = "drop") |>
    mutate(strand = if_else(strand == "fw", "5p", "3p")) |>
    pivot_wider(
      names_from = strand,
      values_from = dmg_ntot, 
      names_prefix = "dmg_ntot_"
    )
  r1 |>
    transmute(
      library_id = library_id,
      contig_id = contig_id,
      dmg_f_5p = fwf0,
      dmg_f_3p = bwf0,
      dmg_amp = A,
      dmg_decay = q,
      dmg_noise = c,
      dmg_Zfit = Zfit
    ) |>
    left_join(dmg_totals, by = "contig_id") |>
    select(library_id, contig_id, dmg_f_5p:dmg_Zfit, dmg_ntot_5p, dmg_ntot_3p)
})


key_cols <- c("sample_id", "library_id", "contig_id")
row_keys <- list(r_cov, r_entropy, r_nm, r_stats) |>
  discard(~ !nrow(.x)) |>
  map(~ select(.x, any_of(key_cols))) |>
  bind_rows() |>
  distinct()

res <- row_keys |>
  left_join(r_cov, by = key_cols) |>
  left_join(r_entropy, by = key_cols) |>
  left_join(r_nm, by = key_cols) |>
  left_join(r_stats, by = key_cols) |>
  left_join(metadmg_tbl, by = c("library_id", "contig_id")) |>
  select(
    sample_id, library_id, contig_id, contig_l,
    n_reads, dp_avg:cov_pos_rel_entropy_1000,
    ani, edit_dist_mode:edit_dist_decay_end, everything()
  ) |>
  mutate(across(where(is.numeric), ~ round(.x, 5)))

write_tsv(res, file = args$out_file, na = "NaN")
