#!/usr/bin/env Rscript
# Copyright 2026 Martin Sikora <martin.sikora@sund.ku.dk>
#
#  This file is free software: you may copy, redistribute and/or modify it
#  under the terms of the GNU General Public License as published by the
#  Free Software Foundation, either version 2 of the License, or (at your
#  option) any later version.
#
#  This file is distributed in the hope that it will be useful, but
#  WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#  General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <http://www.gnu.org/licenses/>.

# plot metadmg and read length distributions at library and sample level

## ---------------------------------------------------------
## libraries

suppressPackageStartupMessages({
  library(argparse)
  library(tidyverse)
  library(viridis)
  library(Rsamtools)
  library(ggrepel)
  library(cowplot)
  library(scales)
})


## --------------------------------------------------
## command line argument setup and parsing

parser <- ArgumentParser()

parser$add_argument("-b", "--bam_file",
  action = "store",
  dest = "bam_file",
  help = "Input BAM file"
)

parser$add_argument("-s", "--sample_id",
  action = "store",
  dest = "sample_id",
  help = "Sample ID to process"
)

parser$add_argument("-p", "--prefix",
  action = "store",
  dest = "prefix",
  help = "Prefix for run"
)

parser$add_argument("-r", "--ref",
  action = "store",
  dest = "ref",
  help = "Reference name"
)

parser$add_argument("-o", "--out_dir",
  action = "store",
  dest = "out_dir",
  help = "Output directory"
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

args <- parser$parse_args()

th <- theme_bw() +
  theme(
    panel.grid.major = element_line(
      linetype = "dotted",
      linewidth = 0.25
    ),
    panel.grid.minor = element_blank(),
    strip.background = element_blank(),
    strip.text.y = element_text(angle = 0, hjust = 0)
  )


## ---------------------------------------------------------
## alignment data

hdr <- scanBamHeader(args$bam_file)

rg_idx <- which(names(hdr[[1]]$text) == "@RG")
rg_info <- map_dfr(hdr[[1]]$text[rg_idx], ~ {
  tibble(
    rg_id = gsub("ID:", "", .x[grepl("^ID", .x)]),
    library_id = gsub("LB:", "", .x[grepl("^LB", .x)]),
    sample_id = gsub("SM:", "", .x[grepl("^SM", .x)])
  )
})

# Filter to only read groups with library_id
valid_rg <- rg_info |>
  filter(!is.na(library_id)) |>
  pull(rg_id)

param <- ScanBamParam(
  what = "qwidth",
  tag = c("RG", "NM"),
  flag = scanBamFlag(
    isDuplicate = FALSE,
    isSupplementaryAlignment = FALSE,
    isSecondaryAlignment = FALSE
  ),
  mapqFilter = args$mq
)

# get read length data
bam_data <- scanBam(args$bam_file, param = param)[[1]]

bam_tbl <- if (length(bam_data[["qwidth"]]) > 0) {
  rg_id <- bam_data$tag$RG
  keep <- rg_id %in% valid_rg

  if (any(keep)) {
    tibble(
      qwidth = as.integer(bam_data[["qwidth"]][keep]),
      nm = as.integer(bam_data$tag$NM[keep]),
      rg_id = rg_id[keep]
    ) |>
      left_join(rg_info, by = "rg_id")
  } else {
    tibble()
  }
} else {
  tibble()
}

## ---------------------------------------------------------
## metadamage data

infiles <- list.files(
  path = file.path(args$out_dir, args$sample_id, "metadamage"),
  pattern = "global.dfit.gz$",
  full.names = TRUE,
  recursive = TRUE
)

dmg_tbl <- map_dfr(infiles, ~ {
  r <- read_tsv(.x,
    col_types = cols(.default = col_double()),
    show_col_types = FALSE
  )
  run_level <- str_split_i(.x, "/", 4)
  library_id <- ifelse(
    run_level == "library_level",
    str_split_i(.x, "/", 5),
    "all"
  )

  r1 <- r |>
    select(
      starts_with("fwdx"), starts_with("fwf")
    ) |>
    pivot_longer(
      cols = starts_with(c("fwdx", "fwf")),
      names_to = "variable",
      values_to = "value"
    ) |>
    separate(
      col = variable,
      into = c("variable", "position"),
      sep = "(?<=\\D)(?=\\d)",
      remove = FALSE,
      extra = "merge",
      fill = "right"
    ) |>
    pivot_wider(names_from = variable, values_from = value) |>
    mutate(position = as.integer(position)) |>
    arrange(position) |>
    mutate(
      library_id = library_id,
    )
  r1
})

## ---------------------------------------------------------
## summary stats

stats_file <- file.path(
  args$out_dir, args$sample_id, "tables",
  paste0(args$sample_id, ".", args$prefix, ".", args$ref, ".filter.bam_stats.tsv")
)

stats_tbl <- read_tsv(stats_file) |>
  filter(
    contig_id == "genome"
  )


## ---------------------------------------------------------
## plot

pdf(
  file = file.path(
    args$out_dir, args$sample_id, "plots",
    paste0(args$sample_id, ".", args$prefix, ".", args$ref, ".filter.summary.pdf")
  ),
  width = 12,
  height = 3.5,
)

walk(unique(stats_tbl$library_id), ~ {
  p <- list()

  ## summry stats and alignent stats for library
  d_stats <- stats_tbl |>
    filter(library_id == .x) |>
    select(
      n_reads,
      dp_avg,
      ani,
      cov_pos_rel_entropy_1000,
      read_l_avg,
      dmg_amp,
      dmg_Zfit
    )

  if (.x == "all") {
    d_bam <- bam_tbl
  } else {
    d_bam <- bam_tbl |>
      filter(library_id == .x)
  }

  ## read edit distance distribution
  d_nm <- d_bam |>
    count(nm) |>
    mutate(
      p = n / sum(n),
    )
  lab <- paste0(
    "Total reads: ", comma(d_stats$n_reads),
    "\nAvg read depth: ", round(d_stats$dp_avg, 3),
    "\nANI: ", round(d_stats$ani, 3),
    "\nCoverage entropy: ", round(d_stats$cov_pos_rel_entropy_1000, 3)
  )
  d_lab <- tibble(
    x = 10,
    y = max(d_nm$p),
    label = lab
  )
  p[[1]] <- d_nm |>
    ggplot(aes(
      x = nm,
      y = p
    )) +
    geom_line(
      linewidth = 0.25,
    ) +
    geom_point(
      aes(size = n),
      shape = 22,
      color = "navy",
      fill = "navy",
    ) +
    geom_text(
      aes(
        x = x,
        y = y,
        label = label
      ),
      hjust = 1,
      vjust = 1,
      size = 3,
      data = d_lab,
    ) +
    scale_x_continuous(
      limits = c(0, 10),
      breaks = 0:10
    ) +
    scale_size_continuous(
      range = c(0.5, 3),
    ) +
    xlab("Read edit distance") +
    ylab("Proportion of reads") +
    ggtitle(paste0(
      "Sample: ", args$sample_id,
      "\nLibrary: ", .x,
      "\nEdit distance"
    )) +
    th +
    theme(
      strip.text.y = element_blank(),
      plot.title = element_text(size = 10),
    )

  ## read length distribution
  d_rl <- d_bam |>
    count(qwidth)

  lab <- paste0(
    "Avg read length: ", round(d_stats$read_l_avg, 2), " bp"
  )
  d_lab <- tibble(
    x = max(d_rl$qwidth),
    y = max(d_rl$n),
    label = lab
  )
  p[[2]] <- d_rl |>
    ggplot(aes(
      x = qwidth,
      y = n
    )) +
    geom_col(fill = "navy") +
    geom_text(
      aes(
        x = x,
        y = y,
        label = label
      ),
      hjust = 1,
      vjust = 1,
      size = 3,
      data = d_lab,
    ) +
    scale_x_continuous(
      name = "Read length (bp)",
      breaks = seq(0, 300, by = 10)
    ) +
    scale_y_continuous(
      name = "Read count",
      labels = comma,
    ) +
    ggtitle("\n\nRead length") +
    th +
    theme(
      strip.text.y = element_blank(),
      plot.title = element_text(size = 10),
    )

  ## damage
  d_dmg <- dmg_tbl |>
    filter(library_id == .x)

  lab <- paste0(
    "Dfit amplitude: ", round(d_stats$dmg_amp, 3),
    "\nDfit Z: ", round(d_stats$dmg_Zfit, 3)
  )
  d_lab <- tibble(
    x = 15,
    y = 0.3,
    label = lab
  )

  p[[3]] <- d_dmg |>
    ggplot(aes(
      x = position,
    )) +
    geom_hline(
      yintercept = 0,
      linewidth = 0.25,
      linetype = "dashed"
    ) +
    geom_ribbon(
      aes(
        ymin = fwdx - fwdxConf,
        ymax = fwdx + fwdxConf
      ),
      color = NA,
      fill = "grey",
      alpha = 0.5,
    ) +
    geom_line(
      aes(y = fwf),
      color = "firebrick1",
      linewidth = 0.25,
    ) +
    geom_point(
      aes(y = fwf),
      color = "firebrick1",
      size = 1,
      shape = 4
    ) +
    geom_line(aes(y = fwdx),
      linewidth = 0.25,
    ) +
    geom_text(
      aes(
        x = x,
        y = y,
        label = label
      ),
      hjust = 1,
      vjust = 1,
      size = 3,
      data = d_lab,
    ) +
    xlab("Position from read end") +
    ylab("Damage frequency") +
    scale_x_continuous(
      breaks = 0:15
    ) +
    coord_cartesian(
      ylim = c(0, 0.3),
      xlim = c(0, 15)
    ) +
    ggtitle("\n\nDamage") +
    th +
    theme(
      strip.text.y = element_blank(),
      plot.title = element_text(size = 10)
    )

  ## assemble and save
  pl_row <- plot_grid(
    p[[1]] + theme(legend.position = "none"),
    p[[2]] + theme(legend.position = "none"),
    p[[3]] + theme(legend.position = "none"),
    rel_widths = c(1, 1, 1), nrow = 1
  )
  print(pl_row)
})
