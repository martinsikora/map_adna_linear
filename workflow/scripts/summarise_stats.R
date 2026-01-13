#!/usr/bin/env Rscript
# Copyright 2025 Martin Sikora <martin.sikora@sund.ku.dk>
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

# summary statistics for libraries across samples

## ---------------------------------------------------------
## libraries

suppressPackageStartupMessages({
  library(argparse)
  library(tidyverse)
  library(ggrepel)
  library(viridis)
  library(scales)
})


## --------------------------------------------------
## command line argument setup and parsing

parser <- ArgumentParser()

parser$add_argument("files",
  nargs = "+",
  help = "Sample bam stats files to process"
)

parser$add_argument("-p", "--prefix_out_plots",
  action = "store",
  dest = "prefix_out_plots",
  help = "Output file name prefix for plots"
)

parser$add_argument("-t", "--prefix_out_tables",
  action = "store",
  dest = "prefix_out_tables",
  help = "Output file name prefix for tables"
)
args <- parser$parse_args()


## ---------------------------------------------------------
## helpers

th <- theme_bw() +
  theme(
    panel.grid.major = element_line(
      linetype = "dotted",
      linewidth = 0.25
    ),
    panel.grid.minor = element_blank(),
    strip.background = element_blank()
  )


## ---------------------------------------------------------
## bam stats

bam_stats <- map_dfr(args$files, read_tsv)

## read depth vs breadth
d <- bam_stats |>
  filter(
    library_id == "all",
    contig_id == "genome"
  )

pdf(paste0(args$prefix_out_plots, ".bam_stats_coverage.pdf"),
  width = 6.5,
  height = 4.5
)
p <- ggplot(d, aes(
  x = dp_avg,
  y = coverage_p
))
p +
  geom_hline(
    yintercept = 1,
    size = 0.25,
    linetype = "dashed"
  ) +
  geom_vline(
    xintercept = 1,
    size = 0.25,
    linetype = "dashed"
  ) +
  geom_point(
    size = 2,
    shape = 21,
    fill = "grey80",
    color = "white"
  ) +
  geom_function(
    fun = function(x) 1 - exp(-x),
    size = 0.5,
    colour = "grey30"
  ) +
  geom_point(
    aes(
      fill = cov_pos_rel_entropy_1000,
    ),
    size = 2,
    stroke = 0.25,
    shape = 21,
  ) +
  geom_text_repel(
    aes(label = sample_id),
    segment.size = 0.25,
    min.segment.length = 0,
    segment.color = "grey",
    size = 1.5,
  ) +
  annotation_logticks(sides = "bl") +
  scale_fill_viridis(
    name = "Relative entropy\nread start positions",
    limits = c(min(d$cov_pos_rel_entropy_1000, na.rm = TRUE), 1)
  ) +
  scale_x_continuous(
    name = "Average read depth",
    trans = "log10",
    labels = trans_format("log10", math_format(10^.x))
  ) +
  scale_y_continuous(
    name = "Breadth of coverage",
    limits = c(5e-4, 1),
    trans = "log10",
    labels = trans_format("log10", math_format(10^.x))
  ) +
  th
dev.off()

## write table

o <- bam_stats |>
  filter(library_id == "all") |>
  select(-library_id)

write_tsv(o, file = paste0(args$prefix_out_tables, ".bam_stats.tsv"))


## ---------------------------------------------------------
## lib stats

f <- gsub("bam_stats", "lib_stats", args$files)
lib_stats <- map_dfr(f, read_tsv)

## duplicate rate
w <- nrow(lib_stats) %/% 10 + 4

pdf(paste0(args$prefix_out_plots, ".lib_stats_p_dup.pdf"),
  width = w,
  height = 5
)
p <- ggplot(lib_stats, aes(
  x = library_id,
  y = prop_dup
))
p +
  geom_hline(
    yintercept = 1,
    linetype = "dotted"
  ) +
  geom_col(aes(y = 1),
    fill = "gainsboro"
  ) +
  geom_col(fill = "grey20") +
  facet_grid(. ~ sample_id,
    scales = "free_x",
    space = "free_x"
  ) +
  xlab("Library ID") +
  ylab("Duplicate rate") +
  th +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      size = 6,
      vjust = 0.5
    ),
    strip.text.x = element_text(
      angle = 90,
      hjust = 0,
      vjust = 0.5
    ),
    panel.spacing = unit(0, "lines")
  )
dev.off()

## total reads vs duplicate rate
sample_ids <- sort(unique(lib_stats$sample_id))
x_lim <- c(1, max(lib_stats$num_reads_total, na.rm = TRUE) / 1e6)

pdf(paste0(args$prefix_out_plots, ".lib_stats_n_reads_p_dup.pdf"),
  width = 6,
  height = 5
)

p <- ggplot(lib_stats, aes(
  x = num_reads_total / 1e6,
  y = prop_dup
))
walk(sample_ids, ~ {
  d1 <- lib_stats |>
    filter(sample_id == .x)
  print(p +
    geom_hline(
      yintercept = 1,
      linetype = "dotted"
    ) +
    geom_point(
      shape = 19,
      color = "grey"
    ) +
    geom_point(
      shape = 21,
      size = 2,
      color = "black",
      fill = "royalblue",
      data = d1
    ) +
    geom_text_repel(aes(label = library_id),
      size = 2,
      segment.size = 0.25,
      min.segment.length = 0,
      data = d1
    ) +
    scale_x_continuous(
      trans = "log10",
      limits = x_lim
    ) +
    annotation_logticks(sides = "b") +
    xlab("N Reads (million)") +
    ylab("Duplicate rate") +
    th +
    ggtitle(.x))
})
dev.off()

## unique mapped reads vs dup rate
x_lim <- c(1, max(lib_stats$num_reads_final, na.rm = TRUE))

pdf(paste0(args$prefix_out_plots, ".lib_stats_n_map_uniq_p_dup.pdf"),
  width = 6,
  height = 5
)

p <- ggplot(lib_stats, aes(
  x = num_reads_final,
  y = prop_dup
))

walk(sample_ids, ~ {
  d1 <- lib_stats |>
    filter(sample_id == .x)
  print(p +
    geom_hline(
      yintercept = 1,
      linetype = "dotted"
    ) +
    geom_point(
      shape = 19,
      color = "grey"
    ) +
    geom_point(
      shape = 21,
      size = 2,
      color = "black",
      fill = "royalblue",
      data = d1
    ) +
    geom_text_repel(aes(label = library_id),
      size = 2,
      segment.size = 0.25,
      min.segment.length = 0,
      data = d1
    ) +
    scale_x_continuous(
      trans = "log10",
      limits = x_lim
    ) +
    annotation_logticks(sides = "b") +
    xlab("N unique reads ") +
    ylab("Duplicate rate") +
    th +
    ggtitle(.x))
})
dev.off()


## total reads vs unique mapped reads
x_lim <- c(
  1,
  max(lib_stats$num_reads_total, na.rm = TRUE) / 1e6
)
y_lim <- c(
  1,
  max(lib_stats$num_reads_final, na.rm = TRUE)
)

pdf(paste0(args$prefix_out_plots, ".lib_stats_n_reads_n_map_uniq.pdf"),
  width = 4,
  height = 5
)

p <- ggplot(lib_stats, aes(
  x = num_reads_total / 1e6,
  y = num_reads_final
))

walk(sample_ids, ~ {
  d1 <- lib_stats |>
    filter(sample_id == .x) |>
    mutate(lab = paste0(
      library_id, "\np final = ",
      formatC(prop_final, format = "g", digits = 2)
    ))
  print(p +
    geom_abline(
      slope = 1,
      intercept = seq(-10, 10),
      size = 0.25,
      linetype = "dashed",
      color = "grey60"
    ) +
    geom_abline(
      slope = 1,
      intercept = 0,
      size = 0.25,
      linetype = "solid",
      color = "grey60"
    ) +
    geom_point(
      shape = 19,
      color = "grey"
    ) +
    geom_point(
      shape = 21,
      size = 2,
      color = "black",
      fill = "royalblue",
      data = d1
    ) +
    geom_text_repel(aes(label = lab),
      size = 1.5,
      hjust = 0,
      segment.size = 0.25,
      min.segment.length = 0,
      data = d1,
      max.overlaps = Inf
    ) +
    scale_y_continuous(
      trans = "log10",
      labels = comma,
      limits = y_lim
    ) +
    scale_x_continuous(
      trans = "log10",
      labels = comma,
      limits = x_lim
    ) +
    coord_equal() +
    annotation_logticks(sides = "bl") +
    xlab("N Reads (million)") +
    ylab("N unique reads") +
    th +
    ggtitle(.x))
})
dev.off()

## write final lib stats table with libraries ranked

o <- lib_stats |>
  arrange(sample_id, desc(prop_final)) |>
  group_by(sample_id) |>
  mutate(
    lib_rank = row_number()
  ) |>
  ungroup()

write_tsv(o, file = paste0(args$prefix_out_tables, ".lib_stats.tsv"))
