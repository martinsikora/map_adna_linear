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

# summary statistics per library file for mapped reads of a sample

## ---------------------------------------------------------
## libraries

suppressPackageStartupMessages({
  library(argparse)
  library(tidyverse)
})


## --------------------------------------------------
## command line argument setup and parsing

parser <- ArgumentParser()

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

args <- parser$parse_args()


## ---------------------------------------------------------
## fq stats

infiles <- list.files(
  path = file.path(args$out_dir, args$sample_id, "logs"),
  pattern = paste0(args$sample_id, ".*\\.map_init.log$"),
  full.names = TRUE
)

r_fq <- map_dfr(infiles, ~ {
  r1 <- readLines(.x)

  num_reads_paired <- gsub(
    " +([0-9]+) .*", "\\1",
    r1[grepl("were paired", r1)]
  ) |>
    as.numeric()

 num_reads_unpaired <- gsub(
   " +([0-9]+) .*", "\\1",
   r1[grepl("were unpaired", r1)]
 ) |>
   as.numeric()
 
  num_reads_total <- num_reads_paired + num_reads_unpaired

  r2 <- gsub(".*\\/", "", .x) |>
    strsplit("\\.") |>
    unlist()

  tibble(
    sample_id = r2[1],
    library_id = r2[2],
    rg_id = r2[3],
    num_reads_total = num_reads_total,
  )
})


## ---------------------------------------------------------
## markdup stats

infiles <- list.files(
  path = file.path(args$out_dir, args$sample_id, "logs"),
  pattern = paste0(args$sample_id, ".*\\.mrkdup.metrics.txt"),
  full.names = TRUE
)

r_markdup <- map_dfr(infiles, ~ {
  r1 <- read_tsv(.x, skip = 6, n_max = 1)
  r2 <- gsub(".*\\/", "", .x) |>
    strsplit("\\.") |>
    unlist()
  tibble(
    sample_id = r2[1],
    library_id = r1$LIBRARY,
    num_reads_mapped = r1$UNPAIRED_READS_EXAMINED + r1$READ_PAIRS_EXAMINED,
    num_reads_dup = r1$UNPAIRED_READ_DUPLICATES + r1$READ_PAIR_DUPLICATES,
  )
})


## ---------------------------------------------------------
## bam stats

infile <- file.path(
  args$out_dir, args$sample_id, "tables",
  paste0(
    args$sample_id, ".",
    args$prefix, ".", args$ref, ".filter.bam_stats.tsv"
  )
)

r_bam <- read_tsv(infile) |>
  filter(
    library_id != "all",
    contig_id == "genome"
  ) |>
  select(sample_id, library_id, n_reads) |>
  rename(num_reads_final = n_reads)


## ---------------------------------------------------------
## final result table

res <- r_fq |>
  left_join(r_bam) |>
  group_by(sample_id, library_id) |>
  summarise(
    num_reads_total = sum(num_reads_total, na.rm = TRUE),
    num_reads_final = sum(num_reads_final, na.rm = TRUE),
    .groups = "drop"
  ) |>
  left_join(r_markdup) |>
  mutate(
    prop_mapped = num_reads_mapped / num_reads_total,
    prop_dup = num_reads_dup / num_reads_mapped,
    prop_final = num_reads_final / num_reads_total
  ) |>
  select(
    sample_id, library_id, num_reads_total, num_reads_mapped,
    prop_mapped, num_reads_dup, prop_dup, num_reads_final, prop_final
  ) |>
  arrange(sample_id, library_id)

outfile <- file.path(
  args$out_dir, args$sample_id, "tables",
  paste0(
    args$sample_id, ".",
    args$prefix, ".", args$ref, ".filter.lib_stats.tsv"
  )
)
write_tsv(res, file = outfile, na = "NaN")
