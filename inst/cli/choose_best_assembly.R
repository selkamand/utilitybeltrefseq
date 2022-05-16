#!/usr/bin/env Rscript

library(docopt)
library(utilitybeltrefseq)
'Download Best Assembly

Usage:
  choose_best_assembly.R [--all --verbose] <taxid_of_interest>
  choose_best_assembly.R -h | --help

Options:



Example:
  # Download Mycoplasma genitalium reference genome
  choose_best_assembly 2097
' -> doc

arguments <- docopt(doc)

#print(arguments)

choose_best_assembly(taxid_of_interest = arguments$taxid_of_interest,return_accession_only = !arguments$verbose,  break_ties_based_on_newest_sequence_added = !arguments$all)


