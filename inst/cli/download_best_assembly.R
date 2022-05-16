#!/usr/bin/env Rscript

library(docopt)
library(utilitybeltrefseq)
'Download Best Assembly

Usage:
  download_best_assembly <taxid_of_interest>
  download_best_assembly -h | --help

Example:
  # Download Mycoplasma genitalium reference genome
  download_best_assembly 2097
' -> doc

arguments <- docopt(doc)

print(arguments)

download_best_assembly(taxid_of_interest = arguments$taxid_of_interest)


