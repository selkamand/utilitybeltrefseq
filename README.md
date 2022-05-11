
<!-- README.md is generated from README.Rmd. Please edit that file -->

# utilitybeltrefseq

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)

<!-- badges: end -->

The goal of utilitybeltrefseq is to make it easy to find and download
the best refseq assemblies available for any given species. All you need
to know is the species [NCBI taxonomy
ID](https://www.ncbi.nlm.nih.gov/taxonomy)

## Installation

You can install the development version of utilitybeltrefseq from
[GitHub](https://github.com/selkamand/utilitybeltrefseq) with:

``` r
# install.packages("devtools")
devtools::install_github("selkamand/utilitybeltrefseq")
```

## Example

Say you want to find the best reference genome for *Escherichia coli*.

First we load the library then download a list of available assemblies
from refseq.

``` r
library(utilitybeltrefseq)
```

``` r
update_refseq_data_cache()
```

For best results, run the above command every couple of months to stay
up to date with whats currently in refseq.

Next, we need to find the NCBI taxonomy ID of our species of interest.
We can find that [here](https://www.ncbi.nlm.nih.gov/taxonomy) (taxid:
562). We could also have used the **taxize** package if we wanted to
stay within R.

Now that we have the species level taxid, we can run
`choose_best_assembly`

``` r
choose_best_assembly(taxid_of_interest = 562)
#> Multiple (2) best hits with score = 101400. We will just return the the most recently added assembly with this quality
#> Chosen Assembly:
#>   assembly_accession      GCF_000008865.2 
#>   bioproject      PRJNA57781 
#>   biosample   SAMN01911278 
#>   wgs_master       
#>   refseq_category     reference genome 
#>   taxid   386585 
#>   species_taxid   562 
#>   organism_name   Escherichia coli O157:H7 str. Sakai 
#>   infraspecific_name      strain=Sakai substr. RIMD 0509952 
#>   isolate      
#>   version_status      latest 
#>   assembly_level      Complete Genome 
#>   release_type    Major 
#>   genome_rep      Full 
#>   seq_rel_date    2018/06/08 
#>   asm_name    ASM886v2 
#>   submitter   GIRC 
#>   gbrs_paired_asm     GCA_000008865.2 
#>   paired_asm_comp     identical 
#>   ftp_path    https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/008/865/GCF_000008865.2_ASM886v2 
#>   excluded_from_refseq     
#>   relation_to_type_material    
#>   asm_not_live_date   na 
#>   ref_score   101400 
#>   best_ref    TRUE
#> [1] "GCF_000008865.2"
```

Note that we get a warning that there are multiple ‘best’ assemblies. If
there are multiple strains in a species with complete assemblies - it
his hard to know which to return. We choose the most recently added
assembly. To return all of the high quality assemblies you can run:

``` r
 choose_best_assembly(
   taxid_of_interest = 562, 
   break_ties_based_on_newest_sequence_added = FALSE,
   return_accession_only = TRUE # set this to false to get more info about each assembly
   )
#> Multiple (2) best hits with score = 101400. Please choose one manually or if you haven't already - add an intraspecific filter and try again
#> [1] "GCF_000005845.2" "GCF_000008865.2"
```

Keep an eye on the intraspecific name column. Often can use this to
choose which strain you’re after. See `?choose_best_assembly` for
details.

Once we have the assembly accession we’re interested in, we can download
it:

``` r
 download_assembly(
   target_assembly_accession = GCF_000005845.2
   )
```

If we just want to quickly download a ref-genome of our species of
interest we can just run

``` r
download_best_assembly(taxid_of_interest = 562)
```
