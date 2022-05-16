

#' Cache refseq contents
#'
#' Download a file that lists metadata for all refseq entries (including ftp path).
#' We will use this later to download assemblies.
#'
#' @return refseq dataframe describing available assemblies
#'
#' @examples
#' \dontrun{
#' load_refseq_data_frame_from_ftp()
#' }
load_refseq_data_frame_from_ftp <- function(){
  refseq_data_frame = utils::read.csv("https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt", sep = "\t", header=TRUE, skip = 1, check.names = FALSE)
  names(refseq_data_frame) <- sub(pattern = "# ?", replacement = "", x = names(refseq_data_frame))
  # assembly_accession

  required_columns = c("species_taxid", "refseq_category", "organism_name", "assembly_level", "assembly_accession", "version_status")
  utilitybeltassertions::assert_names_include(refseq_data_frame, expected_names = required_columns)
  assertthat::assert_that(nrow(refseq_data_frame) > 0, msg = utilitybeltassertions::fmterror("Refseq Dataframe has 0 rows. Somethings probably wrong with the ftp download link. Please ensure you're using the latest version of this package, and if the problem persists, open a github issue so package maintaner can fix this"))

  # Guess best reference
  refseq_data_frame = refseq_data_frame %>%
    dplyr::group_by(species_taxid) %>%
    dplyr::mutate(
      ref_score = 0 +
        (refseq_category=="reference genome") * 10^5 +
        (refseq_category=="representative genome") * 10^4 +
        (genome_rep=="Full") * 10^3 +
        (assembly_level == "Complete Genome") * 4*10^2 +
      (assembly_level == "Chromosome") * 3*10^2 +
      (assembly_level == "Scaffold") * 2*10^2 +
      (assembly_level == "Contig") * 1*10^2,
      best_ref = ref_score == max(ref_score)
    ) %>%
    dplyr::ungroup()

  return(refseq_data_frame)
}

cache_dir_location <- function(){
  return("~/.utilitybeltrefseq")
}

cache_file_location <- function(){
  return(paste0(cache_dir_location(),"/assembly_summary_refseq.txt"))
}

#' Write refseq dataframe
#'
#' @param refseq_data_frame dataframe produced by [load_refseq_data_frame_from_ftp]
#'
#' @return NULL
#'
write_refseq_dataframe <- function(refseq_data_frame){
  cache_dir_location=cache_dir_location()
  if(!dir.exists(cache_dir_location)){
    message("Creating folder to store refseq cache: ", cache_dir_location)
    dir.create(cache_dir_location())
  }
  else{
    message("Using existing folder to store refseq cache: ", cache_dir_location)
  }

  outfile = cache_file_location()

  message("Writing refseq table to file: [", outfile, "]")
    utils::write.table(
      x = refseq_data_frame,
      file = outfile,
      row.names = FALSE,
      sep = "\t"
    )
}

#' Update refseq data cache
#'
#' Download a file that lists metadata for all refseq entries (including ftp path).
#' We will use this later to download assemblies.
#'
#' @return NULL. Run for side-effects
#' @export
#'
#' @examples
#' \dontrun{
#' update_refseq_data_cache()
#' }
update_refseq_data_cache <- function(){
  message("Updating the refseq data cache usually take several minutes ...")
  write_refseq_dataframe(load_refseq_data_frame_from_ftp())
  message("Done")
}

#' Load data from cache
#'
#' Loads a dataframe describing reference genomes and their ftp links
#'
#' @return data.frame
#' @export
#'
#' @examples
#' refseq_data_frame = load_refseq_data_frame_from_cache()
load_refseq_data_frame_from_cache <- function(){
  expected_filepath = cache_file_location()
  assertthat::assert_that(file.exists(expected_filepath), msg = paste0("Could not find a cached refseq dataframe at [",expected_filepath,"]. Please run update_refseq_data_cache"))

  #browser()
  refseq_data_frame = data.table::fread(expected_filepath, sep = "\t", check.names = FALSE)
  refseq_data_frame
  required_columns = c("taxid", "species_taxid", "refseq_category", "organism_name", "assembly_level", "assembly_accession", "version_status", "ref_score", "seq_rel_date")
  #browser()
  utilitybeltassertions::assert_names_include(refseq_data_frame, expected_names = required_columns)
  assertthat::assert_that(nrow(refseq_data_frame) > 0)
  return(refseq_data_frame)
}

#' Delete the refseq Cache
#'
#' Remove existing cache. Good practice to run before uninstalling package.
#'
#' @export
#'
delete_refseq_data_cache <- function(){
  if( dir.exists(cache_dir_location())){
    confirmed = utils::askYesNo("Are you sure you want to delete the folder: ", dir.exists())

    if(confirmed){
      unlink(cache_dir_location(),recursive=TRUE)
      message("Cache deleted")
    }
    else
      message("Cache not deleted")
  }
  else{
   message("Nothing to delete. No cache at: ", cache_dir_location())
  }
}


#' Best Assembly
#'
#' @param taxid_of_interest taxonomyID (int)
#' @param intraspecific_filter e.g. 'strain=Raji'. For options look at the infraspecific_name column of the refseq_data_frame (string)
#' @param break_ties_based_on_newest_sequence_added if we cant decide which is the best ref to use because they;re so similar - just pick the most recently uploaded sequence
#' @param refseq_data_frame from \strong{load_refseq_data_frame_from_cache()}
#' @param return_accession_only Return just the best assembly accession. If false will return a dataframe with accession ID plus lots of other info (bool)
#'
#' @return if return_accession_only=TRUE then best assembly accession (string). Otherwise returns dataframe with accession ID of best assembly plus lots of other info (bool)
#' @export
#'
choose_best_assembly <- function(taxid_of_interest, intraspecific_filter = NA, break_ties_based_on_newest_sequence_added = TRUE, return_accession_only = TRUE, refseq_data_frame = load_refseq_data_frame_from_cache()){

  # Trigger error if refseq_data_frame doesnt exist
  invisible(utils::head(refseq_data_frame, n=1))
  #browser()

  # is species in refseq database
  assertthat::assert_that(taxid_of_interest %in% refseq_data_frame[["species_taxid"]],
                          msg = utilitybeltassertions::fmterror("Can't find taxid [", taxid_of_interest, "]", " in refseq database. Are you sure this is a species-level taxid with a refseq assembly?"))


  best_assemblies = refseq_data_frame %>%
    dplyr::filter(species_taxid == taxid_of_interest) %>%
    dplyr::slice_max(ref_score, with_ties = TRUE)

  if(!is.na(intraspecific_filter)){
    message("running through intraspecific filter: ", intraspecific_filter)
    best_assemblies <- best_assemblies %>%
      dplyr::filter(infraspecific_name == intraspecific_filter)
    assertthat::assert_that(nrow(best_assemblies) > 0, msg = utilitybeltassertions::fmterror("The addition of intraspecific species filter: [",intraspecific_filter,"] does not contain any of the species level best assemblies"))

  assertthat::assert_that(
    best_assemblies[["ref_score"]] > 0,
    msg = utilitybeltassertions::fmterror("Couldn't find any decent assemblies to use: perhaps try to manually find a good refseq assembly accession?"))
    }

  score_of_tophit = max(best_assemblies[["ref_score"]])
  if (nrow(best_assemblies) > 1){
    num_of_tied_assemblies = nrow(best_assemblies)

    if(break_ties_based_on_newest_sequence_added){
      best_assemblies_final=dplyr::slice_max(best_assemblies, seq_rel_date)
      message("Multiple (", num_of_tied_assemblies,")", " best hits with score = ",score_of_tophit,". We will just return the the most recently added assembly with this quality\n")
      prettyprint_single_row_df(best_assemblies_final, title = "Chosen Assembly: ")
      if(return_accession_only) return(best_assemblies_final[["assembly_accession"]]) else return(dplyr::tibble(best_assemblies_final))
    }
    else{
      message("Multiple (", num_of_tied_assemblies,")", " best hits with score = ",score_of_tophit, ". Please choose one manually or if you haven't already - add an intraspecific filter and try again")
      best_assemblies_final=dplyr::tibble(best_assemblies)
      if(return_accession_only) return(best_assemblies_final[["assembly_accession"]]) else return(dplyr::tibble(best_assemblies_final))
    }

  }

  if(return_accession_only) {prettyprint_single_row_df(best_assemblies, title = "Chosen Assembly: "); return(best_assemblies[["assembly_accession"]]) } else return(dplyr::tibble(best_assemblies))
}

#' Pretty Print a Single Refseq Entry
#'
#' @param title text printed at the top of the prettyprint summary (string)
#' @param single_row_of_tabular_data any dataframe with 1 row. Doesn't require any particular column_names / row typing. (data.frame)
#'
#' @return NULL
prettyprint_single_row_df <- function(single_row_of_tabular_data, title = NULL){
  assertthat::assert_that(nrow(single_row_of_tabular_data) == 1, msg = utilitybeltassertions::fmterror("Function `prettyprint_assembly_entry` was designed to work with single row dataframes. Input dataframe has ", nrow(single_row_of_tabular_data), " entries"))
  names = names(single_row_of_tabular_data)
  values = as.character(single_row_of_tabular_data)

  prettyprint_format = paste("\t", utilitybeltassertions::fmtbold(names), "\t",values, "\n")

  if(!is.null(title)){
    message(utilitybeltassertions::fmtbold(title))
  }

  message(prettyprint_format)

  return(invisible(NULL))
}

#' Download reference genome
#'
#' Actually download the genome using aria2c
#'
#' @inheritParams choose_best_assembly
#' @param target_assembly_accession accession of target assembly
#' @param output_folder directory to save assembly fasta file into (string)
#'
#' @return filepath of downloaded file
#' @export
download_assembly <- function(target_assembly_accession, output_folder = getwd(), refseq_data_frame = load_refseq_data_frame_from_cache()){

  #browser()
  assertthat::assert_that(target_assembly_accession %in% refseq_data_frame[["assembly_accession"]])

  refseq_data_frame = refseq_data_frame %>%
    dplyr::filter(assembly_accession %in% target_assembly_accession)

  assertthat::assert_that(nrow(refseq_data_frame) == 1)

  sequence_file_name = basename(refseq_data_frame[["ftp_path"]])
  ftp_root = refseq_data_frame[["ftp_path"]]
  assembly_fasta_path = paste0(ftp_root, "/", sequence_file_name, "_genomic.fna.gz")
  #browser()

  assertthat::assert_that(dir.exists(output_folder), msg = utilitybeltassertions::fmterror("Could not find folder: ", output_folder))
  full_dest_filepath = paste0(output_folder, "/", basename(assembly_fasta_path))

  message("downloading assembly:\n\t[", target_assembly_accession, "]\n\nfrom the RefSeq ftp link \n\t[", assembly_fasta_path, "]\nto:\t", full_dest_filepath)

  utils::download.file(assembly_fasta_path, destfile =  full_dest_filepath)

  return(full_dest_filepath)
}


download_assembly_aria2c <- function(refseq_data_frame = load_refseq_data_frame_from_cache(), target_assembly_accession){

  utilitybeltassertions::assert_program_exists_in_path("aria2c")
  #browser()
  assertthat::assert_that(target_assembly_accession %in% refseq_data_frame[["assembly_accession"]])

  refseq_data_frame = refseq_data_frame %>%
    dplyr::filter(assembly_accession %in% target_assembly_accession)

  assertthat::assert_that(nrow(refseq_data_frame) == 1)

  sequence_file_name = basename(refseq_data_frame[["ftp_path"]])
  ftp_root = refseq_data_frame[["ftp_path"]]
  assembly_fasta_path = paste0(ftp_root, "/", sequence_file_name, "_genomic.fna.gz")
#browser()
  message("downloading assembly:\n\t[", target_assembly_accession, "]\n\nfrom the RefSeq ftp link \n\t[", assembly_fasta_path, "]")

  res=system(paste0("aria2c ", assembly_fasta_path))
  if (res != 0)
    message("Download failed")
  else
    message("download complete")

  #sys.calls(
}

#' Download_best_assembly
#'
#' @inheritDotParams download_assembly
#' @inheritParams choose_best_assembly
#' @inherit download_assembly return
#' @export
#'
download_best_assembly <- function(taxid_of_interest, ...){
  best_asssembly=target_assembly_accession = choose_best_assembly(taxid_of_interest = taxid_of_interest, return_accession_only = TRUE)
  message("Attempting to download the best assembly for taxid ", taxid_of_interest, " (", best_asssembly, ")")
  res=download_assembly(target_assembly_accession = best_asssembly, ...)
  return(res)
}

# download_annotations <- function(refseq_data_frame = load_refseq_data_frame_from_cache(), target_assembly_accession, outfile_dir = "."){
#
#   utilitybeltassertions::assert_program_exists_in_path("aria2c")
#   #browser()
#   assertthat::assert_that(target_assembly_accession %in% refseq_data_frame[["assembly_accession"]])
#
#   refseq_data_frame = refseq_data_frame %>%
#     dplyr::filter(assembly_accession %in% target_assembly_accession)
#
#   assertthat::assert_that(nrow(refseq_data_frame) == 1)
#
#   sequence_file_name = basename(refseq_data_frame[["ftp_path"]])
#   ftp_root = refseq_data_frame[["ftp_path"]]
#   assembly_fasta_path = paste0(ftp_root, "/", sequence_file_name, "_genomic.gtf.gz")
#   #browser()
#   message("downloading assembly:\n\t[", target_assembly_accession, "]\n\nfrom the RefSeq ftp link \n\t[", assembly_fasta_path, "]")
#
#   res=system(paste0("aria2c ", assembly_fasta_path))
#   if (res != 0)
#    message("Download failed")
#   else
#     message("download complete")
# }

#' List commandline scripts
#'
#' @export
#'
#' @examples
#' cli()
cli <- function(){
  executable_paths = dir(system.file("cli/", package = "utilitybeltrefseq"), pattern = ".R$", full.names = TRUE)
  message(
    "Ensure scripts are executable by running the following in your terminal:
    chmod +x ", paste0(executable_paths, collapse=" "),
    "\n\nThen run any of the following commands scripts to learn about usage: \n\t", paste0(executable_paths, collapse="\n")
    )
}

