


# RefSeq Metadata Cache ------------------------------------------------------------


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

cache_dir_location <- function(verbose = TRUE){
  cache_dir_path_environment_var=Sys.getenv("UTILITYBELTREFSEQ_CACHE_LOCATION")
  if(cache_dir_path_environment_var != ""){
    if(verbose) message("Basing cache directory on environment variable (UTILITYBELTREFSEQ_CACHE_LOCATION):\n", cache_dir_path_environment_var)
    return(cache_dir_path_environment_var) # Allow environmental variable $UTILITYBELTREFSEQ_CACHE_LOCATION to overrwrite default cache location (~/.utilitybeltrefseq)
  }
  else{
    return("~/.utilitybeltrefseq")
  }
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



#' Pretty Print a Single Refseq Entry
#'
#' @param title text printed at the top of the prettyprint summary (string)
#' @param single_row_of_tabular_data any dataframe with 1 row. Doesn't require any particular column_names / row typing. (data.frame)
#' @param write_output_to_file should we write a simplified output to file (boolean)
#' @param file_title required only if write_output_to_file == TRUE. Title of file to write output to (string)
#' @param verbose verbosity level (boolean)
#'
#' @return NULL
prettyprint_single_row_df <- function(single_row_of_tabular_data, title = NULL, write_output_to_file = FALSE, file_title = NULL, verbose = TRUE){
  assertthat::assert_that(nrow(single_row_of_tabular_data) == 1, msg = utilitybeltassertions::fmterror("Function `prettyprint_assembly_entry` was designed to work with single row dataframes. Input dataframe has ", nrow(single_row_of_tabular_data), " entries"))
  names = names(single_row_of_tabular_data)
  values = as.character(single_row_of_tabular_data)

  prettyprint_format = paste("\t", utilitybeltassertions::fmtbold(names), "\t",values, "\n")

  if(!is.null(title)){
    message(utilitybeltassertions::fmtbold(title))
  }

  message(prettyprint_format)

  if (write_output_to_file){
    assertthat::assert_that(!is.null(file_title), msg = "[prettyprint_single_row_df] please specify the name of the text file you want to write output to using the `file_title` argument")
    if(verbose) message("writing to file: [", file_title, "]")
    utils::write.table(data.frame("Properties" = names, "Values" = values), file = file_title, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
  }

  return(invisible(NULL))
}



# Choose Best Assembly ----------------------------------------------------



#' Best Assembly
#'
#' @param taxid_of_interest taxonomyID (int)
#' @param intraspecific_filter e.g. 'strain=Raji'. For options look at the infraspecific_name column of the refseq_data_frame (string)
#' @param break_ties_based_on_newest_sequence_added if we cant decide which is the best ref to use because they;re so similar - just pick the most recently uploaded sequence
#' @param refseq_data_frame from \strong{load_refseq_data_frame_from_cache()}
#' @param return_accession_only Return just the best assembly accession. If false will return a dataframe with accession ID plus lots of other info (bool)
#' @param log_chosen_assembly should we log details of chosen assembly to a file in addition to printing the results
#' @param logtitle title of logfile (best_assembly_for_taxid_`taxid`.txt")
#'
#' @return if return_accession_only=TRUE then best assembly accession (string). Otherwise returns dataframe with accession ID of best assembly plus lots of other info (bool)
#' @export
#'
choose_best_assembly <- function(taxid_of_interest, intraspecific_filter = NA, break_ties_based_on_newest_sequence_added = TRUE, return_accession_only = TRUE, refseq_data_frame = load_refseq_data_frame_from_cache(), log_chosen_assembly = TRUE, logtitle = paste0("best_assembly_for_taxid_", taxid_of_interest, ".txt")){

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
      prettyprint_single_row_df(best_assemblies_final, title = "Chosen Assembly: ", write_output_to_file = log_chosen_assembly, file_title = logtitle)
      #if(log_chosen_assembly){ write(taxid_of_interest}
      if(return_accession_only) return(best_assemblies_final[["assembly_accession"]]) else return(dplyr::tibble(best_assemblies_final))
    }
    else{
      message("Multiple (", num_of_tied_assemblies,")", " best hits with score = ",score_of_tophit, ". Please choose one manually or if you haven't already - add an intraspecific filter and try again")
      best_assemblies_final=dplyr::tibble(best_assemblies)
      if(return_accession_only) return(best_assemblies_final[["assembly_accession"]]) else return(dplyr::tibble(best_assemblies_final))
    }

  }

  if(return_accession_only) {prettyprint_single_row_df(best_assemblies, title = "Chosen Assembly: ", write_output_to_file = log_chosen_assembly, file_title = logtitle); return(best_assemblies[["assembly_accession"]]) } else return(dplyr::tibble(best_assemblies))
}

# Download Refseq Files ---------------------------------------------------

#' download refseq assembly files
#'
#' @param target_assembly_accession accession of assembly you're interested in. For example GCF_000371765.1 (string)
#' @param assembly_file_type type of file to download from refseq repo (string)
#' @param output_folder directory to output file into (string)
#' @param refseq_data_frame do not reccomend changing this. Automatically loads from \strong{load_refseq_data_frame_from_cache()}
#'
#' @return path to downloaded file
#' @export
#'
#' @examples
#' \dontrun{
#' download_assembly_files(
#'   target_assembly_accession = "GCF_000371765.1",
#'   assembly_file_type = "genome"
#' )
#' }
#'
download_assembly_files <- function(target_assembly_accession, assembly_file_type = c("genome", "stats", "genome_annotation"), output_folder = getwd(), refseq_data_frame = load_refseq_data_frame_from_cache()){
  assembly_file_type <- rlang::arg_match(assembly_file_type)
  assembly_file_extension = dplyr::case_when(
    assembly_file_type == "genome" ~ "_genomic.fna.gz",
    assembly_file_type == "stats" ~ "_assembly_stats.txt",
    assembly_file_type == "annotation" ~ "_genomic.gff.gz",
    TRUE ~ NA_character_
    )

  assertthat::assert_that(!is.na(assembly_file_extension), msg = paste0("Failed to recognise assembly_file_type == ", assembly_file_type))
  assertthat::assert_that(target_assembly_accession %in% refseq_data_frame[["assembly_accession"]])

  message("Seaching for ", assembly_file_type, " of assembly: ", target_assembly_accession)

  refseq_data_frame = dplyr::filter(
      .data = refseq_data_frame,
      assembly_accession %in% target_assembly_accession
      )

  assertthat::assert_that(nrow(refseq_data_frame) != 0, msg = paste0("Could not find assembly [", target_assembly_accession, "] in the local cache of refseq assemblies. If you're sure the assembly accession code is correct, try updating the local cache by running `update_refseq_data_cache()`"))
  assertthat::assert_that(nrow(refseq_data_frame) <= 1, msg = paste0("Found multiple  entries for assembly [", target_assembly_accession, "] in the local cache. This is unexpected ... please send assembly code used and this error message to package maintainer"))

  sequence_file_name = basename(refseq_data_frame[["ftp_path"]])
  ftp_root = refseq_data_frame[["ftp_path"]]
  assembly_file_path = paste0(ftp_root, "/", sequence_file_name, assembly_file_extension)

  assertthat::assert_that(dir.exists(output_folder), msg = utilitybeltassertions::fmterror("Could not find folder: ", output_folder))
  full_dest_filepath = paste0(output_folder, "/", basename(assembly_file_path))

  message("downloading assembly file [type = ",assembly_file_type ,"]:\n\t[", target_assembly_accession, "]\n\nfrom the RefSeq ftp link \n\t[", assembly_file_path, "]\nto:\t", full_dest_filepath)

  utils::download.file(assembly_file_path, destfile =  full_dest_filepath)

  return(full_dest_filepath)
}


#' Download reference genome
#'
#' Download the genome assembly
#'
#' @inheritParams download_assembly_files
#' @inherit download_assembly_files
#' @export
refseq_download_genome_assembly <- function(target_assembly_accession, output_folder = getwd(), refseq_data_frame = load_refseq_data_frame_from_cache()){
  download_assembly_files(target_assembly_accession = target_assembly_accession, assembly_file_type = "genome", output_folder = output_folder, refseq_data_frame = refseq_data_frame)
}

#' download_assembly_stats
#'
#' @inheritParams download_assembly_files
#' @inherit download_assembly_files
#' @export
#' @examples
#' refseq_download_genome_assembly_stats("GCF_000371765.1")
refseq_download_genome_assembly_stats <- function(target_assembly_accession, output_folder = getwd(), refseq_data_frame = load_refseq_data_frame_from_cache()){
  download_assembly_files(target_assembly_accession = target_assembly_accession, assembly_file_type = "stats", output_folder = output_folder, refseq_data_frame = refseq_data_frame)
}

#' Download_best_assembly
#'
#' @inheritDotParams refseq_download_genome_assembly
#' @inheritParams choose_best_assembly
#' @inherit refseq_download_genome_assembly return
#' @export
#'
download_best_assembly <- function(taxid_of_interest, download_stats = TRUE, ...){
  best_asssembly=target_assembly_accession = choose_best_assembly(taxid_of_interest = taxid_of_interest, return_accession_only = TRUE, log_chosen_assembly = TRUE)
  message("Attempting to download the best assembly for taxid ", taxid_of_interest, " (", best_asssembly, ")")
  res=refseq_download_genome_assembly(target_assembly_accession = best_asssembly, ...)

  if(download_stats){
    stats = refseq_download_genome_assembly_stats(target_assembly_accession = best_asssembly, ...)
  }

  return(res)
}


# Dependent on CLI tools---------------------------------------------------------------------



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


#' Prep a refseq assembly for alignment
#'
#' @param path_to_assembly_fasta path to assembly fasta
#'
#' @return prefix of bwa-indexed (character)
#' @export
#'
prep_for_alignment <- function(path_to_assembly_fasta){
  utilitybeltassertions::assert_filenames_have_valid_extensions(
    filenames = path_to_assembly_fasta,
    valid_extensions = c("fa.gz", "fasta.gz", "fa", "fasta", "fna", "fna.gz"),
    ignore_case = TRUE
    )

  utilitybeltassertions::assert_program_exists_in_path(program_names = "bwa-mem2")

  indexing_command = paste0("bwa-mem2 index ", path_to_assembly_fasta)

  message("Running: ", indexing_command)
  return_code = system(indexing_command)

  if(return_code != 0){
   stop(paste0("bwa-mem2 indexing of reference genome has failed. exit code ", return_code))
  }

  message(utilitybeltassertions::fmtsuccess("\nReference Genome Indexing Complete (bwa2-mem index)"))
  return(path_to_assembly_fasta)
}

#' Align Reads
#'
#' @param bwa_index_prefix prefix of indexed reference genome (e.g. output of prep_for_alignment()) (string)
#' @param outfile_prefix prefix used for output alignment files
#' @param forward_reads path to fastq file containing forward reads (string)
#' @param reverse_reads path to fastq file containing reverse reads (string)
#' @param bwa_options supply any other options to bwa-mem2 mem command. For example '-k 50' will change min seed length to 50. (string)
#' @param threads how many cup threads to use (int)
#'
#' @return outfile path
#' @export
#'
#' @examples
#' \dontrun{
#'
#' prep_for_alignment()
#' }
align_reads <- function(bwa_index_prefix, forward_reads, reverse_reads, outfile_prefix, bwa_options="", threads=parallel::detectCores()){

  utilitybeltassertions::assert_program_exists_in_path(program_names = c("bwa-mem2", "samtools"))

  utilitybeltassertions::assert_filenames_have_valid_extensions(
    filenames = c(forward_reads, reverse_reads),
    valid_extensions = c("fastq", "fastq.gz", "fq", "fq.gz"),
    ignore_case = TRUE
    )

  utilitybeltassertions::assert_files_exist(filepaths = c(
    paste0(bwa_index_prefix, ".0123"),
    paste0(bwa_index_prefix, ".amb"),
    paste0(bwa_index_prefix, ".ann"),
    paste0(bwa_index_prefix, ".pac"),
    paste0(bwa_index_prefix, ".bwt.2bit.64"),
    paste0(bwa_index_prefix, ".bwt.8bit.32")
    ), supplementary_error_message = "Please run prep_for_alignment() on your refseq assembly to generate these files")


  message("Running all processes using: ", threads, " cores")

  outfile_name = paste0(outfile_prefix, ".sorted.bam")

  if(bwa_options != ""){
    bwa_options = paste0(bwa_options, " ")
  }

  alignment_command = paste0(
    "bwa-mem2 mem -t ", threads, " ", bwa_options, bwa_index_prefix, " ", forward_reads, " ", reverse_reads, " | samtools sort -O bam -@ ", threads, " -o ", outfile_name
    )

  message("Running: ", alignment_command)

  exit_code = system(alignment_command)

  #exit code 6 doesnt effect actual alignment results - just the description of time taken after. Fixed in latest version of bwa-mem but these are not readily installable on mac
  assertthat::assert_that(exit_code %in% c(0, 6), msg = utilitybeltassertions::fmterror("Alignment failed. Exit code [", exit_code,"]"))

  outfile_path = normalizePath(outfile_name)

  message(utilitybeltassertions::fmtsuccess("\n\nAlignment Complete. Result at: ", outfile_path))

  return(outfile_path)
}



bam_index_create_command <- function(path_to_sorted_bam, threads){
  indexing_command = "samtools index -@ <threads> <path_to_sorted_bam>" %>%
    sub(pattern = "<threads>", replacement = threads) %>%
    sub(pattern = "<path_to_sorted_bam>", replacement = path_to_sorted_bam)
  return(indexing_command)
}

#' Bam Index
#'
#' Index a bam file
#'
#' Produce a bai index for a sorted bam. Requires samtools to be installed.
#'
#' @param path_to_sorted_bam path to a coordinate sorted bam file (string)
#' @param threads number of cores to use for multithreading (int)
#'
#' @return path to sorted bam (string)
#' @export
#'
#' @examples
#' \dontrun{
#' bam_index("path/to/sorted/bam")
#' }
bam_index <- function(path_to_sorted_bam, threads = parallel::detectCores()){
  utilitybeltassertions::assert_program_exists_in_path(program_names = c("samtools"))
  utilitybeltassertions::assert_files_exist(path_to_sorted_bam)


  indexing_command = bam_index_create_command(path_to_sorted_bam, threads)
  message("Running: ", indexing_command)
  exit_code = system(indexing_command)

  assertthat::assert_that(exit_code == 0, msg = utilitybeltassertions::fmterror("Indexing of bam file [",path_to_sorted_bam,"] failed. Exit code [", exit_code, "]"))
  message(utilitybeltassertions::fmtsuccess("\n\nIndexing Complete"))

  return(path_to_sorted_bam)
}




