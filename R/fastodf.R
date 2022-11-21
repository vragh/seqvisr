#' Read in a FASTA file as a data.frame.
#'
#' @description
#' Utility function that uses \link[seqinr]{read.fasta} to load a FASTA file as a
#' \link[base]{data.frame}.
#'
#' @usage fastodf(path = NULL, seqtype = c("DNA", "AA"), incl_filepath = TRUE)
#'
#' @param path (character string, mandatory) the path to the input FASTA file.
#'
#' @param seqtype (character string, optional) the type of sequence being read in;
#' DNA ("DNA") or amino acid ("AA"). Defaults to "DNA".
#'
#' @param incl_filepath (logical, optional) should the path to the file being read
#' be included in the \code{data.frame}. Irrespective of whether the user sets this
#' to TRUE or FALSE, a column (\code{filename}) will be included in the output
#' \code{data.frame} to keep the number of columns consistent.
#'
#' @returns a \link[base]{data.frame} with the following columns: \code{seqname},
#' \code{seq}, and \code{filename}. If \code{incl_filepath} is set to TRUE, then
#' \code{filename} will include the full path to the input file. It will be set to
#' \code{NA} otherwise.
#'
#' @examples
#' \dontrun{
#' #Input data
#' inpath <- system.file("extdata", "cdsearchr_testdata.fasta", package = "seqvisr", mustWork = TRUE)
#' #Reading in some sample amino acid sequences
#' fastodf(path = inpath, seqtype = "AA")
#' }
#'
#' @importFrom seqinr read.fasta
#'
#' @export
#'
fastodf <- function(path = NULL, seqtype = c("DNA", "AA"), incl_filepath = TRUE){

  #require(seqinr)
  df <- seqinr::read.fasta(file = path, as.string = TRUE, whole.header = TRUE, seqtype = seqtype)
  df <- base::data.frame(seqname = base::names(base::unlist(df)),
                       seq = base::unlist(df),
                       stringsAsFactors = FALSE)
  if(incl_filepath){
    df$filename <- base::rep_len(path, length.out = base::nrow(df))
  } else{
    df$filename <- NA
  }

  base::row.names(df) <- NULL

  return(df)

}


#' Read in all FASTA files in a directory as a data.frame.
#'
#' @description
#' Wrapper around \link{fastodf} that reads in all FASTA files in a
#' directory into a single \link[base]{data.frame}. Uses \link{list_files}
#' internally also. At present, \code{fasdirf} requires that all files being read in
#' are of the same type (DNA or amino acid).
#'
#' @usage fasdirdf(path = NULL, pat = NULL, seqtype = c("DNA", "AA"), incl_filepath = TRUE)
#'
#' @param path (character string, mandatory) the path to the directory containing
#' the input FASTA files.
#'
#' @param pat (character string, optional) a regex string as used by \link[base]{list.files}
#' to specify which file/directory names should be returned.
#'
#' @param seqtype (character string, optional) the type of sequence being read in;
#' DNA ("DNA") or amino acid ("AA"). Defaults to "DNA".
#'
#' @param incl_filepath (logical, optional) should the path to the file being read
#' be included in the \code{data.frame}. Irrespective of whether the user sets this
#' to TRUE or FALSE, a column (\code{filename}) will be included in the output
#' \code{data.frame} to keep the number of columns consistent.
#'
#' @returns a \link[base]{data.frame} with the following columns: \code{seqname},
#' \code{seq}, and \code{filename}. If \code{incl_filepath} is set to TRUE, then
#' \code{filename} will include the full path to the input file. It will be set to
#' \code{NA} otherwise.
#'
#' @examples
#' \dontrun{
#' #Input data
#' inpath <- dirname(system.file("extdata", "cdsearchr_testdata.fasta",
#'                               package = "seqvisr", mustWork = TRUE))
#' #Reading in some sample amino acid sequences
#' fasdirdf(path = inpath, seqtype = "AA", pat = "*_testdata.fasta")
#' }
#'
#' @importFrom seqinr read.fasta
#'
#' @export
#'
fasdirdf <- function(path = NULL, pat = NULL, seqtype = c("DNA", "AA"), incl_filepath = TRUE){

  if(is.null(path)){ path <- getwd() }
  if(is.null(pat)){ pat <- "*" }

  return(do.call("rbind", lapply(list_files(path = path, pattern = pat, full.names = TRUE, incl_dirs = FALSE), fastodf, seqtype = seqtype, incl_filepath = incl_filepath)))
}
