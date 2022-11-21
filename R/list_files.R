#' Wrapper around \code{base::list.files()} that adds an option to exclude directories
#' automatically.
#'
#' @description
#' \link[base]{list.files} returns all objects in the search path including directories.
#' \code{list_files} is a wrapper that adds an option to exclude directories. In specific,
#' it allows the user to override the default behavior of \code{list.files()} including
#' subdirectory names when \code{recursive} is set to FALSE. This is controlled by the
#' \code{incl_dirs} argument which filters our directories when set to TRUE.
#'
#' @usage list_files(path = ".", pattern = NULL, all.files = FALSE, full.names = TRUE,
#' recursive = FALSE, ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE,
#' incl_dirs = FALSE)
#'
#' @inheritParams base::list.files
#'
#' @param incl_dirs (binary, optional) set to FALSE to remove all directories from
#' the search results irrespective of whether or not \code{recursive=TRUE}.
#'
#' @examples
#' \dontrun{
#' #Input data
#' inpath <- dirname(system.file("extdata", "cdsearchr_testdata.fasta",
#'                               package = "seqvisr", mustWork = TRUE))
#' #Listing files
#' list_files(inpath)
#' }
#'
#' @seealso \link[base]{list.files} for documentation regarding all other arguments.
#'
list_files <- function(path = ".", pattern = NULL, all.files = FALSE, full.names = TRUE,
                       recursive = FALSE, ignore.case = FALSE, include.dirs = FALSE, no.. = FALSE,
                       incl_dirs = FALSE){
  #path <- paste0(mypath, "/", refpath)
  if(path == ".") { path = base::getwd() }

  #Include directories if recursive is set.
  if(incl_dirs & recursive){include.dirs = TRUE}

  files <- base::list.files(path = path, pattern = pattern, all.files = all.files, full.names = TRUE,
                      recursive = recursive, ignore.case = ignore.case, include.dirs = include.dirs,
                      no.. = no..)

  if(!incl_dirs){
    files <- files[!base::dir.exists(files)]
  }

  if(!full.names){
    return(base::basename(files))
  } else{
    return(files)
  }

}
