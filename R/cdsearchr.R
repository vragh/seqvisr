#' Submit a set of queries to NCBI's CD-SEARCH via R.
#'
#' @description
#' Core internal \code{\link{cdsearchr}} function that prepares and submit queries to the
#' CD-SEARCH webserver located at a specified URL. Uses httr for POST\'ing to
#' the server.
#'
#' @usage cdssubmit(queries = NA, cdsurl = NULL,
#' db = c("cdd", "pfam", "smart", "tigrfam", "cog", "kog"),
#' smode = c("auto", "prec", "live"), evalue = 0.01, useid1 = TRUE, compbasedadj = 1,
#' biascompfilter = TRUE, tdata = c("hits", "aligns", "feats"), alnfmt = NA,
#' dmode = c("rep", "std", "full"), qdefl = TRUE, cddefl = TRUE, maxhit = 500,
#' nseqs = NULL)
#'
#' @param queries (character string, mandatory) path to a FASTA file containing the
#' query protein sequences.
#'
#' @param cdsurl (character string, mandatory) the URL at which the remote CD-SEARCH
#' server is accessible. Should be inherited from [`cdsearchr`]. (Set to
#' NULL by default.)
#'
#' @param db (character string, optional) controls which databases CD-SEARCH should
#' search the queries against. Please refer to "database selection" under the URL
#' https://www.ncbi.nlm.nih.gov/Structure/cdd/cdd_help.shtml#BatchRPSBSearchMode
#' for particulars on the databases. This parameter only has an effect if `smode` is
#' set to "live" (see below). (Set to "cdd" by default.)
#'
#' @param smode (character string, optional) controls which search mode CD-SEARCH
#' should use. "auto" will check the queries first against a set of precalculated
#' results (by checking query identifiers; really only works if these are sequences
#' in NCBI already), and if that fails, it performs a "live" search against the
#' CD-SEARCH database. "prec" would return results only for queries that have a
#' result in the precalculated database. "live" will search every query anew against
#' its databases even if precalculated results exist for that query. (Set to "auto"
#' by default.)
#'
#' @param evalue (numeric, optional) expect value (statistical significance threshold)
#' used for filtering and reporting annotation matches. (Set to 0.01 by default.)
#'
#' @param useid1 (binary, optional) controls whether queries should also be
#' searched against archived sequence identifiers if the query's identifier (if
#' it happens to be an NCBI identifier) does not match anything in the current
#' Entrez Protein database records. (Set to TRUE by default.)
#'
#' @param compbasedadj (integer, optional) should CD-SEARCH use compositionally-
#' corrected scoring? (0 - correction turned off; 1 - correction turned on.) (Set
#' to 1 by default.)
#'
#' @param biascompfilter (binary, optional) should compositionally biased regions
#' of the queries be filtered out? (Set to TRUE by default.)
#'
#' @param tdata (character string, optional) what type of target data should be
#' returned: "hits" (domain hits), "aligns" (domain alignments), or "feats" (domain
#' features). Changing from the default might break functionality as of the current
#' version of `seqvisr`. (Set to "hits" by default.)
#'
#' @param alnfmt (character string, optional) data format to be used for downloading
#' alignment data in the event `tmode` is set to "aligns". This will never be the
#' case for cdsearchr, and this option exists only for the sake of completeness.
#' (Set to `NA` by default.)
#'
#' @param dmode (character string, optional) which data mode must be used for the results.
#' This dictates what set of domains are returned: the highest scoring hit for each region
#' of the sequence ("rep"), the best hits from each database available in CD-SEARCH (so
#' multiple hits per query region are possible; "std"), or all hits ("full"). (Set to
#' "rep" by default.)
#'
#' @param qdefl (binary, optional) should query titles be included in the results? (Set to
#' TRUE by default.)
#'
#' @param cddefl (binary, optional) should domain titles be included in the results? (Set
#' to FALSE by default.)
#'
#' @param maxhit (integer, optional) maximum number of results per query that should be retrieved.
#' Only matters if `smode` is set to "live".
#'
#' @param nseqs (integer, mandatory) number of sequences being submitted to the CD-SEARCH
#' server. This value should be automatically calculated and passed on from [`cdsearchr`].
#' (Set to NULL by default.)
#'
#' @returns
#' A \pkg{httr} response object containing (among other things, but most importantly) the
#' unique `cdsid` referring to the submitted queries and search request which can subsequently
#' be used to query the server.
#'
#' @importFrom rlang .data flatten
#'
#' @import httr
#'
#'
cdssubmit <- function(queries = NA, cdsurl = NULL, db = c("cdd", "pfam", "smart", "tigrfam", "cog", "kog"),
                      smode = c("auto", "prec", "live"), evalue = 0.01, useid1 = TRUE, compbasedadj = 1,
                      biascompfilter = TRUE, tdata = c("hits", "aligns", "feats"), alnfmt = NA,
                      dmode = c("rep", "std", "full"), qdefl = TRUE, cddefl = TRUE, maxhit = 500, nseqs = NULL){

  # #For testing.
  # db <- "cdd"; smode <- "auto"; useid1 <- TRUE; compbasedadj <- 1; biascompfilter <- TRUE; evalue <- 0.01;
  #tdata <- "hits"; alnfmt <- NA; dmode <- "rep"; qdefl <- TRUE; cddefl <- TRUE; queries <- "test1.fasta"

  if(base::is.na(queries)){
    base::stop("Please provide a path to the file containing a set of queries (FASTA sequences or NCBI accessions).\n")
  }

  if(base::is.null(cdsurl)){
    base::stop("No URL to the CD-SEARCH API supplied.\n")
  }

  #Converting binary logicals to text.
  cddefl <- base::ifelse(cddefl, "true", "false")
  qdefl <- base::ifelse(qdefl, "true", "false")
  useid1 <- base::ifelse(useid1, "true", "false")
  biascompfilter <- base::ifelse(biascompfilter, "true", "false")

  #Bundling the CD-Search parameters into a named list.
  cdsparams <- base::list("db" = db, "smode" = smode, "useid1" = useid1, "compbasedadj" = compbasedadj,
                          "filter" = biascompfilter, "evalue" = evalue, "tdata" = tdata, "alnfmt" = alnfmt,
                          "dmode" = dmode, "qdefl" = qdefl, "cddefl" = cddefl, "maxhit" = maxhit)

  #Removing parameters that are set to NA, i.e., not in use (alnfmt basically).
  cdsparams <- cdsparams[-base::which(base::sapply(cdsparams, function(x){base::any(base::is.na(x))}), arr.ind = TRUE)]

  #Constructing the full set of parameters that will be passed to httr::POST.
  urlparams <- rlang::flatten(base::list(cdsparams, base::list(queries = httr::upload_file(path = queries, type = "application/fasta"))))

  base::cat("Submitting ", nseqs, " queries to the CD-SEARCH webserver.\n")

  #Submitting the queries to CD-SEARCH.
  cdsresp1 <- httr::POST(url = cdsurl, body = urlparams)

  return(cdsresp1)

}




#' Check for the status of a submitted search on NCBI's CD-SEARCH.
#'
#' @description
#' Check on the status of a submitted search associated with the given `cdsid`
#' (@seealso [`seqvisr::cdssubmit`]). The values for all the parameters below
#' are inherited from [`cdsearchr`] (@seealso [cdsearchr()] for some details).
#' This function is also responsible for interpreting (and printing to the
#' terminal) the CD-SEARCH-specific error-codes.
#'
#' @usage cdscheck(cdsid_check = NULL, tdata_check = NULL, cdsurl_check = NULL,
#' maxreps = NULL, sleep_duration = NULL)
#'
#' @param cdsid_check (character string, mandatory) the `cdsid` returned by
#' the NCBI CD-SEARCH server upon successful submission of a search request.
#'
#' @param tdata_check (character string, mandatory) the format in which the
#' target data should be retrieved.
#'
#' @param cdsurl_check (character string, mandatory) URL at which the CD-SEARCH
#' server can be accessed.
#'
#' @param maxreps (integer, mandatory) maximum number of attempts to be made to
#' query the server on the status of the search before giving up.
#'
#' @param sleep_duration (integer, mandatory) maximum duration in seconds between
#' successive queries to the server.
#'
#' @returns 1 if the search terminated successfully, otherwise 0.
#'
#' @import httr
#'
#' @importFrom stringr str_extract str_detect
#'
cdscheck <- function(cdsid_check = NULL, tdata_check = NULL, cdsurl_check = NULL, maxreps = NULL, sleep_duration = NULL){

  if(base::is.null(cdsid_check) | base::is.null(cdsurl_check)) { base::stop("No CDS ID or URL provided. Aborting.") }
  if(base::is.null(tdata_check)) { base::warning("No value provided for tdata, using default \"hits\""); tdata_check <- "hits" }

  #Keep checking on the search status until it returns status 0.
  search_success <- 0
  counter <- 0
  if(base::is.null(sleep_duration)) { sleep_duration <- 10 }
  if(base::is.null(maxreps)) { maxreps <- 20 }

  base::cat("Search in progress.\n")
  base::cat("Querying search status in ", sleep_duration, " seconds.\n")
  base::cat("Will query a total of ", maxreps, " times.\n")

  while(search_success == 0 & counter < maxreps){

    #if(counter > 0){
    base::Sys.sleep(sleep_duration)
    #}

    #Query the server with the CDS ID and extracting the status of the search.
    cdsresp2 <- httr::POST(url = cdsurl_check, body = list("cdsid" = cdsid_check, "tdata" = tdata_check))
    cdsresp2_cont <- httr::content(cdsresp2, as = "text", encoding = "utf8")

    cdsstat <- (-1)
    cdsstat <- stringr::str_extract_all(cdsresp2_cont, "status\\s+[\\d]")
    cdsstat <- stringr::str_extract(cdsstat, "\\d")

    #cat("Current status is: ", cdsstat, ".\n")

    if(cdsstat == 0){

      # cat("cdsid ", cdsid_check, " tdata", tdata_check, " url", cdsurl_check, "\n")
      # cat(cdsresp2_cont)
      # View(data.table::fread(text = cdsresp2_cont, skip = 7, sep = "\t"))
      # cat(stringr::str_detect(cdsresp2_cont, "Superfamily$|Definition$"))

      if(stringr::str_detect(cdsresp2_cont, "Superfamily$|Definition$")){
        base::cat("Server returned no annotations. Please check manually with the URL indicated earlier.\n")
        break
      } else{
        base::cat("Remote search on CD-SEARCH completed!! Will retrieve the results now.\n")
        search_success <- 1
        break
      }
    } else if(cdsstat == 1){
      base::message("Error: invalid CD-SEARCH ID.\n")
      break
    } else if(cdsstat == 2){
      base::message("Error: invalid input to CD-SEARCH for retrieval (missing query information or CD-SEARCH ID).")
      break
    } else if(cdsstat == 3){
      base::cat("Search in progress.\n")
    } else if(cdsstat == 4){
      base::message("Error: query manager service error.")
      break
    } else if(cdsstat == 5){
      base::message("Error: the data is corrupted or no longer available. Please re-run/re-submit the search.")
      break
    } else{
      base::message("Error: unknown error. These were the contents of the response from the CD-SEARCH server:\n",
                    cdsresp2_cont, "\n#-----------#")
      break
    }

    counter = counter + 1
    if(counter == maxreps){
      base::message("Exceeded maximum number of attempts (", maxreps, ") to retrieve the search results. Aborting.")
      break
    }

  }

  return(search_success)

}




#' Retrieve annotation results from CD-SEARCH.
#'
#' @description
#' cdsearchr() internal function that retrieves the (HTML-formatted)
#' tabular data containing the annotations for a given `cdsid`
#' (@seealso [cdssubmit()]). All parameters are inherited from
#' [cdsearchr()].
#'
#' @usage cdsretrieve(cdsid, tdata, qdefl, cddefl, dmode, cdsurl)
#'
#' @param cdsid (character string, mandatory) the `cdsid` returned by
#' the NCBI CD-SEARCH server upon successful submission of a search request.
#'
#' @param tdata (character string, mandatory) the format in which the
#' target data should be retrieved.
#'
#' @param qdefl (binary, optional) should query titles be included in the
#' results? (Set to TRUE by default.)
#'
#' @param cddefl (binary, optional) should domain titles be included in the
#' results? (Set to FALSE by default.)
#'
#' @param dmode (character string, optional) which data mode must be used for
#' the results.
#'
#' @param cdsurl (character string, mandatory) URL at which the CD-SEARCH
#' server can be accessed.
#'
#' @returns
#' A [`httr`] response object containing a (HTML-formatted) text table of
#' annotation results from CD-SEARCH.
#'
#' @import httr
#'
cdsretrieve <- function(cdsid, tdata, qdefl, cddefl, dmode, cdsurl){

  cdsresp3 <- httr::POST(url = cdsurl,
                         body = list("cdsid" = cdsid, "tdata" = tdata,
                                     "qdefl" = qdefl, "cddefl" = cddefl,
                                     "dmode" = dmode))
  cdsresp3_cont <- httr::content(cdsresp3, as = "text", encoding = "utf8")

  return(cdsresp3_cont)

}


#' Extract CD-SEARCH results as data.frame.
#'
#' @description
#' [`cdsearchr`] internal function that takes the HTML-formatted results from
#' [`cdsretrieve`] and converts it into a [`base::data.frame`] using
#' [`data.table::fread`].
#'
#' @usage cdsrestodf(content, rowstoskip = 7, sep = "\t")
#'
#' @param content ([`httr`] response object, mandatory) [`httr`] response object
#' containing the annotation results from the CD-SEARCH server.
#'
#' @param rowstoskip (numeric, mandatory) the response contains header rows for
#' the results table. These need to be skipped. This indicates how many such rows
#' there are.
#'
#' @param sep (character, mandatory) separator symbol delineating successive rows
#' of the results table.
#'
#' @returns A [`base::data.frame`].
#'
#' @importFrom data.table fread
cdsrestodf <- function(content, rowstoskip = 7, sep = "\t"){

  df <- data.table::fread(text = content, skip = rowstoskip, sep = sep)

  return(df)

}




#' Access NCBI's (batch) CD-SEARCH from R.
#'
#' @description
#' [`cdsearchr()`] provides an `R` interface for NCBI's CD-SEARCH sequence annotation tool.
#' It takes the path to a FASTA file containing the query protein sequences as input
#' and returns a [`data.frame`] containing the annotation results.
#'
#' @note cdsearchr() respects the CD-SEARCH API's 4000 queries submission limit.
#' Therefore, users should pre-chunk FASTA files into files of 4000 sequences or
#' fewer before passing them on cdsearchr.#'
#'
#' @usage cdsearchr(queries = NA, db = c("cdd", "pfam", "smart", "tigrfam", "cog", "kog"),
#' smode = c("auto", "prec", "live"), useid1 = TRUE, compbasedadj = 1,
#' biascompfilter = TRUE, evalue = 0.01, tdata = c("hits", "aligns", "feats"),
#' alnfmt = NA, dmode = c("rep", "std", "full"), qdefl = TRUE, cddefl = FALSE,
#' maxhit = 500, check_max = 10, check_wait = 20)
#'
#' @param queries (character string, mandatory) path to a FASTA file containing the
#' query protein sequences.
#'
#' @param db (character string, optional) controls which databases CD-SEARCH should
#' search the queries against. Please refer to "database selection" under the URL
#' <https://www.ncbi.nlm.nih.gov/Structure/cdd/cdd_help.shtml#BatchRPSBSearchMode>
#' for particulars on the databases. This parameter only has an effect if `smode` is
#' set to "live" (see below). (Set to "cdd" by default.)
#'
#' @param smode (character string, optional) controls which search mode CD-SEARCH
#' should use. "auto" will check the queries first against a set of precalculated
#' results (by checking query identifiers; really only works if these are sequences
#' in NCBI already), and if that fails, it performs a "live" search against the
#' CD-SEARCH database. "prec" would return results only for queries that have a
#' result in the precalculated database. "live" will search every query anew against
#' its databases even if precalculated results exist for that query. (Set to "auto"
#' by default.)
#'
#' @param useid1 (binary, optional) controls whether queries should also be
#' searched against archived sequence identifiers if the query's identifier (if
#' it happens to be an NCBI identifier) does not match anything in the current
#' Entrez Protein database records. (Set to TRUE by default.)
#'
#' @param compbasedadj (integer, optional) should CD-SEARCH use compositionally-
#' corrected scoring? (0 - correction turned off; 1 - correction turned on.) (Set
#' to 1 by default.)
#'
#' @param biascompfilter (binary, optional) should compositionally biased regions
#' of the queries be filtered out? (Set to TRUE by default.)
#'
#' @param evalue (numeric, optional) expect value (statistical significance threshold)
#' used for filtering and reporting annotation matches. (Set to 0.01 by default.)
#'
#' @param tdata (character string, optional) what type of target data should be
#' returned: "hits" (domain hits), "aligns" (domain alignments), or "feats" (domain
#' features). Changing from the default might break functionality as of the current
#' version of `seqvisr`. (Set to "hits" by default.)
#'
#' @param alnfmt (character string, optional) data format to be used for downloading
#' alignment data in the event `tmode` is set to "aligns". This will never be the
#' case for cdsearchr, and this option exists only for the sake of completeness.
#' (Set to `NA` by default.)
#'
#' @param dmode (character string, optional) which data mode must be used for the results.
#' This dictates what set of domains are returned: the highest scoring hit for each region
#' of the sequence ("rep"), the best hits from each database available in CD-SEARCH (so multiple
#' hits per query region are possible; "std"), or all hits ("full"). (Set to "rep" by default.)
#'
#' @param qdefl (binary, optional) should query titles be included in the results? (Set to TRUE
#' by default.)
#'
#' @param cddefl (binary, optional) should domain titles be included in the results? (Set to FALSE
#' by default.)
#'
#' @param maxhit (integer, optional) maximum number of results per query that should be retrieved.
#' Only matters if `smode` is set to "live".
#'
#' @param check_max (numeric, optional) how many times should cdsearchr() query for results before
#' giving up? (Set to 10 attemps by default.)
#'
#' @param check_wait (numeric, optional) how long -- in seconds -- must cdsearchr() wait between
#' successive requests to the CD-SEARCH API while querying for the results. (Set to 20 by
#' default.)
#'
#' @details
#' cdsearchr() is an `R`-based interface to the NCBI CD-SEARCH application. It uses `httr` internally
#' to submit and retrieve data. Once the queries have been submitted, cdsearchr() will repeatedly
#' query the CD-SEARCH server until it receives a response (success/failure) or the number of attempts
#' exceeds `check_max`. Although `check_wait` has been set to 20 (seconds) it is recommended that the
#' user adjusts this based on the size of the data set (set to a smaller value for smaller data sets
#' and vice versa).
#'
#' @examples
#' \dontrun{
#' inpath <- system.file("extdata", "cdsearchr_testdata.fasta", package = "seqvisr", mustWork = TRUE)
#' cdsearchr(queries = inpath, check_wait = 2)
#' }
#'
#' @import httr
#'
#' @importFrom rlang .data flatten
#'
#' @importFrom stringr str_extract str_detect
#'
#' @importFrom data.table fread
#'
#' @export
#'
cdsearchr <- function(queries = NA, db = c("cdd", "pfam", "smart", "tigrfam", "cog", "kog"),
                      smode = c("auto", "prec", "live"), useid1 = TRUE, compbasedadj = 1,
                      biascompfilter = TRUE, evalue = 0.01, tdata = c("hits", "aligns", "feats"),
                      alnfmt = NA, dmode = c("rep", "std", "full"), qdefl = TRUE, cddefl = FALSE,
                      maxhit = 500, check_max = 10, check_wait = 20){

  #Filling in the "missing" defaults for arguments in case the user hasn't supplied them.
  db <- base::match.arg(db)
  smode <- base::match.arg(smode)
  tdata <- base::match.arg(tdata)
  dmode <- base::match.arg(dmode)

  #Check if input file has <4000 sequences (maximum allowed by the API),
  #If not, ask the user to chunk the input file, and halt execution.
  os_is_win <- stringr::str_detect(base::Sys.info()[1], "Windows")
  if(os_is_win){
    nseqs <- as.numeric(base::system(command = base::paste0("(Select-String -Path \"", queries, "\" -Pattern \"^>\" -AllMatches).Matches.Count"), intern = TRUE))
  } else{
    nseqs <- as.numeric(base::system(command = base::paste0("grep -c \"^>\" ", queries), intern = TRUE))
  }

  if(nseqs > 4000){

    stopmsg <- base::paste0("The CD-SEARCH webserver has a limit of 4000 queries per submission. The query file contains ", nseqs, " queries.\nPlease provide an input file with fewer than 4000 queries.")
    base::message(stopmsg)

  } else{

    #----#
    #Executing main search submission and retrieval functions.
    #----#

    #Base CD-Search URL
    burl <- "https://www.ncbi.nlm.nih.gov/Structure/bwrpsb/bwrpsb.cgi"

    #For testing only.
    #cdsresp <- cdssubmit(queries = queries, cdsurl = cdsurl)
    #For testing.
    #db <- "cdd"; smode <- "auto"; useid1 <- TRUE; compbasedadj <- 1; biascompfilter <- TRUE; evalue <- 0.01; tdata <- "hits"; alnfmt <- NA; dmode <- "rep"; qdefl <- TRUE; cddefl <- TRUE; queries <- "/home/owner/Nextcloud/laptop_rplace/test.fasta"; nseqs <- 4; check_max <- 4; check_wait <- 4;

    base::message("cdsearchr")

    cdsresp <- cdssubmit(queries = queries, db = db, smode = smode, useid1 = useid1,
                         compbasedadj = compbasedadj, biascompfilter = biascompfilter,
                         tdata = tdata, alnfmt = alnfmt, dmode = dmode, qdefl = qdefl,
                         cddefl = cddefl, maxhit = maxhit, cdsurl = burl, nseqs = nseqs)

    cdssubmit_statcode <- httr::status_code(cdsresp)

    if(cdssubmit_statcode == 200){

      cdsid <- stringr::str_extract(httr::content(cdsresp, as = "text", encoding = "utf8"),
                                    "QM3\\-qcdsearch\\-[A-Z0-9]+\\-[A-Z0-9]+(?=\\n)")

      base::cat("Queries submitted successfully.\n")
      base::cat(base::paste0("If cdsearchr fails, the status of the search can be checked manually at the following URL:\nhttps://www.ncbi.nlm.nih.gov/Structure/bwrpsb/bwrpsb.cgi?cdsid=", cdsid), "\n")

      cdsstatcheck <- cdscheck(cdsid_check = cdsid, tdata_check = tdata, cdsurl_check = burl,
                               maxreps = check_max, sleep_duration = check_wait)

      if(cdsstatcheck == 1){
        cdsresp <- cdsretrieve(cdsid = cdsid, tdata = tdata, qdefl = qdefl,
                               cddefl = cddefl, dmode = dmode, cdsurl = burl)

        outdf <- cdsrestodf(cdsresp)

        base::cat("Sequence annotations successfully retrieved.\n")

        return(outdf)
      }
    } else{
      base::message("Something went wrong. The server returned status code ", cdssubmit_statcode)
    }
  }
}


####
