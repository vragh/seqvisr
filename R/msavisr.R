#'Multiple sequence alignment (MSA) visualization and manipulation thereof.
#'
#'Given a fasta formatted MSA, msavisr will attempt to produce a visualization of the MSA
#'with matches, mismatches, gaps, and (optional) regions of interest (ROIs) being highlighted
#'in different colors. The matches, mismatches, and gaps are calculated by comparing all other
#'sequences to a designated reference sequence within the MSA. ROIs are specified manually (see
#'below).
#'
#'
#' @usage msavisr(mymsa = "name_of_msa_file.fasta", myref = "fasta_header_of_reference_sequence",
#' refontop = TRUE, myroi = NULL, hnon = NULL, hmat = NULL, hroi = NULL, wnon = NULL, wmat = NULL,
#' wroi = NULL, #' anon = NULL, amat = NULL, aroi = NULL, basecolors = NULL, cbfcols = FALSE)
#'
#' @param mymsa (character string, mandatory) the name of the fasta-formatted file containing the
#' multiple sequence alignment (MSA). The full path to the file can also be provided here (in which
#' case, the mypath argument must be set to NULL).
#'
#' @param myref (character string, mandatory) the fasta header (not including the ">") of one of the
#' sequences in the MSA which act as the reference sequence.
#'
#' @param mypath (character string, optional) the path to the directory containing the MSA. (Default
#' : looks in the current directory.)
#'
#' @param refontop (boolean, optional) should the reference sequence be placed at the top of the MSA plot
#' (default) or at the bottom? (Set to FALSE for bottom.)
#'
#' @param myroi (list of vectors, optional) the user can manually specify regions of interest (ROI;
#' i.e., positions in the MSA) that they wish to highlight manually via myroi. This must be supplied
#' as a list of vectors wherein each vector is of the format c(seqname, pos, desc), or c(seqname, pos),
#' or c(pos, desc), or c(pos). The seqname (name of a specific sequence in the MSA) and the desc
#' (description of the feature) are optional. pos (i.e., the positions that denote the ROI) can be a
#' single integer (e.g., 100), or a range of integers (1:10). Each such ROI vector can also hold one
#' or more (and combinations) of single integer pos values and integer range pos values (e.g., something
#' like c("Seq1", 1:10, 100, "Helix")). Thus, for example, if a SNP at position 100 in "seq2" and a the
#' region of nucleotides/amino-acids 20:30 representing a domain in all sequences were to be ROIs, they
#' would be supplied like so list(c("seq2", 100, "SNP"), c(20:30, "domain")). Every ROI is highlighted
#' with a separate color (from a color blind friendly palette; the colors are chosen automatically).
#'
#' @param hnon (double, optional) the height of the feature "bar(s)" for all mismatches and gaps in the MSA. (Default: 0.4.)
#'
#' @param hmat (double, optional) the height of the feature "bar(s)" for all matches in the MSA. (Default: 0.4.)
#'
#' @param hroi (double, optional) the height of the feature "bar(s)" for all ROIs in the MSA. (Default: 0.4.)
#'
#' @param wnon (double, optional) the width of the feature "bar(s)" for all mismatches and gaps in the MSA. (Default: 2.0.)
#'
#' @param wmat (double, optional) the width of the feature "bar(s)" for all matches in the MSA. (Default: 1.0.)
#'
#' @param wroi (double, optional) the width of the feature "bar(s)" for all ROIs in the MSA. (Default: 4.0.)
#'
#' @param anon (double, optional) the transparency of the feature "bar(s)" for all mismatches and gaps in the MSA. (Default: 1.0.)
#'
#' @param amat (double, optional) the transparency of the feature "bar(s)" for all matches in the MSA. (Default: 1.0.)
#'
#' @param aroi (double, optional) the transparency of the feature "bar(s)" for all ROIs in the MSA. (Default: 1.0.)
#'
#' @param basecolors (vector of 3 character strings, optional) although the colors for the ROIs are chosen automatically
#' the colors for the matches, mismatches, and gaps can optionally be supplied by the user via this argument. Defaults are
#' c("gray", "black", "white") when no colors are supplied and cbfcols (see below) is set to FALSE.
#'
#' @param cbfcols (boolean, optional) allows for the user to choose whether the coloring scheme used should also be color
#' blind friendly for the matches, mismatches, and gaps (like the colors for the ROIs, if any) in the event the user is not
#' supplying their own colors. (Default: FALSE.)
#'
#' @details
#' msavisr plots the matches as a separate geom_tile() layer, the gaps + mismatches as a geom_tile() layer, and the ROIs as a
#' separate geom_tile() layer (if ROIs are supplied).
#'
#' The user will have to tweak the values for the widths and heights of these layers (via the arguments outlined above) to achieve
#' the desired visualization "effects". In general, it is advisable to set the widths of the mismatches + gaps (and/or ROIs, if any;
#' so wnon and wroi respectively) larger than that of the matches (wmat). The heights can be increased if necessary. Altering the
#' transparency levels does not seem to be very useful. Note: altering the transparency levels does not update the transparency of
#' the colors shown in the legend!!
#'
#' ROIs are especially useful for visualizing features such as single nucleotide polymorphisms (SNPs) in nucleotide MSAs and other
#' such isolated features that might normally become "buried" within the bulk of the sequence. This can be easily achieved by
#' indicating the SNPs position as an ROI and jacking up its width (wroi) and/or height (hroi) values.
#'
#' @return A ggplot2 object is returned to the parent environment for plotting and/or further downstream processing/manipulation.
#'
#' @note
#' The only issue with specifying ROIs in the manner implemented here is that, if for instance, a ROI in "seq2" at position 100
#' and in "seq4" at the same position need to be highlighted, it cannot be supplied like so c("seq2", "seq4", 100), and instead must
#' be supplied as two separate vectors c("seq2", 100), c("seq4", 100). The lowdown: specifying common ROIs in a SUBSET of the MSA can
#' be a tedious process. Unfortunately, as of now, no internal workarounds have been implemented.
#'
#' @examples
#' \dontrun{
#' #Input data
#' testmsa <- system.file("extdata", "testaln_mrna.fasta", package = "seqvisr", mustWork = TRUE)
#'
#' #Basic visualization
#' msavisr(mymsa = testmsa, myref = "Ref0") #No ROIs
#'
#' #Defining an example ROI
#' testroi <- list(c("Ref0", 100:110, "Ref0 Domain1"), c(14, "Pseudouridine"), c(20:30, "Domain2"), c("Seq2", 55, "SNP"))
#'
#' #MSA with ROIs
#' msavisr(mymsa = testmsa, myref = "Ref0", myroi = testroi)
#' }
#'
#' @importFrom utils read.table
#'
#' @importFrom viridis viridis
#'
#' @importFrom tidyr separate_rows
#'
#' @importFrom magrittr %>% %<>%
#'
#' @importFrom stringr str_detect str_replace str_squish
#'
#' @import dplyr ggplot2
#'
#' @export

msavisr <- function(mymsa = NULL, myref = NULL, mypath = NULL, refontop = TRUE, myroi = NULL,
                       hnon = NULL, hmat = NULL, hroi = NULL,
                       wnon = NULL, wmat = NULL, wroi = NULL,
                       anon = NULL, amat = NULL, aroi = NULL,
                       basecolors = NULL, cbfcols = FALSE){

  ##############################################################################################################

  #General description:

  #msavis is an R function to visualize multiple sequence alignments within R.
  #The MSA should be supplied as a fasta file.
  #It is absolutely required that all sequences in the MSA are of the same length!!

  #The function has two mandatory arguments: the name of the MSA file (mymsa, character string),
  #and the name of one of the sequences to be used as the reference against which gaps and mismatches
  #will be calculated (myref, character string).
  #A path (mypath, character string) can also be supplied if necessary.

  #The function accepts both nucleotide and amino acid MSAs (it is actually alphabet agnostic).
  #It produces a ggplot2-based plot (using geom_tile()) of the MSA with matches, mismatches, and gaps highlighted using
  #different colors.
  #By default, the reference sequene is located at the "top" of the MSA visualization.
  #The default coloring scheme is a color blind-friendly combination, but this can be
  #overridden by supplying a three string vector containing the colors for matches, mismatches, and gaps
  #respectively.

  #The function also offers the ability to adjust the widths and heights of the geom_tiles() representing
  #the matches and non-matches (i.e., gaps + mismatches). As the matches and non-matches are plotted as
  #separate layers, adjusting the widths and heights allows for some level of control over highlighting the features
  #in the MSA.

  #Additionally, the use can elect to supply the function with a list of vectors (list(c(), c(), ...))
  #indicating specific regions of interest (ROIs) in the MSA. These ROIs can be either a position in the MSA,
  #or a position in a specific sequence. ROIs are passed as a separate layer with its own height and width parameters
  #permitting this layer to be used for purposes such as highlighting SNPs that might otherwise not be visible in the MSA.


  #Description of the function arguments

  #mymsa (character string, mandatory) the name of the fasta-formatted file containing the multiple sequence alignment (MSA).

  #myref (character string, mandatory) the fasta header (not including the ">") of one of the sequences in the MSA which
  #act as the reference sequence.

  #mypath (character string, optional) the path to the directory containing the MSA.

  #refontop (boolean, optional) should the reference sequence be placed at the top of the MSA plot (default)
  #or at the bottom? (Set to FALSE for bottom.)

  #myroi (list of vectors, optional) the user can manually specify regions of interest (ROI; i.e., positions in the MSA)
  #that they wish to highlight manually via myroi. This must be supplied as a list of vectors wherein each vector
  #is of the format c(seqname, pos, desc), or c(seqname, pos), or c(pos, desc), or c(pos).
  #The seqname (name of a specific sequence in the MSA) and the desc (description of the feature)
  #are optional. pos (i.e., the positions that denote the ROI) can be a single integer (e.g., 100),
  #or a range of integers (1:10). Each such ROI vector can also hold one or more (and combinations) of
  #single integer pos values and integer range pos values (e.g., something like c("Seq1", 1:10, 100, "Helix")).
  #Thus, for example, if a SNP at position 100 in "seq2" and a the region of nucleotides/amino-acids 20:30 representing a domain
  #in all sequences were to be ROIs, they would be supplied like so list(c("seq2", 100, "SNP"), c(20:30, "domain")).
  #Every ROI is highlighted with a separate color (from a color blind friendly palette).

  #hnon (double, optional) the height of the feature "bar(s)" for all mismatches and gaps in the MSA.
  #hmat (double, optional) the height of the feature "bar(s)" for all matches in the MSA.
  #hroi (double, optional) the height of the feature "bar(s)" for all ROIs in the MSA.

  #wnon (double, optional) the width of the feature "bar(s)" for all mismatches and gaps in the MSA.
  #wmat (double, optional) the width of the feature "bar(s)" for all matches in the MSA.
  #wroi (double, optional) the width of the feature "bar(s)" for all ROIs in the MSA.

  #anon (double, optional) the transparency of the feature "bar(s)" for all mismatches and gaps in the MSA.
  #amat (double, optional) the transparency of the feature "bar(s)" for all matches in the MSA.
  #aroi (double, optional) the transparency of the feature "bar(s)" for all ROIs in the MSA.

  #basecolors (vector of 3 character strings, optional) although the colors for the ROIs are chosen automatically
  #the colors for the matches, mismatches, and gaps can optionally be supplied by the user via this argument.
  #Defaults are c("gray", "black", "white") when no colors are supplied and cbfcols (see below) is set to FALSE.

  #cbfcols (boolean, optional) allows for the user to choose whether the coloring scheme used should also be color
  #blind friendly for the matches, mismatches, and gaps (like the colors for the ROIs, if any) in the event the user
  #is not supplying their own colors. (Default: FALSE.)



  ##############################################################################################################


  #require(RColorBrewer)
  #require(viridis)
  #require(magrittr)
  #require(stringr)
  #require(dplyr)
  #require(tidyr)
  #require(ggplot2)

  #Data for testing
  #mymsa <- "/home/owner/Nextcloud/laptop_rplace/google.fasta"
  #  myref <- "KY485192.1"
  #mymsa <- "/home/owner/Nextcloud/laptop_rplace/testaln_mrna.fasta"
  #  myref <- "Ref0"


  #Basic checks
  if(is.null(mymsa)) { stop("Please supply the name of a MSA file! (Must be fasta formatted!)") }
  if(is.null(myref)) { stop("Please indicate which sequence you would like to use as the reference!") }

  #Setting up the path to the object
  if(is.null(mypath)){
    mypath <- mymsa
  } else{
    mypath <- paste0(mypath, "/", mymsa)
  }


  #Reading in the fasta data as a newline-delimited data.frame
  df <- utils::read.table(mypath, sep = "\n")

  #Container data.frame where all the computations will occur
  fasdf <- data.frame(curhead = c(), curseq = c(), stringsAsFactors = FALSE)

  #Identifying those rows that contain the fasta headers
  headerpos <- which(stringr::str_detect(df$V1, "^>"))

  #Using headerpos and a for loop to create a two column data.frame
  #wherein each column in each row represents the fasta header and fasta sequence respectively
  for(i in 1:length(headerpos)){

    if(i == length(headerpos)){
      curpos <- headerpos[i]
      nextpos <- nrow(df)
      curhead <- df$V1[curpos]
      curseq <- paste0(df$V1[(curpos+1):nextpos], collapse = "")
    } else{
      curpos <- headerpos[i]
      nextpos <- headerpos[i+1]
      curhead <- df$V1[curpos]
      curseq <- paste0(df$V1[(curpos+1):(nextpos-1)], collapse = "")
    }

    fasdf <- dplyr::bind_rows(fasdf, data.frame(curhead, curseq, stringsAsFactors = FALSE))

    if(i == length(headerpos)){
      rm(curpos, nextpos, curhead, curseq, headerpos)
    }

  }
  rm(i)


  #Need a column that will contain the length of each sequence
  #Since this is a MSA--all sequences are of the same length, so length of just one of the is sufficient
  curlen <- nchar(fasdf$curseq[1])

  #Moving each character of each sequence into its own row
  fasdf %<>% dplyr::mutate(curseq = stringr::str_squish(curseq))
  fasdf %<>% tidyr::separate_rows(curseq, sep = "")
  #separate_rows() adds an additional column for whatever reason (probably because "" is used as the separator)
  #Removing this; however, the sequence length calculated above is the "correct" one, so removing this, does not affect
  #that (this blank is in addition)
  fasdf %<>% dplyr::filter(curseq != "")

  #Using the curlen variable above, creating a position column representing each character's position
  #in the sequence
  fasdf %<>% dplyr::mutate(curpos = rep(1:curlen, len = nrow(fasdf)))


  #Cleaning the > in the fasta header column
  fasdf %<>% dplyr::mutate(curhead = stringr::str_replace(curhead, "^>", ""))




  #The user can select whether the reference sequence goes on the top or bottom of the plot
  olvls <- fasdf %>% dplyr::filter(curhead != myref) %>% dplyr::distinct(curhead)
  #olvls <- as.character(olvls)
  #refontop <- FALSE
  if(isTRUE(refontop)){
    olvls <- rbind(olvls, myref)
  } else{
    olvls <- rbind(myref, olvls)
  }
  olvls <- olvls$curhead
  #Reordering
  fasdf$curhead <- factor(fasdf$curhead, levels = olvls)
  rm(olvls)


  #To visualize SNPs and gaps, some one sequence needs to be set as the reference
  #This is passed to the function
  #  myref <- "Ref0"
  #  myref <- ">LT734359.1/2709-2956 Human ORFeome Gateway entry vector pENTR223-CLOCK, complete sequence"
  #First taking this reference sequence out separately and replicate()'ing it vector
  #of length (num_sequences) * (num_chars_in_seq); the idea is to cbind() this to the entire data.frame
  #and then mark each character as a perfect match (*), gap (-), or SNP (+) by simply comparing
  #the sequence column (curseq) against this reference column (refcol)
  #That's what's being done below
  myrefseq <- fasdf %>% dplyr::filter (curhead == myref) %>% dplyr::select(curseq)
  refcol <- unlist(replicate(length(unique(fasdf$curhead)), myrefseq$curseq, simplify = FALSE))

  fasdf <- cbind(fasdf, refcol)
  rm(myrefseq, refcol)

  for(i in 1:nrow(fasdf)){

    if(fasdf$curseq[i] == fasdf$refcol[i]){
      fasdf$outcol[i] <- "Match" #"*" Identical
    } else if(fasdf$curseq[i] != fasdf$refcol[i] & fasdf$curseq[i] == "-"){
      fasdf$outcol[i] <- "Gap" #"-" Gap
    } else{
      fasdf$outcol[i] <- "Mismatch" #"+" SNP
    }

  }







  #The user can also manually specify regions of interest (ROI; i.e., positions in the MSA)
  #that they wish to highlight manually via myroi as a list of vectors wherein each vector
  #is of the format c(seqname, pos, desc), or c(seqname, pos), or c(pos, desc), or c(pos).
  #The seqname (name of a specific sequence in the MSA) and the desc (description of the feature)
  #are optional. pos (i.e., the positions that denote the ROI) can be a single integer (e.g., 100),
  #or a range of integers (1:10). Each such ROI vector can also hold one or more (and combinations) of
  #single integer pos values and integer range pos values (e.g., something like c("Seq1", 1:10, 100, "Helix")).

  #If such are specified, they also need to be indicated in fasdf, and additonal legend
  #items and such must be set up, which is what is being done in the if block below along with the
  #concomitant plotting.

  #Values for testing (must have the appropriate MSA loaded)
  #
  #myroi <- list(c("VM_4359", 100, "Hydrophobic"), c("VM_4359", 200:210), c("Novel_4124", 220), 180)
  #myroi <- list(c("Seq2", 100, "Hydrophobic"), c("Seq3", 200:210), c("Seq4", 220), 180, c(50:60, "Helix"))


  #Checking if regions of interest have been supplied and plot accordingly

  #If they have been supplied
  if(!is.null(myroi)){


    #Create an empty column that will indicate whether the particular position
    #is a region of interest or not
    fasdf$roicol <- NA

    #Iterate through the ROIs and indicate them as ROIs in fasdf's outcol
    #Get the items in the current position in the list
    for(i in 1:length(myroi)){


      curset <- myroi[[i]]
      curset

      #Interpreting the structure of the ROIs depends upon whether or not the first and
      #last items in the vector are alphanumeric
      lind <- length(curset)
      hasroiname <- ifelse(stringr::str_detect(curset[1], "[A-Za-z]+"), TRUE, FALSE)
      hasroidesc <- ifelse(stringr::str_detect(curset[lind], "[A-Za-z]+"), TRUE, FALSE)

      #The user can supply specific positions in a specific sequence, or specific positions
      #universal to all sequences
      #To check for this, need to first see if the very first item in curset is a character
      #string or not, and proceed from there

      #If a specific sequence, positions, and description have been provided
      if(hasroiname & hasroidesc){
        roiseq <- curset[1]
        roipos <- curset[-1]
        roipos <- roipos[-length(roipos)]
        roidesc <- curset[length(curset)]

        roiidx <- paste0(roipos[1], "-", roipos[length(roipos)])

        fasdf %<>% dplyr::mutate(roicol = ifelse(curhead == roiseq & curpos %in% roipos,
                                          paste0("ROI ", roiseq, " ", roidesc, " ", roiidx), roicol))
      } else if(hasroiname & !hasroidesc){

        #If a specific sequence and positions have been provided but no description

        roiseq <- curset[1]
        roipos <- curset[-1]

        roiidx <- paste0(roipos[1], "-", roipos[length(roipos)])

        fasdf %<>% dplyr::mutate(roicol = ifelse(curhead == roiseq & curpos %in% roipos,
                                          paste0("ROI ", roiseq, " ", outcol, " ", roiidx), roicol))
      } else if(!hasroiname & hasroidesc){

        #If positions and a descriptions have been provided, but no specific sequence

        roipos <- curset[-length(curset)]
        roidesc <- curset[length(curset)]

        roiidx <- paste0(roipos[1], "-", roipos[length(roipos)])

        fasdf %<>% dplyr::mutate(roicol = ifelse(curpos %in% roipos,
                                          paste0("ROI ", roidesc, " ", roiidx), roicol))
      } else if(!hasroiname & !hasroidesc){
        #If only positions have been provided

        roiidx <- paste0(roipos[1], "-", roipos[length(roipos)])

        fasdf %<>% dplyr::mutate(roicol = ifelse(curpos %in% roipos,
                                          paste0("ROI ", outcol, " ", roiidx), roicol))
      }

    }
    rm(roidesc, roipos, roiseq, lind, i, hasroiname, hasroidesc, curset)

    #For all remaining roicol positions that are just NA, fill them with their equivalent
    #outcol values
    #
    fasdf %<>% dplyr::mutate(roicol = ifelse(is.na(roicol), outcol, roicol))

    #unique(fasdf$roicol)


    #For ordering the legend items
    baselvls <- c("Match", "Mismatch", "Gap")
    roilvls <- unique(fasdf$roicol)[!(unique(fasdf$roicol) %in% baselvls)]
    roilvls <- roilvls[order(roilvls)]
    baselvls <- c(baselvls, roilvls)

    fasdf$roicol <- factor(fasdf$roicol, levels = baselvls)

    #fasdf$outcol <- factor(fasdf$outcol, levels = c("Match", "Mismatch", "Gap"))

    #Setting up the colors for plotting
    #The color palette consists of colors for the basic match, mismatch, gap trio
    #plus whatever are the n categories of ROIs
    #The colors for the ROIs are selected from viridis by passing the length of roilvls
    #from above to it
    #The basecolors may be optionally supplied by the user; if the user sets cbf to TRUE,
    #then the basecolors are also chosen from viridis()
    #basecolors <- NULL
    #basecolors <- c("red", "green", "yellow")
    #cbfcols <- FALSE
    cbPalette <- c()

    if(isTRUE(cbfcols) & is.null(basecolors)){
      cbPalette <- viridis::viridis(length(baselvls))
    } else if(is.null(basecolors) & length(cbPalette) != length(baselvls)){
      basecolors <- c("gray", "black", "white")
      cbPalette <- c(basecolors, viridis::viridis(length(roilvls)))
    } else if(!is.null(basecolors)){
      cbPalette <- c(basecolors, viridis::viridis(length(roilvls)))
    }




    #rm(baselvls, roilvls)


    #For plotting, the matches, nonmatches (gaps + mismatches), and ROIs need to be
    #supplied as separate layers
    matches <- fasdf %>% dplyr::filter(roicol == "Match")
    nonmatches <- fasdf %>% dplyr::filter(roicol != "Match" & !stringr::str_detect(roicol, "^ROI"))
    rois <- fasdf %>% dplyr::filter(stringr::str_detect(roicol, "^ROI"))


    #To make specific features (e.g. a ROI or gaps) stand out, each layer has its own width and height
    #parameters
    #Defaults supplied below can be overridden by user inputs
    #    hmat <- 0.4
    #    hnon <- 0.8
    #    hroi <- 0.8
    #    wmat <- 1.0
    #    wnon <- 4.0
    #    wroi <- 4.0
    #    amat <- 1.0
    #    anon <- 1.0
    #    aroi <- 1.0
    #Defaults
    if(is.null(hmat)) {hmat <- 0.4}
    if(is.null(hnon)) {hnon <- 0.4}
    if(is.null(hroi)) {hroi <- 0.4}

    if(is.null(wmat)) {wmat <- 1.0}
    if(is.null(wnon)) {wnon <- 2.0}
    if(is.null(wroi)) {wroi <- 4.0}

    if(is.null(amat)) {amat <- 1.0}
    if(is.null(anon)) {anon <- 1.0}
    if(is.null(aroi)) {aroi <- 1.0}




    #Plotting
    myplt <- ggplot2::ggplot() +
      ggplot2::geom_tile(data = matches,
                         ggplot2::aes(x = curpos, y = curhead, fill = roicol), width = wmat, height = hmat, alpha = amat) +
      ggplot2::geom_tile(data = nonmatches,
                         ggplot2:: aes(x = curpos, y = curhead, fill = roicol), width = wnon, height = hnon, alpha = anon) +
      ggplot2::geom_tile(data = rois,
                         ggplot2::aes(x = curpos, y = curhead, fill = roicol), width = wroi, height = hroi, alpha = aroi) +
      ggplot2::xlab("Position in sequence") + ggplot2::ylab("") +
      ggplot2::scale_fill_manual(name = "", values = cbPalette, drop = FALSE) +
      ggplot2::theme_classic() +
      ggplot2::theme(legend.key = ggplot2::element_rect(colour = "black", size = 1.0),
            axis.line.y = ggplot2::element_blank(), axis.line.x = ggplot2::element_blank(),
            axis.text.y = ggplot2::element_text(face = "bold"), axis.ticks.y = ggplot2::element_blank())
    myplt


  }


  #If not regions of interest have been supplied
  if(is.null(myroi)){


    #For ordering the legend items
    baselvls <- c("Match", "Mismatch", "Gap")
    fasdf$outcol <- factor(fasdf$outcol, levels = baselvls)

    #Setting up the colors
    #basecolors <- NULL
    #basecolors <- c("red", "green", "yellow")
    #cbfcols <- FALSE
    cbPalette <- c()

    if(isTRUE(cbfcols) & is.null(basecolors)){
      cbPalette <- viridis::viridis(length(baselvls))
    } else if(is.null(basecolors) & length(cbPalette) != length(baselvls)){
      basecolors <- c("gray", "black", "white")
      cbPalette <- basecolors
    } else if(!is.null(basecolors)){
      cbPalette <- basecolors
    }

    #Not ROIs so only two layers: matches and non-matches as layers for ggplot()
    matches <- fasdf %>% dplyr::filter(outcol == "Match")
    nonmatches <- fasdf %>% dplyr::filter(outcol != "Match")


    #Width, height, alpha parameters for each geom_tile layer
    #Defaults supplied below can be overridden by user inputs
    #    hmat <- 0.4
    #    hnon <- 0.8
    #hroi <- 0.8
    #    wmat <- 1.0
    #    wnon <- 4.0
    #wroi <- 4.0
    #    amat <- 1.0
    #    anon <- 1.0
    #aroi <- 1.0
    #Defaults
    if(is.null(hmat)) {hmat <- 0.4}
    if(is.null(hnon)) {hnon <- 0.4}
    if(is.null(hroi)) {hroi <- 0.4}

    if(is.null(wmat)) {wmat <- 1.0}
    if(is.null(wnon)) {wnon <- 2.0}
    if(is.null(wroi)) {wroi <- 4.0}

    if(is.null(amat)) {amat <- 1.0}
    if(is.null(anon)) {anon <- 1.0}
    if(is.null(aroi)) {aroi <- 1.0}




    #Plotting
    myplt <- ggplot2::ggplot() +
      ggplot2::geom_tile(data = matches,
                         ggplot2::aes(x = curpos, y = curhead, fill = outcol), width = wmat, height = hmat, alpha = amat) +
      ggplot2::geom_tile(data = nonmatches,
                         ggplot2::aes(x = curpos, y = curhead, fill = outcol), width = wnon, height = hnon, alpha = anon) +
      ggplot2::xlab("Position in sequence") + ggplot2::ylab("") +
      ggplot2::scale_fill_manual(name = "", values = cbPalette, drop = FALSE) +
      ggplot2::theme_classic() +
      ggplot2::theme(legend.key = ggplot2::element_rect(colour = "black", size = 1.0),
            axis.line.y = ggplot2::element_blank(), axis.line.x = ggplot2::element_blank(),
            axis.text.y = ggplot2::element_text(face = "bold"), axis.ticks.y = ggplot2::element_blank())
    myplt


  }


  return(myplt)

}