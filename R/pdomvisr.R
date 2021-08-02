#' Convert a single sequence's quasi-standard TSV "row" input data.frame format to a position-per-row ggplot()-friendly format.
#'
#' @usage
#' tsvtogginp(inpdf = NULL)
#'
#' @param inpdf (character, mandatory) the name of the input data.frame().
#'
#' @details
#' Function to take a data.frame() indicating the sequence name, length, sequence start offset, and domains' start and end
#' positions (one row for each domain), FOR ONE SEQUENCE, and convert this into (another) data.frame that can be parsed by
#' ggplot2. Namely, the output data.frame has one row for each residue in the sequence, and indicates the positions (important
#' for ggplot2), and the associated annotation (e.g., "ABC domain"). The function also has additional columns indicating data
#' for labeling that should be visualized on the plot or used to create it--this is basically done by taking the domain's name
#' and assigning it to a label column but only for a single residue at the middle of the domain. (If the labels of all residues
#' from each domain are displayed, the labeling will be illegible.)
#'
#' tsvtogginp() is invoked by pdomvisr() through the wrapper script tsvtgginp_multi().
#'
#' @return a data.frame() with one row per residue indicating the sequence name, position (indicated by that row), the description
#' associated with that position (e.g., "ABC domain", or NA for no description), a column indicating which feature the position
#' belongs to (e.g., say there are two domains named "ABC", this numerical column is a way to distinguish between them), a column
#' indicating the layer (offset, base, feature; see ?pdomvisr()), label for the particular feature (assigned to one position within
#' the feature; see ?pdomvisr()), and a column indicating the original positions of the rows (e.g., when a sequence offset has been
#' applied; see ?pdomvisr()).
#'
#' @seealso tsvtogginp_multi(), pdomvisr()
#'
#' @examples
#' \dontrun{
#' #Input data
#' inpath <- system.file("extdata", "pdomvisr_testdata.tsv", package = "seqvisr", mustWork = TRUE)
#' inpdat <- read.table(inpath, header = TRUE)
#' inpdat <- inpdat[2:3, ]
#'
#' #Example run
#' tsvtogginp(inpdat)
#' }
#'
#' @import dplyr
#'
#' @importFrom magrittr %>% %<>%
#'
#' @importFrom tibble tibble
#'
#' @export

tsvtogginp <- function(inpdf = NULL){

  #To avoid the "no visible binding for global variable" and "Undefined global functions or variables:"
  #warnings during devtools::check()
  start_loc <- NULL
  posdesc <- NULL
  fset <- NULL
  pos <- NULL
  mid <- NULL
  prot_acc <- NULL
  orig_pos <- NULL


  if(is.null(inpdf)){ stop("No input data.frame provided!!") }


  #Check if the input data.frame has 5 columns, and rename them as follows.
  #c("prot_acc", "seq_len", "signature_desc", "start_loc", "stop_loc")
  if(ncol(inpdf) != 6){
    stopmsg <- base::strwrap("Please ensure all required input columns are provided!!
       These include must include the following in THIS order:
       Sequence name, sequence length, start offset, domain descriptor, domain start position, domain end position.
       Please supply one domain per row, and use multiple rows per sequence in case of multi-domain sequences!!")
    stop(stopmsg)
  }

  #Renaming columns properly, so that the user can supply pretty much whatever column
  #names they like.
  names(inpdf) <- c("prot_acc", "seq_len", "offset", "signature_desc", "start_loc", "stop_loc")



  #Pulling the annotations out separately.
  annots <- inpdf[ , -c(1:3)]
  #There may be more than one feature with the same name (e.g., PAS).
  #So using the fset column to distinguish them. This will be carried
  #over to the output data.frame.
  #Arrange annots by start_loc so that the fset assignment is sequential
  #in the following data.frames (i.e., domain at 200 gets fset 1, domain at
  #400 gets fset2).
  annots %<>% dplyr::arrange(start_loc)
  annots$fset <- base::seq(base::nrow(annots))

  #Need the sequence length as a numeric, and the sequence name to
  #start creating the output data.frame.
  slen <- as.numeric(unique(inpdf[ , 2]))
  sname <- as.character(unique(inpdf[ , 1]))

  #These steps below basically create the output data.frame
  #by taking a single-column data.frame consisting of the sequence
  #name, and replicating it slen times. Then a position column
  #(basically the row number) is added, yielding a data.frame with
  #one row for each residue.
  tmp <- base::data.frame(prot_acc = sname, stringsAsFactors = FALSE)
  outdf <- base::do.call("rbind", base::replicate(slen, tmp, simplify = FALSE))
  outdf$pos <- base::seq(base::nrow(outdf))


  #Assigning the annotations. This will first be initialized as an
  #empty column.
  outdf$posdesc <- NA
  #There may be more than one feature with the same name (e.g., PAS).
  #So I will distinguish these with the fset column as used in the
  #annot data.frame above. Also initialized to NA at first.
  outdf$fset <- NA

  #I also want to have layered control over the features.
  #So creating a layer column, and initializing it as base.
  #All rows corresponding to a domain or feature will have this
  #overwritten with the string "feature".
  outdf$layer <- "base"

  #Assigning the annotation using a for loop.
  #The idea is to loop through every row (i.e., residue) in outdf,
  #and check whether it lies within the range of residues identifed
  #as a domain or feature in the annots data.frame. If it does, then
  #it will be assigned the annotation (e.g., "PAS") corresponding to
  #the appropriate row from the annots data.frame.
  for(i in 1:base::nrow(outdf)){
    for(j in 1:base::nrow(annots)){
      if(is.na(outdf$posdesc[i])){

        #The as.numeric() transformation of the start and stop locations is important as
        #in some edge cases, the values in those cells are being treated as characters.
        if(outdf$pos[i] >= as.numeric(annots$start_loc[j]) & outdf$pos[i] <= as.numeric(annots$stop_loc[j])){
          outdf$posdesc[i] <- annots$signature_desc[j]
          outdf$layer[i] <- "feature"
          outdf$fset[i] <- annots$fset[j]


        }
      }
    }
  }


  #For each feature on the protein, I need to assign the descriptor
  #string to just one column belonging to that feature. That one
  #column's position is calculated as the start of the domain/feature
  #+ half the length of the feature (i.e., the mid-point of the feature).
  #This will be done through the label column.
  outdf$label <- NA

  outdf %<>%
    dplyr::group_by(layer, posdesc, fset) %>%
    dplyr::arrange(pos, .by_group = TRUE) %>%
    dplyr::mutate(mid = base::min(pos) + base::round( (base::max(pos)-base::min(pos))/2 ,0)) %>%
    dplyr::mutate(label = base::ifelse(pos == mid, posdesc, NA)) %>%
    dplyr::ungroup() %>%
    dplyr::select(-mid)



  #The user can also supply a sequence start offset value.
  #This is basically to position sequences relative to other sequences,
  #or to highlight the fact that a sequence is--for example--a partial
  #fragement.
  #For this, what's being done is that as many rows as the offset
  #indicates are being rbind'd to the head of outdf. Subsequently,
  #the position counters are all re-incremented accordingly.
  poffset <- as.numeric(base::unique(inpdf[ , 3]))
  if(poffset != 0){

    #Need a single row of outdf to replicate.
    offdf <- outdf[1, ]
    #Setting all the column values except the first (sequence identifier)
    #to NA.
    offdf %<>% dplyr::mutate(across(.cols = -prot_acc, ~ NA))
    #Setting the layer column to "offset"
    offdf$layer <- "offset"
    #Setting the label to ""
    offdf$label <- NA

    #Replicating this single row the requisite number of times.
    offdf <- base::do.call("rbind", base::replicate(poffset, offdf, simplify = FALSE))

    #Updating the position counter for offdf.
    #This will count up from a negative number to 0.
    offdf$pos <- 1 - base::seq(base::nrow(offdf))

    #Binding offdf to the head of outdf.
    outdf <- dplyr::bind_rows(offdf, outdf)

    #Updating the position counter.
    outdf$orig_pos <- outdf$pos
    outdf %<>% dplyr::arrange(orig_pos)
    outdf$pos <- base::seq(base::nrow(outdf))


  } else{
    #No offset supplied, just creating that orig_pos column.
    outdf$orig_pos <- outdf$pos
  }



  #Returning the output data.frame as a tibble.
  #This is because the parent function that will be calling this
  #uses dplyr::group_map, and that seems to work best with a tibble
  #as an input.
  return(tibble::tibble(outdf))
}





#' Wrapper script around tsvtogginp() to perform conversion of a data.frame() of multiple sequences and their domain annotations.
#'
#' @usage
#' tsvtogginp_multi(inpdf = NULL)
#'
#' @param inpdf (character, mandatory) the name of the input data.frame().
#'
#' @details
#' Basically applies tsvtogginp() on multiple sequences worth of rows using dplyr::group_map(). tsvtogginp_multi() is invoked by
#' pdomvisr().
#'
#' @return a data.frame() with one row per residue indicating the sequence name, position (indicated by that row), the description
#' associated with that position (e.g., "ABC domain", or NA for no description), a column indicating which feature the position
#' belongs to (e.g., say there are two domains named "ABC", this numerical column is a way to distinguish between them), a column
#' indicating the layer (offset, base, feature; see ?pdomvisr()), label for the particular feature (assigned to one position within
#' the feature; see ?pdomvisr()), and a column indicating the original positions of the rows (e.g., when a sequence offset has been
#' applied; see ?pdomvisr()).
#'
#' @seealso tsvtogginp(), pdomvisr()
#'
#' @examples
#' \dontrun{
#' #Input data
#' inpath <- system.file("extdata", "pdomvisr_testdata.tsv", package = "seqvisr", mustWork = TRUE)
#' inpdat <- read.table(inpath, header = TRUE)
#'
#' #Example run
#' tsvtogginp_multi(inpdat)
#' }
#'
#' @import dplyr
#'
#' @importFrom magrittr %>% %<>%
#'
#' @importFrom tibble tibble
#'
#' @export

tsvtogginp_multi <- function(inpdf = NULL){

  #To avoid the "no visible binding for global variable" and "Undefined global functions or variables:"
  #warnings during devtools::check()
  prot_acc <- NULL

  if(is.null(inpdf)){ stop("No input data.frame provided!!") }


  #Check if the input data.frame has 5 columns, and rename them as follows.
  #c("prot_acc", "seq_len", "signature_desc", "start_loc", "stop_loc")
  if(ncol(inpdf) != 6){
    stopmsg <- base::strwrap("Please ensure all required input columns are provided!!
       These include must include the following in THIS order:
       Sequence name, sequence length, start offset, domain descriptor, domain start position, domain end position.
       Please supply one domain per row, and use multiple rows per sequence in case of multi-domain sequences!!")
    stop(stopmsg)
  }

  #Renaming columns properly, so that the user can supply pretty much whatever column
  #names they like.
  base::names(inpdf) <- c("prot_acc", "seq_len", "offset", "signature_desc", "start_loc", "stop_loc")

  #Applying tsvtogginp() to every unique sequence's set of rows in the input data.frame.
  outdf <- inpdf %>%
    dplyr::group_by(prot_acc) %>%
    dplyr::group_map(~ tsvtogginp(.x), .keep = TRUE) %>%
    dplyr::bind_rows() %>%
    dplyr::ungroup()

  return(tibble::tibble(outdf))

}




#' Protein domain structure visualization in R.
#'
#' @description
#' pdomvisr() is a simple function to plot a diagram of the domain structure of one or more
#' sequences. pdomvisr() uses ggplot2::ggplot() internally. The only mandatory input is a table
#' with the following information (in this particular order): sequence name, sequence length,
#' sequence offset, domain description, domain start coordinate, domain end coordinate.
#'
#'
#' @usage pdomvisr(inpdat = NULL, mypath = NULL,
#' nbreaks = NULL, cbfcols = FALSE,
#' legend = TRUE, show_offsets = TRUE,
#' label_size = "auto", featcols = NULL,
#' hbase = 0.2, hfeat = 2.4*hbase, hoff = 0.8*hbase,
#' alpbase = 1.0, alpfeat = 1.0, alpoff = 0.05,
#' fillbase = "black", filloff = "white",
#' colorbase = "black", coloroff = "gray",
#' nudge_x = 0.0, nudge_y = 0.5)
#'
#'
#' @param inpdat (character string or name of an object, mandatory) the character string may be the
#' name of a file or the full path to it (in which case mypath should be set to NULL). The file must
#' be a table containing the information necessary to plot the domain structure diagram. Alternatively,
#' inpdat can also be supplied the name of an object in R's environment (e.g., a data.frame containing
#' the requisite data). This is useful when the input data needs to be pre-processed in R first. The
#' path/filename option is more suitable when pdomvisr is being called only for plotting. In this case,
#' data.table's fread() is used to read the data into the function first. Please see the 'Details' section
#' for information on how the input data must be formatted.
#'
#' @param mypath (character string, optional) in the event that inpdat is supplied the name of a file, the
#' path to where this file is located can be supplied through mypath.
#'
#' @param nbreaks (numeric, optional) controls the number of X-axis ticks in the plotted domain structure
#' diagram. If the user does not supply a number, this is automatically calculated to produce a tick every
#' 100 residues (based on the length of the longest sequence included in the plot).
#'
#' @param cbfcols (boolean, optional) controls whether colorblind-friendly colors will be used for the
#' domains. (Set to FALSE by default.)
#'
#' @param legend (boolean, optional) controls whether the legend associating the domain colors to the domain
#' descriptions should be plotted along with the main plot. (Set to TRUE by default.)
#'
#' @param show_offsets (boolean, optional) controls whether sequence offsets should be plotted or hidden.
#' A sequence offset is a whole number (supplied as a part of the input table) indicating how far off from
#' the actual first residue of the sequence the first residue indicated in the input data is. This is relevant
#' when plotting partial sequences for instance (e.g., an internal fragment). (Set to TRUE by default.)
#'
#' @param label_size (character or numeric, optional) controls whether the domain descriptions are displayed as
#' labels on the domains. Also controls the size of the text if the labels are displayed. The size can be
#' controlled by supplying a positive integer > 0. Supplying 0 prevents the labels from being displayed. Passing
#' "auto" leaves the size estimation to R. If set to "repel", then the labels are drawn offset from the features
#' and connected to them by straight lines. If set to "repel", the arguments nudge_x and nudge_y (see below) can
#' be adjusted by the user to vary the positioning of the labels. Note: this argument's values do not affect the
#' legend. (Set to "auto" by default.)
#'
#' @param featcols (character vector, optional) controls the colors assigned to the features/domains. As many
#' colors as there are unique features in the data set must be supplied.
#'
#' @param hbase (numeric, optional) controls the height of the tiles corresponding to the non-domain portions
#' of the sequence. (Set to 0.2 by default.)
#'
#' @param hfeat (numeric, optional) controls the height of the tiles corresponding to the DOMAIN portions of the
#' sequence. Under default settings, this scales automatically with hbase. (Set to 2.4 * hbase by default.)
#'
#' @param hoff (numeric, optional) controls the height of the tiles representing the sequence offset. Under
#' default settings, this scales automatically with hbase. (Set to 0.8 * hbase by default.)
#'
#' @param alpbase (numeric, optional) controls the alpha level of the tiles representing the non-domain portions
#' of the sequence. (Set to 1.0 by default.)
#'
#' @param alpfeat (numeric, optional) controls the alpha level of the tiles representing the DOMAIN portions of
#' the sequence(s). (Set to 1.0 by default.)
#'
#' @param alpoff (numeric, optional) controls the alpha level of the tiles representing the sequence offset. (Set
#' to 0.05 by default.)
#'
#' @param fillbase (character, optional) fill color for the tiles representing the non-domain portions of the
#' sequence. Any value accepted by ggplot2's "fill" is accepted here, as this just passes the value on
#' to that particular argument. (Set to "black" by default.)
#'
#' @param filloff (character, optional) fill color for the tiles representing the offset sequence. Any value
#' accepted by ggplot2's "fill" is accepted here, as this just passes the value on to that particular argument.
#' (Set to "white" by default.)
#'
#' @param colorbase (character, optional) line color for the tiles representing the non-domain portions of the
#' sequence. Any value accepted by ggplot2's "color" is accepted here, as this just passes the value on
#' to that particular argument. (Set to "black" by default.)
#'
#' @param coloroff (character, optional) line color for the tiles representing the offset sequence. Any value
#' accepted by ggplot2's "color" is accepted here, as this just passes the value on to that particular argument.
#' (Set to "gray" by default.)
#'
#' @param nudge_x (numeric, optional) if label_size is set to "repel", adjusting this value changes the horizontal
#' starting position of the label (this is the same parameter as ggrepel::geom_text_repel()'s nudge_x; so see that
#' function's help page for more details).
#'
#' @param nudge_y (numeric, optional) if label_size is set to "repel", adjusting this value changes the vertical
#' starting position of the label (this is the same parameter as ggrepel::geom_text_repel()'s nudge_y; so see that
#' function's help page for more details).
#'
#' @details
#' pdomvisr() plots a domain structure diagram given the coordinates of domains in one or more sequences.
#'
#' The only mandatory input is a table with the following information (in this particular order): sequence name,
#' sequence length, sequence offset, domain description, domain start coordinate, domain end coordinate. Most column
#' names are self-explanatory, and should be readily produced by most domain annotation tools (or should be producible
#' by hand). The only column that will have to be created by the user is the sequence offset column. This is a special
#' feature that is strictly optional. The objective of this column is to indicate how far away from the actual start
#' of the sequence the "indicated" start of the sequence is. For example, if a sequence is listed as being 200 residues
#' long, and has an offset value of 10, this implies that the sequence is actually 210 residues long, but only the last
#' 200 residues "exist". The offset column is useful for visually demonstrating that a sequence is partial. If none of
#' the sequences are partial and/or in no need of an offset, the values in this column should be set to 0. It is mandatory
#' that this column exists (with a default value of 0 in all rows) even if it is not in use!!
#'
#' Each row in the input table should correspond to the coordinates for a particular domain in a particular sequence.
#' Therefore, if a sequence contains more than one domain, it will have to be represented by as many rows as there are
#' domains in it. The sequence name, sequence length, and sequence offset columns will (unfortunately) have to be repeated
#' in all such rows, despite being redundant in this manner.
#'
#' The input data can be any tabular file that data.table's fread() can parse. Alternatively, the user can also supply
#' the name of an R object containing the data. The R object in question can be a data.frame, data.table::data.table,
#' or tibble::tibble. This is useful in a situation where the data has had to have been munged in order to prepare it
#' for plotting with pdomvisr(). This input (file name, file name + path, or object name) is the only mandatory argument
#' required by pdomvisr(). The tabular input format was chosen as it is tool and platform agnostic, and most prominent
#' annotation tools (e.g., Hmmer3 and InterProScan) are capable of producing outputs in this format (or produce outputs
#' coercible into this format).
#'
#' pdomvisr() uses ggplot2::ggplot() internally to draw the domain structure diagram. In specific, it uses ggplot2::
#' geom_tile() to render each position in the sequence as its own tile. pdomvisr() uses an internal function
#' (tsvtogginp_multi()) to transform the input data into a "long" style data.frame consisting of one row per position
#' per sequence. Each position is identified as belonging to one of three (ggplot2) layers: offset, base, or feature.
#' The offset layer, as the name suggests, represents all positions that form the "offset sequence". The base layer
#' includes all positions that are neither a part of the offset sequence nor a part of a domain in the sequence. The
#' domain layer includes all residues that belong to a domain. Additional columns exist to indicate labeling (and other
#' logistics) for the domains. As there is no upper or lower bound on the length of a domain as far as pdomvisr() is
#' concerned, even single residues (e.g., an active site) can be annotated by adding a row carrying the same value for the
#' start and end positions for the "domain".
#'
#' pdomvisr() plots the diagram on the current device (the plot pane in RStudio, for example) and also returns the ggplot2
#' object itself to the parent environment from which the function was called. The user therefore has complete control over
#' the output.
#'
#' @return A ggplot2 object is returned to the parent environment for plotting and/or further downstream processing/manipulation.
#'
#' @note
#' In some cases it may be necessary to include a sequence that has no annotated domains. Most domain annotation tools would
#' not include a row for such sequences. But as long as the user adds a row to the input table manually with the sequence
#' description column set to NA, and the start and end positions of the domain set to 0, a sequence with no annotated domains
#' will still show up in the plot.
#'
#' @examples
#' \dontrun{
#' #Input data
#' inpath <- system.file("extdata", "pdomvisr_testdata.tsv", package = "seqvisr", mustWork = TRUE)
#'
#' #Default function call with colorblind-friendly colors.
#' pdomvisr(inpdat = inpath, cbfcols = TRUE)
#' }
#'
#' @importFrom data.table fread
#'
#' @importFrom tibble tibble
#'
#' @importFrom viridis viridis
#'
#' @importFrom magrittr %>% %<>%
#'
#' @importFrom scales pretty_breaks
#'
#' @importFrom ggrepel geom_text_repel
#'
#' @import dplyr ggplot2
#'
#' @export

pdomvisr <- function(inpdat = NULL, mypath = NULL,
                     nbreaks = NULL, cbfcols = FALSE,
                     legend = TRUE, show_offsets = TRUE,
                     label_size = "auto", featcols = NULL,
                     hbase = 0.2, hfeat = 2.4*hbase, hoff = 0.8*hbase,
                     alpbase = 1.0, alpfeat = 1.0, alpoff = 0.05,
                     fillbase = "black", filloff = "white",
                     colorbase = "black", coloroff = "gray",
                     nudge_x = 0.0, nudge_y = 0.5){

  #To avoid the "no visible binding for global variable" and "Undefined global functions or variables:"
  #warnings during devtools::check()
  pos <- NULL
  prot_acc <- NULL
  posdesc <- NULL
  label <- NULL

  #Stop if input is empty.
  if(is.null(inpdat)){
    stop("The domain annotation table is a mandatory input!! Please provide this via the inpdat argument!!")
  }

  #Checking to see if what is supplied to inpdat is a path or an object.
  if(is.data.frame(inpdat) | data.table::is.data.table(inpdat) | tibble::is_tibble(inpdat)){
    #Data is already loaded into R, do nothing.
    inpdf <- inpdat
  } else{
    #Setting up the path to the object
    if(!is.null(mypath)){
      inpdat <- paste0(mypath, "/", inpdat)
    }

    #Reading in the data using data.table::fread.
    inpdf <- data.table::fread(inpdat)
  }


  #Converting data for plotting with tsvtogginp_multi.
  inpdf <- tsvtogginp_multi(inpdf)



  #The tsvtogginp() setup (accessed through tsvtogginp_multi)
  #classifies the input into 3 layers: offset, base, and feature.
  #These layers will be plotted separately to maximize control.
  lbase <- inpdf %>% dplyr::filter(layer == "base")
  lfeat <- inpdf %>% dplyr::filter(layer == "feature")
  loff <- inpdf %>% dplyr::filter(layer == "offset")


  #Heights are set relative to the base layer's height (for now).
  #if(is.null(hbase)) { hbase <- 0.2 }
  #if(is.null(hfeat)) { hfeat <- 2.4 * hbase }
  #if(is.null(hoff)) { hoff <- 0.8 * hbase }
  #hbase <- 0.2
  #hfeat <- 2.4 * hbase
  #hoff <- 0.8 * hbase


  #Alphas can be set by the user to control the transparency of the various layers.
  #Not exposing these to the user for now.
  #alpbase <- 1.0
  #alpfeat <- 1.0
  #alpoff <- 0.05


  #Fills and line colors for the various layers.
  #fillbase <- "black"
  #ffeat is set by the posdesc column
  #ffeat <- "yellow"
  #filloff <- "white"

  #colorbase <- "black"
  #colorfeat is set automatically.
  #colorfeat <- "yellow"
  #coloroff <- "gray"


  #Setting the number of X-axis breaks.
  #This is exposed to the user for customization.
  #Set to produce an X-axis tick every 100 residues (approximately) based on the
  #highest value of pos in the dataset (i.e., longest sequence).
  #nbreaks <- NULL
  if(is.null(nbreaks)){
    nbreaks <- base::round(base::max(inpdf$pos)/100)
  }



  #Plotting the basic plot first.
  #Need to account for the fact that a set of featureless sequences
  #might get plotted. So enclosing the main plot code in an if-else statement.
  if(nrow(lfeat) != 0){

    myplt <- ggplot2::ggplot() +
      ggplot2::geom_tile(data = lbase, ggplot2::aes(x = pos, y = prot_acc, fill = posdesc, height = hbase),
                         alpha = alpbase, fill = fillbase, color = colorbase) +
      ggplot2::geom_tile(data = lfeat, ggplot2::aes(x = pos, y = prot_acc, fill = posdesc, height = hfeat),
                         alpha = alpfeat) +
      ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n = nbreaks)) +
      ggplot2::theme_classic() +
      ggplot2::labs(x = "Position", y = "Sequence", fill = "Features")


    #VERY IMPORTANT - CONTROL OF GEOM_TEXT SIZE, AND WHETHER IT IS ACTIVATED AT ALL
    #OR NOT.
    #Only relevant if the lfeat layer exists at all, hence why it is inside this if statement.
    #label_size <- "auto"
    #label_size <- "0.5"
    #label_size <- "repel"; nudge_x <- 0.5; nudge_y <- 0.5
    if(label_size == "auto"){
      myplt <- myplt +
        ggplot2::geom_text(data = lfeat, ggplot2::aes(x = pos, y = prot_acc, label = label),
                           na.rm = TRUE)
    } else if(is.numeric(label_size)){
      if(label_size != 0){
        myplt <- myplt +
          ggplot2::geom_text(data = lfeat, ggplot2::aes(x = pos, y = prot_acc, label = label),
                             size = label_size, na.rm = TRUE)
      }
    } else if(label_size == "repel"){
      myplt <- myplt +
        ggrepel::geom_text_repel(data = lfeat, ggplot2::aes(x = pos, y = prot_acc, label = label),
                                 max.iter = 1.0, min.segment.length = 0,
                                 nudge_y = nudge_y, nudge_x = nudge_x,
                                 na.rm = TRUE, direction = "both")
    }



  } else{

    #Else just plot the base layer.
    myplt <- ggplot2::ggplot() +
      ggplot2::geom_tile(data = lbase, ggplot2::aes(x = pos, y = prot_acc, fill = posdesc, height = hbase),
                         alpha = alpbase, fill = fillbase, color = colorbase) +
      ggplot2::scale_x_continuous(breaks = scales::pretty_breaks(n = nbreaks)) +
      ggplot2::theme_classic() +
      ggplot2::labs(x = "Position", y = "Sequence", fill = "Features")

  }


  #Should sequence offsets be plotted if supplied?
  #Meaningful if and only if the offset layer has anything in it
  #to begin with.
  #show_offsets = TRUE
  if(show_offsets && nrow(loff) != 0){
    myplt <- myplt +
      ggplot2::geom_tile(data = loff, ggplot2::aes(x = pos, y = prot_acc, height = hoff),
                         alpha = alpoff, fill = filloff, color = coloroff)
  }

  #Use colorblind-friendly colors?
  #cbfcols <- TRUE
  if(cbfcols){
    #colpal <- sample(unique(viridis::turbo(n = 100)), length(unique(lfeat$posdesc)))
    colpal <- viridis::turbo(length(unique(lfeat$posdesc)))

    myplt <- myplt +
      ggplot2::scale_fill_manual(name = "", values = colpal, drop = TRUE)
  }

  #The user can also supply their own colors for the features.
  #featcols <- NULL
  if(!is.null(featcols)){
    if(length(featcols) != length(unique(lfeat$posdesc))){
      stop("Please supply as many colors as there are unique features in your annotation data!!")
    }
    myplt <- myplt +
      ggplot2::scale_fill_manual(name = "", values = featcols, drop = TRUE)
  }

  #Show the legend?
  #legend = TRUE
  if(!legend){
    myplt <- myplt + ggplot2::theme(legend.position = "none")
  }

  #Plotting to the plot pane.
  myplt

  #Returning object to parent environment so the user can manipulate it further if necessary.
  return(myplt)

}


