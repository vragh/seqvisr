#' Convert a single sequence's quasi-standard TSV "row" input data.frame format to a position-per-row ggplot()-friendly format.
#'
#' @usage
#' tsvtogginp(inpdf = NULL)
#'
#' @param inpdf (character, mandatory) the name of the input data.frame().
#'
#' @details
#' Function to take a data.frame() indicating the sequence name, length, sequence start offset, feature height, feature color,
#' and feature start and end positions (one row for each feature), FOR ONE SEQUENCE, and convert this into (another) data.frame
#' that can be parsed by ggplot2. Namely, the output data.frame has one row for each residue in the sequence, and indicates the
#' positions (important for ggplot2), and the associated annotation (e.g., "ABC domain"). The function also has additional columns
#' indicating data for labeling that should be visualized on the plot or used to create it--this is basically done by taking the
#' feature's name and assigning it to a label column but only for a single residue at the middle of the feature. (If the labels of
#' all residues from each feature are displayed, the labeling will be illegible.)
#'
#' tsvtogginp() is invoked by pdomvisr() through the wrapper script tsvtgginp_multi().
#'
#' @return a data.frame() with one row per residue indicating the sequence name, position (indicated by that row), the description
#' associated with that position (e.g., "ABC domain", or NA for no description), a column indicating which feature the position
#' belongs to (e.g., say there are two domains named "ABC", this numerical column is a way to distinguish between them), a column
#' indicating the layer (offset, base, feature; see ?pdomvisr()), a numerical column indicating the height for the feature's tiles,
#' an alphanumeric column indicating the colors assigned to the item's tiles, labels for the particular feature (assigned to one
#' position within the feature; see ?pdomvisr()), and a column indicating the original positions of the rows (e.g., when a sequence
#' offset has been applied; see ?pdomvisr()). Column names are: prot_acc, pos, posdesc, fset, layer, height, color, label, orig_pos.
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
#' @importFrom rlang .data
#'
#' @export

tsvtogginp <- function(inpdf = NULL){

  #To avoid the "no visible binding for global variable" and "Undefined global functions or variables:"
  #warnings during devtools::check()
  #start_loc <- NULL
  #posdesc <- NULL
  #fset <- NULL
  #pos <- NULL
  #mid <- NULL
  #prot_acc <- NULL
  #orig_pos <- NULL


  if(base::is.null(inpdf)){ base::stop("No input data.frame provided!!") }


  #Check if the input data.frame has 5 columns, and rename them as follows.
  #c("prot_acc", "seq_len", "signature_desc", "start_loc", "stop_loc")
  if(base::ncol(inpdf) != 8){
    stopmsg <- base::strwrap("Please ensure all required input columns are provided!!
       These include must include the following in THIS order:
       Sequence name, sequence length, start offset, feature height, feature color, feature descriptor, feature start position, feature end position.
       Please supply one feature per row, and use multiple rows per sequence in case of multi-feature sequences!!")
    base::stop(stopmsg)
  }

  #Renaming columns properly, so that the user can supply pretty much whatever column
  #names they like.
  base::names(inpdf) <- c("prot_acc", "seq_len", "offset", "hfeat",
                          "colorfeat", "signature_desc", "start_loc", "stop_loc")



  #Pulling the annotations out separately.
  annots <- inpdf[ , -c(1:3)]
  #There may be more than one feature with the same name (e.g., PAS).
  #So using the fset column to distinguish them. This will be carried
  #over to the output data.frame.
  #Arrange annots by start_loc so that the fset assignment is sequential
  #in the following data.frames (i.e., feature at 200 gets fset 1, feature at
  #400 gets fset2).
  annots %<>% dplyr::arrange(.data$start_loc)
  annots$fset <- base::seq(base::nrow(annots))

  #Need the sequence length as a numeric, and the sequence name to
  #start creating the output data.frame.
  slen <- as.numeric(unique(inpdf[ , 2]))
  #sname <- as.character(unique(inpdf[ , 1]))
  #Do not as.character() for sname, or it will convert factors
  #into integers. Not helpful when users are passing that kind of
  #data in for plotting where they're already ordered things a certain
  #way.
  sname <- unique(inpdf[ , 1])

  #These steps below basically create the output data.frame
  #by taking a single-column data.frame consisting of the sequence
  #name, and replicating it slen times. Then a position column
  #(basically the row number) is added, yielding a data.frame with
  #one row for each residue.
  tmp <- base::data.frame(prot_acc = sname, stringsAsFactors = FALSE)
  outdf <- base::do.call(base::rbind, base::replicate(slen, tmp, simplify = FALSE))
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
  #All rows corresponding to a feature or feature will have this
  #overwritten with the string "feature".
  outdf$layer <- "base"


  #The base layer's height and color are set by the parent function, so
  #it just gets set to NA here.
  outdf$height <- NA
  outdf$colorfill <- NA


  #outdf could be looped through directly, and have positions
  #annotated within that correspond to features. However, that
  #is only a clean solution when overlapping annotations do not exist.
  #The more accommodating solution would be to actually have a separate
  #set of rows (as in outdf) for every unique feature in the sequence.
  #These concatenated with outdf (and a similar data.frame for the
  #optional sequence offset) is devoid of problems with overlapping
  #annotations, as each unique feature is represented by its own set
  #of tiles, and successive, overlapping sets of features are essentially
  #then tiles superimposed upon one another.

  #For this, the main structure is a loop that goes through each unique
  #feature (identified by fset) in the annots data.frame.
  for(i in 1:base::length(base::unique(annots$fset))){

    #Take one feature and its columns. Each feature only has one row
    #in the annotations table.
    curdat <- annots %>% dplyr::filter(.data$fset == base::unique(annots$fset)[i])

    #To create the outdf-equivalent positions data.frame,
    #need to calculate the length of the feature (i.e., number of
    #rows in that data.frame) first.

    #The as.numeric() transformation of the start and stop locations is important as
    #in some edge cases, the values in those cells are being treated as characters.
    curlen <- as.numeric(curdat$stop_loc) - as.numeric(curdat$start_loc)


    #If a feature is a single residue, it will a length of zero,
    #which is impractical here, as it still need a data.frame (albeit
    #with just one row). So setting curlen to 1 for such cases.
    if(curlen == 0){ curlen <- 1 }


    #Creating the outdf-equivalent called curdf. This is actually just
    #reusing the outdf creation code from above, but modified to take curlen
    #instead of slen in replicate().
    #See https://stackoverflow.com/a/48704973/9494044 for this do.call's base::rbind
    #no quotation weirdness.
    curdf <- base::do.call(base::rbind, base::replicate(curlen, tmp, simplify = FALSE))
    curdf$pos <- as.numeric(curdat$start_loc) + base::seq(base::nrow(curdf))

    #Adding in the position-wise descriptors, fset, and layer identifiers.
    #The layer identifier is "feature" in this case.
    curdf$posdesc <- curdat$signature_desc
    curdf$fset <- curdat$fset
    curdf$layer <- "feature"

    #Adding in the feature draw height values.
    curdf$height <- curdat$hfeat
    curdf$colorfill <- curdat$colorfeat

    #Collecting the data.frames for successive features in a data.frame using
    #bind_rows().
    if(i == 1){
      featdf <- curdf
    } else{
      featdf <- dplyr::bind_rows(featdf, curdf)
    }

  }
  rm(curdat, curdf, curlen)

  #Row-binding the features position-data.frame with outdf which has
  #the position-data.frame for the entire sequence background.
  #This is conditional. No point in adding that in for zero-length
  #"features" (which only exist to present a sequence w/o any
  #features).
  if(base::nrow(featdf) == 1){
    if(!base::is.na(featdf$posdesc)){
      outdf <- dplyr::bind_rows(outdf, featdf)
    }
  } else{
    outdf <- dplyr::bind_rows(outdf, featdf)
  }
  base::rm(featdf)



  #For each feature on the protein, I need to assign the descriptor
  #string to just one column belonging to that feature. That one
  #column's position is calculated as the start of the domain/feature
  #+ half the length of the feature (i.e., the mid-point of the feature).
  #This will be done through the label column.
  outdf$label <- NA

  outdf %<>%
    dplyr::group_by(.data$layer, .data$posdesc, .data$fset) %>%
    dplyr::arrange(.data$pos, .by_group = TRUE) %>%
    dplyr::mutate(mid = base::min(.data$pos) + base::round( (base::max(.data$pos)-base::min(.data$pos))/2 ,0)) %>%
    dplyr::mutate(label = base::ifelse(.data$pos == .data$mid, .data$posdesc, NA)) %>%
    dplyr::ungroup() %>%
    dplyr::select(-.data$mid)



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
    offdf %<>% dplyr::mutate(across(.cols = -.data$prot_acc, ~ NA))
    #Setting the layer column to "offset"
    offdf$layer <- "offset"
    #Setting the label to ""
    offdf$label <- NA

    #Offset heights and colors, like the base layer heights and colors, are
    #also set by the calling function. So they just get set to NA here.
    offdf$height <- NA
    offdf$colorfill <- NA

    #Replicating this single row the requisite number of times.
    offdf <- base::do.call(base::rbind, base::replicate(poffset, offdf, simplify = FALSE))

    #Updating the position counter for offdf.
    #This will count up from a negative number to 0.
    offdf$orig_pos <- rev(1 - base::seq(base::nrow(offdf)))
    offdf$pos <- base::seq(base::nrow(offdf))

    #Updating outdf's position counter column.
    outdf$orig_pos <- outdf$pos
    outdf$pos <- outdf$pos + poffset

    #Binding offdf to the head of outdf.
    outdf <- dplyr::bind_rows(offdf, outdf)

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





#' Wrapper script around tsvtogginp() to perform conversion of a data.frame() of multiple sequences and their feature annotations.
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

  if(base::is.null(inpdf)){ base::stop("No input data.frame provided!!") }


  #Check if the input data.frame has 5 columns, and rename them as follows.
  #c("prot_acc", "seq_len", "signature_desc", "start_loc", "stop_loc")
  if(base::ncol(inpdf) != 8){
    stopmsg <- base::strwrap("Please ensure all required input columns are provided!!
       These include must include the following in THIS order:
       Sequence name, sequence length, start offset, feature height, feature color, feature descriptor, feature start position, feature end position.
       Please supply one feature per row, and use multiple rows per sequence in case of multi-feature sequences!!")
    base::stop(stopmsg)
  }

  #Renaming columns properly, so that the user can supply pretty much whatever column
  #names they like.
  base::names(inpdf) <- c("prot_acc", "seq_len", "offset", "hfeat",
                          "colorfeat", "signature_desc", "start_loc", "stop_loc")

  #Applying tsvtogginp() to every unique sequence's set of rows in the input data.frame.
  outdf <- inpdf %>%
    dplyr::group_by(.data$prot_acc) %>%
    dplyr::group_map(~ tsvtogginp(.x), .keep = TRUE) %>%
    dplyr::bind_rows() %>%
    dplyr::ungroup()

  return(tibble::tibble(outdf))

}




#' Protein domain structure visualization in R.
#'
#' @description
#' pdomvisr() is a simple function to plot a diagram of the domain/feature structure of one or more
#' sequences. pdomvisr() uses ggplot2::ggplot() internally. The only mandatory input is a table
#' with the following information (in this particular order): sequence name, sequence length,
#' sequence offset, feature height, feature color, feature description, feature start coordinate,
#' feature end coordinate.
#'
#'
#' @usage pdomvisr(inpdat = NULL, mypath = NULL,
#' xlabel = "Position", ylabel = "Sequence",
#' leglabel = "Features", nbreaks = NULL,
#' hide_y_axis = FALSE, legend = TRUE,
#' show_offsets = TRUE, label_size = "auto",
#' hbase = 0.2, hoff = 0.8*hbase, alpbase = 1.0,
#' alpfeat = 1.0, alpoff = 0.05,
#' fillbase = "gray80", filloff = "gray60",
#' colorbase = "gray80", coloroff = "gray60",
#' nudge_x = 0.0, nudge_y = 0.5)
#'
#'
#' @param inpdat (character string or name of an object, mandatory) the character string may be the
#' name of a file or the full path to it (in which case mypath should be set to NULL). The file must
#' be a table containing the information necessary to plot the domain/feature structure diagram. Alternatively,
#' inpdat can also be supplied the name of an object in R's environment (e.g., a data.frame containing
#' the requisite data). This is useful when the input data needs to be pre-processed in R first. The
#' path/filename option is more suitable when pdomvisr() is being called only for plotting. In this case,
#' data.table's fread() is used to read the data into the function first. Please see the 'Details' section
#' for information on how the input data must be formatted.
#'
#' @param mypath (character string, optional) in the event that inpdat is supplied the name of a file, the
#' path to where this file is located can be supplied through mypath.
#'
#' @param xlabel (character string, optional) sets the label for the X-axis. Set to "" to disable. (Set to
#' "Position" by default.)
#'
#' @param ylabel (character string, optional) sets the label for the Y-axis. Set to "" to disable. (Set to
#' "Sequence" by default.)
#'
#' @param leglabel (character string, optional) sets the label for the legend. Set to "" to disable. (Set to
#' "Features" by default.)
#'
#' @param nbreaks (numeric, optional) controls the number of X-axis ticks in the plotted domain structure
#' diagram. If the user does not supply a number, this is automatically calculated to produce a tick every
#' 100 residues (based on the length of the longest sequence included in the plot).
#'
#' @param hide_y_axis (boolean, optional) controls whether the Y-axis grid line and ticks must be visible or
#' not. (Set to TRUE by default.)
#'
#' @param legend (boolean, optional) controls whether the legend associating the feature colors to the feature
#' descriptions should be plotted along with the main plot. (Set to TRUE by default.)
#'
#' @param show_offsets (boolean, optional) controls whether sequence offsets should be plotted or hidden.
#' A sequence offset is a whole number (supplied as a part of the input table) indicating how far off from
#' the actual first residue of the sequence the first residue indicated in the input data is. This is relevant
#' when plotting partial sequences for instance (e.g., an internal fragment). (Set to TRUE by default.)
#'
#' @param label_size (character or numeric, optional) controls whether the feature descriptions are displayed as
#' labels on the features. Also controls the size of the text if the labels are displayed. The size can be
#' controlled by supplying a positive integer > 0. Supplying 0 prevents the labels from being displayed. Passing
#' "auto" leaves the size estimation to R. If set to "repel", then the labels are drawn offset from the features
#' and connected to them by straight lines. If set to "repel", the arguments nudge_x and nudge_y (see below) can
#' be adjusted by the user to vary the positioning of the labels. Note: this argument's values do not affect the
#' legend. (Set to "auto" by default.)
#'
#' @param hbase (numeric, optional) controls the height of the tiles corresponding to the non-feature portions
#' of the sequence. (Set to 0.2 by default.)
#'
#' @param hoff (numeric, optional) controls the height of the tiles representing the sequence offset. Under
#' default settings, this scales automatically with hbase. (Set to 0.8 * hbase by default.)
#'
#' @param alpbase (numeric, optional) controls the alpha level of the tiles representing the non-feature portions
#' of the sequence. (Set to 1.0 by default.)
#'
#' @param alpfeat (numeric, optional) controls the alpha level of the tiles representing the features of
#' the sequence(s). (Set to 1.0 by default.)
#'
#' @param alpoff (numeric, optional) controls the alpha level of the tiles representing the sequence offset. (Set
#' to 0.05 by default.)
#'
#' @param fillbase (character, optional) fill color for the tiles representing the non-feature portions of the
#' sequence. Any value accepted by ggplot2's "fill" is accepted here, as this just passes the value on
#' to that particular argument. (Set to "black" by default.)
#'
#' @param filloff (character, optional) fill color for the tiles representing the offset sequence. Any value
#' accepted by ggplot2's "fill" is accepted here, as this just passes the value on to that particular argument.
#' (Set to "white" by default.)
#'
#' @param colorbase (character, optional) line color for the tiles representing the non-feature portions of the
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
#' pdomvisr() plots a domain/feature structure diagram given the coordinates of features in one or more sequences.
#'
#' The only mandatory input is a table with the following information (in this particular order): sequence name,
#' sequence length, sequence offset, feature height, feature color, feature description, feature start coordinate,
#' feature end coordinate. Most column names are self-explanatory, and should be readily produced by most feature
#' annotation tools (or should be producible by hand). The sequence offset, feature height, and feature color columns
#' must be typically defined by the user (no annotation tool produces these). (More on these columns later.)
#'
#' Each row in the input table should correspond to the coordinates for a particular feature in a particular sequence.
#' Therefore, if a sequence contains more than one feature, it will have to be represented by as many rows as there are
#' features in it. The sequence name, sequence length, and sequence offset columns will (unfortunately) have to be repeated
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
#' includes all positions that are neither a part of the offset sequence nor a part of a feature in the sequence. The
#' feature layer includes all residues that belong to a feature. Additional columns exist to indicate labeling (and other
#' logistics) for the features. As there is no upper or lower bound on the length of a feature as far as pdomvisr() is
#' concerned, even single residues (e.g., an active site) can be annotated by adding a row carrying the same value for the
#' start and end positions for the "feature".
#'
#' pdomvisr() plots the diagram on the current device (the plot pane in RStudio, for example) and also returns the ggplot2
#' object itself to the parent environment from which the function was called. The user therefore has complete control over
#' the output.
#'
#' About the user defined columns:
#'
#' #' The offset column is strictly optional. The objective of this column is to indicate how far away from the actual start
#' of the sequence the "indicated" start of the sequence is. For example, if a sequence is listed as being 200 residues long,
#' and has an offset value of 10, this implies that the sequence is actually 210 residues long, but only the last 200 residues
#' "exist". The offset column is useful for visually demonstrating that a sequence is partial. If none of the sequences are
#' partial and/or in no need of an offset, the values in this column should be set to 0. It is mandatory that this column exists
#' (with a default value of 0 in all rows) even if it is not in use!!
#'
#' The feature height column is mandatory and must be numeric. It defines the height of the tiles corresponding to a particular feature.
#' Unless the non-feature portions of the sequences have been assigned heights greater than 1.0, this column need not contain values greater
#' than 1.0 either (a value of 0.4 should suffice). Each unique feature, can of course be assigned its own height value. pdomvisr() permits
#' inclusion of "dummy" rows where only the sequence identifier, sequence length (and optionally, offset) are set; these rows are to visually
#' represent sequences with no domains in them. For such cases, the height should be set to 0. Finally, there is no reason that different domains
#' cannot have different heights, and therefore every row in the input data.frame can have a different height assigned to it.
#'
#' The feature color column is analogous to the feature height column, but assigns the colors for the tiles instead. Everything discussed above
#' also applies for this column. For "dummy" rows, the feature color must be set to NA.
#'
#' @return A ggplot2 object is returned to the parent environment for plotting and/or further downstream processing/manipulation.
#'
#' @note
#' In some cases it may be necessary to include a sequence that has no annotated features. Most feature annotation tools would
#' not include a row for such sequences. But as long as the user adds a row to the input table manually with the sequence
#' description column set to NA, and the start and end positions of the feature set to 0, a sequence with no annotated features
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
#' @importFrom rlang .data
#'
#' @import dplyr ggplot2
#'
#' @export

pdomvisr <- function(inpdat = NULL, mypath = NULL,
                     xlabel = "Position", ylabel = "Sequence",
                     leglabel = "Features", nbreaks = NULL,
                     hide_y_axis = FALSE, legend = TRUE,
                     show_offsets = TRUE, label_size = "auto",
                     hbase = 0.2, hoff = 0.8*hbase, alpbase = 1.0,
                     alpfeat = 1.0, alpoff = 0.05,
                     fillbase = "gray80", filloff = "gray60",
                     colorbase = "gray80", coloroff = "gray60",
                     nudge_x = 0.0, nudge_y = 0.5){

  #To avoid the "no visible binding for global variable" and "Undefined global functions or variables:"
  #warnings during devtools::check()
  #pos <- NULL
  #prot_acc <- NULL
  #posdesc <- NULL
  #label <- NULL
  #Longer term solution is this:
  #https://community.rstudio.com/t/how-to-solve-no-visible-binding-for-global-variable-note/28887/2

  #Stop if input is empty.
  if(base::is.null(inpdat)){
    base::stop("The feature annotation table is a mandatory input!! Please provide this via the inpdat argument!!")
  }

  #Checking to see if what is supplied to inpdat is a path or an object.
  if(base::is.data.frame(inpdat) | data.table::is.data.table(inpdat) | tibble::is_tibble(inpdat)){
    #Data is already loaded into R, do nothing.
    inpdf <- inpdat
  } else{
    #Setting up the path to the object
    if(!base::is.null(mypath)){
      inpdat <- base::paste0(mypath, "/", inpdat)
    }

    #Reading in the data using data.table::fread.
    inpdf <- data.table::fread(inpdat)
  }


  #Converting data for plotting with tsvtogginp_multi.
  inpdf <- tsvtogginp_multi(inpdf)



  #The tsvtogginp() setup (accessed through tsvtogginp_multi)
  #classifies the input into 3 layers: offset, base, and feature.
  #These layers will be plotted separately to maximize control.
  lbase <- inpdf %>% dplyr::filter(.data$layer == "base")
  lfeat <- inpdf %>% dplyr::filter(.data$layer == "feature")
  loff <- inpdf %>% dplyr::filter(.data$layer == "offset")


  #Heights are set relative to the base layer's height (for now).
  #if(is.null(hbase)) { hbase <- 0.2 }
  #if(is.null(hfeat)) { hfeat <- 2.4 * hbase }
  #if(is.null(hoff)) { hoff <- 0.8 * hbase }
  #  hbase <- 0.2
  #hfeat <- 2.4 * hbase
  #  hoff <- 0.8 * hbase

  #Alphas can be set by the user to control the transparency of the various layers.
  #Not exposing these to the user for now.
  #  alpbase <- 1.0
  #  alpfeat <- 1.0
  #  alpoff <- 0.05

  #Fills and line colors for the various layers.
  #  fillbase <- "gray"
  #ffeat is set by the posdesc column
  #ffeat <- "yellow"
  #  filloff <- "white"

  #  colorbase <- "gray"
  #colorfeat is set automatically.
  #colorfeat <- "yellow"
  #  coloroff <- "gray"


  #Setting the number of X-axis breaks.
  #This is exposed to the user for customization.
  #Set to produce an X-axis tick every 100 residues (approximately) based on the
  #highest value of pos in the dataset (i.e., longest sequence).
  #nbreaks <- NULL
  if(is.null(nbreaks)){
    nbreaks <- base::round(base::max(inpdf$pos)/100)
  }
  #Based on this https://stackoverflow.com/a/51019485/9494044
  #Pretty breaks will start with 0.
  pretty_br <- base::pretty(inpdf$pos, n = nbreaks)
  #So setting first break to 1 manually.
  pretty_br[1] <- 1



  #Plotting the basic plot first.
  #Need to account for the fact that a set of featureless sequences
  #might get plotted. So enclosing the main plot code in an if-else statement.
  if(base::nrow(lfeat) != 0){

    #Need to arrange lfeat in such a manner that the smallest (potentially overlapping)
    #features get plotted last.
    lfeat %<>%
      dplyr::group_by(.data$prot_acc, .data$fset) %>%
      dplyr::mutate(featlen = n()) %>%
      dplyr::arrange(.data$featlen, .by_group = TRUE) %>%
      dplyr::ungroup() %>%
      dplyr::select(-.data$featlen)

    #Feature color and fill are set directly via user inputs from the input data.frame.
    #The feature heights, defined by the user, are set here.
    myplt <- ggplot2::ggplot() +
      ggplot2::geom_tile(data = lbase,
                         ggplot2::aes(x = .data$pos, y = .data$prot_acc, fill = .data$posdesc, height = as.numeric(hbase)),
                         alpha = alpbase, fill = fillbase, color = colorbase) +
      ggplot2::geom_tile(data = lfeat,
                         ggplot2::aes(x = .data$pos, y = .data$prot_acc, fill = .data$posdesc, height = as.numeric(.data$height)),
                         alpha = alpfeat) +
      #ggplot2::scale_x_continuous(limits = c(1, NA), breaks = scales::extended_breaks(n = nbreaks)) +
      ggplot2::scale_x_continuous(breaks = pretty_br) +
      ggplot2::theme_classic() +
      ggplot2::labs(x = xlabel, y = ylabel, fill = "Features")


    #VERY IMPORTANT - CONTROL OF GEOM_TEXT SIZE, AND WHETHER IT IS ACTIVATED AT ALL
    #OR NOT.
    #Only relevant if the lfeat layer exists at all, hence why it is inside this if statement.
    #label_size <- "auto"
    #label_size <- "0.5"
    #label_size <- "repel"; nudge_x <- 0.5; nudge_y <- 0.5
    if(label_size == "auto"){
      myplt <- myplt +
        ggplot2::geom_text(data = lfeat, ggplot2::aes(x = .data$pos, y = .data$prot_acc, label = .data$label),
                           na.rm = TRUE)
    } else if(is.numeric(label_size)){
      if(label_size != 0){
        myplt <- myplt +
          ggplot2::geom_text(data = lfeat, ggplot2::aes(x = .data$pos, y = .data$prot_acc, label = .data$label),
                             size = label_size, na.rm = TRUE)
      }
    } else if(label_size == "repel"){
      myplt <- myplt +
        ggrepel::geom_text_repel(data = lfeat, ggplot2::aes(x = .data$pos, y = .data$prot_acc, label = .data$label),
                                 force_pull = 0,
                                 max.time = 1, max.iter = 1e5, min.segment.length = 0,
                                 nudge_y = nudge_y, nudge_x = nudge_x,
                                 na.rm = TRUE, direction = "both",
                                 #segment.curvature = -1e-20,
                                 arrow = ggplot2::arrow(length = ggplot2::unit(0.015, "npc")))
    }

  } else{

    #Else just plot the base layer.
    myplt <- ggplot2::ggplot() +
      ggplot2::geom_tile(data = lbase, ggplot2::aes(x = .data$pos, y = .data$prot_acc, fill = .data$posdesc, height = as.numeric(hbase)),
                         alpha = alpbase, fill = fillbase, color = colorbase) +
      #ggplot2::scale_x_continuous(limits = c(1, NA), breaks = scales::extended_breaks(n = nbreaks)) +
      ggplot2::scale_x_continuous(breaks = pretty_br) +
      ggplot2::theme_classic() +
      ggplot2::labs(x = "Position", y = "Sequence", fill = "Features")

  }


  #Should sequence offsets be plotted if supplied?
  #Meaningful if and only if the offset layer has anything in it
  #to begin with.
  #show_offsets = TRUE
  if(show_offsets && base::nrow(loff) != 0){
    myplt <- myplt +
      ggplot2::geom_tile(data = loff, ggplot2::aes(x = .data$pos, y = .data$prot_acc, height = hoff),
                         alpha = alpoff, fill = filloff, color = coloroff)
  }

  #Use colorblind-friendly colors?
  #cbfcols <- TRUE
  #if(cbfcols){
  #  #colpal <- sample(unique(viridis::turbo(n = 100)), length(unique(lfeat$posdesc)))
  #  colpal <- viridis::turbo(length(unique(lfeat$posdesc)))
  #
  #  myplt <- myplt +
  #    ggplot2::scale_fill_manual(name = "", values = colpal, drop = TRUE)
  #}

  #The user can also supply their own colors for the features.
  #featcols <- NULL
  #if(!is.null(featcols)){
  #  if(length(featcols) != length(unique(lfeat$posdesc))){
  #    stop("Please supply as many colors as there are unique features in your annotation data!!")
  #  }
  #  myplt <- myplt +
  #    ggplot2::scale_fill_manual(name = "", values = featcols, drop = TRUE)
  #}

  #Colors are supplied by the user via the input data.frame, and this is a column in lfeat, loff, and lbase.
  #The colors are stored in the column colorfill. Only feature rows have colors, base and offset rows just have
  #NA since those colors are set universally.
  #
  #inplabs <- lfeat %>% distinct(colorfill, .keep_all = TRUE) %>% ungroup() %>% select(c(colorfill, posdesc))
  inpcols <- lfeat %>%
    dplyr::group_by(.data$prot_acc, .data$fset) %>%
    dplyr::distinct(.data$colorfill, .keep_all = TRUE) %>%
    dplyr::ungroup() %>%
    dplyr::select(c(.data$posdesc, .data$colorfill)) %>%
    tibble::deframe()
    #select(colorfill)

  myplt <- myplt +
    ggplot2::scale_fill_manual(name = "Features", values = inpcols,
                               labels = ggplot2::waiver(), guide = ggplot2::guide_legend())

  #Show the legend?
  #legend = TRUE
  if(!legend){
    myplt <- myplt + ggplot2::theme(legend.position = "none")
  }

  #Hide Y-axis and ticks?
  if(hide_y_axis){
    myplt <- myplt +
      ggplot2::theme(axis.line.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank())
  }


  #Plotting to the plot pane.
  myplt

  #Returning object to parent environment so the user can manipulate it further if necessary.
  return(myplt)

}
