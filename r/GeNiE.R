#' Populate a GeNiE model with species data
#'
#' Given a GeNiE model file this function will attempt to populate it
#' with information in a dataframe with rows corresponding to species
#' and columns corresponding to nodes in the model with an appropriate
#' value (or NA) in each cell
#' @param fn filename of the GeNiE .xdsl file
#' @param spp a dataframe corresponding to rows of species and columns
#' of nodes in the GeNiE model, identifying the value of each node. The
#' dataframe must have a column called "sp_name" to identify each species.
#' @param fn.out output filename -- defaults to the input file name
#' @param NA.syn synonyms for NA's or missign values. Defaults to c("none")
#' @param nodelist a character vector of node names to attempt to populate.
#' Defaults to all nodes deterumined by column names. "sp_name" and
#' an optional column named "sp_code" are ignored.
#' @return No value returned, but creates a new .xdsl file containing the
#' newly populated GeNiE model
#' @export
populate.genie <- function(fn, spp, fn.out=fn, NA.syn=c("none"),
                           nodelist=NULL) {
    if (is.null(nodelist)) {
        nodelist=names(spp)[!(names(spp) %in% c("sp_name", "sp_code"))]
    }
    ## set the row names of spp to the spp_name column
    ## Right now this assumes no duplicates in spp_name
    row.names(spp) <- spp$spp_name

    d <- xmlInternalTreeParse(fn)
    # find the cases node
    dcases <- getNodeSet(d, "//cases")[[1]]

    splist <- spp$sp_name
    for (s in splist) {
        ## create a new case node for the species
        newspp <- newXMLNode("case", attrs=list(name=s), parent=dcases)

        ## now add a new evidence node for each variable
        for (i in nodelist) {
            v <- spp[s, i]
            if (is.na(v) | (v %in% NA.syn))
                next
            newXMLNode("evidence", attrs=list(node=i, state=v), parent=newspp)
        }
    }
    saveXML(d, fn.out)
    return(NULL)
}

#' Load a GeNiE model from an .xdsl file and return a catnet object
#'
#' Given a filename this function loads it as XML, extracts the nodes,
#' cpt's, and parents and returns a \code{\link{catNetwork}} object.
#' @param fn the filename to read in
#' @return a \code{\link{catNetwork}} object representation of the GeNiE network
#' given by the filename
#' @export
genie2catnet <- function(fn) {
    d <- xmlInternalTreeParse(fn)
    ## for exploring the file
    dnodes <- getNodeSet(d, "//node")
    dcpts <- getNodeSet(d, "//cpt")
    dpars <- getNodeSet(d, "//parents")

    ## pull out the nodes and node definition objects
    nodes <- xpathSApply(d, "//node", xmlGetAttr, "id")
    cpts <- xpathSApply(d, "//cpt", xmlGetAttr, "id")
    if (length(setdiff(nodes, cpts)) > 0) {
        stop("Node names and definitions do not match up!")
    }

    ## create a list of category levels for each model node
    catlist <- list(NULL)
    z <- 1
    for (i in nodes) {
        path <- paste0("//cpt[@id='", i, "'][last()]/state[@id]")
        states <- xpathSApply(d, path, xmlGetAttr, "id")
        catlist[[z]] <- states
        z <- z + 1
    }

    ## create a list of parent nodes for each node
    parlist <- list(NULL)
    z <- 1
    for (i in nodes) {
        path <- paste0("//cpt[@id='", i, "'][last()]/parents[last()]")
        parents.node <- getNodeSet(d, path)
        if (length(parents.node) == 0) {
            parlist[[z]] <- NULL
        } else {
            ## in case there's more than one parents definition (shouldn't be)
            ## take the last one
            parentstring <- xmlValue(parents.node[[length(parents.node)]])
            parlist[[z]] <- unlist(strsplit(parentstring, " "))
        }
        z <- z + 1
    }

    ## If the final element of a list is NULL it isn't counted for the
    ## length of the list. This will cause problems when we use cnNew
    ## to create the network, so force the list length to match
    ## the length of the nodes list -- any new nodes forced will have
    ## a NULL value, which is what we want.
    length(parlist) <- length(nodes)

    ## create a nested list of probabilities for each node given
    ## the probabilities for each parent node
    problist <- list(NULL)
    z <- 1
    for (i in nodes) {
        path <- paste0("//cpt[@id='", i, "'][last()]/probabilities[last()]")
        probs.node <- getNodeSet(d, path)
        if (length(probs.node) == 0) {
            problist[[z]] <- NULL  # should never happen
        } else {
            ## in case there's more than one probabilities definition
            ## (which should not happen) take the last one
            probstring <- xmlValue(probs.node[[length(probs.node)]])
            p <- as.numeric(unlist(strsplit(probstring, " ")))
            ## id number of this node
            nid <- which(nodes %in% i)
            ## get the number of levels in each of the node's parents
            levels <- parent.parms(nid, nodes, parlist, catlist)$levels
            if (is.null(levels)) {  ## this is a top node
                problist[[z]] <- p
            } else {
                ## add on the number of levels for this node
                levels <- c(levels, length(catlist[[nid]]))
                problist[[z]] <- parseprobs(levels, p)
            }
            z <- z + 1
        }
    }

    ## If the final element of a list is NULL it isn't counted for the
    ## length of the list. This will cause problems when we use cnNew
    ## to create the network, so force the list length to match
    ## the length of the nodes list -- any new nodes forced will have
    ## a NULL value, which is what we want.
    length(problist) <- length(nodes)

    return(cnNew(nodes, catlist, parlist, problist))
}
