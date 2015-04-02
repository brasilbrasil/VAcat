#' Get the relative weight of a parent node's value
#'
#' Calculate the relative weight of a parent node on a child node.
#' Takes a vector of possible values of the parent node and returns
#' a scaled weighting on the child node, usually in (0, 1) or (-1, 0, 1).
#'
#' @param obs the observed level of the parent node.
#' @param values a vector of possible values the parent node might have,
#' in order from most positive to most negative.
#' @param weight a vector of weights to assign each possible value. If
#' NULL it will assign weights of 1, 0 to values with two levels,
#' or 1, 0, -1 to values with three levels. For more than three levels
#' weights needs to be explicitly assigned.
#' @return a relative weight for a parent node value. Using default
#' parameters either 1, 0, or -1.
#' @export
get.relative.weight <- function(obs, values, wt=NULL) {
    if (!(obs %in% values))
        stop("Observed value ", obs, " not in list of supplied values.")
    nv <- length(values)
    if (nv != length(unique(values)))
        stop("Supplied list of possible values are not unique.")
    if (is.null(wt)) {
        if (nv == 2)
            wt <- c(1, 0)
        else if (nv == 3)
            wt <- c(1, 0, -1)
        else
            wt <- rep(0, nv + 1)
    }
    if (length(wt) != nv)
        stop("No default weights assigned for values of ",
             "length ", nv, ".")
    return(wt[which(values %in% obs)])
}

#' Calculate a set of parent node weights
#'
#' Calculates a set of output state probabilities given a dataframe
#' with rows for each parent of a node and columns
#' describing how to calculate the weights for each parent, a set of
#' observed values for each parent, and a base probability of the
#' positive (mathematically) outcome.
#'
#' @param df a dataframe with rows for each parent node and columns
#' as described in notes below
#' @param levels a vector of observed states for the parent nodes
#' @param base the base (starting) probability of the positive outcome
#' @param minw the minimum probability of the positive outcome. If set to
#' NULL no minimum is checked. Defaults to 0.1.
#' @param maxw the maximum probability of the postive outcome. If set to
#' NULL no maximum is checked. Defaults to 0.9.
#' @return a probability for the node to have the "Favorable" value
#' as determined by labels used outside of this function.
#' @note the dataframe df must have a row for each parent node and columns
#' parent containing the name of the parent node, func containing the name
#' of the function used to compute the weights. "func" can be "general"
#' meaning the generic linear model approach is used, or the name of the
#' weighting function to be used (retrieved from the local
#' environment with the \code{\link{get}} function)
#' which is then passed the df dataframe and
#' levels vector.
#'
#' The weight column is a vector of weights for the parent nodes when
#' using the general generic linear approach.
#'
#' There are also columns "hilow1", "hilow2", and "hilow3"
#' listing possible values of the parent node in order from most positive
#' to least positive effect on the output.
#'
#' The parent node values in levels must be in the same order as the rows
#' of df.
#' @export
get.favorable.weight <- function(df, levels, base=0.1,
                                 minw=0.1, maxw=0.9) {
    if (length(levels) != nrow(df))
        stop("get.node.weights: number of levels does not match the ",
             "rows in df.")
    if (df$func[1] == "general") {
        ## make a list of possible parent values from hilow columns
        pvs <- NULL
        for (i in 1:nrow(df)) {
            foo <- df[i, c("hilow1", "hilow2", "hilow3")]
            if (foo[1, 3] == "" | is.na(foo[1, 3]))
                foo <- foo[, -3]
            pvs[[i]] <- as.vector(t(foo))
        }
        ## check to make sure the values in levels are legitimate
        for (i in 1:nrow(df)) {
            if (!(levels[i] %in% pvs[[i]])) {
                cat("\ni = ", i, "level =", levels[i], "\n")
                stop("get.node.weights: given a level not in one of",
                     " the hilow columns.")
            }
        }
        ## calculate the weight of the favorable outcome
        outwt <- base
        for (i in 1:nrow(df)) {
            outwt <- outwt + df$weight[i] *
                get.relative.weight(levels[i], pvs[[i]])
        }
        ## truncate to minimum or maximum if set
        if (!is.null(minw))
            outwt <- max(minw, outwt)
        if (!is.null(maxw))
            outwt <- min(outwt, maxw)
    } else {
        wtfun <- get(df$func[1])
        if (is.null(wtfun))
            stop("get.favorable.weight: weighting function ", df$func[1],
                 " not found in the environment.")
        outwt <- wtfun(df, levels)
    }
    return(outwt)
}

#' Build cpts for a node given a description of the parents
#'
#' Constructs a dataframe with columns describing all possible values
#' of the parent nodes and probabilities that their child node will
#' have a certain categorical value. Only works for computing nodes
#' with two possible categorical values.
#'
#' @param df a dataframe with rows for each parent node and columns
#' as described in notes below
#' @return a dataframe with columns for each of the parent nodes and the
#' two output probabilities, and a row for each possible combination of
#' the parent values. The positive output probability is named with
#' the value in the "positive" column of the dataframe.
#' @note the dataframe df must have a row for each parent node and columns
#' "parent" containing the name of the parent node, func containing the name
#' of the function used to compute the weights as described in
#' \code{\link{get.favorable.weight}}. If the function is not "general"
#' it must also have columns for the base, weight, and hilow as described
#' in \code{\link{get.favorable.weight}}.
#' @export
build.cpt <- function(df) {
    ## build a list of vectors of possible values of the parent nodes
    pvs <- NULL
    for (i in 1:nrow(df)) {
        if (df$func[i] == "general") {
            foo <- df[i, c("hilow1", "hilow2", "hilow3")]
            if (foo[1, 3] == "" | is.na(foo[1, 3]))
                foo <- foo[, -3]
            pvs[[i]] <- as.vector(t(foo))
        } else {
            pvs[[i]] <- unlist(strsplit(df$hilow1[i], ";"))
        }
    }
    ## create a dataframe of all possible parent combinations, with
    ## columns for the probabilities of node values
    parent.matrix <- expand.grid(pvs, KEEP.OUT.ATTRS=FALSE,
                                 stringsAsFactors=FALSE)
    names(parent.matrix) <- df$parent
    parent.matrix$positive <- rep(0.5, nrow(parent.matrix))
    parent.matrix$negative <- parent.matrix$positive
    posname <- df$positive[1]
    names(parent.matrix)[c(nrow(df)+1, nrow(df)+2)] <-
        c(posname, paste0("not", posname))
    ## calculate the positive weight for each row of the parent matrix
    for (i in 1:nrow(parent.matrix)) {
        obslevels <- as.vector(t(parent.matrix[i, 1:(ncol(parent.matrix)-2)]))
        parent.matrix[i, posname] <-
            get.favorable.weight(df, obslevels, base=df$base[1])
    }
    ## the contrapostive weight is simply its complement
    parent.matrix[, paste0("not", posname)] <- 1 - parent.matrix[, posname]
    return(parent.matrix)
}


#' Set the cpt for a node in a catNet model
#'
#' Uses a dataframe of parent node value combinations and value probabilities
#' (for instance, as generated by \code{\link{build.cpt}})
#' to set the cpt for a node in a catNet model.
#'
#' @param N a catNet model
#' @param nname the name of a node to alter in model N
#' @param cptdf a dataframe consisting of rows for all possible combinations
#' of node nname's parents, and probabilities for node nname's value in each
#' combination. More details below.
#' @return a catNet model the same as N, with the probabilities for node
#' nname updated to match those in cptdf.
#' @note the dataframe cptdf must have column names for each of the parent
#' nodes, plus two columns for the probabilities of node states for each
#' combination. The second column from the right (column ncol(cptdf) - 1)
#' must have a name matching one of the defined values for node nname.
#' @export
set.cpt <- function(N, nname, cptdf) {
    require(catnet)
    nodelist <- cnNodes(N)
    if (!(nname %in% nodelist))
        stop("Node ", nname, " was not found in the model.")
    nodenum <- which(nodelist %in% nname)
    ## set the notpositive cpt column to the proper value
    posname <- names(cptdf)[ncol(cptdf)-1]
    nodevals <- N@categories[[nodenum]]
    if (length(nodevals) != 2)
        stop("Trying to set the cpt of a node with other than 2 values.")
    names(cptdf)[ncol(cptdf)] <- nodevals[(nodevals[1] == posname) + 1]
    ## make sure the parent node names are legitimate
    parentnodes <- unlist(cnParents(N, nname))
    parentnames <- names(cptdf)[1:(ncol(cptdf) - 2)]
    if (any(!(parentnames %in% parentnodes)))
        stop("One or more column names of cptdf not model node parents.")
    nparents <- length(parentnames)
    ## Re-order the cpt dataframe so columns are in the right order
    cptdf <- cptdf[, c(parentnodes, nodevals)]
    ## Re-order the rows so they match the order of the model
    ordermat <- matrix(as.numeric(NA), nrow=nrow(cptdf), ncol=nparents)
    for (i in 1:nparents) {
        colnodenum <- which(nodelist %in% names(cptdf)[i])
        colvals <- N@categories[[colnodenum]]
        for (j in 1:length(colvals)) {
            ordermat[, i][cptdf[, i] == colvals[j]] <- j
        }
    }
    if (any(is.na(ordermat)))
        stop("Found node values in the given cpt that don't match ",
             "node values in the model.")
    ordermat <- as.data.frame(ordermat)
    if (ncol(ordermat) == 1) {
        o <- order(ordermat[, 1])
    } else {
        o <- do.call(order, ordermat[, 1:nparents])
    }
    cptdf <- cptdf[o, ]

    ## get the current cpt structure for this node
    plist <- N@probabilities[[nodenum]]
###    if (length(plist) != length(parentnames))
###        stop("The number of parent nodes in the given cpt does not ",
###             "match the number of parents in the model.")
    ## recursively parse the list, setting the cpt probabilities at
    ## the lowest level
    traverse.list <- function(pl, df) {
        if (ncol(df) == 3) {
            if (length(pl) != nrow(df))
                stop("Error traversing list, number of parent node ",
                     "combinations don't match.")
            for (i in 1:nrow(df)) {
                pl[[i]] <- as.vector(t(df[i, 2:3]))
            }
        } else {
            ## loop over the possible values of this parent node
            for (i in 1:length(pl)) {
                subdf <- df[df[, 1]==unique(df[, 1])[i], 2:ncol(df)]
                pl[[i]] <- traverse.list(pl[[i]], subdf)
            }
        }
        return(pl)
    }
    ## set the cpt for the node, and return the altered model
    N@probabilities[[nodenum]] <- traverse.list(plist, cptdf)
    return(N)
}

#' Set or get a vector of node values
#'
#'
#' Uses a dataframe of parent node value combinations and value probabilities
#' (for instance, as generated by \code{\link{build.cpt}})
#' to set the cpt for a node in a catNet model.
#'
#' @param N a catNet model
#' @param nname the name of a node to alter in model N
#' @param value the desired value of the node. Stops on an error if the
#' value is not one specified in the catNet model. If NULL it returns a
#' vector of legal values.
#' @return a named vector, whose name is the specified nodename (nname). If
#' value is NULL the vector contains all legal values.
#' @export
set.node.levels <- function(N, nname, value=NULL) {
    require(catnet)
    nodelist <- cnNodes(N)
    if (!(nname %in% nodelist))
        stop("Node ", nname, " was not found in the model.")
    nodenum <- which(nodelist %in% nname)
    nodevalues <- N@categories[[nodenum]]
    if (is.null(value)) {
        outvals <- nodevalues
    } else {
        if (!(value %in% nodevalues))
            stop("Node value not in list specified by model.")
        outvals <- value
    }
    outvals <- list(foo=outvals)
    names(outvals) <- nname
    return(outvals)
}

#' Summarize the posterior output nodes of a model
#'
#' Given a model and a set of data, calculate the distribution of
#' values for a set of unknown, modeled nodes.
#'
#' @param N a catNet model
#' @param data a dataframe of data for the model. Rows represent
#' observations and columns are nodes in the model. NA values in a
#' column mean the node is to be estimated.
#' @param targets a character vector of node names to be summarized
#' @return a dataframe with columns for the node name, the value of
#' the node, and the proportion of results with that value.
#' @export
summarize.output.nodes <- function(N, data, targets) {
    require(catnet)
    mnodes <- cnNodes(N)
    if (length(mnodes) != ncol(data))
        stop("Number of model nodes and dataframe columns do not match.")
    if (length(names(data)) != length(unique(names(data))))
        stop("Column names in dataframe are not unique.") ## not possible
    if (any(!(names(data) %in% mnodes)))
        stop("Column names in dataframe do not match model nodes.")
    if (any(!(targets %in% names(data))))
        stop("Not all targets found in node names.")
    if (any(!is.na(data[, targets])))
        stop("Target nodes must be missing (NA) in the dataframe.")
    ## fit the model given the data; extract target nodes
    data <- data[, mnodes]
    newdata <- cnPredict(N, data)[, targets]
    ## if just one column, force it into a dataframe
    if (is.null(dim(newdata))) {
        newdata <- data.frame(as.character(newdata), stringsAsFactors=FALSE)
        names(newdata) <- targets
    }
    ## placeholder object to hold outputs
    outdf <- NULL
    for (i in 1:ncol(newdata)) {
        nodename <- names(newdata)[i]
        ## get the possible levels of this node
        nodenum <- which(mnodes %in% nodename)
        nodelevels <- N@categories[[nodenum]]
        foo <- data.frame(node=rep(nodename, length(nodelevels)),
                          value=nodelevels,
                          p=rep(0, length(nodelevels)),
                          stringsAsFactors=FALSE)
        bar <- table(newdata[, i]) / nrow(newdata)
        foo$p <- bar[foo$value]
        outdf <- rbind(outdf, foo)
        rm(foo, bar)
    }
    return(outdf)
}
