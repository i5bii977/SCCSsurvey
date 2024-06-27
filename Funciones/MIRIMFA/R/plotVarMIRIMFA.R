#' @export
#' @import polycor
plotVarMIRIMFA <- function (object,
                               comp = c(1, 2),
                               varNames = FALSE,
                               cex = 3,
                               pch = 19,
                               col = NULL,
                               alpha = 0.7,
                               spty = TRUE,
                               cutoff = 0,
                               radIn = 0.5,
                               overlap = TRUE,
                               ncols = 2,
                               legTitle = "Variables",
                               graph = TRUE) {

  ## Checking general input arguments & initialization of variables####
  ## object
  if (!any(class(object) != "MIRIMFA")) {
    stop("'object' must be an object of class 'MIRIMFA'.", call. = FALSE)
  }

  RefData <- split.default(object[["res.imputeMFA"]],
                           rep(seq_along(object[["call"]][["group"]]),
                               object[["call"]][["group"]]))
  imputations <- object[["res.MI"]]
  for (i in seq_along(imputations)) {
    imputations[[i]] <- split.default(imputations[[i]],
                                      rep(seq_along(object[["call"]][["group"]]),
                                          object[["call"]][["group"]]))
  }
  comprConf <- object[["compromise"]]
  dtab <- seq_along(RefData)
  if (is.null(object[["call"]][["name.group"]])) {
    nmTables <- paste0("G", dtab)
  } else {
    nmTables <- object[["call"]][["name.group"]]
  }
  nbTables <- length(RefData)

  ## comp
  ncomp <- object[["call"]][["ncp"]]

  if (length(comp) != 2) {
    stop("'comp' must be a vector of length equal to 2", call. = FALSE)
  } else {
    if (any(!is.finite(comp)))
      stop("the elements of 'comp' must be positive integers",
           call. = FALSE)

    if (!is.numeric(comp) || any(comp < 1))
      stop("the elements of 'comp' must be positive integers",
           call. = FALSE)

    if (any(comp > ncomp))
      stop("the elements of 'comp' must be smaller or equal than ",
           ncomp, ".", call. = FALSE)
  }

  comp <- round(comp)

  ## col

  ##- internal function for character color checking -----------#
  ##------------------------------------------------------------#
  isColor <- function(x) { vapply(x, function(x) {
    tryCatch(is.matrix(col2rgb(x)), error=function(e) FALSE) },
    TRUE)
  }
  ##------------------------------------------------------------#

  ##- brewer palette 'Set1' ------------------------------------#
  brewerPal <- c( "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
                  "#FF7F00", "#FFFF33", "#A65628", "#F781BF" )
  ##------------------------------------------------------------#


  if (is.null(col)) {
    col <- rainbow(length(colnames(object[["res.imputeMFA"]])))
    names(col) <- colnames(object[["res.imputeMFA"]])
  } else {
    if (length(col) != length(colnames(object[["res.imputeMFA"]]))) {
      stop("'col' must be a color names vector of length ",
           length(colnames(object[["res.imputeMFA"]])), ".", call. = FALSE)
    } else {
      if (any(!isColor(col))) {
        stop("'col' must be a character vector of recognized colors.",
             call. = FALSE)
      }
    }
  }

  ## varNames
  flagVarNames <- rep(FALSE, nbTables)

  if (is.logical(varNames)) {
    if (varNames) {
      varNames <- nmTables
      flagVarNames <- rep(TRUE, nbTables)
    }
  } else {
    if (is.character(varNames)) {
      if (!any(varNames %in% nmTables)) {
        stop("One element of 'varNames' does not match with the ",
             "names of the data tables. \n'varNames' should be one of: ",
             toString(nmTables), call. = FALSE)
      } else {
        id <- (nmTables %in% varNames)
        flagVarNames[id] <- TRUE
      }
    } else {
      stop("Incorrect value for 'varNames'", call. = FALSE)
    }
  }

  ## cex
  if (length(cex) == 1) {
    if (is.null(cex) || !is.numeric(cex) || !is.finite(cex) || (cex < 0)) {
      stop("'cex' must be a positive number", call. = FALSE)
    }
    cex <- rep(cex, nbTables)
  } else {
    if (length(cex) > nbTables) {
      stop("'cex' must be a numeric vector of either length 1",
           " or length ", nbTables, call. = FALSE)
    }
    if (any(is.null(cex)) | any(!is.numeric(cex)) | any(!is.finite(cex)) |
        any(cex < 0)) {
      stop("'cex' must be a numeric vector of either length 1",
           " or length ", nbTables, call. = FALSE)
    }
  }

  ## pch
  if (length(pch) == 1) {
    pch <- rep(pch, nbTables)
  } else {
    if (length(pch) > nbTables) {
      stop("'pch' must be a vector of either length 1 or length ",
           nbTables, call. = FALSE)
    }
  }

  ## alpha
  if (length(alpha) != 1) {
    stop("'alpha' must be a single value", call. = FALSE)
  }

  if (is.null(alpha) || !is.numeric(alpha) || !is.finite(alpha)) {
    stop("'alpha' transparency value must be numeric", call. = FALSE)
  }

  if ((alpha < 0) || (alpha > 1)) {
    stop("'alpha' transparency value ", alpha,
         ", outside the interval [0,1]", call. = FALSE)
  }

  ## spty
  if (length(spty) != 1 || !is.logical(spty)) {
    stop("'spty' must be either TRUE or FALSE", call. = FALSE)
  }

  ## cutoff
  if (length(cutoff) != 1) {
    stop("'cutoff' must be a single value", call. = FALSE)
  }

  if (is.null(cutoff) || !is.numeric(cutoff) || !is.finite(cutoff)) {
    stop("'cutoff' value must be numeric", call. = FALSE)
  }

  if ((cutoff > 1) || (cutoff < 0)) {
    stop("'cutoff' value ", cutoff, ", outside the interval [0,1]",
         call. = FALSE)
  }

  ## radIn
  if (length(radIn) != 1) {
    stop("'radIn' must be a single value", call. = FALSE)
  }

  if (is.null(radIn) || !is.numeric(radIn) || !is.finite(radIn)) {
    stop("'radIn' value must be numeric", call. = FALSE)
  }

  if ((radIn > 1) || (radIn < 0)) {
    stop("'radIn' value ", radIn, ", outside the interval [0,1]",
         call. = FALSE)
  }

  ## overlap
  if (length(overlap) != 1 || !is.logical(overlap)) {
    stop("'overlap' must be a logical value", call. = FALSE)
  }

  ## ncols
  if (length(ncols) != 1) {
    stop("'ncols' must be of length 1", call. = FALSE)
  }

  if (is.null(ncols) || !is.numeric(ncols) || !is.finite(ncols) ||
      ncols < 1) {
    stop("invalid number of columns 'ncols' for facetting data tables",
         call. = FALSE)
  }

  ## legTitle
  legTitle <- as.graphicsAnnot(legTitle)

  ## graph
  if(!is.logical(graph)){
    stop("'graph' must be either TRUE or FALSE", call. = FALSE)
  }


  ## variables scatter plot ##################################


  coordRef <- lapply(dtab, function(x) { suppressWarnings({
    hetcor(cbind(RefData[[x]], comprConf))$correlations[1:ncol(RefData[[x]]),
                                                        (ncol(RefData[[x]])+1):
                                                          (ncol(RefData[[x]])+ncol(comprConf))]
  })
  })
  coordImp <- list()
  for (i in seq_along(imputations)) {
    coordImp[[i]] <- lapply(dtab, function(x) { suppressWarnings({
      hetcor(cbind(imputations[[i]][[x]], comprConf))$correlations[1:ncol(imputations[[i]][[x]]),
                                                                   (ncol(imputations[[i]][[x]])+1):
                                                                     (ncol(imputations[[i]][[x]])+ncol(comprConf))]
    })
    })
  }

  idx <- vector("list", nbTables)

  if (cutoff > 0) {
    for (i in seq_len(nbTables - 1)) {
      for (j in 2:nbTables) {
        mat <- coordRef[[i]] %*% t(coordRef[[j]])
        tmp <- (abs(mat) >= cutoff)

        if (all(!tmp)) {
          stop("'cutoff' must be < ", round(max(abs(mat)), 2),
               call. = FALSE)
        }

        tmp <- which(tmp, arr.ind = TRUE)
        idx[[i]] <- rbind(idx[[i]], tmp[, 1, drop = FALSE])
        idx[[j]] <- rbind(idx[[j]], tmp[, 2, drop = FALSE])
      }
    }

    idx <- lapply(idx, function(x) unique(x))
    namevar <- NULL

    for (i in seq_len(nbTables)) {
      coordRef[[i]] <- matrix(coordRef[[i]][idx[[i]], ], ncol = length(comp))
      for (j in 1:length(coordImp)){
        coordImp[[j]][[i]] <- matrix(coordImp[[j]][[i]][idx[[i]], ], ncol = length(comp))
      }
      namevar <- c(namevar, colnames(RefData[[i]])[idx[[i]]])
    }
  }

  ## data frame for ggplot
  combined_coord <- list()
  for(i in seq_along(coordImp)){
    combined_coord[[i]] <- do.call(rbind, coordImp[[i]])
  }
  combined_coord <- as.data.frame(do.call(rbind, combined_coord))
  combined_coord$source <- "Individual imputations"
  temp_df <- as.data.frame(do.call(rbind, coordRef))
  temp_df$source <- "Compromise"
  df <- rbind(temp_df, combined_coord)

  nbVar <- vapply(coordRef, nrow, 1L)
  #  df <- do.call(rbind, coordRef)

  #  if (cutoff == 0) { namevar <- rownames(df) }
  if (cutoff == 0) { namevar <- rownames(do.call(rbind, coordRef)) }
  rownames(df) <- NULL

  df <- data.frame(df, namevar, rep(rep(nmTables, nbVar),length(coordImp)+1))
  names(df) <- c("x", "y", "Source", "Var", "Tables")
  df$Tables <- factor(df$Tables, levels = nmTables)
  df$Source <- factor(df$Source, levels = c("Individual imputations","Compromise"))

  x <- y <- Tables <- Var <- NULL
  circle <- data.frame(x = cos(seq(0, 2 * pi, l = 1000)),
                       y = sin(seq(0, 2 * pi, l = 1000)))

  circleInt <- data.frame(x = radIn * cos(seq(0, 2 * pi, l = 1000)),
                          y = radIn * sin(seq(0, 2 * pi, l = 1000)))

  ## carried out ggplot
  g <- ggplot(df, aes(x = x, y = y)) + theme_bw() +
    scale_x_continuous(limits = c(-1, 1)) +
    scale_y_continuous(limits = c(-1, 1)) +
    geom_path(data = circle, aes(x = x, y = y), color = "black") +
    geom_path(data = circleInt, aes(x = x, y = y), color = "black") +
    geom_vline(aes(xintercept = 0), linetype = 2, colour = "gray40") +
    geom_hline(aes(yintercept = 0), linetype = 2, colour = "gray40")

  if (spty) { g <- g + coord_equal(ratio = 1) }

  for (i in seq_len(nbTables)) {
    id <- (df$Tables == nmTables[i])

    if (flagVarNames[i]) {
      g <- g + geom_text(data = df[id, ],
                         aes(x = x, y = y, colour = Var, label = Var),
                         size = cex[i], alpha = alpha,
                         show.legend = FALSE) +
        geom_point(data = df[id, ][1, ],
                   aes(x = 0, y = 0, colour = Var, shape = Source),
                   size = 0, alpha = alpha)
    } else {
      g <- g + geom_point(data = df[id, ],
                          aes(x = x, y = y, colour = Var, shape = Source),
                          size = cex[i], alpha = alpha)
    }
  }

  g <- g + scale_colour_manual(name = legTitle, values = col,
                               breaks = unique(namevar))  +
    scale_shape_manual(values = c("Compromise" = 1, "Individual imputations" = 3)) +
    #      guides(colour = guide_legend(override.aes = list(shape = pch)))+#,
    #                                                       size = cex))) +
    labs(x = paste0('Comp ', comp[1]), y = paste0('Comp ', comp[2])) +
    theme(legend.text = element_text(size = 12),
          legend.title = element_text(size = 12))

  if (!overlap) {
    g <- g + facet_wrap(~ Tables, ncol = ncols, as.table = TRUE) +
      theme(strip.text = element_text(size = 12))
  }

  if (graph) print(g)
  return(invisible(list(df = df, ggp = g)))
  #  return(list(id=id, nbTables=nbTables, nmTables=nmTables, df = df))
}
