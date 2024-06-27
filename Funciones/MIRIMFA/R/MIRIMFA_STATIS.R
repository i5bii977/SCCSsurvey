#' @export
#' @importFrom FactoMineR MFA
#' @importFrom missMDA imputeMFA
MIRIMFA_STATIS <- function(X,
                           folder_path,
                           prefix_imputations = "RIMFA_number_"){
  STATIS <- function (Ktab, nf, threshold) {

    nlig <- nrow(Ktab[[1]])
    ntab <- length(Ktab)
    result <- sep <- list()

    sep <- lapply(Ktab, function(wk) { wk <- as.matrix(wk)
    wk <- t(wk) * sqrt(1/ncol(wk))
    wk <- t(wk) %*% wk
    wk })

    sep <- matrix(unlist(sep), nlig * nlig, ntab)

    RV <- t(sep) %*% sep
    ak <- sqrt(diag(RV))
    RV <- sweep(RV, 1, ak, "/")
    RV <- sweep(RV, 2, ak, "/")

    eig1 <- eigen(RV, symmetric=TRUE)

    if (any(eig1$vectors[, 1] < 0)) {
      eig1$vectors[, 1] <- -eig1$vectors[, 1]
    }

    tabw <- eig1$vectors[, 1]

    sep <- t(t(sep)/ak)
    Cro <- rowSums(sweep(sep, 2, tabw, "*"))
    Cro <- matrix(unlist(Cro), nlig, nlig)
    eig1 <- eigen(Cro, symmetric=TRUE)
    result$Cro <- Cro
    rm(Cro)
    eig <- eig1$values
    rank <- sum((eig/eig[1]) > threshold)

    if (nf <= 0) { nf <- 2 }

    if (nf > rank) { nf <- rank }

    wref <- eig1$vectors[, seq_len(nf)]
    rm(eig1)
    w <- data.frame(t(t(wref) * sqrt(eig[seq_len(nf)])))
    names(w) <- paste("C", seq_len(nf), sep="")
    result$Cli <- w

    return(result)
  }


  create_total_disjunctive_table <- function(X, group, type) {
    tab_disj <- matrix(0, nrow(X), 0)  # Initialize total disjunctive table

    for (g in 1:length(group)) {
      # Extract subset of columns for the current group
      if (g == 1)
        aux_base <- X[, 1:group[1], drop = FALSE]
      else
        aux_base <- X[, (cumsum(group)[g - 1] + 1):cumsum(group)[g], drop = FALSE]

      # Construct disjunctive table based on variable type
      if (type[g] == "n") {
        # For qualitative variables
        tab_disj_group <- tab.disjonctif.prop(aux_base)
      } else {
        # For other variable types (s or c), directly include them in the table
        tab_disj_group <- as.matrix(aux_base)
      }

      # Append the disjunctive table of the current group to the total table
      tab_disj <- cbind(tab_disj, tab_disj_group)
    }

    return(tab_disj)
  }
  ## IMPUTATION ITSELF #################
  imputeMissingMFA <- function(incompData, U, comp,
                               maxiter, threshold) {
    id.na <- is.na(incompData)
    real.na <- matrix(FALSE, nrow = nrow(incompData), ncol = ncol(incompData))

    # Identify missing values for dummy variables
    for (i in 1:ncol(incompData)) {
      if (grepl("\\.\\d+", colnames(incompData)[i])) {  # Check if variable name indicates it's a dummy variable
        # For dummy variables, values other than 0 or 1 are considered missing
        real.na[, i] <- incompData[, i] != 0 & incompData[, i] != 1
      } else {
        # For non-dummy variables, use is.na to identify missing values
        real.na[, i] <- is.na(incompData[, i])
      }
    }

    colnames(real.na) <- colnames(incompData)
    rownames(real.na) <- rownames(incompData)

    Xold <- incompData
    Xold[id.na] <- 0
    U <- as.matrix(U[, comp])
    Vold <- matrix(999, nrow=ncol(incompData), ncol=length(comp))

    iter <- 0

    repeat {
      Xold <- scale(Xold)
      center <- attr(Xold, "scaled:center")
      sigma <- attr(Xold, "scaled:scale")

      V <- t(solve(crossprod(U)) %*% t(U) %*% Xold)
      Xnew <- U %*% t(V)
      Xnew <- sweep(Xnew, 2, sigma, "*")
      Xnew <- sweep(Xnew, 2, center, "+")
      Xnew[!real.na] <- incompData[!real.na]

      if (sqrt(sum((V - Vold)^2)) < threshold) break

      Xold <- Xnew

      if (iter == maxiter) {
        warning(paste("maximum number of iterations reached in data imputation with",
                      sqrt(sum((V - Vold)^2)),
                      "/",
                      threshold,
                      "=",
                      round(sqrt(sum((V - Vold)^2))/threshold,2)),
                call.=FALSE)
        break
      }

      Vold <- V
      iter <- iter + 1
    }

    imputed <- Xnew[real.na]
    original <- incompData[!real.na]

    return(list(imputed = imputed, original = original, total = Xnew))

  }
  ## FIND CATEGORIES #######################################
  find.category <- function(X, tabdisj) {
    nbdummy <- rep(1, ncol(X))
    is.quali <- which(!unlist(lapply(X, is.numeric)))
    nbdummy[is.quali] <- unlist(lapply(X[, is.quali,
                                         drop = FALSE], nlevels))
    vec = c(0, cumsum(nbdummy))
    Xres <- X
    for (i in 1:ncol(X)) {
      if (i %in% is.quali)
        Xres[, i] <- as.factor(levels(X[, i])[apply(tabdisj[,
                                                            (vec[i] + 1):vec[i + 1]],
                                                    1,
                                                    which.max)])
      else Xres[, i] <- tabdisj[, vec[i] + 1]
    }
    return(Xres)
  }

  # EXTRACTING FILES #############################################

  all_files <- list.files(path = folder_path, full.names = TRUE)

  # CHECKING ######################################################
  ## ncp ##########################################################

  # Load ncps from MFAs
  ncps <- lapply(all_files, function(file) {
    res <- readRDS(file)
    return(res$call$ncp)
  })

  # Check if all ncps are the same
  if (length(unique(ncps)) > 1) {
    stop("The ncp values obtained from the MFA files are not consistent.")
  }
  ncp <- as.numeric(unique(ncps))

  ## Groups ############################################################

  # Load groups from MFAs
  groups <- lapply(all_files, function(file) {
    res <- readRDS(file)
    return(res$call$group)
  })

  # Check if all groups are the same
  if (length(unique(groups)) > 1) {
    stop("The groups obtained from the MFA files are not consistent.")
  }
  group <- groups[[1]]

  ## Types ############################################################

  # Load types from MFAs
  types <- lapply(all_files, function(file) {
    res <- readRDS(file)
    return(res$call$type)
  })

  # Check if all groups are the same
  if (length(unique(types)) > 1) {
    stop("The types obtained from the MFA files are not consistent.")
  }
  type <- types[[1]]

  ## Row weights ####################################################

  # Load types from MFAs
  row.ws <- lapply(all_files, function(file) {
    res <- readRDS(file)
    return(res$call$row.w)
  })

  # Check if all groups are the same
  if (length(unique(row.ws)) > 1) {
    stop("The row weight values obtained from the MFA files are not consistent.")
  }
  row.w <- row.ws[[1]]

  ## Row weights ####################################################

  # Load types from MFAs
  coeff.ridges <- lapply(all_files, function(file) {
    res <- readRDS(file)
    return(res$call$coeff.ridge)
  })

  # Check if all groups are the same
  if (length(unique(coeff.ridges)) > 1) {
    stop("The regularization coefficient values obtained from the MFA files are not consistent.")
  }
  coeff.ridge <- coeff.ridges[[1]]

  ## Thresholds ####################################################

  # Load types from MFAs
  thresholds <- lapply(all_files, function(file) {
    res <- readRDS(file)
    return(res$call$threshold)
  })

  # Check if all groups are the same
  if (length(unique(thresholds)) > 1) {
    stop("The threshold values obtained from the MFA files are not consistent.")
  }
  threshold <- thresholds[[1]]

  ## Supplementary individuals #######################################

  # Load types from MFAs
  ind.sups <- lapply(all_files, function(file) {
    res <- readRDS(file)
    return(res$call$ind.sup)
  })

  # Check if all groups are the same
  if (length(unique(ind.sups)) > 1) {
    stop("The indexes of the supplemmentary individuals obtained from the MFA files are not consistent.")
  }
  ind.sup <- ind.sups[[1]]

  ## Supplementary variables ######################################

  # Load types from MFAs
  num.group.sups <- lapply(all_files, function(file) {
    res <- readRDS(file)
    return(res$call$num.group.sup)
  })

  # Check if all groups are the same
  if (length(unique(num.group.sups)) > 1) {
    stop("The indexes of the supplemmentary groups of variables obtained from the MFA files are not consistent.")
  }

  num.group.sup <- num.group.sups[[1]]

  ## Maximum iterations ##########################################

  # Load types from MFAs
  maxiters <- lapply(all_files, function(file) {
    res <- readRDS(file)
    return(res$call$maxiter)
  })

  # Check if all groups are the same
  if (length(unique(maxiters)) > 1) {
    stop("The maximum number of iterations obtained from the MFA files are not consistent.")
  }

  maxiter <- maxiters[[1]]

  ## Name groups ######################################

  # Load types from MFAs
  name.groups <- lapply(all_files, function(file) {
    res <- readRDS(file)
    return(res$call$name.group)
  })

  # Check if all groups are the same
  if (length(unique(name.groups)) > 1) {
    stop("The name of the groups obtained from the MFA files are not consistent.")
  }

  name.group <- name.groups[[1]]

  # Assemble arguments ##########################################

  # Load imputations
  imputations <- lapply(all_files, function(file) {
    res <- readRDS(file)
    return(res$imputation$completeObs)
  })

  # Load tab.disj
  tab.disj <- lapply(all_files, function(file) {
    res <- readRDS(file)
    return(res$imputation$tab.disj)
  })

  # Load U from MFAs
  U <- lapply(all_files, function(file) {
    res <- readRDS(file)
    return(res$MFA)
  })



  # Calculate the compromise space (STATIS method)
  tmp <- STATIS(U, nf = ncp, threshold = threshold)

  ##- data imputation --------------------------------------#
  ##--------------------------------------------------------#
  incompData <- create_total_disjunctive_table(X,group,type)

  impD <- imputeMissingMFA(incompData,
                           tmp$Cli,
                           comp=seq_len(ncp),
                           maxiter=maxiter,
                           threshold=threshold)
  result <- find.category(X,impD$total)


  # Determine the dimensions of your data frames
  num_rows <- nrow(tab.disj[[1]])
  num_cols <- ncol(tab.disj[[1]])
  num_slices <- length(tab.disj)

  # Create an empty three-dimensional array
  array_3d <- array(NA, dim = c(num_rows, num_cols, num_slices))

  # Populate the array with the values from each data frame
  for (i in 1:num_slices) {
    array_3d[,,i] <- as.matrix(tab.disj[[i]])
  }
  tab.disj <- array_3d

  object <- list(res.MI = imputations,
                 res.imputeMFA = result,
                 call = list(X = X,
                             group = group,
                             ncp = ncp,
                             type = type,
                             M = length(all_files),
                             row.w = row.w,
                             coeff.ridge = coeff.ridge,
                             threshold = threshold,
                             ind.sup = ind.sup,
                             num.group.sup = num.group.sup,
                             maxiter = maxiter,
                             name.group = name.group,
                             tab.disj = tab.disj),
                 compromise = tmp$Cli,
                 res.tab.disj = impD$total)
  class(object) <- c("MIRIMFA", "list")
  return(object)
  #  return(threshold)

  #  return(list(incompData =incompData, Cli = tmp$Cli, comp=seq_len(ncp), maxiter=maxiter, threshold=threshold))
}
