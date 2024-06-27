#' @export
#' @importFrom FactoMineR MFA
#' @importFrom missMDA imputeMFA
#' @import FactoMineR
#' @import missMDA

MIRIMFA <- function (X,
                     group,
                     ncp = 2,
                     type = rep("s", length(group)),
                     M = 5,
                     row.w = NULL,
                     coeff.ridge = 1,
                     threshold = 1e-06,
                     ind.sup = NULL,
                     num.group.sup = NULL,
                     maxiter = 1000,
#                     folder_path = here("imputationsRMFA")),
                     name.group = NULL)
{
  # Define internal functions ##############################
  ## STATIS for the compromise space  ######################
  STATIS <- function (Ktab, nf, threshold=1e-07) {

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
  ## TOTAL DISJUNCTIVE TABLE  ##############################
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
  ## IMPUTATION WHILE FactoMineR still has bugs ############
  imputeMFAseedSolved <- function (X,
                                   group,
                                   ncp = 2,
                                   type = rep("s", length(group)),
                                   method = c("Regularized", "EM"),
                                   row.w = NULL,
                                   coeff.ridge = 1,
                                   threshold = 1e-06,
                                   ind.sup = NULL,
                                   num.group.sup = NULL,
                                   seed = NULL,
                                   maxiter = 1000, ...)
  {
    impute <- function(X, group, ncp = 2, type = rep("s", length(group)),
                       method = NULL, threshold = 1e-06, ind.sup = NULL, num.group.sup = NULL,
                       seed = NULL, maxiter = 1000, row.w = NULL, coeff.ridge = 1,
                       ...) {
      moy.p <- function(V, poids) {
        res <- sum(V * poids, na.rm = TRUE)/sum(poids[!is.na(V)])
      }
      ec <- function(V, poids) {
        res <- sqrt(sum(V^2 * poids, na.rm = TRUE)/sum(poids[!is.na(V)]))
      }
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
                                                                (vec[i] + 1):vec[i + 1]], 1, which.max)])
          else Xres[, i] <- tabdisj[, vec[i] + 1]
        }
        return(Xres)
      }
      method <- match.arg(method, c("Regularized", "regularized",
                                    "EM", "em"), several.ok = T)[1]
      method <- tolower(method)
      if (!is.null(seed))
        set.seed(seed)
      X <- as.data.frame(X)
      X <- droplevels(X)
      if ("n" %in% type) {
        niveau = NULL
        for (j in 1:ncol(X)) {
          if (!is.numeric(X[, j]))
            niveau = c(niveau, levels(X[, j]))
        }
        for (j in 1:ncol(X)) {
          if (!is.numeric(X[, j])) {
            if (sum(niveau %in% levels(X[, j])) != nlevels(X[,
                                                             j]))
              levels(X[, j]) = paste(colnames(X)[j], levels(X[,
                                                              j]), sep = "_")
          }
        }
      }
      group.mod = group
      Xhat <- matrix(0, nrow(X), 0)
      Xhat2 <- matrix(0, nrow(X), 0)
      MM <- vector(mode = "list", length = length(group))
      ET <- vector(mode = "list", length = length(group))
      tab.disj.comp <- vector(mode = "list", length = length(group))
      ponderation <- rep(1, length(group))
      if (ncp == 0) {
        result <- list()
        if (sum(unlist(lapply(X, is.numeric))) > 0)
          ind.quanti <- which(unlist(lapply(X, is.numeric)))
        else ind.quanti <- NULL
        if (sum(lapply(X, class) == "factor") > 0)
          ind.quali <- which((lapply(X, class)) == "factor")
        else ind.quali <- NULL
        result$completeObs <- X
        if (!is.null(ind.quali))
          result$completeObs[, ind.quali] <- find.category(X[,
                                                             ind.quali, drop = F], tab.disjonctif.prop(X[,
                                                                                                         ind.quali, drop = F], row.w = row.w))
        if (!is.null(ind.quanti)) {
          tab.disj <- X[, ind.quanti, drop = F]
          Moy <- matrix(colMeans(tab.disj, na.rm = T),
                        nrow = nrow(tab.disj), ncol = ncol(tab.disj),
                        byrow = T)
          tab.disj[is.na(tab.disj)] <- Moy[is.na(tab.disj)]
          result$completeObs[, ind.quanti] <- tab.disj
        }
        nbdummy <- rep(1, ncol(X))
        is.quali <- which(!unlist(lapply(X, is.numeric)))
        if (length(is.quali) != 0) {
          nbdummy[is.quali] <- unlist(lapply(X[, is.quali,
                                               drop = FALSE], nlevels))
          tabdisj <- matrix(NA, nrow(X), ncol = sum(nbdummy))
          tabdisj[, cumsum(nbdummy)[which(nbdummy == 1)]] <- as.matrix(result$completeObs[,
                                                                                          ind.quanti, drop = F])
          auxQuali <- tab.disjonctif.prop(X[, ind.quali,
                                            drop = F], row.w = row.w)
          tabdisj[, -cumsum(nbdummy)[which(nbdummy == 1)]] <- auxQuali
          rownames(tabdisj) <- rownames(X)
          colnames(tabdisj) <- paste0("v", 1:ncol(tabdisj))
          colnames(tabdisj)[cumsum(nbdummy)[which(nbdummy ==
                                                    1)]] <- colnames(result$completeObs[, ind.quanti,
                                                                                        drop = F])
          colnames(tabdisj)[-cumsum(nbdummy)[which(nbdummy ==
                                                     1)]] <- colnames(auxQuali)
          result$tab.disj <- tabdisj
        }
        for (g in 1:length(group)) {
          if (g == 1)
            aux.base <- X[, 1:group[1], drop = FALSE]
          else aux.base <- X[, (cumsum(group)[g - 1] +
                                  1):cumsum(group)[g], drop = FALSE]
          if (type[g] == "n") {
            tab.disj = tab.disjonctif.prop(aux.base, seed,
                                           row.w = row.w)
            tab.disj.comp[[g]] = tab.disj
            group.mod[g] <- ncol(tab.disj)
          }
        }
        ind.var <- vector(mode = "list", length = length(group))
        ind.var[[1]] <- 1:group.mod[1]
        for (g in 2:length(group)) ind.var[[g]] <- (cumsum(group.mod)[g -
                                                                        1] + 1):cumsum(group.mod)[g]
        result$call$group.mod <- group.mod
        ind.var <- vector(mode = "list", length = length(group))
        ind.var[[1]] <- 1:result$call$group.mod[1]
        for (g in 2:length(group)) ind.var[[g]] <- (cumsum(result$call$group.mod)[g -
                                                                                    1] + 1):cumsum(result$call$group.mod)[g]
        result$call$ind.var = ind.var
        return(result)
      }
      for (g in 1:length(group)) {
        if (g == 1)
          aux.base <- X[, 1:group[1], drop = FALSE]
        else aux.base <- X[, (cumsum(group)[g - 1] + 1):cumsum(group)[g],
                           drop = FALSE]
        if (type[g] == "s") {
          Xhat2 <- cbind.data.frame(Xhat2, aux.base)
          MM[[g]] <- apply(as.data.frame(aux.base), 2,
                           moy.p, row.w)
          aux.base <- t(t(as.matrix(aux.base)) - MM[[g]])
          ET[[g]] <- apply(as.data.frame(aux.base), 2,
                           ec, row.w)
          aux.base <- t(t(as.matrix(aux.base))/ET[[g]])
          missing <- which(is.na(as.matrix(aux.base)))
          if (any(is.na(aux.base)))
            aux.base[missing] <- 0
          ponderation[g] <- FactoMineR::svd.triplet(aux.base,
                                                    ncp = 1, row.w = row.w)$vs[1]
          Xhat <- cbind.data.frame(Xhat, aux.base/ponderation[g])
          if (!is.null(seed) & (length(missing) != 0)){
            Xhat <- as.matrix(Xhat)
            Xhat[missing] <- rnorm(length(missing))
            Xhat <- as.data.frame(Xhat)}
        }
        if (type[g] == "c") {
          Xhat2 <- cbind.data.frame(Xhat2, aux.base)
          MM[[g]] <- apply(as.data.frame(aux.base), 2,
                           moy.p, row.w)
          aux.base <- t(t(as.matrix(aux.base)) - MM[[g]])
          missing <- which(is.na(as.matrix(aux.base)))
          if (any(is.na(aux.base)))
            aux.base[missing] <- 0
          ponderation[g] = FactoMineR::svd.triplet(aux.base,
                                                   ncp = 1, row.w = row.w)$vs[1]
          Xhat <- cbind.data.frame(Xhat, aux.base/ponderation[g])
          if (!is.null(seed) & (length(missing) != 0)){
            Xhat <- as.matrix(Xhat)
            Xhat[missing] <- rnorm(length(missing))
            Xhat <- as.data.frame(Xhat)}
        }
        if (type[g] == "n") {
          tab.disj = tab.disjonctif.prop(aux.base, seed,
                                         row.w = row.w)
          tab.disj.comp[[g]] = tab.disj
          group.mod[g] <- ncol(tab.disj)
          MM[[g]] = apply(tab.disj, 2, moy.p, row.w)/ncol(aux.base)
          Z = t(t(tab.disj)/apply(tab.disj, 2, moy.p, row.w))
          Z = t(t(Z) - apply(Z, 2, moy.p, row.w))
          Zscale = t(t(Z) * sqrt(MM[[g]]))
          ponderation[g] <- FactoMineR::svd.triplet(Zscale,
                                                    row.w = row.w)$vs[1]
          Xhat <- cbind.data.frame(Xhat, Zscale/ponderation[g])
          Xhat2 <- cbind.data.frame(Xhat2, as.data.frame(tab.disjonctif(aux.base)))
        }
      }
      ind.var <- vector(mode = "list", length = length(group))
      ind.var[[1]] <- 1:group.mod[1]
      for (g in 2:length(group)) ind.var[[g]] <- (cumsum(group.mod)[g -
                                                                      1] + 1):cumsum(group.mod)[g]
      fittedX <- Xhat <- as.matrix(Xhat)
      if (ncp >= min(nrow(Xhat) - 2, ncol(Xhat) - 1))
        stop("ncp is too large")
      ncp <- min(ncp, ncol(X) - 1, nrow(X) - 2)
      missing <- which(is.na(as.matrix(Xhat2)))
      nb.iter <- 1
      old <- Inf
      nrX <- nrow(Xhat)
      ncX <- sum(group.mod) - sum(group[type == "n"])
      if (length(num.group.sup) > 0)
        ncX <- ncX - (sum(group.mod[num.group.sup]) - sum(group[num.group.sup][type[num.group.sup] ==
                                                                                 "n"]))
      while (nb.iter > 0) {
        for (g in 1:length(group)) {
          if (g == 1)
            aux.base <- Xhat[, 1:group.mod[1], drop = FALSE]
          else aux.base <- Xhat[, (cumsum(group.mod)[g -
                                                       1] + 1):cumsum(group.mod)[g], drop = FALSE]
          aux.base <- aux.base * ponderation[g]
          if (type[g] == "s") {
            aux.base <- t((t(aux.base) * ET[[g]]) + MM[[g]])
            MM[[g]] <- apply(aux.base, 2, moy.p, row.w)
            aux.base <- t(t(aux.base) - MM[[g]])
            ET[[g]] <- apply(aux.base, 2, ec, row.w)
            aux.base <- t(t(aux.base)/ET[[g]])
            ponderation[g] = FactoMineR::svd.triplet(aux.base,
                                                     ncp = 1, row.w = row.w)$vs[1]
          }
          if (type[g] == "c") {
            aux.base <- t(t(aux.base) + MM[[g]])
            MM[[g]] <- apply(aux.base, 2, moy.p, row.w)
            aux.base <- t(t(aux.base) - MM[[g]])
            ponderation[g] = FactoMineR::svd.triplet(aux.base,
                                                     ncp = 1, row.w = row.w)$vs[1]
          }
          if (type[g] == "n") {
            tab.disj = t(t(aux.base)/sqrt(MM[[g]])) + matrix(1,
                                                             nrow(aux.base), ncol(aux.base))
            tab.disj = t(t(tab.disj) * apply(tab.disj.comp[[g]],
                                             2, moy.p, row.w))
            tab.disj.comp[[g]] = tab.disj
            MM[[g]] = apply(tab.disj, 2, moy.p, row.w)/ncol(aux.base)
            if (any(MM[[g]] < 0)) {
              stop(paste("The algorithm fails to converge. Choose a number of components (ncp) less or equal than ",
                         ncp - 1, " or a number of iterations (maxiter) less or equal than ",
                         maxiter - 1, sep = ""))
            }
            Z = t(t(tab.disj)/apply(tab.disj, 2, moy.p,
                                    row.w))
            Z = t(t(Z) - apply(Z, 2, moy.p, row.w))
            aux.base = t(t(Z) * sqrt(MM[[g]]))
            ponderation[g] <- FactoMineR::svd.triplet(aux.base,
                                                      row.w = row.w, ncp = 1)$vs[1]
          }
          if (g == 1)
            Xhat[, 1:group.mod[1]] <- aux.base/ponderation[g]
          else Xhat[, (cumsum(group.mod)[g - 1] + 1):cumsum(group.mod)[g]] <- aux.base/ponderation[g]
        }
        if (!is.null(num.group.sup)) {
          for (g in num.group.sup) {
            if (g == 1)
              Xhat[, 1:group.mod[1]] <- Xhat[, 1:group.mod[1]] *
                1e-08
            else Xhat[, (cumsum(group.mod)[g - 1] + 1):cumsum(group.mod)[g]] <- Xhat[,
                                                                                     (cumsum(group.mod)[g - 1] + 1):cumsum(group.mod)[g]] *
                1e-08
          }
        }
        svd.res <- FactoMineR::svd.triplet(Xhat, row.w = row.w,
                                           ncp = ncp)
        if (length(num.group.sup) > 0)
          sigma2 <- nrX * ncX/min(ncX, nrX - 1) * sum((svd.res$vs[-c(1:ncp)]^2)/((nrX -
                                                                                    1) * ncX - (nrX - 1) * ncp - ncX * ncp + ncp^2))
        else sigma2 <- nrX * ncX/min(ncX, nrX - 1) * sum((svd.res$vs[-c(1:ncp)]^2)/((nrX -
                                                                                       1) * ncX - (nrX - 1) * ncp - ncX * ncp + ncp^2))
        sigma2 <- min(sigma2 * coeff.ridge, svd.res$vs[ncp +
                                                         1]^2)
        if (method == "em")
          sigma2 <- 0
        lambda.shrinked = (svd.res$vs[1:ncp]^2 - sigma2)/svd.res$vs[1:ncp]
        if (ncp == 1)
          fittedX = tcrossprod((svd.res$U[, 1, drop = FALSE] *
                                  row.w) * lambda.shrinked, svd.res$V[, 1, drop = FALSE])
        else fittedX = tcrossprod(t(t(svd.res$U[, 1:ncp] *
                                        row.w) * lambda.shrinked), svd.res$V[, 1:ncp])
        fittedX <- fittedX/row.w
        diff <- Xhat - fittedX
        diff[missing] <- 0
        objective <- sum(diff^2 * row.w)
        criterion <- abs(1 - objective/old)
        old <- objective
        nb.iter <- nb.iter + 1
        Xhat[missing] <- fittedX[missing]
        if (!is.null(num.group.sup)) {
          for (g in num.group.sup) {
            if (g == 1)
              Xhat[, 1:group.mod[1]] <- Xhat[, 1:group.mod[1]] *
                1e+08
            else Xhat[, (cumsum(group.mod)[g - 1] + 1):cumsum(group.mod)[g]] <- Xhat[,
                                                                                     (cumsum(group.mod)[g - 1] + 1):cumsum(group.mod)[g]] *
                1e+08
          }
        }
        if (!is.nan(criterion)) {
          if ((criterion < threshold) && (nb.iter > 5))
            nb.iter <- 0
          if ((objective < threshold) && (nb.iter > 5))
            nb.iter <- 0
        }
        if (nb.iter > maxiter) {
          nb.iter <- 0
          warning(paste("Stopped after ", maxiter, " iterations"))
        }
      }
      for (g in 1:length(group)) {
        if (g == 1)
          aux.base <- Xhat[, 1:group.mod[1], drop = FALSE]
        else aux.base <- Xhat[, (cumsum(group.mod)[g - 1] +
                                   1):cumsum(group.mod)[g], drop = FALSE]
        aux.base <- aux.base * ponderation[g]
        if (type[g] == "s") {
          aux.base <- t(t(aux.base) * ET[[g]])
          aux.base <- t(t(aux.base) + MM[[g]])
        }
        if (type[g] == "c")
          aux.base <- sweep(aux.base, 2, MM[[g]], FUN = "+")
        if (type[g] == "n") {
          tab.disj = t(t(aux.base)/sqrt(MM[[g]])) + matrix(1,
                                                           nrow(aux.base), ncol(aux.base))
          aux.base = t(t(tab.disj) * apply(tab.disj.comp[[g]],
                                           2, moy.p, row.w))
        }
        if (g == 1)
          Xhat[, 1:group.mod[1]] <- aux.base
        else Xhat[, (cumsum(group.mod)[g - 1] + 1):cumsum(group.mod)[g]] <- aux.base
      }
      completeObs <- as.matrix(Xhat2)
      completeObs[missing] <- Xhat[missing]
      result <- list()
      result$tab.disj <- completeObs
      result$completeObs <- find.category(X, completeObs)
      result$call$group.mod = group.mod
      result$call$ind.var = ind.var
      return(result)
    }
    obj = Inf
    method <- tolower(method)
    if (is.null(row.w))
      row.w = rep(1, nrow(X))/nrow(X)
    if (length(ind.sup) > 0)
      row.w[ind.sup] <- row.w[ind.sup] * 1e-08
    if (!any(is.na(X)))
      stop("no missing values in X, this function is not useful. Perform MFA on X.")
    res.impute <- impute(X, group = group, ncp = ncp, type = type,
                         method = method, threshold = threshold, seed = seed,
                         maxiter = maxiter, row.w = row.w, ind.sup = ind.sup,
                         num.group.sup = num.group.sup, coeff.ridge = coeff.ridge)
    return(res.impute)
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
                                                            (vec[i] + 1):vec[i + 1]], 1, which.max)])
      else Xres[, i] <- tabdisj[, vec[i] + 1]
    }
    return(Xres)
  }

  # Checking variables #####################################
  X <- as.data.frame(X)
  if(M < 2 | is.null(M) | !is.numeric(M)){
    stop("'M' must be a number equal or larger than two.")
  } else {
    M <- round(M)
  }
  is.quali <- which(!unlist(lapply(X, is.numeric)))
  X[, is.quali] <- lapply(X[, is.quali, drop = FALSE], as.factor)
  X <- droplevels(X)
  OneCat <- sapply(X, nlevels) == 1
  if (sum(OneCat) > 0) {
    if (sum(OneCat) == 1)
      warning("The following variable is constant and has been suppressed form the analysis: ",
              names(which(OneCat > 0)))
    if (sum(OneCat) > 1)
      warning("The following variable are constant and have been suppressed form the analysis: ",
              paste(names(which(OneCat > 0)), collapse = ", "))
    X <- X[, -which(OneCat)]
    if (ncol(X) <= 1)
      stop("No sufficient variables have 2 categories or more")
  }
  # Perform multiple imputations
  U <- list()
  imputations <- list()
  tab.disj <- list()
  successful_imputations <- 0
  current_seed <- 1
  while (successful_imputations < M) {
    seed <- current_seed
    tryCatch({
      result <- imputeMFAseedSolved(X,
                                    group,
                                    ncp = ncp,
                                    type = type,
                                    method = "Regularized",
                                    row.w = row.w,
                                    coeff.ridge = coeff.ridge,
                                    threshold = threshold,
                                    ind.sup = ind.sup,
                                    num.group.sup = num.group.sup,
                                    seed = seed,
                                    maxiter = maxiter)
      mfa <- FactoMineR::MFA(result$completeObs,
                 group,
                 type = type,
                 ncp = ncp,
                 name.group = name.group,
                 graph = F)
      tab.disj[[successful_imputations + 1]] <- result$tab.disj
      imputations[[successful_imputations + 1]] <- result$completeObs
      U[[successful_imputations + 1]] <- mfa[["ind"]][["coord"]]
      successful_imputations <- successful_imputations + 1
      #      if (!file.exists(folder_path)) {
      #        dir.create(folder_path, recursive = TRUE)
      #      }
      #      saveRDS(result,
      #              here(folder_path,
      #                   paste0("imputationRIMFA_number_",
      #                          successful_imputations,".rds")))
      #      saveRDS(mfa,
      #              here(folder_path,
      #                   paste0("RIMFA_number_",
      #                          successful_imputations,".rds")))
    }, error = function(e) {
      cat("Imputation failed for seed:", seed, "- Error:", conditionMessage(e), "\n")
      # Increment current_seed to try a new seed in the next iteration
    })
    current_seed <- current_seed + 1
  }

  ##- calculation of the compromise space (STATIS method) --#
  ##--------------------------------------------------------#
  tmp <- STATIS(U, nf=ncp)
  colnames(tmp$Cli) <- paste0("comp ", seq_len(ncp))
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
                             M = M,
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
}
