#' @export
#' @importFrom FactoMineR MFA
#' @importFrom missMDA imputeMFA

MIRIMFAextraimp <- function (X,
                             group,
                             ncp = 2,
                             type = rep("s", length(group)),
                             manual_threshold = 1,
                             row.w = NULL,
                             coeff.ridge = 1,
                             threshold = 1e-06,
                             ind.sup = NULL,
                             num.group.sup = NULL,
                             maxiter = 1000,
                             folder_path = here("imputationsRIMFA"),
                             name.group = NULL,
                             start = 0,
                             initial_seed = 1) {
  ## Imputation while FactoMineR still has bugs ############

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

  # Checking variables #############################
  X <- as.data.frame(X)
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
  successful_imputations <- start
  current_seed <- initial_seed
  while (successful_imputations < manual_threshold) {
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
      mfa <- MFA(result$completeObs,
                 group,
                 type = type,
                 ncp = ncp,
                 name.group = name.group,
                 graph = F)
      call <- list(X = X,
                   group = group,
                   ncp = ncp,
                   type = type,
                   method = "Regularized",
                   row.w = row.w,
                   coeff.ridge = coeff.ridge,
                   threshold = threshold,
                   ind.sup = ind.sup,
                   num.group.sup = num.group.sup,
                   seed = seed,
                   maxiter = maxiter,
                   name.group = name.group)
      RIMFA <- list(imputation = result, MFA = mfa, call = call)
      successful_imputations <- successful_imputations + 1
      if (!file.exists(folder_path)) {
        dir.create(folder_path, recursive = TRUE)
      }
      saveRDS(RIMFA,
              file.path(here(folder_path,
                             paste0("RIMFA_number_",
                                    successful_imputations,".rds"))))
    }, error = function(e) {
      cat("Imputation failed for seed:", seed, "- Error:", conditionMessage(e), "\n")
      # Increment current_seed to try a new seed in the next iteration
    })
    current_seed <- current_seed + 1
  }
  return(list(final_seed = current_seed, final_index = successful_imputations))
}
