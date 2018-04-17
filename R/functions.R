####### 1. NBLDA FUNCTIONS #######

#' @title Control parameters for trained NBLDA model.
#'
#' @description Define control parameters to be used within \code{\link{trainNBLDA}} function.
#'
#' @param folds A positive integer. The number of folds for k-fold model validation.
#' @param repeats A positive integer. This is the number of repeats for k-fold model validation. If NULL, 0 or negative, it is set to 1.
#' @param foldIdx a list with indices of hold-out samples for each fold. It should be a list where folds are nested within repeats.
#' If NULL, \code{folds} and \code{repeats} are used to define hold-out samples.
#' @param rhos A vector of tuning parameters that control the amount of soft thresholding performed. If NULL, it is automatically generated within
#' \code{\link{trainNBLDA}} using \code{tuneLength}, i.e. the length of grid search. See details.
#' @param beta A smoothing term. A Gamma(beta,beta) prior is used to fit the Poisson model. Recommendation is to just leave it at 1, the default value.
#' See Witten (2011) and Dong et. al. (2016) for details.
#' @param prior A vector of length equal to the number of classes indicating the prior class probabilities. If NULL, all classes are assumed to be
#' equally distributed.
#' @param transform a logical. If TRUE, count data is transformed using power transformation. If \code{alpha} is not specified the power transformation
#' parameter is automatically calculated using goodness-of-fit test. See Witten (2011) for details.
#' @param apha a numeric value within [0, 1] to be used for power transformation.
#' @param truephi a vector of length equal to the number of variables representing the true overdispersion parameters for each variable. If a single
#' value is given, it is replicated for all variables. If a vector of length unequal to the number of variables is given, the first element of this
#' vector is used and replicated for all variables. If NULL, estimated overdispersions are used in the classifier. See details.
#' @param target BURADA KALDIM...
#' @param phi.epsilon a positive value for controlling the number of features whose dispersions are shrinked towards 0. See details.
#'
#' @return give information about returned objects here.
#'
#' @details \code{rhos} is used to control the level of sparsity, i.e. the number of variables (or features) used in classifier. If a variable has
#' no contribution to discrimination function, it should be removed from the model. By setting rhos within the interval [0, Inf], it is possible
#' control the amount of variables that is removed from the model. As the upper bound of rhos decreases towards 0, fewer variables are removed.
#' If rhos = 0, all variables are included in classifier.
#'
#' \code{truephi} controls how Poisson model differs from Negative Binomial model. If overdispersion is zero, Megative Binomial model converges to
#' Poisson model. Hence, the results from \code{\link{trainNBLDA}} is identical to PLDA results from \code{\link[PoiClaClu]{Classify}} when truephi = 0.
#'
#' \code{phi.epsilon} is a value used to shring estimated overdispersions towards 0. Poisson model assumes that there is no overdispersion in the
#' observed counts. However, this is not a valid assumption in highly overdispersed count data. \code{NBLDA} performs a shrinkage on estimated
#' overdispersions. Although the amount of shrinkage is dependent on several parameters such as \code{delta}, \code{target} and \code{truephi}, some
#' of the shrinked overdispersions might be very close to 0. By defining a threshold value for shrinked overdispersions, it is possible to shrink
#' very small overdispersions towards 0. If estiamted overdispersion is below \code{phi.epsilon}, it is shrinked to 0.
#'
#' @author Turcosa Developing Team
#'
#' @references Witten, DM (2011). Classification and clustering of sequencing data using a Poisson model.
#' Ann. Appl. Stat. 5(4), 2493--2518. doi:10.1214/11-AOAS493.
#'
#' Dong, K., Zhao, H., Tong, T., & Wan, X. (2016). NBLDA: negative binomial linear discriminant analysis for RNA-Seq data.
#' BMC Bioinformatics, 17(1), 369. http://doi.org/10.1186/s12859-016-1208-1.
#'
#' @keywords keywords_here
#'
#' @seealso <add seealso info here>
#'
#' @examples
#' # Enter examples here.
#' 1L
#'
#' @name control-NBLDA
#' @rdname control
#'
#' @export
control <- function(folds = 5, repeats = 2, foldIdx = NULL, rhos = NULL, beta = 1,
                    prior = NULL, transform = FALSE, alpha = NULL, truephi = NULL, target = 0,
                    phi.epsilon = 0.15, normalize.target = FALSE, delta = NULL, return.selected.gene.names = FALSE,
                    multicore = FALSE, ...){

  if (repeats <= 0 | is.null(repeats)){
    repeats <- 1
  }

  list(folds = folds, repeats = repeats, foldIdx = foldIdx, rhos = rhos, phi.epsilon = phi.epsilon, beta = beta,
       prior = prior, transform = transform, alpha = alpha, truephi = truephi,
       target = target, normalize.target = normalize.target, delta = delta,
       return.selected.gene.names = return.selected.gene.names, multicore = multicore, ...)
}


nblda <- function(x, y, xte = NULL, rhos = 0, beta = 1, type = c("mle", "deseq", "quantile", "tmm"),
                  prior = NULL, transform = FALSE, alpha = NULL, truephi = NULL, phi.epsilon = 0.05,
                  target = NULL, normalize.target = TRUE, delta = NULL, return.selected.gene.names = FALSE, ...){

  type <- match.arg(type)

  if (!is.null(target)){
    if (target < 0){
      warning("'target' should be non-negative. It is set to 0.")
      target <- 0
    }
  }

  if (is.null(xte)){
    xte <- x
    warning("Since no xte was provided, testing was performed on training data set.")
  }

  if (is.null(rhos)){
    warning("'rhos' cannot be NULL. Setting rho = 0.")
  }

  # if (length(rhos) > 1){
  #   warning("Multiple values of rho is given. Only the first element of rho is used.")
  # }

  if (transform){
    if (is.null(alpha)){
      ### Train ve test seti buradaki "alpha" değerine göre transform ediliyor.
      ### Alpha değeri MLE altında size factor hesabı kullanılarak yapılıyor.
      ### TO DO ####
      ### Further Research: Farklı size factor kestirimleri altında power transformations'lar kullanılabilir.
      alpha <- FindBestTransform(x)
    }

    if (alpha <= 0 || alpha > 1) stop("alpha must be between 0 and 1")
    x <- x^alpha
    xte <- xte^alpha
  }

  if (length(unique(y)) < 2){
    stop("Response variable must have at least two classes.")
  }

  if (is.null(prior)){
    prior <- rep(1/length(unique(y)), length(unique(y)))
  }

  null.out <- NullModel(x, type = type)
  ns <- null.out$n
  nste <- NullModelTest(null.out, xte)$nste

  uniq <- sort(unique(y))

  ## Dispersion estimates.
  if (!is.null(truephi)){
    if (length(truephi) >= 2 & (length(truephi) != ncol(ns))){
      truephi <- truephi[1]
      disperhat <- rep(truephi, ncol(ns))
      warning("The length of \"truephi\" should be the same as number of features. Only the first element is used and replicated for all features.")

    } else if (length(truephi) == 1){
      disperhat <- rep(truephi, ncol(ns))

    } else if (length(truephi) >= 2 & (length(truephi) == ncol(ns))){
      disperhat <- truephi
    }
    disperhat.list = list(initial = disperhat, adj = disperhat, cmp = NULL, delta = NULL, target = NULL)

  } else {
    ####### TO DO #######
    # tt değerinin 0'a eşit olması durumu ile analiz yapılıyor ama normalde üstteki değere shrink ediliyor. Yu et al çalışmasındaki gibi.
    # Bunun nedeni araştırılacak. 0 Değerini almak problem yaratır mı diye incelenecek.
    if (is.null(target)){
      if (normalize.target){
        tt <- getT(t(x), sizeFactors = null.out$rawsizestr, verbose = FALSE, propForSigma = c(0, 1))$target
      } else {
        tt <- getT(t(x), sizeFactors = rep(1, nrow(x)), verbose = FALSE, propForSigma = c(0, 1))$target
      }

    } else {
      tt <- target
    }

    ### Moment estimation of gene-wise dispersions.
    rM = rowMeans(t(x))
    rV = sSeq:::rowVars(t(x))

    disp = (rV - rM)/rM^2   ## Dispersion estimates using method-of-moments.
    # disp0 = numeric()

    ## Negative dispersions are set to 0.
    disp0 <- sapply(disp, function(x){
      max(0, x)
    })

    # disperhat = getAdjustDisp(disp0, shrinkTarget = tt, verbose = FALSE)$adj
    disperhat.list = getShrinkedDispersions(obs = disp0, shrinkTarget = tt, verbose = FALSE, delta = delta)
    disperhat <- disperhat.list$adj
    disperhat <- pmax(rep(0, length(disperhat)), disperhat) ## Negative dispersions are set to 0.
  }

  phihat <- as.numeric(disperhat)

  ### Shrink dispersion estimates to 0 using phi.epsilon.
  if (any(phihat <= phi.epsilon)){
    phihat[phihat <= phi.epsilon] <- 0
    disperhat <- disperhat.list$adj <- phihat
  }

  # Replace Infinity with zero dispersions.
  inv.phihat <- (1/phihat)
  if (any(inv.phihat == Inf | inv.phihat == -Inf)){
    id.inf <- which(inv.phihat == Inf | inv.phihat == -Inf)
    inv.phihat[id.inf] <- 0
  } else {
    id.inf <- NULL
  }

  if (length(rhos) == 1){
    ds <- GetD(ns = ns, x = x, y = y, rho = rhos, beta = beta)
    discriminant <- matrix(NA, nrow = nrow(xte), ncol = length(uniq))

    # NBLDA
    for (k in 1:length(uniq)){
      dstar = ds[k, ]
      part2 = 1 + t(nste) * dstar * phihat
      part1 = dstar/part2

      part3 <- inv.phihat * log(part2)
      if (!is.null(id.inf)){
        part3.limits <- t(nste[ ,id.inf]) * dstar[id.inf]
        part3[id.inf, ] <- part3.limits  # Replace zero dispersed genes with limit values.
      }

      discriminant[ ,k] <- rowSums(xte * t(log(part1))) - colSums(part3) + log(prior[k])
    }

    save <- list(list(ns = ns, nste = nste, ds = ds, discriminant = discriminant, ytehat = uniq[apply(discriminant, 1, which.max)],
                      alpha = alpha, rho = rhos, x = x, y = y, xte = xte, type = type, prior = prior, dispersions = disperhat.list,
                      transform = transform, trainParameters = null.out))
  } else {
    save <- lapply(rhos, function(u){
      ds <- GetD(ns = ns, x = x, y = y, rho = u, beta = beta)
      discriminant <- matrix(NA, nrow = nrow(xte), ncol = length(uniq))

      # NBLDA
      for (k in 1:length(uniq)){
        dstar = ds[k, ]
        part2 = 1 + t(nste) * dstar * phihat
        part1 = dstar/part2

        part3 <- inv.phihat * log(part2)
        if (!is.null(id.inf)){
          part3.limits <- t(nste[ ,id.inf])*dstar[id.inf]
          part3[id.inf, ] <- part3.limits  # Replace zero dispersed genes with limit values.
        }

        discriminant[ ,k] <- rowSums(xte * t(log(part1))) - colSums(part3) + log(prior[k])
      }

      save.i <- list(ns = ns, nste = nste, ds = ds, discriminant = discriminant, ytehat = uniq[apply(discriminant, 1, which.max)],
                     alpha = alpha, rho = u, x = x, y = y, xte = xte, type = type, prior = prior, dispersions = disperhat.list,
                     transform = transform, trainParameters = null.out)
      save.i
    })
  }

  # class(save) <- "poicla"
  return(save)
}



#' @importFrom sSeq getT getAdjustDisp rowVars
#'
#' @export
trainNBLDA <- function(x, y, type = c("mle", "deseq", "quantile", "tmm"), tuneLength = 10, metric = c("accuracy", "error"),
                       train.control = control(), ...){

  type <- match.arg(type)
  metric <- match.arg(metric)
  y <- as.factor(y)
  tc <- train.control

  if (tc$transform){
    if (is.null(tc$alpha)){
      alpha <- FindBestTransform(x)
    }
    if (alpha <= 0 || alpha > 1) stop("alpha must be between 0 and 1")
    x <- x^alpha
    tc$alpha <- alpha
  } else {
    alpha <- NULL
  }

  # 1. Define sets of model parameter values to evaluate
  if (is.null(tc$rhos)){
    ns <- NullModel(x, type = type)$n
    uniq <- sort(unique(y))
    maxrho <- rep(NA, length(uniq))

    for (k in 1:length(uniq)){
      a <- colSums(x[y == uniq[k],]) + tc$beta
      b <- colSums(ns[y == uniq[k],]) + tc$beta
      maxrho[k] <- max(abs(a/b - 1) * sqrt(b), na.rm = TRUE)
    }

    tc$rhos <- rhos <- seq(0, max(maxrho, na.rm = TRUE) * (2/3), len = tuneLength)
  }

  # 2. For each parameter set do
  if (is.null(tc$foldIdx)){
    # Each element of the returned list is the indices of test samples in this fold.
    allFolds <- lapply(1:tc$repeats, function(x){
      tmp <- balanced.folds(y, nfolds = tc$folds)
      names(tmp) <- paste("Fold.", 1:tc$folds, sep = "")
      tmp
    })
    names(allFolds) <- paste("Repeat.", 1:tc$repeats, sep = "")

    tc$foldIdx <- allFolds
    repeats <- tc$repeats
    folds <- tc$folds

  } else {
    allFolds <- tc$foldIdx
    repeats <- length(allFolds)
    folds <- length(allFolds[[1]])
  }

  tuningRes <- tuningResults <- NULL
  tuneGrid <- expand.grid(rho = rhos, epsilon = tc$epsilon)

  # for each fold within each repat, fit a plda model.
  # tuningResults <- foreach(iii = 1:nrow(tuneGrid), .combine = "c", .inorder = TRUE, .verbose = FALSE,
  #                          .export = c("tuneGrid", "allFolds")) %dopar% {
  #
  #
  # }

  foldRes <- lapply(allFolds, function(ii){
    errors <- acc <- nonzeros <- matrix(NA, nrow = length(rhos), ncol = length(ii))

    for (jj in 1:length(ii)){
      idx <- ii[[jj]]
      # out <- plda(x = x[-idx, ], y = y[-idx], xte = x[idx, ], rhos = rhos, beta = beta,
      #             type = "quantile", prior = prior, transform = FALSE)

      plist <- tc
      plist$transform <- FALSE
      plist <- c(list(x = x[-idx, ], y = y[-idx], xte = x[idx, ], type = type), plist)

      out <- do.call("nblda", plist)
      # out <- nblda(x = x[-idx, ], y = y[-idx], xte = x[idx, ], rhos = rhos, beta = tc$beta, type = type,
      #              prior = tc$prior, transform = FALSE, alpha = alpha, truephi = tc$true.dispersions,
      #              target = tc$target, normalize.target = tc$normalize.target, delta = tc$delta,
      #              return.selected.gene.names = tc$return.selected.gene.names)

      err.i <- nonzero.i <- acc.i <- NULL
      for (i in 1:length(out)){
        acc.i <- c(acc.i, sum(out[[i]][["ytehat"]] == y[idx]) / length(y[idx]))
        err.i <- c(err.i, sum(out[[i]][["ytehat"]] != y[idx]))
        nonzero.i <- c(nonzero.i, sum(colSums(out[[i]][["ds"]] != 1) != 0 | out[[i]]$dispersions$adj != 0))
      }

      acc[ ,jj] <- acc.i
      errors[ ,jj] <- err.i
      nonzeros[ ,jj] <- nonzero.i
    }

    list(errors = rowMeans(errors), accuracies = rowMeans(acc), nonzeros = rowMeans(nonzeros))
    # list(errors = errors, nonzeros = nonzeros)
  })

  err.sum <- nonzero.sum <- acc.sum <- 0
  for (r in 1:length(foldRes)){
    acc.sum <- acc.sum + foldRes[[r]]$accuracies
    err.sum <- err.sum + foldRes[[r]]$errors
    nonzero.sum <- nonzero.sum + foldRes[[r]]$nonzeros
  }

  tuningRes <- cbind(rhos, round(acc.sum/repeats, 10), round(err.sum/repeats, 10), round(nonzero.sum/repeats, 10))
  colnames(tuningRes) <- c("rho", "accuracy", "errors", "nonzero")

  if (metric == "accuracy"){
    best.rho <- tuningRes[, "rho"][max(which(tuningRes[, "accuracy"] == max(tuningRes[ ,"accuracy"])))]
  } else if (metric == "error"){
    best.rho <- tuningRes[, "rho"][max(which(tuningRes[, "errors"] == min(tuningRes[ ,"errors"])))]
  }

  # Set transform to FALSE since data is already transformed above.
  final.model <- nblda(x = x, y = y, xte = NULL, rhos = best.rho, beta = tc$beta, type = type, prior = tc$prior,
                       transform = FALSE, alpha = alpha, return.selected.gene.names = tc$return.selected.gene.names,
                       truephi = tc$truephi, target = tc$target, normalize.target = tc$normalize.target,
                       delta = tc$delta, phi.epsilon = tc$phi.epsilon)

  final.model[[1]]$transform <- tc$transform
  res <- list(tuning.results = tuningRes, best.rho = best.rho, final.model = final.model[[1]])
  res
}


predict.nblda <- function(trained.model, xte, return = c("predictions", "everything"), ...){
  ## Args:
  ##  trained.model: an object returned form trainNBLDA.
  ##  xte: a data.frame or matrix for test set. Samples in the rows and features in the columns.
  ##  return: should return everything or only predictions?

  return <- match.arg(return)

  fitted.model <- trained.model$final.model
  type <- fitted.model$type

  if (fitted.model$transform){
    alpha <- fitted.model$alpha
    if (alpha <= 0 || alpha > 1){
      stop("alpha must be between 0 and 1")
    }
    xte <- xte^alpha
  }

  ### prior trained model'den alınacak.
  prior <- fitted.model$prior

  # null.out <- NullModel(x, type = type)  ### trained modelden alınacak
  # ns <- null.out$n
  nste <- NullModelTest(fitted.model$trainParameters, xte = xte)$nste

  uniq <- sort(unique(fitted.model$y))
  # ds <- GetD(ns, x, y, rho, beta)   ### trained modelden alınacak
  ds <- fitted.model$ds
  phihat <- fitted.model$dispersions$adj
  discriminant <- matrix(NA, nrow = nrow(xte), ncol = length(uniq))

  # Replace Infinity with zero dispersions.
  inv.phihat <- (1/phihat)
  if (any(inv.phihat == Inf | inv.phihat == -Inf)){
    id.inf <- which(inv.phihat == Inf | inv.phihat == -Inf)
    inv.phihat[id.inf] <- 0
  } else {
    id.inf <- NULL
  }

  for (k in 1:length(uniq)){
    dstar = ds[k, ]
    part2 = 1 + t(nste) * dstar * phihat
    part1 = dstar/part2

    part3 <- inv.phihat * log(part2)
    if (!is.null(id.inf)){
      part3.limits <- t(nste[ ,id.inf]) * dstar[id.inf]
      part3[id.inf, ] <- part3.limits  # Replace zero dispersed genes with limit values.
    }

    discriminant[ ,k] <- rowSums(xte * t(log(part1))) - colSums(part3) + log(prior[k])
  }

  pred <- uniq[apply(discriminant, 1, which.max)]

  if (return == "predictions"){
    return(pred)
  } else if (return == "everything"){
    list(xte = xte, nste = nste, ds = ds, discriminant = discriminant, prior = prior,
         ytehat = pred, alpha = alpha, type = type, dispersions = phihat)
  }
}


getShrinkedDispersions <- function(obs, shrinkTarget = NULL, verbose = FALSE, delta = NULL){
  # This function returns the shrinked dispersion estimates.
  # Args:
  #   obs: a numeric vector. dispersion estimates for each feature.
  #   shrinkTarket: a numeric value. initial dispersion estimates are shrinked towards this value. If NULL, target value
  #                 is estimated from initial dispersion estimates.
  #   delta: a numeric value. This is the weight that is used for shrinkage algorithm. If 0, no shrinkage is
  #          performed on intiial values. If equals 1, initial values are forced to shrinked to target value.
  #          If NULL, weight are automatically estimated from initial disperson estimates.

  if (!is.null(delta)){
    if (delta < 0 | delta > 1){
      stop(warning("'delta' should be within [0, 1]. Calculation stopped..."))
    }
  }

  obs[is.na(obs)] = 0
  upBound = shrinkTarget

  if (verbose) {
    print(paste("Shrink initial dispersion estimates toward ", shrinkTarget, ".", sep = ""))
  }

  S.mat = var(obs, na.rm = T)
  cmp = data.frame(mean = mean(obs, na.rm = T), sigmasq.part = S.mat)
  mean.mat = rep(upBound, length(obs))
  dif.mat = obs - mean.mat
  dif2.mat = sum(dif.mat^2)

  if (is.null(delta)){
    delta = ((length(obs) - 2) * S.mat/(dif2.mat))   ## 1 - delta
    delta <- pmax(0, delta)
  }

  disp.shrinked <- delta * mean.mat + (1 - delta) * obs

  return(list(initial = obs, adj = disp.shrinked, cmp = cmp, delta = delta, target = shrinkTarget))
}

######## 2. PLDA FUNCTIONS ##########
## type seçenekleri "TMM", "deseq" ve "none" olacak.
plda <- function(x, y, xte = NULL, rhos = 0, beta = 1, type = c("mle", "deseq", "quantile", "tmm"),
                 prior = NULL, transform = TRUE, alpha = NULL, return.selected.gene.names = FALSE){

  type <- match.arg(type)

  if (is.null(xte)){
    xte <- x
    warning("Since no xte was provided, testing was performed on training data set.")
  }

  if (is.null(rhos)){
    warning("'rhos' cannot be NULL. Setting rho = 0.")
  }

  # if (length(rhos) > 1){
  #   warning("Multiple values of rho is given. Only the first element of rho is used.")
  # }

  if (transform){
    if (is.null(alpha)){
      ### Train ve test seti buradaki "alpha" değerine göre transform ediliyor.
      ### Alpha değeri MLE altında size factor hesabı kullanılarak yapılıyor.
      ### TO DO ####
      ### Further Research: Farklı size factor kestirimleri altında power transformations'lar kullanılabilir.
      alpha <- FindBestTransform(x)
    }
    if (alpha <= 0 || alpha > 1) stop("alpha must be between 0 and 1")
    x <- x^alpha
    xte <- xte^alpha
  }

  if (length(unique(y)) < 2){
    stop("Response variable must have at least two classes.")
  }

  if (is.null(prior)){
    prior <- rep(1/length(unique(y)), length(unique(y)))
  }

  null.out <- NullModel(x, type = type)
  ns <- null.out$n
  nste <- NullModelTest(null.out, xte)$nste

  uniq <- sort(unique(y))

  if (length(rhos) == 1){
    ds <- GetD(ns = ns, x = x, y = y, rho = rhos, beta = beta)
    discriminant <- matrix(NA, nrow = nrow(xte), ncol = length(uniq))

    for (k in 1:length(uniq)){
      discriminant[ ,k] <- rowSums(scale(xte, center = FALSE, scale = (1/log(ds[k, ])))) - rowSums(scale(nste, center = FALSE, scale = (1/ds[k, ]))) + log(prior[k])
    }

    save <- list(list(ns = ns, nste = nste, ds = ds, discriminant = discriminant, ytehat = uniq[apply(discriminant, 1, which.max)],
                      alpha = alpha, rho = rhos, x = x, y = y, xte = xte, type = type, prior = prior,
                      transform = transform, trainParameters = null.out))
  } else {
    save <- lapply(rhos, function(rho.i){
      ds <- GetD(ns = ns, x = x, y = y, rho = rho.i, beta = beta)
      discriminant <- matrix(NA, nrow = nrow(xte), ncol = length(uniq))

      for (k in 1:length(uniq)){
        discriminant[ ,k] <- rowSums(scale(xte, center = FALSE, scale = (1/log(ds[k, ])))) - rowSums(scale(nste, center = FALSE, scale = (1/ds[k, ]))) + log(prior[k])
      }

      save.i <- list(ns = ns, nste = nste, ds = ds, discriminant = discriminant, ytehat = uniq[apply(discriminant, 1, which.max)],
                     alpha = alpha, rho = rho.i, x = x, y = y, xte = xte, type = type, prior = prior,
                     transform = transform, trainParameters = null.out)
      save.i
    })
  }

  # class(save) <- "poicla"
  return(save)
}

#### Training algorithm (PLDA)#####
# 1. Define sets of model parameter values to evaluate
# 2. For each parameter set do
#   2.1 For each resampling iteration do
#     2.1.1 Hold-out specific samples
#     2.1.2 [Optional] Pre-process the data
#     2.1.3 Fit the model on the remainder
#     2.1.4 Predict hold-out samples
#   End
#   2.2 Calculate the average performance across hold-out predictions
# End
# 3. Determine the optimal parameter set
# Fit the final model to all the training data using the optimal parameter set

trainPLDA <- function(x, y, rhos = NULL, beta = 1, type = c("mle", "deseq", "quantile", "tmm"),
                      prior = NULL, transform = TRUE, alpha = NULL, return.selected.gene.names = FALSE,
                      nfolds = 5, repeats = 10, foldIdx = NULL, tuneLength = 30, metric = c("accuracy", "error"), ...){

  type <- match.arg(type)
  metric <- match.arg(metric)
  y <- as.factor(y)

  if (transform){
    if (is.null(alpha)){
      alpha <- FindBestTransform(x)
    }
    if (alpha <= 0 || alpha > 1) stop("alpha must be between 0 and 1")
    x <- x^alpha
  }

  # 1. Define sets of model parameter values to evaluate
  if (is.null(rhos)){
    ns <- NullModel(x, type = type)$n
    uniq <- sort(unique(y))
    maxrho <- rep(NA, length(uniq))
    for (k in 1:length(uniq)){
      a <- colSums(x[y == uniq[k],]) + beta
      b <- colSums(ns[y == uniq[k],]) + beta
      maxrho[k] <- max(abs(a/b - 1) * sqrt(b), na.rm = TRUE)
    }
    rhos <- seq(0, max(maxrho, na.rm = TRUE) * (2/3), len = tuneLength)
  }

  # 2. For each parameter set do
  if (is.null(foldIdx)){
    # Each element of the returned list is the indices of test samples in this fold.
    allFolds <- lapply(1:repeats, function(x){
      tmp <- balanced.folds(y, nfolds = nfolds)
      names(tmp) <- paste("Fold.", 1:nfolds, sep = "")
      tmp
    })
    names(allFolds) <- paste("Repeat.", 1:repeats, sep = "")
  } else {
    allFolds <- foldIdx
    repeats <- length(allFolds)
    nfolds <- length(allFolds[[1]])
  }

  tuningRes <- NULL
  # for each fold within each repat, fit a plda model.
  foldRes <- lapply(allFolds, function(ii){
    errors <- acc <- nonzeros <- matrix(NA, nrow = length(rhos), ncol = length(ii))

    for (jj in 1:length(ii)){
      idx <- ii[[jj]]
      # out <- plda(x = x[-idx, ], y = y[-idx], xte = x[idx, ], rhos = rhos, beta = beta,
      #             type = "quantile", prior = prior, transform = FALSE)

      out <- plda(x = x[-idx, ], y = y[-idx], xte = x[idx, ], rhos = rhos, beta = beta,
                  type = type, prior = prior, transform = FALSE)

      err.i <- nonzero.i <- acc.i <- NULL
      for (i in 1:length(out)){
        acc.i <- c(acc.i, sum(out[[i]][["ytehat"]] == y[idx]) / length(y[idx]))
        err.i <- c(err.i, sum(out[[i]][["ytehat"]] != y[idx]))
        nonzero.i <- c(nonzero.i, sum(colSums(out[[i]][["ds"]] != 1) != 0))
      }

      acc[ ,jj] <- acc.i
      errors[ ,jj] <- err.i
      nonzeros[ ,jj] <- nonzero.i
    }

    list(errors = rowMeans(errors), accuracies = rowMeans(acc), nonzeros = rowMeans(nonzeros))
    # list(errors = errors, nonzeros = nonzeros)
  })

  err.sum <- nonzero.sum <- acc.sum <- 0
  for (r in 1:length(foldRes)){
    acc.sum <- acc.sum + foldRes[[r]]$accuracies
    err.sum <- err.sum + foldRes[[r]]$errors
    nonzero.sum <- nonzero.sum + foldRes[[r]]$nonzeros
  }

  tuningRes <- cbind(rhos, round(acc.sum/repeats, 10), round(err.sum/repeats, 10), round(nonzero.sum/repeats, 10))
  colnames(tuningRes) <- c("rho", "accuracy", "errors", "nonzero")

  if (metric == "accuracy"){
    best.rho <- tuningRes[, "rho"][max(which(tuningRes[, "accuracy"] == max(tuningRes[ ,"accuracy"])))]
  } else if (metric == "error"){
    best.rho <- tuningRes[, "rho"][max(which(tuningRes[, "errors"] == min(tuningRes[ ,"errors"])))]
  }

  # Set transform to FALSE since data is already transformed above.
  final.model <- plda(x = x, y = y, xte = NULL, rhos = best.rho, beta = beta, type = type, prior = prior,
                      transform = FALSE, alpha = alpha, return.selected.gene.names = return.selected.gene.names)

  final.model[[1]]$transform <- transform
  res <- list(tuning.results = tuningRes, best.rho = best.rho, final.model = final.model[[1]])
  res
}


predict.plda <- function(trained.model, xte, return = c("predictions", "everything"), ...){
  ## Args:
  ##  x: raw count matrix which is used while training PLDA classifier. Samples in the rows and genes in the columns
  ##  y: a vector of classes for train set samples.
  ##  xte: a raw count matrix for test samples.
  ##  rho: a numeric value. This is the best value of tuning parameter.
  ##  beta: a numeric value. This is the beta parameter from trained model.
  ##  type: normalization type. Should be same as in trained model.
  ##  prior: prior probabilities of each class. Should be same as in trained model.
  ##  transform: TRUE/FALSE.
  ##  alpha: a numeri value between 0 and 1. It is used to perform power transformation on raw data.
  ##  null.out: train set parameters.
  ##  ds: offset parameters for each gene. gene expression parameters.

  return <- match.arg(return)

  fitted.model <- trained.model$final.model
  type <- fitted.model$type

  if (fitted.model$transform){
    alpha <- fitted.model$alpha
    if (alpha <= 0 || alpha > 1){
      stop("alpha must be between 0 and 1")
    }
    xte <- xte^alpha
  }

  ### prior trained model'den alınacak.
  prior <- fitted.model$prior

  # null.out <- NullModel(x, type = type)  ### trained modelden alınacak
  # ns <- null.out$n
  nste <- NullModelTest(fitted.model$trainParameters, xte = xte)$nste

  uniq <- sort(unique(fitted.model$y))
  # ds <- GetD(ns, x, y, rho, beta)   ### trained modelden alınacak
  ds <- fitted.model$ds
  discriminant <- matrix(NA, nrow = nrow(xte), ncol = length(uniq))
  for (k in 1:length(uniq)){
    discriminant[ ,k] <- rowSums(scale(xte, center = FALSE, scale = (1/log(ds[k, ])))) - rowSums(scale(nste, center = FALSE, scale = (1/ds[k, ]))) + log(prior[k])
  }
  pred <- uniq[apply(discriminant, 1, which.max)]

  if (return == "predictions"){
    return(pred)
  } else if (return == "everything"){
    list(xte = xte, nste = nste, ds = ds, discriminant = discriminant, prior = prior,
         ytehat = pred, alpha = alpha, type = type)
  }
}


## Find best value of trainsformation parameter "alpha" for power transformation.
FindBestTransform <- function(x){
  alphas <- seq(.01, 1, len = 50)
  gof <- rep(NA, length(alphas))
  for (alpha in alphas){
    gof[alphas == alpha] <- GoodnessOfFit(x^alpha, type = "mle")
  }
  return(alphas[which.min(abs(gof - (nrow(x) - 1)*(ncol(x) - 1)))])
}

NullModel <- function(x, type = c("mle", "deseq", "quantile", "none", "tmm")){
  # x MUST be a n times p matrix - i.e. observations on the rows and features on the columns
  type <- match.arg(type)
  rowsum <- rowSums(x)
  colsum <- colSums(x)

  if (type == "mle"){
    sizes <- rowSums(x)/sum(x)  # size factors for each sample
    mle <- outer(sizes, colsum, "*")
    return(list(n = mle, sizes = sizes, lib.size = rowsum, gene.length = colsum, type = type))

  } else if (type == "quantile"){ ### BU KISIM DÜZENLENECEK.
    # This is quantile normalization idea of Bullard et al 2010 -- quantile-normalize using 75th quantile of observed counts for each sample, excluding zero-counts
    sample.qts <- apply(x, 1, quantile, .75)
    sample.qts.pmax <- pmax(sample.qts, 1) # Don't wait to standardize by 0... min allowed is 1
    sizes <- sample.qts.pmax/sum(sample.qts.pmax)   # size factors for each sample
    fit <- outer(sizes, colsum, "*")
    return(list(n = fit, sizes = sizes, lib.size = rowsum, gene.length = colsum, sample.qts.tr = sample.qts,
                rawsizetr = sample.qts.pmax, type = type))

  } else if (type == "deseq"){ #Trying to implement idea from Anders and Huber Genome Biology 2010 paper.
    # I stole the next 3 lines from the DESeq bioconductor package.... it was in the function estimateSizeFactorsForMatrix
    counts <- t(x)
    geomeans <- exp(rowMeans(log(counts)))   # Geometric means for each features.
    sizes <- apply(counts, 2, function(cnts) median((cnts/geomeans)[geomeans > 0]))
    rawsizestr <- sizes  ## This part will be used to normalize test samples using train sample information.
    sizes <- sizes/sum(sizes)  # size factors for each sample
    fit <- outer(sizes, colsum, "*")   # s.i * g.j  eşitliği.
    return(list(n = fit, sizes = sizes, geomeans = geomeans, rawsizestr = rawsizestr,
                lib.size = rowsum, gene.length = colsum, type = type))

  } else if (type == "tmm"){
    counts <- t(x)
    null.out <- calcTMMnormFactors(object = counts)

    lib.size <- colSums(counts)  ## Library sizes
    gene.length <- rowSums(counts)   # Gene lengths

    nf <- null.out$norm.factors   ## Normalization factors
    sizes <- nf / sum(nf)
    fit <- outer(sizes, gene.length, "*")

    return(list(n = fit, refSample = null.out$refSample, geomean = null.out$geomean,
                sizes = sizes, normalization.factors = nf, lib.size = lib.size, type = type))

    # counts <- t(x)
    # countsDGE <- DGEList(counts = counts, genes = colnames(x))
    # countsDGE.normalized <- calcNormFactors(countsDGE, method = "TMM")   ## RLE: DESeq mantigi ile normalize ediyor.
    # countsDGE.normalized <- estimateCommonDisp(countsDGE.normalized)
    #
    # # TMM Normalized Counts
    # # From: http://grokbase.com/t/r/bioconductor/127r11kp23/bioc-edger-and-libsize-normalization
    # # n.normalized <- 1e6 * (x / (nf * lsize))   ## same results with cpm(..., log = FALSE)
    # fit <- t(countsDGE.normalized$pseudo.counts)
    #
    # # This function is used to find reference sample.
    # # Codes are copied from edgeR and modified here.
    # # rawCounts are count matrix with genes in the row and samples in the column
    # findRefSample <- function (rawCounts, lib.size, p = 0.75){
    #   y <- t(t(rawCounts) / lib.size)
    #   f75 <- apply(y, 2, function(x){
    #     quantile(x, p = p)
    #   })
    #   refSample <- which.min(abs(f75 - mean(f75)))
    #   return(refSample)
    # }
    #
    # nf <- countsDGE.normalized$samples$norm.factors   ## Normalization factors
    # lsize <- countsDGE.normalized$samples$lib.size  ## Library sizes
    #
    # refSampleID <- findRefSample(counts, lib.size = colSums(counts), p = 0.75)
    #
    # return(list(n = fit, refSampleID = refSampleID, refSample = x[refSampleID, ],
    #             normalization.factors = nf, lib.size = lsize, type = type))

  } else if (type == "none"){
    fit <- x
    lsize <- rowSums(x)

    return(list(n = fit, lib.size = lsize, type = type))
  }
}

NullModelTest <- function(null.out, xte = NULL){
  # xte MUST be a n times p matrix - i.e. observations on the rows and features on the columns

  type <- null.out$type

  if (is.null(xte)){
    stop("'xte' cannot be NULL. Which data set should be normalized?")
  }

  if (type == "mle"){
    sizeste <- rowSums(xte)/sum(null.out$lib.size)
    nste <- outer(sizeste, null.out$gene.length, "*")
    rawsize.xte <- NULL

  } else if (type == "quantile"){   ### BU KISIM DÜZENLENECEK.
    sizeste <- pmax(1, apply(xte, 1, quantile, .75)) # don't want to scale by anything less than 1...
    rawsize.xte <- sizeste
    sizeste <- sizeste/sum(null.out$sample.qts.tr)
    nste <- outer(sizeste, null.out$gene.length, "*")

  } else if (type == "deseq"){
    countste <- t(xte)
    geomeans <- null.out$geomeans
    geom.xte  = countste / geomeans
    rawsize.xte = apply(geom.xte, 2, function(x) median(x[geomeans > 0], na.rm = TRUE))

    sizeste <- rawsize.xte / sum(null.out$rawsizestr)  #### PoiClaClu paketi rawsize toplamları olarak train setin değerini kullanıyor. Bu değerlerin toplamı
    #### 1'e eşit çıkmıyor. Eşitliğin 1 olabilmesi için test setinin rawsize değerlerine bölmek gerekiyor.
    #### Burada bir problem olabilir mi???
    # sizete <- rawsize.xte / sum(rawsize.xte)

    nste <- outer(sizeste, null.out$gene.length, "*")

  } else if (type == "tmm"){
    referenceSample <- null.out$refSample

    ## Check if the feature names of reference sample are in the same order with rawCounts.
    if (identical(colnames(xte), names(referenceSample))){
      xte <- rbind(xte, referenceSample)
    } else {
      stop(warning("Reference sample either does not have same features or the features are not in the same order as test set. Calculation stops.."))
    }

    tmp.out <- calcTMMnormFactors(t(xte), refColumn = nrow(xte), geomean = null.out$geomean, test.sample = TRUE)
    nf <- tmp.out$norm.factors   ## Normalization factors
    # lsize <- countsDGE.nf$samples$lib.size  ## Library sizes

    xte <- xte[-nrow(xte), ]   ## Remove reference sample from test data.

    # sizeste <- nf / sum(nf)
    sizeste <- nf / sum(null.out$normalization.factors)   ## sum of normalization factors from training set. "deseq" ile benzer mantıkla.
    rawsize.xte <- nf

    nste <- outer(sizeste, colSums(xte), "*")
    lib.size <- rowSums(xte)  ## Library sizes
    gene.length <- colSums(xte)

  } else if (type == "none"){


  }

  return(list(nste = nste, sizeste = sizeste, type = type, rawsizete = rawsize.xte))
}


# Poisson Dissimilarities for Clustering.
# PoissonDistance <- function(x, beta = 1, type = c("mle", "deseq", "quantile"), transform = TRUE, alpha = NULL, perfeature = FALSE){
#   type <- match.arg(type)
#   if (!transform && !is.null(alpha)) stop("You have asked for NO transformation but have entered alpha.")
#   if (transform && !is.null(alpha)){
#     if (alpha > 0 && alpha <= 1) x <- x^alpha
#     if (alpha <= 0 || alpha>1) stop("alpha must be between 0 and 1")
#   }
#
#   if (transform && is.null(alpha)){
#     alpha <- FindBestTransform(x)
#     x <- x^alpha
#   }
#   dd <- matrix(0, nrow=nrow(x), ncol=nrow(x))
#   ddd <- NULL
#
#   if (perfeature) ddd <- array(0, dim=c(nrow(x), nrow(x), ncol(x)))
#
#   for (i in 2:nrow(dd)){
#     xi <- x[i,]
#
#     for (j in 1:(i-1)){
#       xj <- x[j,]
#       n <- NullModel(x[c(i,j),],type=type)$n
#       ni <- n[1,]
#       nj <- n[2,]
#       di <- (xi+beta)/(ni+beta)
#       dj <- (xj+beta)/(nj+beta)
#       dd[i,j] <- sum(ni+nj-ni*di-nj*dj+xi*log(di)+xj*log(dj))
#
#       if (perfeature) ddd[j,i,] <- ddd[i,j,] <- ni+nj-ni*di-nj*dj+xi*log(di)+xj*log(dj)
#     }
#   }
#
#   save <- list(dd=as.dist(dd+t(dd)), alpha=alpha, x=x, ddd=ddd, alpha=alpha, type=type)
#   class(save) <- "poidist"
#   return(save)
# }

# print.poidist <- function(x,...){
#   if(!is.null(x$alpha)) cat("Value of alpha used to transform data: ", x$alpha, fill = TRUE)
#   if(is.null(x$alpha)) cat("No transformation performed.", fill = TRUE)
#   cat("This type of normalization was performed:", x$type, fill = TRUE)
#   cat("Dissimilarity computed for ", nrow(x$x), " observations.", fill = TRUE)
# }




######### 3. HELPER FUNCTIONS   ##########
# 02/07/11 -- Helper functions.  (Copied from PoiClaClu source files.)
GoodnessOfFit <- function(x, type){
  ns <- NullModel(x, type = type)$n
  return(sum(((x - ns)^2)/ns, na.rm=TRUE))
}

# ColorDendrogram <- function(hc, y, main = "", branchlength = 0.7, labels = NULL, xlab = NULL, sub = NULL, ylab = "", cex.main = NULL){
#   if (is.null(labels))
#     labels <- rep("", length(y))
#   plot(hc, hang = 0.2, main = main, labels = labels, xlab = xlab,sub = sub, ylab = ylab, cex.main = cex.main)
#   y <- y[hc$ord]
#   if (is.numeric(y)) {
#     y <- y + 1
#     y[y == 7] <- "orange"
#   }
#   for (i in 1:length(hc$ord)) {
#     o = hc$merge[, 1] == -hc$ord[i] | hc$merge[, 2] == -hc$ord[i]
#     segments(i, hc$he[o] - branchlength, i, hc$he[o], col = y[i])
#   }
# }

permute.rows <- function(x){
  dd <- dim(x)
  n <- dd[1]
  p <- dd[2]
  mm <- runif(length(x)) + rep(seq(n) * 10, rep(p, n))
  matrix(t(x)[order(mm)], n, p, byrow = T)
}

#' @importFrom caret createFolds
foldIndex <- function(data = NULL, n = NULL, nFolds = 5, repeats = 2){
  if (!is.null(data)){
    n = nrow(data)
  }

  indIn <- indOut <- list()

  for (j in 1:repeats){
    tmpIn = createFolds(y = 1:n, k = nFolds, list = TRUE, returnTrain = TRUE)
    tmpOut = lapply(tmpIn, function(x)c(1:n)[-x])

    indIn = c(indIn, tmpIn)
    indOut = c(indOut, tmpOut)
  }

  nms = paste(rep(paste("Fold", 1:nFolds, sep = ""), repeats),
              rep(paste(".Rep", 1:repeats, sep = ""), c(rep(nFolds, repeats))), sep = "")
  names(indIn) <- names(indOut) <- nms
  return(list(indexIn = indIn, indexOut = indOut))
}

# Create folds and determine the indices of samples in each fold.
balanced.folds <- function(y, nfolds = min(min(table(y)), 10)){
  totals <- table(y)
  fmax <- max(totals)
  nfolds <- min(nfolds, fmax)
  # makes no sense to have more folds than the max class size
  folds <- as.list(seq(nfolds))
  yids <- split(seq(y), y)
  # nice way to get the ids in a list, split by class
  ###Make a big matrix, with enough rows to get in all the folds per class
  bigmat <- matrix(NA, ceiling(fmax/nfolds) * nfolds, length(totals))
  for(i in seq(totals)) {
    bigmat[seq(totals[i]), i] <- sample(yids[[i]])
  }
  smallmat <- matrix(bigmat, nrow = nfolds) # reshape the matrix
  ### Now do a clever sort to mix up the NAs
  smallmat <- permute.rows(t(smallmat)) ### Now a clever unlisting
  x <- apply(smallmat, 2, function(x) x[!is.na(x)])
  if(is.matrix(x)){
    xlist <- list()
    for(i in 1:ncol(x)){
      xlist[[i]] <- x[,i]
    }
    return(xlist)
  }
  return(x)
}

## Standart error of a given vector of parameter estimation.
se <- function(vec){
  return(sd(vec)/sqrt(length(vec)))
}

Soft <- function(x, a){
  return(sign(x)*pmax(abs(x) - a, 0))
}

# Find d.kj estimates
GetD <- function(ns, x, y, rho = NULL, beta){
  if (is.null(rho)){
    rho <- 0
    warning("'rho' is not given. It is set to rho = 0.")
  }

  if (missing(beta) || is.null(beta)){
    beta <- 1
  }

  uniq <- sort(unique(y))
  ds <- matrix(1, nrow = length(uniq), ncol = ncol(x))
  for (k in 1:length(uniq)){
    a <- colSums(x[y == uniq[k], ]) + beta
    b <- colSums(ns[y == uniq[k], ]) + beta
    ds[k, ] <- 1 + Soft(a/b - 1, rho/sqrt(b))
  }
  return(ds)
}


###### 4. TMM Normalization functions ######
.calcFactorQuantile_modified <- function (data, lib.size, p = 0.75){
  y <- t(t(data)/lib.size)
  f <- apply(y, 2, function(x) quantile(x, p = p))
  f
}

.calcFactorWeighted_modified <- function (obs, ref, libsize.obs = NULL, libsize.ref = NULL, logratioTrim = 0.3,
                                          sumTrim = 0.05, doWeighting = TRUE, Acutoff = -1e+10){
  obs <- as.numeric(obs)
  ref <- as.numeric(ref)
  if (is.null(libsize.obs))
    nO <- sum(obs)
  else nO <- libsize.obs
  if (is.null(libsize.ref))
    nR <- sum(ref)
  else nR <- libsize.ref
  logR <- log2((obs/nO)/(ref/nR))
  absE <- (log2(obs/nO) + log2(ref/nR))/2
  v <- (nO - obs)/nO/obs + (nR - ref)/nR/ref
  fin <- is.finite(logR) & is.finite(absE) & (absE > Acutoff)
  logR <- logR[fin]
  absE <- absE[fin]
  v <- v[fin]
  if (max(abs(logR)) < 1e-06)
    return(1)
  n <- length(logR)
  loL <- floor(n * logratioTrim) + 1
  hiL <- n + 1 - loL
  loS <- floor(n * sumTrim) + 1
  hiS <- n + 1 - loS
  keep <- (rank(logR) >= loL & rank(logR) <= hiL) & (rank(absE) >=
                                                       loS & rank(absE) <= hiS)
  if (doWeighting)
    f <- sum(logR[keep]/v[keep], na.rm = TRUE)/sum(1/v[keep],
                                                   na.rm = TRUE)
  else f <- mean(logR[keep], na.rm = TRUE)
  if (is.na(f))
    f <- 0
  2^f
}

calcTMMnormFactors <- function(object, lib.size = NULL, refColumn = NULL, logratioTrim = 0.3, geomean = NULL,
                               sumTrim = 0.05, doWeighting = TRUE, Acutoff = -1e+10, p = 0.75, test.sample = FALSE, ...){
  # Args:
  #   test.sample: logical. TRUE ise test örneklemi için hesaplama yapıyor.
  #   geomean: train setinde hesaplanan f değerlerinin geometrik ortalaması.
  #   refColumn: referans sample'in kolon numarası. Train setinden alınan sample. object içerisine eklenmiş olarak gönderilmesi gerekiyor.

  x <- as.matrix(object)
  if (any(is.na(x))) stop("NA counts not permitted")
  if (is.null(lib.size)) lib.size <- colSums(x)
  if (any(is.na(lib.size))) stop("NA lib.sizes not permitted")

  allzero <- .rowSums(x > 0, nrow(x), ncol(x)) == 0

  if (any(allzero)) x <- x[!allzero, , drop = FALSE]
  if (nrow(x) == 0 || ncol(x) == 1) method = "none"

  f75 <- .calcFactorQuantile_modified(data = x, lib.size = lib.size, p = 0.75)
  if (is.null(refColumn)) refColumn <- which.min(abs(f75 - mean(f75)))
  if (length(refColumn) == 0 | refColumn < 1 | refColumn > ncol(x)) refColumn <- 1
  f <- rep(NA, ncol(x))
  for (i in 1:ncol(x)) f[i] <- .calcFactorWeighted_modified(obs = x[ ,i], ref = x[ ,refColumn], libsize.obs = lib.size[i],
                                                            libsize.ref = lib.size[refColumn], logratioTrim = logratioTrim,
                                                            sumTrim = sumTrim, doWeighting = doWeighting, Acutoff = Acutoff)

  if (is.null(geomean)){
    geomean <- exp(mean(log(f)))
  }
  f <- f/geomean

  if (test.sample){
    f <- f[-refColumn]
  }
  return(list(norm.factors = f, geomean = geomean, refSample = x[ ,refColumn], type = "tmm"))
}

##### PLDA ve NBLDA kodlarının uyarlanması.

## Generate count data
generateCountData <- function (n, p, K, param, sdsignal = 1, DE = 0.3, allZero.rm = TRUE, tag.samples = FALSE) {
  # Args:
  #   n: number of samples
  #   p: number of features
  #   K: number of classes
  #   sdsignal: The extent to which the classes are different. If this equals zero then there
  #             are no class differences and if this is large then the classes are very different.
  #             This parameter is used to siulate dkj parameters.
  #   DE: probability of a gene being differentially expressed between classes.
  #   allZero.rm: a logical. If TRUE, features with only zero counts will be removed. Otherwise, it will be kept
  #               in the model as a presence of excess zeros. Might be suitable for zero-inflated scenarios.

  if (n < 4 * K){
    stop("We require n to be at least 4*K.")
  }
  if (DE < 0 || DE > 1){
    stop("DE should be a value within [0, 1].")
  }
  q0 <- rexp(p, rate = 1/25)   # Gene Total
  isDE <- runif(p) < DE
  classk <- matrix(NA, nrow = K, ncol = p)
  for (k in 1:K) {
    lfc <- rnorm(p, sd = sdsignal)   # dkj parametreleri.
    classk[k, ] <- ifelse(isDE, q0 * exp(lfc), q0)   ## Hangi genlerin DE olduğu bilgisi. Bunları da generate edilen veri ile döndürebiliriz.
    ##
  }
  truesf <- runif(n) * 2 + 0.2
  truesfte <- runif(n) * 2 + 0.2
  conds <- sample(c(rep(1:K, 4), sample(1:K, n - 4 * K, replace = TRUE)))
  condste <- sample(c(rep(1:K, 4), sample(1:K, n - 4 * K, replace = TRUE)))
  x <- xte <- matrix(NA, nrow = n, ncol = p)
  for (i in 1:n) {
    for (k in 1:K) {
      if (conds[i] == k)
        x[i, ] <- rnbinom(p, mu = truesf[i] * classk[k, ], size = param)
      if (condste[i] == k)
        xte[i, ] <- rnbinom(p, mu = truesfte[i] * classk[k, ], size = param)
    }
  }
  if (allZero.rm){
    rm <- apply(x, 2, sum) == 0   # Removes all zero columns.
  } else {
    rm <- logical(ncol(x))
  }

  # x <- as.data.frame(x)
  # xte <- as.data.frame(xte)
  colnames(x) <- colnames(xte) <- paste("G", 1:p, sep = "")

  if (tag.samples){
    rownames(x) <- rownames(xte) <- paste("S", 1:n, sep = "")
  }

  return(structure(list(x = t(x[, !rm]), xte = t(xte[, !rm]), y = conds, yte = condste,
                        truesf = truesf, truesfte = truesfte), class = "count.data"))
}


