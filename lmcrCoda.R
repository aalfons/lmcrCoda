# ---------------------------------------
# Authors: Andreas Alfons
#          Erasmus Universiteit Rotterdam
#
#          Nikola Stefelova
#          Palacky University Olomouc
# ---------------------------------------


# Function lmcrCoda() implements a cellwise and rowwise robust estimator of
# linear regression with compositional covariates.  Details on the estimator
# can be found in:
#
# N. Stefelova, A. Alfons, J. Palarea-Albaladejo, P.Filzmoser and K. Hron.
# Robust regression with compositional covariates including cellwise outliers.
# Under review.


# load packages
library("cellWise")
library("robCompositions")
library("robustHD")


# function to compute pivot coordinates
pivotCoord <- function(x, pivotvar = 1) {
  # initializations
  D <- ncol(x)
  order <- seq_len(D)
  if (pivotvar != 1) order <- c(pivotvar, order[-pivotvar])
  # compute pivot coordinates
  sapply(seq_len(D-1), function(j) {
    denominator <- order[j]
    numerator <- order[seq.int(j+1, D)]
    xNumerator <- x[, numerator, drop = FALSE]
    -sqrt((D-j) / (D-j+1)) * log(apply(xNumerator, 1, gm) / x[, denominator])
  })
}

# function to backtransform from pivot coordinates
pivotCoordInv <- function(x) {
  x <- -x
  y <- matrix(0, nrow=nrow(x), ncol=ncol(x)+1)
  D <- ncol(x)+1
  y[,1] <- -sqrt((D-1)/D)*x[,1]
  for (i in 2:ncol(y)){
    for (j in 1:(i-1)){
      y[,i]=y[,i]+x[,j]/sqrt((D-j+1)*(D-j))
    }
  }
  for (i in 2:(ncol(y)-1)){
    y[,i]=y[,i] - x[,i] * sqrt((D-i)/(D-i+1))
  }
  yexp=exp(y)
  x.back=yexp/apply(yexp,1,sum) # * rowSums(derOriginaldaten)
  # if(is.data.frame(x)) x.back <- data.frame(x.back)
  return(x.back)
}

# function for computing geometric mean
gm <- function (x) {
  if(!is.numeric(x)) stop("x has to be a vector of class numeric")
  if (any(na.omit(x == 0))) 0
  else exp(mean(log(unclass(x)[is.finite(x) & x > 0])))
}


# X ............. matrix of compositional explanatory variables
# R ............. matrix of additional real-valued covariates
# tau ........... specifies the quantile to be used for flagging cellwise
#                 outliers
# numDiscrete ... a variable with numDiscrete or fewer values will be
#                 considered discrete and ignored for flagging cellwise
#                 outliers
# pOutLR ........ compositional parts are flagged as cellwise outliers if the
#                 percentage of outlying logratios involving that part are at
#                 least pOutLR
# pOutRow ....... observations are flagged as rowwise outliers if the
#                 percentage of outlying cells is at least pOutRow (or if they
#                 are flagged by DDC as rowwise outliers in the data with
#                 logratios)
ddcCoda <- function(X, R, tau = 0.99, numDiscrete = 3,
                    pOutLR = 0.5, pOutRow = 0.75) {
  ## preparing compositional data
  X <- as.matrix(X)
  n <- nrow(X)
  D <- ncol(X)
  namesX <- colnames(X)
  if (is.null(namesX)) {
    namesX <- paste0("X", seq_len(D))
    colnames(X) <- namesX
  }
  ## preparing additional real-valued variables
  R <- as.matrix(R)
  p <- ncol(R)
  namesR <- colnames(R)
  if (is.null(namesR)) {
    namesR <- paste0("R", seq_len(p))
    colnames(R) <- namesR
  }

  # compute pairwise logratios
  combinations <- combn(namesX, 2, simplify = FALSE)
  Xlr <- lapply(combinations, function(j) log(X[, j[1]] / X[, j[2]]))
  Xlr <- do.call(cbind, Xlr)
  colnames(Xlr) <- sapply(combinations, function(j) paste(j, collapse = "/"))

  # detect outlying cells in pairwise logratios
  Xddc <- data.frame(Xlr, R, check.names = FALSE)
  ddc <- DDC(Xddc, DDCpars = list(tolProb = tau, numDiscrete = numDiscrete,
                                  silent = TRUE))

  # convert index of outlying cells to array indices
  indcells <- arrayInd(ddc$indcells, dim(ddc$remX))
  indcells <- data.frame(row = indcells[, 1],
                         col = colnames(ddc$remX)[indcells[, 2]],
                         stringsAsFactors = FALSE)
  # indices outlying cells in the real-valued variables
  whichR <- which(indcells$col %in% namesR)
  tmp <- indcells[whichR, , drop = FALSE]
  indicesR <- split(tmp$row, factor(tmp$col, levels = namesR))
  # indices of outlying cells in the pairwise logratios
  tmp <- indcells[-whichR, , drop = FALSE]
  indicesXlr <- split(tmp$row, tmp$col)

  # find outlying cells in compositional parts (the percentage of the
  # outlying logratios involving that part must be at least pOutLR)
  indicesX <- sapply(namesX, function(part, namesXlr) {
    whichLR <- grep(part, namesXlr, fixed = TRUE)
    indicesLR <- as.factor(unlist(indicesXlr[whichLR], use.names = FALSE))
    indices <- as.integer(levels(indicesLR))
    thresholdLR <- pOutLR * length(whichLR)
    indices[tabulate(indicesLR) >= thresholdLR]
  }, namesXlr = names(indicesXlr), simplify = FALSE)

  # find outlying rows (observations that are flagged by DDC in the data matrix
  # with pairwise logratios, or observations where the percentage of cellwise
  # outliers in the data matrix with compostions is at least pOutParts)
  indicesXR <- c(indicesX, indicesR, recursive = TRUE, use.names = FALSE)
  indicesXR <- as.factor(indicesXR)
  thresholdRow <- pOutRow * (D + p)
  indicesRow <- as.integer(levels(indicesXR))
  indicesRow <- indicesRow[tabulate(indicesXR) >= thresholdRow]
  indicesRow <- union(as.integer(ddc$indrows), indicesRow)

  # flag only outlying cells that are not part of outlying rows
  indicesX <- lapply(indicesX, setdiff, indicesRow)
  indicesR <- lapply(indicesR, setdiff, indicesRow)

  # return indices of cellwise and rowwise outliers
  list(indicesX = indicesX, indicesR = indicesR, indicesRow = indicesRow)

}


# function to set missing values in a vector
setNA <- function(x, i) {
  x[i] <- NA
  x
}


# function to (multiply) impute missing values
# code is based on function impCoda() from package robCompositions
# X ... matrix or data frame of compositional data with missing values
# R ... matrix or data frame of real-valued variables with missing values
# m ... number of imputations
# see documentation of impCoda() of package robCompositions for other input
# arguments
miCoda <- function(X, R, m = NULL, maxit = 10, eps = 0.5, method = "lmrob",
                   control = lmrob.control(), corrlim = 0.5, init = "KNN",
                   k = 5) {

  ## preparing compositional data
  X <- as.matrix(X)
  n <- nrow(X)
  D <- ncol(X)
  namesX <- colnames(X)
  if (is.null(namesX)) {
    namesX <- paste0("X", seq_len(D))
    colnames(X) <- namesX
  }
  ## preparing additional real-valued variables
  R <- as.matrix(R)
  p <- ncol(R)
  namesR <- colnames(R)
  if (is.null(namesR)) {
    namesR <- paste0("R", seq_len(p))
    colnames(R) <- namesR
  }
  ## further initializations
  method <- match.arg(method)
  init <- match.arg(init)
  if (init == "KNN") {
    if (k > n/4) warning("'k' might be too large")
  }

  ## back-up copy of original data
  xOrig <- cbind(X, R)

  ## index of missings / non-missings
  missingX <- is.na(X)
  observedX <- !missingX
  missingR <- is.na(R)
  observedR <- !missingR

  ## order of the columns of the data according to the amount of missings
  ## (compositional variables are handled first, then real-valued variables)
  nMissingX <- colSums(missingX)
  indX <- sort.int(nMissingX, index.return = TRUE, decreasing = TRUE)$ix
  nMissingR <- colSums(missingR)
  indR <- sort.int(nMissingR, index.return = TRUE, decreasing = TRUE)$ix
  ## skip columns that don't have any missing values
  indX <- setdiff(indX, which(nMissingX == 0))
  indR <- setdiff(indR, which(nMissingR == 0))

  ## first step - replace all NAs with values with 'nearest neighbour' algorithm
  if (init == "KNN") {
    if (length(indX) > 0) {
      X <- impKNNa(X, k = k, metric = "Aitchison", normknn = TRUE)$xImp
    }
    if (length(indR) > 0) {
      Z <- pivotCoord(X)
      RZ <- cbind(R, Z)
      R <- impKNNa(RZ, metric = "Euclidean")$xImp[, seq_len(p), drop = FALSE]
    }
    x <- cbind(X, R)
  }

  ## iterative model-based imputation

  it <- 0
  criteria <- Inf
  errorX <- rep.int(0, D)
  errorR <- rep.int(0, p)
  seqD <- seq_len(D)


  # screen predictor variables (i.e., having robust correlation >= corrlim)
  select <- !is.null(corrlim) && is.finite(corrlim)
  if (select) {
    correlatedX <- lapply(seqD, function(j) {
      observed <- observedX[, j]
      Z <- pivotCoord(X[observed, , drop = FALSE], pivotvar = j)
      predictors <- cbind(Z[, -1, drop = FALSE], R[observed, , drop = FALSE])
      cor <- unname(apply(predictors, 2, corHuber, Z[, 1], fallback = TRUE))
      which(abs(cor) >= corrlim)
    })
    if (length(indR) == 0) Z <- pivotCoord(X)  # otherwise already computed above
    correlatedR <- lapply(seq_len(p), function(j) {
      observed <- observedR[, j]
      predictors <- cbind(R[observed, -j, drop = FALSE],
                          Z[observed, , drop = FALSE])
      cor <- unname(apply(predictors, 2, corHuber, R[observed, j],
                          fallback = TRUE))
      which(abs(cor) >= corrlim)
    })
  }


  ###########################################
  ###  start the iteration

  while (it <= maxit & criteria >= eps) {

    xold <- x
    it <- it + 1

    # loop over compositional parts with missing values
    # (sorted according to amount of NAs)
    for (j in indX) {
      # compute pivot coordinates with current part giving the first coordinate
      Z <- as.matrix(pivotCoord(X, pivotvar = j))
      # fit a regression model for the first pivot coordinate with the other
      # coordinates and all real-valued variables as predictors (use only
      # observations with an observed value in the current compositional part)
      observed <- observedX[, j]
      if (select) {
        predictors <- cbind(Z[observed, -1, drop = FALSE],
                            R[observed, , drop = FALSE])
        predictors <- cbind("(Intercept)" = rep.int(1, n - nMissingX[j]),
                            predictors[, correlatedX[[j]], drop = FALSE])
      } else {
        predictors <- cbind("(Intercept)" = rep.int(1, n - nMissingX[j]),
                            Z[observed, -1, drop = FALSE],
                            R[observed, , drop = FALSE])
      }
      fit <- lmrob.fit(predictors, Z[observed, 1], control = control)
      # predict the value of the pivot coordinate for observations with a
      # missing value in the current compositional part
      missing <- missingX[, j]
      if (select) {
        predictorsNA <- cbind(Z[missing, -1, drop = FALSE],
                              R[missing, , drop = FALSE])
        predictorsNA <- cbind("(Intercept)" = rep.int(1, nMissingX[j]),
                              predictorsNA[, correlatedX[[j]], drop = FALSE])
      } else {
        predictorsNA <- cbind("(Intercept)" = rep.int(1, nMissingX[j]),
                              Z[missing, -1, drop = FALSE],
                              R[missing, , drop = FALSE])
      }
      Z[missing, 1] <- predictorsNA %*% coef(fit)
      # transform back to compositions (ensure correct order of parts)
      order <- c(j, seqD[-j])
      X[, order] <- as.matrix(pivotCoordInv(Z))  # TODO: only transform changed observations (missing)
      # keep track of residual scale for final imputation with draw from
      # estimated conditional distribution
      errorX[j] <- fit$scale
    }
    # loop over real-valued variables with missing values
    # (sorted according to amount of NAs)
    # We don't actually need to compute pivot coordinates again.  It doesn't
    # matter which coordinate system we use, we can use the one from the last
    # iteration of the loop over the compositional parts.
    for (j in indR) {
      # fit a regression model for the current real-valued variable with the
      # other real-valued varialbes and all pivot coordinates as predictors
      # (use only observations with an observed value in the current variable)
      observed <- observedR[, j]
      if (select) {
        predictors <- cbind(R[observed, -j, drop = FALSE],
                            Z[observed, , drop = FALSE])
        predictors <- cbind("(Intercept)" = rep.int(1, n - nMissingR[j]),
                            predictors[, correlatedR[[j]], drop = FALSE])
      } else {
        predictors <- cbind("(Intercept)" = rep.int(1, n - nMissingR[j]),
                            R[observed, -j, drop = FALSE],
                            Z[observed, , drop = FALSE])
      }
      fit <- lmrob.fit(predictors, R[observed, j], control = control)
      # predict the value of the current real-valued variable for observations
      # with a missing value
      missing <- missingR[, j]
      if (select) {
        predictorsNA <- cbind(R[missing, -j, drop = FALSE],
                              Z[missing, , drop = FALSE])
        predictorsNA <- cbind("(Intercept)" = rep.int(1, nMissingR[j]),
                              predictorsNA[, correlatedR[[j]], drop = FALSE])
      } else {
        predictorsNA <- cbind("(Intercept)" = rep.int(1, nMissingR[j]),
                              R[missing, -j, drop = FALSE],
                              Z[missing, , drop = FALSE])
      }
      R[missing, j] <- predictorsNA %*% coef(fit)
      # keep track of residual scale for final imputation with draw from
      # estimated conditional distribution
      errorR[j] <- fit$scale
    }

    # check convergence
    x <- cbind(X, R)
    criteria <- sum(((xold - x)/x)^2, na.rm = TRUE)
    # FIXME: the above is a hack to allow for dummy variables in R (or other
    #        variables which include 0, but this needs to be properly avoided)

  }

  ## add an error to the imputed values (draw from the conditional distribution)
  # default number of imputations: percentage of observations with missings
  if (is.null(m)) {
    pMissing <- sum((rowSums(missingX) + rowSums(missingR)) > 0) / n
    m <- max(2, round(100 * pMissing))
  }
  # multiple imputations with added errors
  xList <- replicate(m, {

    # loop over compositional parts (still sorted according to amount of NAs,
    # but that doesn't matter here)
    XMI <- X
    for (j in indX) {
      Z <- pivotCoord(X, pivotvar = j)
      missing <- missingX[, j]
      Z[missing, 1] <- Z[missing, 1] +
        rnorm(nMissingX[j], mean = 0, sd = errorX[j] * sqrt(1 + nMissingX[j]/n))
      XMI[, j] <- pivotCoordInv(Z)[, 1]
    }
    # loop over real-valued variables (still sorted according to amount of NAs,
    # but that doesn't matter here)
    for (j in indR) {
      missing <- missingR[, j]
      R[missing, j] <- R[missing, j] +
        rnorm(nMissingR[j], mean = 0, sd = errorR[j] * sqrt(1 + nMissingR[j]/n))
    }
    # combine compositional data and real-valued variables
    cbind(XMI, R)

  }, simplify = FALSE)

  ## return multiply imputed data
  res <- list(xOrig = xOrig, xImp = x, xMI = xList, criteria = criteria,
              iter = it, maxit = maxit, w = sum(nMissingX) + sum(nMissingR),
              wind = cbind(missingX, missingR))
  class(res) <- "imp"
  invisible(res)

}


# Pool point estimates and standard errors from multiple imputation
#
# Input:
# fitList  a list as returned by lmrob.fit() containting the regression
#          models fitted to the imputed data sets
#
# Output:
# A matrix that contains the pooled coefficients in the first column, the
# pooled standard errors in the second column, the t-statistic in the third
# column, the estimated degrees of freedom in the fourth column, and the
# p-value in the fifth column (see slide 50 of Lecture 5 for an example)
pool <- function(fitList) {
  # extract number of imputations
  m <- length(fitList)

  # compute average estimates
  coefReps <- sapply(fitList, coef)
  means <- rowMeans(coefReps)

  # compute variance between imputations
  between <- apply(coefReps, 1, var)
  # compute variances within each imputation and average
  vcovReps <- lapply(fitList, vcov)
  varReps <- sapply(vcovReps, diag)
  within <- rowMeans(varReps)
  # compute standard errors by combining within and between variances
  variances <- within + (m+1)/m * between
  se <- sqrt(variances)
  # degrees of freedom if information were complete
  dfComplete <- df.residual(fitList[[1]])
  # degrees of freedom if n -> Inf (Barnard & Rubin 1999)
  gamma <- (m+1)/m * between / variances
  dfInf <- (m-1) / gamma^2
  # estimated degrees of freedom for observed data
  lambda <- (dfComplete + 1) / (dfComplete + 3)
  dfObserved <- lambda * dfComplete * (1 - gamma)
  # compute degrees of freedom for small sample test
  df <- dfInf * dfObserved / (dfInf + dfObserved)  # simpler notation

  # perform t-tests and combine results
  t <- means / se
  pValue <- pValueT(t, df = df)
  coefficients <- cbind(means, se, t, df, pValue)
  tn <- c("Estimate", "Std. Error", "t value", "df", "Pr(>|t|)")
  colnames(coefficients) <- tn
  coefficients
}


## utility functions

# fit linear models to imputed data sets
fitLmrob <- function(data, namesX, namesR, nameY, pivotvar = 1,
                     fun = lmrob, ...) {
  # extract explanatory variables and response
  X <- data[, namesX, drop = FALSE]
  R <- data[, namesR, drop = FALSE]
  y <- data[, nameY, drop = FALSE]
  # prepare predictor matrix
  n <- nrow(data)
  D <- ncol(X)
  Z <- pivotCoord(X, pivotvar = pivotvar)
  colnames(Z) <- paste0("Z", seq_len(D-1))
  data <- data.frame(Z, R, y)
  # fit regression model
  formula <- as.formula(paste(nameY, "~ ."))
  fun(formula, data = data, ...)
}

# hack to get MM standard errors after fitting weighted least squares
lmrobHack <- function(formula, data, fitMM) {
  matchedCall <- match.call()
  # prepare model frame
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  if(is.empty.model(mt)) stop("empty model")
  # extract response and candidate predictors from model frame
  y <- model.response(mf, "numeric")
  x <- model.matrix(mt, mf)
  # compute WLS with the weights of the initial S-estimator
  initS <- fitMM$init.S
  wS <- weights(initS, type = "robustness")
  init <- lm.wfit(x, y, w = wS)
  class(init) <- class(initS)
  # compute WLS with the weights of the final MM-estimator
  wMM <- weights(fitMM, type = "robustness")
  fit <- lm.wfit(x, y, w = wMM)
  # copy some information to make sure the summary() method works correctly
  fit$scale <- fitMM$scale
  fit$converged <- fitMM$converged
  fit$rweights <- wMM
  fit$control <- fitMM$control
  fit$init.S <- init
  fit$call <- matchedCall
  fit$terms <- mt
  fit$model <- mf
  fit$x <- x
  class(fit) <- class(fitMM)
  # estimate the covariance matrix of the covariates
  fit$cov <- vcov(fit)
  # return object
  fit
}

# function to compute p-value based on normal distribution
pValueZ <- function(z, alternative = c("twosided", "less", "greater")) {
  # initializations
  alternative <- match.arg(alternative)
  # compute p-value
  switch(alternative, twosided = 2*pnorm(abs(z), lower.tail = FALSE),
         less = pnorm(z), greater = pnorm(z, lower.tail = FALSE))
}

# function to compute p-value based on t-distribution
pValueT <- function(t, df, alternative = c("twosided", "less", "greater")) {
  # initializations
  alternative <- match.arg(alternative)
  # compute p-value
  switch(alternative, twosided = 2*pt(abs(t), df = df, lower.tail = FALSE),
         less = pt(t, df = df), greater = pt(t, df = df, lower.tail = FALSE))
}


# X ............... matrix of compositional explanatory variables
# y ............... real-valued response variable
# R ............... matrix of additional real-valued covariates
# Dummies ......... matrix of additional dummy variables
#                   (currently only one dummy variable allowed)
# interpretable ... a logical indicating whether to fit models using all
#                   coordinate systems for an interpretable model (TRUE),
#                   or only one coordinate system for prediction (FALSE)
# single .......... a logical indicating whether to add results for single
#                   imputation.  If FALSE, only results for multiple imputation
#                   are returned
# see ddcCoda() and miCoda() for other input arguments
lmcrCoda <- function(X, y, R = NULL, Dummies = NULL,
                     # arguments for the bivariate filter
                     tau = 0.99, numDiscrete = 3, pOutLR = 0.5,
                     pOutRow = 0.75,
                     # arguments for multiple imputation
                     m = NULL, maxit = 10, eps = 0.5, method = "lmrob",
                     control = lmrob.control(), corrlim = NULL,
                     init = "KNN", k = 5,
                     # which results to return
                     interpretable = TRUE, single = FALSE,
                     # for ideal filter: indices of outlying cells/row
                     indicesX = NULL, indicesY = NULL, indicesR = NULL,
                     indicesRow = NULL) {

  ## initializations
  n <- length(y)  # number of observations
  # preparing compositional data
  X <- as.matrix(X)
  D <- ncol(X)
  namesX <- colnames(X)
  if (is.null(namesX)) {
    namesX <- paste0("X", seq_len(D))
    colnames(X) <- namesX
  }
  # preparing additional covariates
  haveR <- !is.null(R)
  if (haveR) {
    R <- as.matrix(R)
    p <- ncol(R)
    namesR <- colnames(R)
    if (is.null(namesR)) {
      namesR <- paste0("R", seq_len(p))
      colnames(R) <- namesR
    }
  } else {
    p <- 0
    namesR <- character()
  }
  # set default threshold for variable screening in imputation step
  if (is.null(corrlim)) corrlim <- if (D + p + 1 > 10) 0.5 else 0.2
  # combine real-valued covariates with response variable
  R <- cbind(R, y = y)
  # preparing additional dummy variables
  haveDummies <- !is.null(Dummies)
  if (haveDummies) {
    Dummies <- as.matrix(Dummies)
    q <- ncol(Dummies)
    if (q > 1) stop("currently only implemented for one dummy variable")
    namesDummies <- colnames(Dummies)
    if (is.null(namesDummies)) {
      namesDummies <- paste0("D", seq_len(q))
      colnames(Dummies) <- namesDummies
    }
  } else {
    q <- 0
    namesDummies <- character()
  }
  # other initializations
  interpretable <- isTRUE(interpretable)
  single <- isTRUE(single)

  ## stage 1: filter cellwise outliers

  useDDC <- is.null(indicesX) || is.null(indicesY) ||
    (haveR && is.null(indicesR)) || is.null(indicesRow)
  if (useDDC) {

    # apply modified DDC algorithm

    if (haveDummies) {

      # indices for subsets
      which0 <- which(Dummies == 0)
      which1 <- which(Dummies == 1)

      # apply modified DDC algorithm separately per subgroup
      ddc0 <- ddcCoda(X[which0, ], R[which0, ], tau = tau,
                      numDiscrete = numDiscrete, pOutLR = pOutLR,
                      pOutRow = pOutRow)
      ddc1 <- ddcCoda(X[which1, ], R[which1, ], tau = tau,
                      numDiscrete = numDiscrete, pOutLR = pOutLR,
                      pOutRow = pOutRow)

      # extract indices
      indicesX <- mapply(function(i0, i1) c(which0[i0], which1[i1]),
                         ddc0$indicesX, ddc1$indicesX, SIMPLIFY = FALSE)
      indicesR <- mapply(function(i0, i1) c(which0[i0], which1[i1]),
                         ddc0$indicesR, ddc1$indicesR, SIMPLIFY = FALSE)
      indicesRow <- c(which0[ddc0$indicesRow], which1[ddc1$indicesRow])

    } else {

      # apply modified DDC algorithm
      ddc <- ddcCoda(X, R, tau = tau, numDiscrete = numDiscrete,
                     pOutLR = pOutLR, pOutRow = pOutRow)

      # extract indices
      indicesX <- ddc$indicesX
      indicesR <- ddc$indicesR
      indicesRow <- ddc$indicesRow

    }

    # compute number of cellwise outliers
    nOut <- sum(sapply(indicesX, length)) + sum(sapply(indicesR, length))

    # combine indices of outliers into one list to be returned with results
    indices <- list(indicesX = indicesX, indicesY = indicesR$y,
                    indicesR = indicesR[namesR], indicesRow = indicesRow)

  } else {

    # combine indices for real-valued covariates and response
    indicesR <- cbind(indicesR, indicesY)
    # compute number of cellwise outliers
    nOut <- sum(indicesX) + sum(indicesR)

  }

  ## if any cells have been set to missing, do all remaining steps, otherwise
  ## perform MM-regression on the original data
  if (nOut > 0) {

    ## filter outlying cells (set to NA)
    if (useDDC) {
      XNA <- mapply(setNA, as.data.frame(X), indicesX)  # returns matrix
      RNA <- mapply(setNA, as.data.frame(R), indicesR)  # returns matrix
    } else {
      XNA <- setNA(X, indicesX)
      RNA <- setNA(R, indicesR)
    }
    if (haveDummies) {
      # combine real-valued variables and dummy variables
      RNA <- cbind(RNA, Dummies)
      namesR <- c(namesR, namesDummies)
      haveR <- TRUE
    }

    ## stage 2: (multiple) imputation
    imputed <- miCoda(XNA, RNA, m = m, maxit = maxit, eps = eps,
                      method = method, control = control, corrlim = corrlim,
                      init = init, k = k)

    ## stages 3 and 4 depend on whether we just want to fit a model in one set
    ## of pivot coordinates (for prediction), or if we want to fit models for
    ## different pivot coordinates (for interpretation)

    if (interpretable) {

      ## stage 3: robust compositional regression
      ## (all sets of pivot coordinates)

      # fit regression models to each imputed data set
      # we need to fully iterate the MM-estimator for the first set of pivot
      # coordinates
      first <- lapply(imputed$xMI, fitLmrob, namesX = namesX, namesR = namesR,
                      nameY = "y", pivotvar = 1, fun = lmrob, control = control)
      # the others can be used with a hack that uses the exact same weights
      # for the other sets pivot coordinates
      other <- lapply(seq(2, D), function(j) {
        mapply(function(XR, fit) {
          fitLmrob(XR, namesX = namesX, namesR = namesR, nameY = "y",
                   pivotvar = j, fun = lmrobHack, fitMM = fit)
        }, XR = imputed$xMI, fit = first, SIMPLIFY = FALSE, USE.NAMES = FALSE)
      })
      # combine model fits for different sets of coordinates
      fitList <- c(list(first), other)
      names(fitList) <- namesX

      ## stage 4: aggregation of regressions from multiple imputed data sets
      coefMI <- lapply(fitList, pool)
      scalesMI <- sapply(first, "[[", "scale")
      scaleMI <- mean(scalesMI)

      ## add results from single imputation
      if (single) {
        # we need to fully iterate the MM-estimator for the first set of pivot
        # coordinates
        first <- fitLmrob(imputed$xImp, namesX = namesX, namesR = namesR,
                          nameY = "y", pivotvar = 1, fun = lmrob,
                          control = control)
        # the others can be used with a hack that uses the exact same weights
        # for the other sets pivot coordinates
        other <- lapply(seq(2, D), function(j) {
          fitLmrob(imputed$xImp, namesX = namesX, namesR = namesR, nameY = "y",
                   pivotvar = j, fun = lmrobHack, fitMM = first)
        })
        # combine model fits for different sets of coordinates
        fitList <- c(list(first), other)
        names(fitList) <- namesX
        # perform t-tests and combine results
        coefSI <- lapply(fitList, function(fit) {
          # extract coefficients and standard errors
          coefficients <- coef(fit)
          se <- sqrt(diag(vcov(fit)))
          df <- df.residual(fit)
          # perform t-tests and combine results
          t <- coefficients / se
          pValue <- pValueT(t, df = df)
          coefmat <- cbind(coefficients, se, t, pValue)
          tn <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
          colnames(coefmat) <- tn
          # return coefficient matrix
          coefmat
        })
        # extract residual scale
        scaleSI <- first$scale
      } else {
        coefSI <- NULL
        scaleSI <- NULL
      }

    } else {

      ## stage 3: robust compositional regression
      ## (only one set of pivot coordinates)

      # fit regression models to each imputed data set
      # we don't need the formula here to make lmrobHack() work, which
      # saves some computation time
      fitList <- lapply(imputed$xMI, function(data) {
        # extract explanatory variables and response
        X <- data[, namesX, drop = FALSE]
        R <- if (haveR) data[, namesR, drop = FALSE]
        y <- data[, "y"]
        # prepare predictor matrix
        n <- nrow(data)
        Z <- as.matrix(pivotCoord(X, pivotvar = 1))
        colnames(Z) <- paste0("Z", seq_len(D-1))
        predictors <- cbind("(Intercept)" = rep.int(1, n), Z, R)
        # fit regression model
        lmrob.fit(predictors, y, control = control)
      })

      ## stage 4: aggregation of regressions from multiple imputed data sets
      coefMI <- pool(fitList)
      scalesMI <- sapply(fitList, "[[", "scale")
      scaleMI <- mean(scalesMI)

      ## add results from single imputation
      if (single) {
        # extract explanatory variables and response
        XImp <- imputed$xImp[, namesX, drop = FALSE]
        RImp <- if (haveR) imputed$xImp[, namesR, drop = FALSE]
        yImp <- imputed$xImp[, "y"]
        # prepare predictor matrix
        ZImp <- as.matrix(pivotCoord(XImp, pivotvar = 1))
        colnames(ZImp) <- paste0("Z", seq_len(D-1))
        predictors <- cbind("(Intercept)" = rep.int(1, n), ZImp, RImp)
        # fit regression model
        fit <- lmrob.fit(predictors, yImp, control = control)
        coefficients <- coef(fit)
        se <- sqrt(diag(vcov(fit)))
        df <- df.residual(fit)
        # perform t-tests and combine results
        t <- coefficients / se
        pValue <- pValueT(t, df = df)
        coefSI <- cbind(coefficients, se, t, pValue)
        tn <- c("Estimate", "Std. Error", "t value", "Pr(>|t|)")
        colnames(coefSI) <- tn
        # extract residual scale
        scaleSI <- fit$scale
      } else {
        coefSI <- NULL
        scaleSI <- NULL
      }

    }

  } else {

    ## no cellwise outliers, perform MM-regression on original data
    # real-valued covariates and dummy variables have not been combined
    # in this branch, they have to be added to predictors separately

    if (interpretable) {

      # apply fully iterated MM-estimator on first coordinate system
      Z <- pivotCoord(X, pivotvar = 1)
      colnames(Z) <- paste0("Z", seq_len(D-1))
      covariates <- if (haveR) R[, namesR]
      if (haveDummies) covariates <- cbind(covariates, Dummies)
      data <- as.data.frame(cbind(Z, covariates, y = y))  # works also with NULL
      first <- lmrob(y ~ ., data = data, control = control)
      # apply WLS with the same weights to other coordinate systems
      other <- lapply(seq(2, D), function(pivotvar, x, y, covariates, fit) {
        z <- pivotCoord(x, pivotvar = pivotvar)
        colnames(z) <- paste0("Z", seq_len(D-1))
        data <- as.data.frame(cbind(z, covariates, y = y))# works also with NULL
        lmrobHack(y ~ ., data = data, fitMM = fit)
      }, x = X, covariates = covariates, y = y, fit = first)
      # combine model fits
      fitList <- c(list(first), other)
      names(fitList) <- colnames(X)
      # perform t-tests and combine results
      coefMI <- lapply(fitList, function(fit) {
        # extract coefficients and standard errors
        coefficients <- coef(fit)
        se <- sqrt(diag(vcov(fit)))
        df <- df.residual(fit)
        # perform t-tests
        t <- coefficients / se
        pValue <- pValueT(t, df = df)
        # mimic results from multiple imputation with extra column for
        # degrees of freedom
        coefmat <- cbind(coefficients, se, t, df, pValue)
        tn <- c("Estimate", "Std. Error", "t value", "df", "Pr(>|t|)")
        colnames(coefmat) <- tn
        # return coefficient matrix
        coefmat
      })
      # extract residual scale
      scaleMI <- first$scale

      ## add results from single imputation
      # (which are the same here but without the column for degrees of freedom)
      if (single) {
        coefSI <- lapply(coefMI, function(cm) cm[, -4, drop = FALSE])
        scaleSI <- scaleMI
      } else {
        coefSI <- NULL
        scaleSI <- NULL
      }

    } else {

      # prepare predictor matrix using only one coordinate system
      Z <- as.matrix(pivotCoord(X, pivotvar = 1))
      colnames(Z) <- paste0("Z", seq_len(D-1))
      predictors <- cbind("(Intercept)" = rep.int(1, n), Z)
      if (haveR) predictors <- cbind(predictors, R[, namesR])
      if (haveDummies) predictors <- cbind(predictors, Dummies)
      # fit regression model
      fit <- lmrob.fit(predictors, y, control = control)
      coefficients <- coef(fit)
      se <- sqrt(diag(vcov(fit)))
      df <- df.residual(fit)
      # perform t-tests and combine results
      t <- coefficients / se
      pValue <- pValueT(t, df = df)
      # mimic results from multiple imputation with extra column for
      # degrees of freedom
      coefMI <- cbind(coefficients, se, t, df, pValue)
      tn <- c("Estimate", "Std. Error", "t value", "df", "Pr(>|t|)")
      colnames(coefMI) <- tn
      # extract residual scale
      scaleMI <- fit$scale

      ## add results from single imputation
      # (which are the same here but without the column for degrees of freedom)
      if (single) {
        coefSI <- coefMI[, -4, drop = FALSE]
        scaleSI <- scaleMI
      } else {
        coefSI <- NULL
        scaleSI <- NULL
      }

    }

  }

  ## return object
  out <- list(SI = coefSI, scaleSI = scaleSI, MI = coefMI, scaleMI = scaleMI,
              interpretable = interpretable)
  if (useDDC) out <- c(out, indices)  # return indices of outliers
  class(out) <- "lmcrCoda"
  out

}


## print method
print.lmcrCoda <- function(x, which = c("MI", "SI"), ...) {
  which <- match.arg(which)
  object <- switch(which, MI = x$MI, SI = x$SI)
  if (is.null(object)) {
    object <- x$MI
    warning("results from single imputation not available; ",
            "returning those from multiple imputation")
  }
  if(x$interpretable) {
    # extract first pivot coordinate from each coefficient matrix
    # and stack them together
    coefmat <- mapply(function(fit, name) {
      coef <- fit[2, , drop = FALSE]
      rownames(coef) <- paste("Z1", name, sep = "_")
      coef
    }, object, names(object), SIMPLIFY = FALSE, USE.NAMES = FALSE)
    coefmat <- do.call(rbind, coefmat)
    # intercept and coefficients of non-compositional variables are identical
    # across coordinate systems, we can extract from any one of them
    intercept <- object[[1]][1, , drop = FALSE]
    D <- nrow(coefmat)
    p <- nrow(object[[1]]) - D
    others <- object[[1]][D + seq_len(p), , drop = FALSE]
    # put everything together
    coefmat <- rbind(intercept, coefmat, others)
  } else {
    # only one pivot coordinate system used
    coefmat <- object
  }
  # print results as usual in R
  printCoefmat(coefmat, ...)
}


## functions for OLS and MM-regression

# OLS estimator for different coordinate systems
lmCoda <- function(X, y, covariates = NULL) {
  # initializations
  D <- ncol(X)
  # apply OLS with different coordinate systems
  fitList <- lapply(seq_len(D), function(pivotvar = 1, x, y, covariates) {
    # prepare data
    z <- pivotCoord(x, pivotvar = pivotvar)
    colnames(z) <- paste0("Z", seq_len(D-1))
    data <- as.data.frame(cbind(z, covariates, y = y))  # also works with NULL
    # fit regression model
    lm(y ~ ., data = data)
  }, x = X, covariates = covariates, y = y)
  names(fitList) <- colnames(X)
  # compute summaries
  out <- lapply(fitList, function(fit) summary(fit)$coefficients)
  class(out) <- "lmCoda"
  out
}

# MM estimator for different coordinate systems
lmrobCoda <- function(X, y, covariates = NULL, control = lmrob.control()) {
  # initializations
  D <- ncol(X)
  # apply fully iterated MM-estimator on first coordinate system
  Z <- pivotCoord(X, pivotvar = 1)
  colnames(Z) <- paste0("Z", seq_len(D-1))
  data <- as.data.frame(cbind(Z, covariates, y = y))  # also works with NULL
  first <- lmrob(y ~ ., data = data, control = control)
  # apply WLS with the same weights to other coordinate systems
  other <- lapply(seq(2, D), function(pivotvar, x, y, covariates, fit) {
    z <- pivotCoord(x, pivotvar = pivotvar)
    colnames(z) <- paste0("Z", seq_len(D-1))
    data <- as.data.frame(cbind(z, covariates, y = y))  # also works with NULL
    lmrobHack(y ~ ., data = data, fitMM = fit)
  }, x = X, covariates = covariates, y = y, fit = first)
  # combine model fits
  fitList <- c(list(first), other)
  names(fitList) <- colnames(X)
  # compute summaries
  out <- lapply(fitList, function(fit) summary(fit)$coefficients)
  class(out) <- c("lmrobCoda", "lmCoda")
  out
}

## print method
print.lmCoda <- function(x, ...) {
  # extract first pivot coordinate from each coefficient matrix
  # and stack them together
  coefmat <- mapply(function(fit, name) {
    coef <- fit[2, , drop = FALSE]
    rownames(coef) <- paste("Z1", name, sep = "_")
    coef
  }, x, names(x), SIMPLIFY = FALSE, USE.NAMES = FALSE)
  coefmat <- do.call(rbind, coefmat)
  # intercept and coefficients of non-compositional variables are identical
  # across coordinate systems, we can extract from any one of them
  intercept <- x[[1]][1, , drop = FALSE]
  D <- nrow(coefmat)
  p <- nrow(x[[1]]) - D
  others <- x[[1]][D + seq_len(p), , drop = FALSE]
  # put everything together
  coefmat <- rbind(intercept, coefmat, others)
  # print results as usual in R
  printCoefmat(coefmat, ...)
}
