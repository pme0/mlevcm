#' @title Train a Machine Learning model
#'
#' @description Train a machine learning model.
#'
#' @param X A numeric functional predictor matrix of size \code{N*P}, where \code{N} is
#' the number of observations (spectra) and \code{P} is the number of points
#' (wavelengths) at which spectra are measured.
#'
#' @param y A numeric or factor response vector of size \code{N}, where \code{N} is the number of
#' observations (spectra).
#'
#' @param Z A numeric non-functional predictor matrix of size \code{N*S}, where \code{N}
#' is the number of observations (spectra) and \code{S} is the number of non-functional predictors
#' after bining/one-hot-encoding has taken place.
#'
#' @param task Regression (\code{"regr"}) for continuous response problems or Classification
#' (\code{"clas"}) for categorical response problems.
#'
#' @param model Currently Linear Model (\code{"lm"}) for Regression or Generalised Linear Model
#' (\code{"glm"}) for Classification.
#'
#' @param reduction Partial Least Squares (\code{"pls"}), Principal Component Analysis
#' (\code{"pca"}), or no further dimension reduction (\code{"n"}).
#'
#' @param smooth Whether to use the original spectra (\code{FALSE}) or smoothed spectra
#' (\code{TRUE}).
#'
#' @param balanced Whether the dataset should be balanced (\code{TRUE}) or not (\code{FALSE}).
#' If \code{TRUE}, observations are discarded so that the number of observations for each
#' level of the response variable is approximately the same. It applies to both Classification
#' and Regression with integer-valued response.
#'
#' @param reps Number of randomisations of the training/validation/testing sets to average over
#' in cross-validation.
#'
#' @param Q_vec Vector of numbers of PCA/PLS components to be tried in cross-validation.
#' The dedault is a vector of evenly spaced values (approximately, due to rounding) between
#' \code{2} and \code{min(80, N_train-1)}, where \code{N_train} is the number of observations
#' in the training subset.
#'
#' @param Q_len Length of \code{Q_vec} when the later is not supplied. The default is 30.
#'
#' @param Q_opt Optimal number of PCA/PLS components. If this is supplied, cross-validation for
#' the number of PCA/PLS components is bypassed.
#'
#' @param tau_Q_opt Threshold for choosing the optimal parameters.
#' If \code{tau_Q_opt=0}, then \code{Q} is the value that minimises RMSD (for Regression)
#' or maximises AUC (for Classification).
#' If \code{tau_Q_opt>0}, then \code{Q} is the smallest value which gives a RMSD/AUC within
#' a margin \code{tau_Q_opt} of the optimal RMSD/AUC.
#'
#' @param lam_cv_type Cross-validation strategy to be used when choosing the penalty
#' parameter \code{lambda}: ordinary cross-validation (\code{ocv}), generalised
#' cross-validation (\code{gcv}) or no cross-validation/penalisation (\code{n}).
#'
#' @param lam_vec Vector of penalty parameters to be tried in cross-validation. The
#' default is a set of 10 values between 0.001 and 20 on an exponential scale.
#'
#' @param split_size Either a vector of length 3 specifying the proportion of observations to be
#' assigned to the training, validation and testing subsets, in this order; or a scalar
#' specifying the proportion of observations to be assigned to the training subset, in which
#' case the validation and testing subsets are assigned a proportion \code{(1-split_size)/2}
#' of observations each.
#'
#' @param intercept Whether to include a model intercept (\code{TRUE}) or not (\code{FALSE}).
#'
#' @param weighted ...
#'
#' @param weights ...
#'
#' @param Bspline_dim The dimension of the cubic B-spline functional representation system.
#'
#' @param t_range ...
#'
#' @param verbose Whether to print a progress bar (\code{TRUE}) or not (\code{FALSE}).
#'
#' @param ll A list whose named elements are the parameters for this function.Provide either
#' the function parameters as usual, or this list, but not both.
#'
#' @details
#'
#' @return An object of class \code{fdaModel}, which is a list containing the results of the model.
#'
#' @references
#' P.M. Esperança, Thomas S. Churcher (2019). "Machine learning based
#' epidemiological vector control monitoring using functional data analysis
#' techniques for near-infrared spectral data".
#' arXiv.
#'
#' P.T. Reiss, R.T. Ogden (2007). "Functional Principal Component Regression
#' and Functional Partial Least Squares".
#' Journal of the American Statistical Association, 102(479), 984-996
#'
#' @seealso
#' \code{\link{fdaML_predict}}
#'
#' @examples
#'
#' @export
#'
fdaML_train <- function(X, y, Z=NULL, task, model=NULL, reduction, intercept=TRUE, smooth=T, balanced=F, reps=100, Q_vec=NULL, Q_len=NULL, Q_opt=NULL, tau_Q_opt=0.0, lam_cv_type="n", lam_vec=NULL, split_size=c(0.5,0.25,0.25), weighted=FALSE, weights=NULL, Bspline_dim=ncol(X), t_range=350:2500, verbose=TRUE, ll=NULL){


  # BASIC ENTRY CHECKS

  if(!is.null(ll)){
    if(nargs() > 1){
      stop("Provide the function parameters either as usual with named arguments or with a named list object 'll', but not both.")
    }
    {
      X <- ll$X
      y <- ll$y
      Z <- ll$Z
      task <- ll$task
      model <- ll$model
      reduction <- ll$reduction
      smooth <- ll$smooth
      intercept <- ll$intercept
      balanced <- ll$balanced
      lam_cv_type <- ll$lam_cv_type
      lam_vec <- ll$lam_vec
      reps <- ll$reps
      Q_len<-ll$Q_len
      Q_vec <- ll$Q_vec
      Q_opt <- ll$Q_opt
      split_size <- ll$split_size
      tau_Q_opt <- ll$tau_Q_opt
      weighted <- ll$weighted
      weights <- ll$weights
      Bspline_dim <- ll$Bspline_dim
      t_range <- ll$t_range
      verbose <- ll$verbose
    }
  }

  if( missing(X) )
    stop("Input 'X' is missing with no default.")

  if( is.list(X) )
    stop("Input 'X' must be a matrix, not a list.")

  if( missing(y) )
    stop("Input 'y' is missing with no default.")

  if( is.list(y) )
    stop("Input 'y' must be a vector, not a matrix/list.")

  if( is.factor(y) ){
    if( any(table(y) <= 1) ){
      stop("Input 'y' (a factor) contains 0 or 1 observations of at least one factor level.") }}

  if( is.list(Z) )
    stop("Input 'Z' must be a numeric matrix, not a list. Categorical variables must be pre-processed using one-hot-encoding.")

  if( missing(task) ){
    stop("Input 'task' is missing with no default.")
  }else if( !(task %in% c("regr","clas")) ){
    stop("Input 'task' must be one of {'regr','clas'}.") }

  if( is.null(model) ){
    if( task == "regr" ){
      model <- "lm"
    }else if( task == "clas" ){
      model <- "glm" }
  }else{
    if( task == "regr" & !(model  %in% c("lm")) ){
      stop("Input 'model' must be one of {'lm'} for regression tasks.")
    }else if( task == "clas" & !(model  %in% c("glm")) ){
      stop("Input 'model' must be one of {'glm'} for classification tasks.") }}

  if( missing(reduction) ){
    stop("Input 'reduction' is missing with no default.")
  }else if( !(reduction %in% c("n", "pca", "pls")) ){
    stop("Input 'reduction' must be one of {'pls','pca','n'}.") }

  if( !is.logical(intercept) )
    stop("Input 'intercept' must be logical (TRUE or FALSE).")

  if( !is.logical(smooth) )
    stop("Input 'smooth' must be logical (TRUE or FALSE).")

  if( !is.logical(balanced) )
    stop("Input 'balanced' must be logical (TRUE or FALSE).")
  if( balanced == TRUE & any(y %% 1 != 0) )
    stop("Input 'balanced' can not be TRUE when the response is not factor or integer valued.")

  if( !missing(reps) ){
    if( !(is.vector(reps) & is.numeric(reps) & length(reps)==1) )
      stop("Input 'reps' must be a numeric scalar.")
    if( reps < 0 | reps %% 1 != 0 )
      stop("Input 'reps' must be a non-negative integer/whole number.") }

  if( !(lam_cv_type %in% c("ocv", "gcv", "n")) )
    stop("Input 'lam_cv_type' must be one of {'ocv','gcv','n'}.")

  if( !is.null(lam_vec) ){
    if( !(is.vector(lam_vec) & is.numeric(lam_vec)) )
      stop("All elements of input 'lam_vec' must be a numeric vector.")
    if( any(lam_vec < 0) )
      stop("All elements of input 'lam_vec' must be non-negative.") }

  if( !missing(split_size) ){
    if( !(is.vector(split_size) & is.numeric(split_size)) ){
      stop("Input 'split_size' must be a numeric vector.") }
    if( !(all(split_size >= 0.05 ) & all(split_size <= 0.95 )) ){
      stop("All elements of input 'split_size' must be between 0.05 and 0.95.") }
    if( length(split_size) == 1 ){
      train_size <- split_size
      valid_size <- test_size <- train_size/2
    }else if( length(split_size) == 3 ){
      if( sum(split_size) == 1 ){
        train_size <- split_size[1]
        valid_size <- split_size[2]
        test_size  <- split_size[3]
      }else{
        stop("Elements of vector input 'split_size' must sum to 1.") }
    }else{
        stop("Input 'split_size' must be a vector of length 1 or 3.") }}

  if( !is.logical(weighted) )
    stop("Input 'weighted' must be logical (TRUE or FALSE).")

  if( !missing(weights) & !is.null(weights) ){
    if( !(is.vector(weights) & is.numeric(weights)) )
      stop("Input 'weights' must be a numeric vector")
    if( any(weights < 0) | any(weights %% 1 != 0) )
      stop("All elements of input 'weights' must be non-negative integer/whole number.") }

  if( !missing(Bspline_dim) ){
    if( !(is.vector(Bspline_dim) & is.numeric(Bspline_dim) & length(Bspline_dim)==1) )
      stop("Input 'Bspline_dim' must be a numeric scalar")
    if( Bspline_dim < 1 | Bspline_dim > ncol(X) )
      stop("Input 'Bspline_dim' must be greater than 1 and smaller than ncol(X).")
    if( Bspline_dim %% 1 != 0 )
      stop("Input 'Bspline_dim' must be an integer/whole number.") }

  if( task == "clas" & lam_cv_type == "gcv" ){
    glmcv <- TRUE
  }else{
    glmcv <- FALSE }

  if( !missing(t_range) ){
    if( !(is.vector(t_range) & is.numeric(t_range)) ){
      stop("Input 't_range' must be a numeric vector.") }
    if( length(t_range) != ncol(X) ){
      stop("Input 't_range' must be a vector of length equal to ncol(X).") }
    if( !all(t_range %% 1 == 0) ){
      stop("All elements of input 't_range' must be integer/whole numbers.") }
  }


  # FORCE UNIFORM/BALANCED RESPONSE BY DISCARDING OBSERVATIONS FROM THE MORE ABUNDANT CLASSES/GROUPS
  lvls <- sort(unique(y)); nlvls <- length(lvls)
  if(balanced){
    min_abundance <- min(sapply(1:nlvls, function(z){ sum(y==lvls[z]) }))
    each_loc <- lapply(1:nlvls, function(z){ which(y==lvls[z]) })
    smplId <- unlist(lapply(1:nlvls, function(z){ sample(each_loc[[z]], min_abundance, replace=F) }))
    X <- X[smplId,]
    y <- y[smplId] }


  # DEFINE TRAINING/VALIDATION/TESTING SUBSETS
  N <- nrow(X);   P <- ncol(X);   S <- ifelse(is.null(Z), 0, ncol(Z))
  splits <- dataSplit(y, train_size, valid_size, test_size, reps, balanced)
  id_train <- splits$id_train;   N_train <- length(id_train[[1]])
  id_valid <- splits$id_valid;   N_valid <- length(id_valid[[1]])
  id_test  <- splits$id_test;    N_test  <- length(id_test[[1]])
  #check: truehist(y[id_test[[1]]], nbins=10, p=F)


  # DEFINE NUMBER OF PCA/PLS COMPONENT GRID
  prioriMaxComponents <- 80
  actualMaxComponents <- min(prioriMaxComponents, N_train-1)
  if(!is.null(Q_opt)){
    if(Q_opt <= min(prioriMaxComponents, N_train-1)){
      Q_vec <- Q_max <- Q_opt
    }else{
      Q_vec <- Q_max <- Q_opt <- min(prioriMaxComponents, N_train-1)
      warning("Input 'Q_opt' is too large and was automatically reduced.")}
  }else{
    if(is.null(Q_vec)){
      if(is.null(Q_len)){
        Q_len <- 30 }
      Q_vec <- floor(seq(2, actualMaxComponents, len=Q_len))
    }else{
      Q_vec <- floor(Q_vec[Q_vec <= (N_train-1)])
      warning(paste0("Input 'Q_vec' contains elements that are too large and which will be ignored. 'Q_vec' is now of length ",length(Q_vec),".")) }}
  Q_vec <- floor(unique(Q_vec))  # rounding can cause doubling
  Q_len <- length(Q_vec)


  # DEFINE REGULARISATION PARAMETER GRID
  if(lam_cv_type != "n"){
    if( is.null(lam_vec) ){
      lam_min <- 0.001;   lam_max <- 20
      lam_vec <- exp(seq(log(lam_min), log(lam_max), len=10)) }}


  # DEFINE MODEL FAMILY
  if(task == "clas"){
    if(length(unique(y)) > 2){
      family <- "multinomial"
    }else if(length(unique(y)) == 2){
      family <- "binomial"
    }else{
      stop("response 'y' must have two or more levels or unique values.") }
  }else if(task == "regr"){
    family <- "gaussian" }


  # USE SMOOTHED SPECTRA
  if(smooth){
    X_original <- X
    X <- fdaSmooth(X)
  }else{
    X_original <- "same as X" }


  # USE FUNCTIONAL REPRESENTATION OF SPECTRA
  if(reduction == "n"){
    breaks <- seq(t_range[1], t_range[length(t_range)], len = P - 2)
  }else if(reduction == "pca" | reduction == "pls"){
    breaks <- seq(t_range[1], t_range[length(t_range)], len = Bspline_dim-2) }
  B    <- bsplineS(t_range, breaks, norder=4, returnMatrix=F)
  B_d2 <- bsplineS(t_range, breaks, norder=4, returnMatrix=F, nderiv=2) # 2nd derivative
  PtP <- t(B_d2) %*% B_d2
  XB <- X %*% t(t(B) / colSums(B))
  # ob <- 1; plot(t_range, X[ob,], type="l"); lines(breaks, XB[ob,-c(1,ncol(XB))], col=2)
  # matplot(t(XB), type="l", col=rgb(0_6,0_6,0_6,0_8)); lines(colMeans(XB[y==1,]), col="red"); lines(colMeans(XB[y==0,]), col="blue")


  # INITIALISATIONS
  y_testpred_opt  <- matrix(0, N_test, reps)
  y_trainpred_opt <- matrix(0, N_train, reps)
  y_validpred_opt <- matrix(0, N_valid, reps)
  resid_test  <- matrix(0, N_test, reps)
  resid_train <- matrix(0, N_train, reps)
  resid_valid <- matrix(0, N_valid, reps)
  err_train <- err_valid <- err_test <- rep(NA,reps)
  perf_cv <- perf_cv_dichotomised <- lam_cv <- messed_up_lambdas <- matrix(NA, reps, Q_len)
  rownames(perf_cv) <- rownames(lam_cv) <- rownames(messed_up_lambdas) <- paste("r",1:reps,sep="")
  colnames(perf_cv) <- colnames(lam_cv) <- colnames(messed_up_lambdas) <- paste("k",1:Q_len,sep="")
  perf_cv_ALT <- array(NA, dim=c(reps, Q_len, nlvls), dimnames=list(paste0("r",1:reps), paste0("k",1:Q_len), paste0("l",lvls)))
  AUC_opt <- matrix(0, reps, 3)
  colnames(AUC_opt) <- c("train","valid","test")
  classProbs <- ROC <- ROC_opt <- ROC_to_plot_binom <- ROC_to_plot_multinom <- mOpt <- D_opt <- U_opt <- cenX <- cenY <- scaX <- scaY <- XBD_train_opt <- XBD_valid_opt <- XBD_test_opt <- penmat_opt <- list()
  if(family == "multinomial"){ for(q in 1:nlvls){ ROC_to_plot_multinom[[q]] <- list() }; names(ROC_to_plot_multinom) <- lvls }
  bias_test <- matrix(0, nrow=reps, ncol=nlvls)
  colnames(bias_test) <- paste0("age",lvls)
  linear_predictor <- matrix(0, N_test, reps)
  alpha <- rep(0, reps)
  beta <- matrix(0, ncol(XB), reps)
  gamma <- if(is.null(Z)){ NULL }else{ matrix(0, S, reps) }
  elnet <- 0  # 0=ridge, 1=lasso
  Z_pvalue <- matrix(0, reps, S)
  deviance_opt <- matrix(0, N_train, reps)


  #===============================================#
  # TRAINING and CROSS-VALIDATION for (Q,lambda)  #
  #===============================================#

  if(verbose){
    cat("","TRAINING and CROSS-VALIDATION:", sep="\n")
    pb = txtProgressBar(min = 0, max = reps, initial = 0, style = 3)
  }

  for(rr in 1:reps){


    # initialisations
    ROC[[rr]] <- list()
    wgts <- GET_weights(weighted, weights, yy = y, id=id_train[[rr]], uniq = lvls)


    # data transformations
    # X (spectra)
    cenX[[rr]] <- colMeans(XB[id_train[[rr]],])
    scaX[[rr]] <- apply(XB[id_train[[rr]],], 2, sd)
    XB_train <- scale(XB[id_train[[rr]],], center = cenX[[rr]], scale=F)
    XB_valid <- scale(XB[id_valid[[rr]],], center = cenX[[rr]], scale=F)
    XB_test  <- scale(XB[id_test[[rr]],],  center = cenX[[rr]], scale=F)
    # y (response)
    if(task == "clas"){
      cenY[[rr]] <- NA
      scaY[[rr]] <- NA
      y_train  <- y[id_train[[rr]]]
      y_valid  <- y[id_valid[[rr]]]
      y_test   <- y[id_test[[rr]]]
    }else{
      cenY[[rr]] <- mean(y[id_train[[rr]]])
      scaY[[rr]] <- sd(y[id_train[[rr]]])
      y_train  <- scale(y[id_train[[rr]]],   center = cenY[[rr]], scale=F)
      y_valid  <- scale(y[id_valid[[rr]]],   center = cenY[[rr]], scale=F)
      y_test   <- scale(y[id_test[[rr]]],    center = cenY[[rr]], scale=F)
    }
    # Z (exogenous variables)
    if(!is.null(Z)){
      Z_train <- Z[id_train[[rr]],,drop=F]
      Z_valid <- Z[id_valid[[rr]],,drop=F]
      Z_test  <- Z[id_test[[rr]] ,,drop=F]
    }


    # dimension reduction
    if(reduction == "pca"){
      svdX <- svd(XB_train)
      D <- svdX$v                   #[loadings] = $rotation in prcomp()  #==#OR#==#  eigX <- eigen(t(XB_train) %*% XB_train); D <- eigX$vectors
      XBD_train <-  XB_train %*% D    #[scores] = XB_train %*% D = $x in prcomp()
    }else if(reduction == "pls"){
      plsX <- plsr(c(y_train) ~ -1 + XB_train, ncomp=N_train-1, method="simpls", validation="none")
      D <- plsX$projection      # [loadings][help(simpls_fit): "the projection matrix used to convert X to scores"]
      XBD_train <- plsX$scores  # = XB_train %*% D  [scores]
      # m1 <- plsr(c(y_train) ~ -1 + XB_train, ncomp=N_train-1, validation="LOO", method="simpls")
      # plot(RMSEP(m1), legendpos = "topright")
      # plot(m1, ncomp = 120, asp = 1, line = TRUE)
    }else if(reduction == "none"){
      stop("*---* need to code D downstream when no reduction is applied *---*")
      D <- diag(ncol(XB_train))
      XBD_train <- XB_train
      XBD_test  <- XB_test
    }


    if(reduction == "none"){
      # WE CAN SKIP THIS STEP AND GO STRAIGHT TO TESTING UNLESS WE NEED CV ON LAMBDA
      #mQ <- taskEstimation(task = task, XX = XBD_train, yy = y_train, optionEst = list(fam=family, wgts=wgts, test=F, elnet=elnet))
    }else{
      for(qq in 1:Q_len){
        Qcomp <- ifelse(is.null(Q_opt), Q_vec[qq], Q_opt)
        D_Q <- D[,1:Qcomp,drop=F]
        #wgts <- if(!is.null(weighted)){ if(isTRUE(weighted)){ FIND_WEIGHTS(as_vector(y_train),bias) }else{ NULL }}
        XBD_trainQ <- XBD_train[,1:Qcomp]
        XBD_validQ <- XB_valid %*% D_Q
        #!!!!!
        if(!is.null(Z)){
          XBD_trainQ <- cbind(Z_train, XBD_trainQ)
          XBD_validQ <- cbind(Z_valid, XBD_validQ)
        }
        #!!!!!
        penmat <- drop(t(D_Q) %*% PtP %*% D_Q)
        optionEst = list(intercept=intercept, S=S, task=task, model=model, fam=family, lam_cv_type=lam_cv_type, wgts=wgts, penmat=penmat, lam_vec=lam_vec, glmcv=glmcv, test=F, elnet=elnet)
        if((model == "glm") & (lam_cv_type != "n")){ optionEst$penmat <- penmat }   #optionEst$D <- D_Q; optionEst$PtP <- PtP}
        lam_cv[rr,qq] <- GET_lambda(XBD_trainQ, XBD_validQ, y_train, y_valid, optionEst)
        optionEst$lam <- lam_cv[rr,qq]
        mQ <- fdaEstimation(XX = XBD_trainQ, yy = y_train, optionEst = optionEst)
        optionPred <- list(intercept=intercept, task=task, model=model, fam=family, lam_cv_type=lam_cv_type, lam=lam_cv[rr,qq], predType="class")
        y_validpred <- fdaPrediction(m = mQ, newX = XBD_validQ, optionPred = optionPred)
        if(task == "clas"){
          if(family == "binomial"){
            perf_cv[rr,qq] <- GET_auc(y_pred = y_validpred, y_true = y_valid, f = family)
          }else if(family == "multinomial"){
            perf_cv[rr,qq] <- GET_auc(y_pred = y_validpred, y_true = y_valid, f = family)
            for(lvls_id in 1:nlvls){
              y_lvl <- as.numeric(y_train == lvls[lvls_id])
              m_lvl <- fdaEstimation(XX=XBD_trainQ, yy=y_lvl, optionEst = list(model=model, fam="binomial", wgts=wgts, penmat = penmat, lam = lam_cv[rr,qq], lam_cv_type=lam_cv_type, glmcv=glmcv, test=T, elnet=elnet))
              optionPred_multinom <- optionPred;  optionPred_multinom$fam <- "binomial"
              y_pred_lvl <- fdaPrediction(m = m_lvl, newX = XBD_validQ, optionPred = optionPred_multinom)
              pred_lvl <- prediction(y_pred_lvl, as.numeric(y_valid == lvls[lvls_id]))
              perf_lvl <- performance(pred_lvl, "auc")
              perf_cv_ALT[rr,qq,paste0("l",lvls[lvls_id])] <- unlist(slot(performance(pred_lvl, "auc"), "y.values"))
            }
          }
        }else if(task == "regr"){
          perf_cv[rr,qq] <- GET_rmsd(y_valid, y_validpred)
        }
        if(family == "multinomial" & lam_cv_type != "n"){ messed_up_lambdas[rr,qq] <- mQ$lambda }
      }#qq
    }#reduction
    if(verbose){
      setTxtProgressBar(pb,rr)
    }
  }#rr

  if((task == "clas") & (family == "multinomial")){
    for(rr in 1:reps){
      for(qq in 1:Q_len){
        perf_cv_dichotomised[rr,qq] <- mean(perf_cv_ALT[rr,qq,])
  }}}


  #===============================================#
  #        TESTING for optimal (Q,lambda)         #
  #===============================================#


  if(is.null(Q_opt)){
    Q_opt_id <- GET_IdxQ_opt(perf_cv, model, threshold=tau_Q_opt)
    Q_opt <- Q_vec[Q_opt_id]
    Q_max_id <- GET_IdxQ_opt(perf_cv, model, threshold=0)
    Q_max <- Q_vec[Q_max_id]
  }else{
    Q_opt_id <- 1
  }
  lam_opt <- mean(lam_cv[,Q_opt_id])
  PCA_variationExplained <- matrix(NA, Q_opt, reps); beta_reduced <- matrix(NA, Q_opt, reps)

  # coefficients' positions
  coef_pos <- matrix(NA,3,2); rownames(coef_pos) <- c("alpha", "gamma", "beta"); colnames(coef_pos) <- c("from", "to")
  if(intercept==T & S != 0){
    coef_pos["alpha",] <- c(1,1)
    coef_pos["gamma",] <- c(2,S+1)
    coef_pos["beta", ] <- c(S+2,1+S+Q_opt)
  }else if(intercept==T & S == 0){
    coef_pos["alpha",] <- c(1,1)
    coef_pos["gamma",] <- c(NA,NA)
    coef_pos["beta", ] <- c(2,1+Q_opt)
  }else if(intercept==F & S != 0){
    coef_pos["alpha",] <- c(NA,NA)
    coef_pos["gamma",] <- c(1,S)
    coef_pos["beta", ] <- c(S+1,S+Q_opt)
  }else if(intercept==F & S == 0){
    coef_pos["alpha",] <- c(NA,NA)
    coef_pos["gamma",] <- c(NA,NA)
    coef_pos["beta", ] <- c(1,Q_opt)
  }

  if(verbose){
    cat("","TESTING:", sep="\n")
    pb = txtProgressBar(min = 0, max = reps, initial = 0, style = 3)
  }

  for(rr in 1:reps){
    ROC_opt[[rr]] <- list()
    XB_train_opt <- scale(XB[id_train[[rr]],], center=cenX[[rr]], scale=F)
    y_train_opt <- if(task == "clas"){ y[id_train[[rr]]] }else{ scale(y[id_train[[rr]]], center=cenY[[rr]], scale=F) }
    if(reduction == "none"){
      ;
    }else if(reduction == "pca"){
      svdX_opt <- svd(XB_train_opt)
      D_opt[[rr]] <- svdX_opt$v[,1:Q_opt]
      XBD_train_opt[[rr]] <- XB_train_opt %*% D_opt[[rr]]
      PCA_variationExplained[,rr] <- svdX_opt$d[1:Q_opt] / sum(svdX_opt$d[1:Q_opt])
    }else if(reduction == "pls"){
      plsX <- plsr(c(y_train_opt) ~ -1 + XB_train_opt, ncomp=N_train-1, method="simpls", validation="none")
      D_opt[[rr]] <- plsX$projection[,1:Q_opt]
      XBD_train_opt[[rr]] <- plsX$scores[,1:Q_opt]
    }
    XBD_valid_opt[[rr]] <- scale(XB[id_valid[[rr]],], cen=cenX[[rr]],scale=F)  %*% D_opt[[rr]]
    XBD_test_opt[[rr]]  <- scale(XB[id_test[[rr]],],  cen=cenX[[rr]],scale=F)  %*% D_opt[[rr]]
    # add exogenous variables
    if(!is.null(Z)){
      XBD_train_opt[[rr]] <- cbind(Z[id_train[[rr]],,drop=F], XBD_train_opt[[rr]])
      XBD_valid_opt[[rr]] <- cbind(Z[id_valid[[rr]],,drop=F], XBD_valid_opt[[rr]])
      XBD_test_opt[[rr]]  <- cbind(Z[id_test[[rr]],, drop=F], XBD_test_opt[[rr]])
    }
    penmat_opt[[rr]] <- t(D_opt[[rr]]) %*% PtP %*% D_opt[[rr]]
    mOpt[[rr]] <- fdaEstimation(XX=XBD_train_opt[[rr]], yy=y_train_opt, optionEst = list(intercept=intercept, S=S, lam_cv_type=lam_cv_type, task=task, model=model, fam=family, wgts=wgts, penmat = penmat_opt[[rr]], lam = lam_opt, glmcv=glmcv, test=T, elnet=elnet))  #coeffs_opt <- solve( StS + lam_cv[rr,which(Q_opt==Q_vec)] * t(D_opt) %*% PtP %*% D_opt ) %*% t(XBD_train_opt[[rr]]) %*% y_train_opt
    if(!is.null(Z)){# p-values for exogenous variables; glmnet() doesn't give std errors by default, thus using glm() here
      Z_pvalue[rr,] <- summary(glm(as.numeric(y_train_opt) ~ 1 + XBD_train_opt[[rr]]))$coefficients[coef_pos["gamma","from"]:coef_pos["gamma","to"], "Pr(>|t|)"]
    }
    optionPred = list(intercept=intercept, lam_cv_type=lam_cv_type, task=task, model=model, fam=family, lam = lam_opt)
    y_trainpred_opt[,rr] <- fdaPrediction(m = mOpt[[rr]], newX = XBD_train_opt[[rr]], optionPred = optionPred)   #XBD_train_opt[[rr]] %*% coeffs_opt
    y_validpred_opt[,rr] <- fdaPrediction(m = mOpt[[rr]], newX = XBD_valid_opt[[rr]], optionPred = optionPred)   #scale(XB[id_valid[[rr]],], cen=cenX[[rr]],scale=F)  %*% D_opt %*% coeffs_opt
    y_testpred_opt[,rr]  <- fdaPrediction(m = mOpt[[rr]], newX = XBD_test_opt[[rr]],  optionPred = optionPred)   #scale(XB[id_test[[rr]],],  cen=cenX[[rr]],scale=F)  %*% D_opt %*% coeffs_opt
    if(task == "regr"){
      # regression coefficients in the space of original predictors
      if(!is.na(coef_pos["alpha",1]))
        alpha[rr] <- mOpt[[rr]]$estimate[1]
      if(!is.na(coef_pos["gamma",1]))
        stop("Gamma not coded for fLM.")
      if(!is.na(coef_pos["beta",1])){
        beta_reduced[,rr] <- mOpt[[rr]]$estimate[coef_pos["beta","from"]:coef_pos["beta","to"]]   # PCA/PLS-reduced space
        beta[,rr] <- c(D_opt[[rr]] %*%  beta_reduced[,rr])                                        # original space
      }
      linear_predictor[,rr] <- alpha[rr] + c(XBD_test_opt[[rr]]  %*% cbind(gamma[,rr], beta_reduced[,rr]))  # 'XBD_test_opt[[rr]]' already includes Z-variables
      # error
      resid_train[,rr] <- y_train_opt - y_trainpred_opt[,rr]
      err_train[rr]    <- GET_rmsd(scale(y[id_train[[rr]]],center=cenY[[rr]],scale=F), y_trainpred_opt[,rr])
      resid_valid[,rr] <- scale(y[id_valid[[rr]]], cen=cenY[[rr]],scale=F) - y_validpred_opt[,rr]
      err_valid[rr]    <- GET_rmsd(scale(y[id_valid[[rr]]],center=cenY[[rr]],scale=F), y_validpred_opt[,rr])
      resid_test[,rr]  <- scale(y[id_test[[rr]]],  cen=cenY[[rr]],scale=F) - y_testpred_opt[,rr]
      err_test[rr]     <- GET_rmsd(scale(y[id_test[[rr]]],center=cenY[[rr]],scale=F), y_testpred_opt[,rr])
      # bias
      for(uu in 1:nlvls){
        id <- y[id_test[[rr]]] == lvls[uu]
        bias_test[rr,uu] <- mean(y_testpred_opt[id,rr]+cenY[[rr]] - lvls[uu])
      }
    }else if(task == "clas"){
      if(family == "binomial"){
        # regression coefficients in the space of original predictors and linear predictor (to plot distributions and cutoff)
        if(model == "glm"){
          if(lam_cv_type == "n"){
            if(!is.na(coef_pos["alpha",1]))
              alpha[rr] <- coefficients(mOpt[[rr]])["(Intercept)",]
            if(!is.na(coef_pos["gamma",1]))
              gamma[,rr] <- as.matrix(coefficients(mOpt[[rr]])[coef_pos["gamma","from"]:coef_pos["gamma","to"]])
            if(!is.na(coef_pos["beta",1])){
              beta_reduced[,rr] <- coefficients(mOpt[[rr]])[coef_pos["beta","from"]:coef_pos["beta","to"],]  # PCA/PLS-reduced space
              beta[,rr] <- c(D_opt[[rr]] %*% beta_reduced[,rr])                                              # original space
            }
          }else{ # when using penalisation, the output from fdaEstimation() is different, so need to distinguish
            if(!is.na(coef_pos["alpha",1]))
              alpha[rr] <- mOpt[[rr]]$estimate[1]
            if(!is.na(coef_pos["gamma",1]))
              gamma[,rr] <- mOpt[[rr]]$estimate[coef_pos["gamma","from"]:coef_pos["gamma","to"]]
            if(!is.na(coef_pos["beta",1])){
              beta_reduced[,rr] <- mOpt[[rr]]$estimate[coef_pos["beta","from"]:coef_pos["beta","to"]]  # PCA/PLS-reduced space
              beta[,rr] <- c(D_opt[[rr]] %*% beta_reduced[,rr])                                        # original space
            }
          }
          linear_predictor[,rr] <- alpha[rr] + c(XBD_test_opt[[rr]]  %*% rbind(gamma[,rr,drop=F], beta_reduced[,rr,drop=F]))  # 'XBD_test_opt[[rr]]' already includes Z-variables
          deviance_opt[,rr] <- (-1)^as.numeric(1 - y_train_opt > y_trainpred_opt[,rr]) * sqrt(-2 * (y_train_opt * log(y_trainpred_opt[,rr]) + (1-y_train_opt) * log(1 - y_trainpred_opt[,rr])))
        }
        # auc/roc
        AUC_opt[rr,] <- c(GET_auc(y_pred = y_trainpred_opt[,rr], y_true = y[id_train[[rr]]], f = family),
                          GET_auc(y_pred = y_validpred_opt[,rr], y_true = y[id_valid[[rr]]], f = family),
                          GET_auc(y_pred = y_testpred_opt[,rr],  y_true = y[id_test[[rr]]],  f = family))
        #ROC_opt[[rr]] <- list(train = performance(prediction(y_trainpred_opt[,rr], y[id_train[[rr]]]), "tpr", "fpr"),
        #                      valid = performance(prediction(y_validpred_opt[,rr], y[id_valid[[rr]]]), "tpr", "fpr"),
        #                      test  = performance(prediction(y_testpred_opt[,rr],  y[id_test[[rr]]]),  "tpr", "fpr"))
        ROC_to_plot_binom$pred_test[[rr]] <- y_testpred_opt[,rr]
        ROC_to_plot_binom$labels[[rr]]    <- y[id_test[[rr]]]
      }else if(family == "multinomial"){
        # get class probabilities
        classProbs[[rr]] <- drop(fdaPrediction(m = mOpt[[rr]], newX = scale(XB[id_test[[rr]],],  cen=cenX[[rr]],scale=F)  %*% D_opt[[rr]], optionPred = list(fam=family, lam_cv_type=lam_cv_type, predType="probs")))
        # regression coefficients in the space of original predictors
        ;
        # auc/roc
        AUC_opt[rr,] <- c(GET_auc(y_pred = y_trainpred_opt[,rr], y_true = y[id_train[[rr]]], f = family),
                          GET_auc(y_pred = y_validpred_opt[,rr], y_true = y[id_valid[[rr]]], f = family),
                          GET_auc(y_pred = y_testpred_opt[,rr],  y_true = y[id_test[[rr]]],  f = family))
        for(lvls_id in 1:nlvls){# see https://stats.stackexchange.com/questions/2151/how-to-plot-roc-curves-in-multiclass-classification and especially https://stats.stackexchange.com/questions/71700/how-to-draw-roc-curve-with-three-response-variable/110550#110550
          y_opt_lvl <- as.numeric(y_train_opt == lvls[lvls_id])
          mOpt_lvl <- fdaEstimation(XX=XBD_train_opt[[rr]], yy=y_opt_lvl, optionEst = list(intercept=intercept, lam_cv_type=lam_cv_type, model=model, fam="binomial", wgts=wgts, penmat = penmat_opt[[rr]], lam = lam_opt, glmcv=glmcv, test=T, elnet=elnet))  #coeffs_opt <- solve( StS + lam_cv[rr,which(Q_opt==Q_vec)] * t(D_opt) %*% PtP %*% D_opt ) %*% t(XBD_train_opt[[rr]]) %*% y_train_opt
          y_testpred_opt_lvl <- fdaPrediction(m = mOpt_lvl, newX = scale(XB[id_test[[rr]],], cen=cenX[[rr]],scale=F) %*% D_opt[[rr]], optionPred = list(fam="binomial", lam_cv_type=lam_cv_type))   #scale(XB[id_test[[rr]],],  cen=cenX[[rr]],scale=F)  %*% D_opt %*% coeffs_opt
          ROC_to_plot_multinom[[as.character(lvls[lvls_id])]]$pred_test[[rr]] <- y_testpred_opt_lvl
          ROC_to_plot_multinom[[as.character(lvls[lvls_id])]]$labels[[rr]]    <- as.numeric(y[id_test[[rr]]] == lvls[lvls_id])
        }
      }#family
    }#task
    if(verbose){
      setTxtProgressBar(pb,rr)
    }
  }#rr
  if(verbose){
    cat("","", sep="\n")
  }


  # OUTPUT
  names(y_testpred_opt) <- names(y_validpred_opt) <- names(y_trainpred_opt) <- NULL;

  ret <- structure( list(perf_cv = perf_cv,   # AUC for 'glm', RMSD for 'lm'
                         perf_cv_dichotomised = perf_cv_dichotomised,
                         lam_vec = lam_vec,
                         lam_cv = lam_cv,
                         lam_opt = lam_opt,
                         m_opt = mOpt,
                         Q_len = Q_len,
                         Q_vec = Q_vec,
                         Q_opt = Q_opt,
                         Q_max = Q_max,
                         classProbs = classProbs,
                         AUC_opt = AUC_opt,
                         ROC_opt = ROC_opt,
                         ROC_to_plot_binom = ROC_to_plot_binom,
                         ROC_to_plot_multinom = ROC_to_plot_multinom,
                         PCA_variationExplained = PCA_variationExplained,
                         alpha = alpha,
                         gamma = gamma,
                         beta_reduced = beta_reduced,
                         beta = beta,
                         beta_roughness = apply(beta, 2, GET_roughness, std=F),
                         beta_roughnessStd = apply(beta, 2, GET_roughness, std=T),
                         beta_avg_roughness = GET_roughness(rowMeans(beta), std=F),
                         beta_avg_roughnessStd = GET_roughness(rowMeans(beta), std=T),
                         coef_pos = coef_pos,
                         id_train = id_train,
                         id_valid = id_valid,
                         id_test = id_test,
                         N_train = N_train,
                         N_valid = N_valid,
                         N_test = N_test,
                         X_original = X_original,
                         X = X,
                         y = y,
                         Z = Z,
                         S = S,
                         Z_pvalue = Z_pvalue,
                         deviance_opt = deviance_opt,
                         breaks = breaks,
                         B = B,
                         B_d2 = B_d2,
                         XB = XB,
                         D_opt = D_opt,
                         penmat_opt = penmat_opt,
                         XBD_train_opt = XBD_train_opt,
                         linear_predictor = linear_predictor,
                         cenX = cenX,
                         cenY = cenY,
                         reps = reps,
                         weights = weights,
                         y_trainpred_opt = y_trainpred_opt,
                         y_validpred_opt = y_validpred_opt,
                         y_testpred_opt  = y_testpred_opt,
                         err_train = err_train,
                         err_valid = err_valid,
                         err_test  = err_test,
                         bias_test = bias_test,
                         resid_train = resid_train,
                         resid_valid = resid_valid,
                         resid_test  = resid_test,
                         lam_cv_type = lam_cv_type,
                         task = task,
                         model = model,
                         family = family,
                         intercept = intercept,
                         t_range = t_range,
                         t_range_mod = seq(t_range[1], t_range[length(t_range)], len = Bspline_dim) # t_range_mod may be different from t_range, depending on the number of breaks used to create the B-spline                         messed_up_lambdas = messed_up_lambdas
                         ), class="fdaModel")

  return(ret)
}

