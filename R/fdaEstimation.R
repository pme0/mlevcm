#' @title Estimate parameters
#'
#' @description Produce parameter estimates for the statistical model.
#'
#' @param XX A numeric predictor matrix.
#'
#' @param yy A numeric response vector.
#'
#' @param optionsEst A list with parameters for estimation.
#'
#' @export
#'
fdaEstimation <- function(XX, yy, optionsEst=NULL){

  if(optionsEst$model == "lm"){
    if(length(optionsEst$wgts)!=1 & !all(optionsEst$wgts==optionsEst$wgts[1])){
      stop("Weighting not included in LossFunction_LM() yet!")
      W <- Matrix(diag(optionsEst$wgts), sparse=T)
      out <- solve( t(XX) %*% W %*% XX + optionsEst$lam * optionsEst$penmat ) %*% t(XX) %*% W %*% yy
    }
    if(optionsEst$intercept){
      XX <- cbind(1,XX) }
    out <- nlm(f = LossFunction_LM, p = rep(0, ncol(XX)), y_vec = yy, XBD_mat = XX, penmat = optionsEst$penmat, lam = optionsEst$lam, l1=ifelse(optionsEst$intercept,2,1), lv=ncol(XX))

    # if required, do cross-validation for lambda via cv.glmnet()
    #out <- glmnet(x = XX, y = as.numeric(yy), intercept = F, family = optionsEst$fam, weights = optionsEst$wgts, alpha = optionsEst$elnet, lambda = optionsEst$lam)
    # run model for optimal cross-validated lambda
    attr(out, "typeEstimation") <- "FDA_lm.fit"

  }else  if(optionsEst$model == "glm"){

    if(optionsEst$lam_cv_type %in% c("ocv", "gcv")){# fPGLM

      if(optionsEst$fam == "multinomial"){
        stop("fPGLM with multinomial response not coded yet in 'fdaEstimation()'.")
      }else{
        if(length(optionsEst$wgts)!=1 & !all(optionsEst$wgts==optionsEst$wgts[1])){
          stop("Weights not yet implemented for fGLM in 'LossFunction_GLM()'.")
        }
        if(optionsEst$intercept){
          XX <- cbind(1,XX) }
        out <- nlm(f = LossFunction_GLM, p = rep(0, ncol(XX)), y_vec = yy, XBD_mat = XX, penmat = optionsEst$penmat, lam = optionsEst$lam, l1 = ifelse(optionsEst$intercept,2,1), lv=ncol(XX))
      }

    }else if(optionsEst$lam_cv_type == "n"){

      if(optionsEst$fam == "multinomial"){
        #
        # https://czep.net/stat/mlelr.pdf
        #
        # if(optionsEst$intercept){
        #   XX <- cbind(1, XX);  colnames(XX) <- c("(Intercept)", paste0("V", 1:(ncol(XX)-1)))
        # }else{
        #   colnames(XX) <- paste0("V", 1:ncol(XX))
        # }
        #out <- glmnet(x = XX, y = yy, intercept = optionsEst$intercept, family = optionsEst$fam, weights = optionsEst$wgts, alpha = 0, lambda = 0)
        # out <- NA; aux <- T; lam <- 0
        # while(aux){
        #   try(out <- glmnet(x = XX, y = yy, intercept = optionsEst$intercept, family = optionsEst$fam, weights = optionsEst$wgts, alpha = 0, lambda = lam))
        #   lam <- nrow(XX) * (0.0001 + 10 * lam)
        #   if(!is.na(out)){ aux <- F }
        # }
        #
        if(optionsEst$intercept){
          out <- multinom(as.factor(yy) ~ .      , data = as.data.frame(XX), weights = optionsEst$wgts, maxit = 1e9, trace=F)
        }else{
          out <- multinom(as.factor(yy) ~ -1 + . , data = as.data.frame(XX), weights = optionsEst$wgts, maxit = 1e9, trace=F)
        }

      }else{
        #out <- glm(factor(yy) ~ XX, family = binomial(link = "logit"))
        out <- glmnet(x = XX, y = yy, intercept = optionsEst$intercept, family = optionsEst$fam, weights = optionsEst$wgts, alpha = 0, lambda = 0, nlambda=1)
      }

    }
    # if required, do cross-validation for lambda via cv.glmnet()
    # if(optionsEst$glmcv){
    #   cv <- cv.glmnet(x = XX, y = yy, intercept = optionsEst$intercept, family = optionsEst$fam, weights = optionsEst$wgts, alpha = optionsEst$elnet) # lambda=exp(seq(log(0.001), log(10), len=10))
    #   lam <- cv$lambda.min
    # }else if( !optionsEst$glmcv | optionsEst$test ){
    #   lam <- optionsEst$lam
    # }
    # run model for optimal cross-validated lambda

    attr(out, "typeEstimation") <- "FDA_glm.fit"
  }
  return(out)
}
