#' @title  ???
#'
#' @description  ???.
#'
#' @param  ???
#'
#' @export
#'
fdaEstimation <- function(XX, yy, optionEst=NULL){

  if(optionEst$model == "lm"){
    if(length(optionEst$wgts)!=1 & !all(optionEst$wgts==optionEst$wgts[1])){
      stop("Weighting not included in LossFunction_LM() yet!")
      W <- Matrix(diag(optionEst$wgts), sparse=T)
      out <- solve( t(XX) %*% W %*% XX + optionEst$lam * optionEst$penmat ) %*% t(XX) %*% W %*% yy
    }
    out <- nlm(f = LossFunction_LM, p = rep(0, ncol(XX) + optionEst$intercept), y_vec = yy, XBD_mat = XX, S=optionEst$S, intercept=optionEst$intercept, penmat = optionEst$penmat, lam = optionEst$lam)

    # if required, do cross-validation for lambda via cv.glmnet()
    #out <- glmnet(x = XX, y = as.numeric(yy), intercept = F, family = optionEst$fam, weights = optionEst$wgts, alpha = optionEst$elnet, lambda = optionEst$lam)
    # run model for optimal cross-validated lambda
    attr(out, "typeEstimation") <- "FDA_lm.fit"

  }else  if(optionEst$model == "glm"){

    if(optionEst$lam_cv_type %in% c("ocv", "gcv")){# fPGLM

      if(optionEst$fam == "multinomial"){
        stop("fPGLM with multinomial response not coded yet in 'fdaEstimation()'.")
      }else{
        if(length(optionEst$wgts)!=1 & !all(optionEst$wgts==optionEst$wgts[1])){
          stop("Weights not yet implemented for fGLM in 'LossFunction_GLM()'.")
        }
        out <- nlm(f = LossFunction_GLM, p = rep(0, ncol(XX) + optionEst$intercept), y_vec = yy, XBD_mat = XX, S=optionEst$S, intercept=optionEst$intercept, penmat = optionEst$penmat, lam = optionEst$lam)
      }

    }else if(optionEst$lam_cv_type == "n"){

      if(optionEst$fam == "multinomial"){
        #
        # https://czep.net/stat/mlelr.pdf
        #
        # if(optionEst$intercept){
        #   XX <- cbind(1, XX);  colnames(XX) <- c("(Intercept)", paste0("V", 1:(ncol(XX)-1)))
        # }else{
        #   colnames(XX) <- paste0("V", 1:ncol(XX))
        # }
        #out <- glmnet(x = XX, y = yy, intercept = optionEst$intercept, family = optionEst$fam, weights = optionEst$wgts, alpha = 0, lambda = 0)
        # out <- NA; aux <- T; lam <- 0
        # while(aux){
        #   try(out <- glmnet(x = XX, y = yy, intercept = optionEst$intercept, family = optionEst$fam, weights = optionEst$wgts, alpha = 0, lambda = lam))
        #   lam <- nrow(XX) * (0.0001 + 10 * lam)
        #   if(!is.na(out)){ aux <- F }
        # }
        #
        if(optionEst$intercept){
          out <- multinom(as.factor(yy) ~ .      , data = as.data.frame(XX), weights = optionEst$wgts, maxit = 1e9, trace=F)
        }else{
          out <- multinom(as.factor(yy) ~ -1 + . , data = as.data.frame(XX), weights = optionEst$wgts, maxit = 1e9, trace=F)
        }

      }else{
        #out <- glm(factor(yy) ~ XX, family = binomial(link = "logit"))
        out <- glmnet(x = XX, y = yy, intercept = optionEst$intercept, family = optionEst$fam, weights = optionEst$wgts, alpha = 0, lambda = 0, nlambda=1)
      }

    }
    # if required, do cross-validation for lambda via cv.glmnet()
    # if(optionEst$glmcv){
    #   cv <- cv.glmnet(x = XX, y = yy, intercept = optionEst$intercept, family = optionEst$fam, weights = optionEst$wgts, alpha = optionEst$elnet) # lambda=exp(seq(log(0.001), log(10), len=10))
    #   lam <- cv$lambda.min
    # }else if( !optionEst$glmcv | optionEst$test ){
    #   lam <- optionEst$lam
    # }
    # run model for optimal cross-validated lambda

    attr(out, "typeEstimation") <- "FDA_glm.fit"
  }
  return(out)
}
