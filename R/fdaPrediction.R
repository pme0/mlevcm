#' @title  ???
#'
#' @description  ???.
#'
#' @param  ???
#'
#' @export
#'
fdaPrediction <- function(m, newX, optionPred){
  ##
  ##  m     :  method object (can be coefficients for lm or more complicated for rf)
  ##  newX  :  matrix of new observation to be predicted

  # LINEAR REGRESSION
  if(attr(m, "typeEstimation") == "FDA_lm.fit"){

    if(optionPred$intercept){
      newX <- cbind(1, newX)
    }
    out <- drop(newX %*% m$estimate)
    attr(out, "typePrediction") <- "FDA_lm.pred"

  # LOGISTIC REGRESSION
  }else if(attr(m, "typeEstimation") == "FDA_glm.fit"){

    if(optionPred$lam_cv_type == "n"){

      if(optionPred$fam == "binomial"){

        out <- predict(object = m, newx = newX, s = optionPred$lam, type = ifelse(optionPred$fam == "binomial", "response", "class"))

      }else if(optionPred$fam == "multinomial"){

        #if(optionPred$intercept){ newX <- cbind(rep(1,nrow(newX)), newX) }
        # if(optionPred$intercept){
        #   newX <- cbind(1, newX);  colnames(newX) <- c("(Intercept)", paste0("V", 1:(ncol(newX)-1)))
        # }else{
        #   colnames(newX) <- paste0("V", 1:ncol(newX))
        # }
        #
        newX <- as.data.frame(newX)
        if(!is.null(optionPred$predType)){
          #out <- predict(object = m, newx = newX, s = 0, type = optionPred$predType)
          out <- predict(m, newdata = newX, type = optionPred$predType)
        }else{
          #out <- predict(object = m, newx = newX, s = optionPred$lam, type = ifelse(optionPred$fam == "binomial", "response", "class"))
          out <- as.numeric(predict(m, newdata = newX, type = "class")) - 1
        }
      }

    }else if(optionPred$lam_cv_type %in% c("ocv","gcv")){

      if(optionPred$intercept){
        newX <- cbind(rep(1,nrow(newX)), newX)
      }
      out <- 1 / (1 + exp(-newX %*% m$estimate))

    }

    attr(out, "typePrediction") <- "FDA_glm.pred"
  }
  return(out)
}
