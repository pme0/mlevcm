#' @title Predict responses
#'
#' @description Produce predictions for new observations.
#'
#' @param m An object of class \code{fdaModel} computed with \code{\link{fdaML_train}}.
#'
#' @param newX A numeric predictor matrix of new observations for which predictions
#' are sought.
#'
#' @param optionsPred A list with parameters for prediction.
#'
#' @export
#'
fdaPrediction <- function(m, newX, optionsPred){

  if(attr(m, "typeEstimation") == "FDA_lm.fit"){

    if(optionsPred$intercept){
      newX <- cbind(1, newX)
    }
    out <- drop(newX %*% m$estimate)
    attr(out, "typePrediction") <- "FDA_lm.pred"

  }else if(attr(m, "typeEstimation") == "FDA_glm.fit"){

    if(optionsPred$lam_cv_type == "n"){
      if(optionsPred$fam == "binomial"){
        out <- predict(object = m, newx = newX, s = optionsPred$lam, type = ifelse(optionsPred$fam == "binomial", "response", "class"))
      }else if(optionsPred$fam == "multinomial"){
        #if(optionsPred$intercept){ newX <- cbind(rep(1,nrow(newX)), newX) }
        # if(optionsPred$intercept){
        #   newX <- cbind(1, newX);  colnames(newX) <- c("(Intercept)", paste0("V", 1:(ncol(newX)-1)))
        # }else{
        #   colnames(newX) <- paste0("V", 1:ncol(newX))
        # }
        #
        newX <- as.data.frame(newX)
        if(!is.null(optionsPred$predType)){
          #out <- predict(object = m, newx = newX, s = 0, type = optionsPred$predType)
          out <- predict(m, newdata = newX, type = optionsPred$predType)
        }else{
          #out <- predict(object = m, newx = newX, s = optionsPred$lam, type = ifelse(optionsPred$fam == "binomial", "response", "class"))
          out <- as.numeric(predict(m, newdata = newX, type = "class")) - 1
        }
      }

    }else if(optionsPred$lam_cv_type %in% c("ocv","gcv")){

      if(optionsPred$intercept){
        newX <- cbind(rep(1,nrow(newX)), newX)
      }
      out <- 1 / (1 + exp(-newX %*% m$estimate))

    }

    attr(out, "typePrediction") <- "FDA_glm.pred"
  }
  return(out)
}
