#' @title  Penalised GLM objective function
#'
#' @description Objective function for the penalised GLM, to be optimised using stats::nlm().
#'
#' @param pars Parameter vector (\code{alpha}, \code{omega}), with \code{alpha} scalar
#' and \code{omega} a \code{Q}-dimensional vector.
#'
#' @param y_vec Response vector.
#'
#' @param XBD_mat Functional predictor matrix after dimension reduction, plus an intercept
#' if appropriate.
#'
#' @param penmat Penalty matrix.
#'
#' @param lam Penalty parameter.
#'
#' @param l1 First column of functional predictor, equal to \code{1} if
#' \code{intercept=FALSE} and \code{2} if \code{intercept=TRUE}.
#'
#' @param lv Last column of functional predictor, equal to\code{col(XBD_mat)}.
#'
#' @details
#'
#' @return The penalised GLM objective function evaluated at \code{pars}.
#'
#' @references
#'
#' @seealso
#'
#' @examples
#'
#' @export
#'
LossFunction_GLM <- function(pars, y_vec, XBD_mat, penmat, lam, l1, lv){
  if(lam == 0){
    return( sum(log(1 + exp((-1)^y_vec * XBD_mat %*% pars))) )
  }else{
    return( sum(log(1 + exp((-1)^y_vec * XBD_mat %*% pars))) + length(y_vec) * lam * (t(pars[l1:lv]) %*% penmat %*% pars[l1:lv]) )
  }
}



#' @title Penalised LM objective function
#'
#' @description Objective function for the penalised LM, to be optimised using stats::nlm()..
#'
#' @param pars Parameter vector (\code{alpha}, \code{omega}), with \code{alpha} scalar
#' and \code{omega} a \code{Q}-dimensional vector.
#'
#' @param y_vec Response vector.
#'
#' @param XBD_mat Functional predictor matrix after dimension reduction, plus an intercept
#' if appropriate.
#'
#' @param penmat Penalty matrix.
#'
#' @param lam Penalty parameter.
#'
#' @param l1 First column of functional predictor (\code{1} if \code{intercept=FALSE},
#' \code{2} if \code{intercept=TRUE}.
#'
#' @param lv Last column of functional predictor (\code{col(XBD_mat)}).
#'
#' @details
#'
#' @return The penalised LM objective function evaluated at \code{pars}.
#'
#' @references
#'
#' @seealso
#'
#' @examples
#'
#' @export
#'
LossFunction_LM <- function(pars, y_vec, XBD_mat, penmat, lam, l1, lv){
  if(lam == 0){
    return( t(y_vec - XBD_mat %*% pars) %*% (y_vec - XBD_mat %*% pars) )
  }else{
    return( t(y_vec - XBD_mat %*% pars) %*% (y_vec - XBD_mat %*% pars) + length(y_vec) * lam * (t(pars[l1:lv]) %*% penmat %*% pars[l1:lv]) )
  }
}



#' @title Compute AUC
#'
#' @description Compute the Area Under the ROC Curve.
#'
#' @param y_pred Predicted response.
#'
#' @param y_true True response.
#'
#' @param fam Model family.
#'
#' @details
#'
#' @return Area Under the ROC Curve (AUC).
#'
#' @references
#'
#' @seealso
#'
#' @examples
#'
#' @export
#'
GET_auc <- function(y_pred, y_true, fam){
  if(fam == "binomial"){
    pred <- prediction(predictions = as.numeric(y_pred), labels = y_true)
    return( as.numeric(performance(pred,"auc")@y.values) )
  }else if(fam =="multinomial"){
    pred <- multiclass.roc(y_true, as.numeric(y_pred))
    return( pred$auc[1] )
  }
}



#' @title Compute RMSD
#'
#' @description Compute the Root Mean Squared Deviation.
#'
#' @param y_true True response.
#'
#' @param y_pred Predicted response.
#'
#' @details
#'
#' @return Root Mean Squared Deviation.
#'
#' @references
#'
#' @seealso
#'
#' @examples
#'
#' @export
#'
GET_rmsd <- function(y_true, y_pred){
  return( sqrt((1/length(y_pred))*sum((y_pred-y_true)^2)) )
}



#' @title Optimal \code{Q} index
#'
#' @description Compute index of best performing number of PCA/PLS components \code{Q}
#' from the cross-validation results.
#'
#' @param perf Performance matrix.
#'
#' @param model Model type.
#'
#' @param threshold Threshold for choosing the optimal parameters.
#' If \code{threshold=0}, then \code{Q} is the value that minimises RMSD (for Regression)
#' or maximises AUC (for Classification).
#' If \code{threshold>0}, then \code{Q} is the smallest value which gives a RMSD/AUC within
#' a margin \code{threshold} of the optimal RMSD/AUC.
#'
#' @details
#'
#' @return The index of the best performing cross-validation parameters.
#'
#' @references
#'
#' @seealso
#'
#' @examples
#'
#' @export
#'
GET_IdxQ_opt <- function(perf, model, threshold){

  if(sum(is.na(perf)) != 0){
    perf <- perf[,1:(min(which(is.na(perf), arr.ind=T)[,2])-1)] }

  if(threshold == 0){

    if(model=="glm"){
      out <- which.max(colMeans(perf, na.rm=T))
    }else if(model=="lm"){
      out <- which.min(colMeans(perf, na.rm=T)) }

  }else{

    if(model=="glm"){
      aux <- colMeans(perf, na.rm=T) - max(colMeans(perf), na.rm=T)
    }else if(model=="lm"){
      aux <- colMeans(perf, na.rm=T) - min(colMeans(perf), na.rm=T) }
    out <- which(abs(aux) < threshold)[1]

  }

  if(is.na(out)){ stop("Something wrong in GET_IdxQ_opt().")}
  return(out)
}



#' @title Compute penalty parameter
#'
#' @description Compute penalty parameter \code{lambda} by either ordinary
#' cross-validation (OCV) or generalised cross-validation (GCV).
#'
#' @param XBD_trainQ Training predictor matrix.
#'
#' @param XBD_validQ Validation predictor matrix.
#'
#' @param y_train Training response vector
#'
#' @param y_valid Validation response matrix.
#'
#' @param optionsEst List of options for estimation.
#'
#' @details
#'
#' @return The value of the penalty parameter \code{lambda}.
#'
#' @references
#'
#' @seealso
#'
#' @examples
#'
#' @export
#'
GET_lambda <- function(XBD_trainQ, XBD_validQ, y_train, y_valid, optionsEst){
  if(optionsEst$lam_cv_type == "n"){
    lam <- 0
  }else if(optionsEst$lam_cv_type=="ocv"){
    lam <- OCV(XX=XBD_trainQ, yy=y_train, newX=XBD_validQ, newY=y_valid, optionsEst=optionsEst)
    # ALTERNATIVELY USE cv.glmnet():
    #cv <- cv.glmnet(x = XBD_trainQ, y = y_train, intercept = F, lambda = lam_vec, family = fam, weights = wgts, alpha = elnet)
    #lam <- cv$lambda.min
  }else if(optionsEst$lam_cv_type=="gcv"){
    stop("Generalised cross validation not fully implemented. Choose lam_cv_type='ocv' instead.")
  }
  if(is.nan(lam)){ lam <- 0 }
  return(lam)
}



#' @title Ordinary Cross-Validation for penalty parameter
#'
#' @description Ordinary Cross-Validation method for optimising the penalty parameter.
#'
#' @param XX Training predictor matrix.
#'
#' @param newX Validation predictor matrix.
#'
#' @param yy Training response vector.
#'
#' @param newY Validation response matrix.
#'
#' @param optionsEst List of options for estimation.
#'
#' @details
#'
#' @return The value of lambda which minimises the Ordinary Cross-Validation criterion,
#' that is, minimises the RMSD (for Regression) or maximises the AUC (for Classification).
#'
#' @references
#'
#' @seealso
#'
#' @examples
#'
#' @export
#'
OCV <- function(XX, yy, newX, newY, optionsEst){
  q <- length(optionsEst$lam_vec)
  res <- rep(0,q)
  optionsEst_ocv <- optionsEst
  for(l in 1:q){
    optionsEst_ocv$lam <- optionsEst_ocv$lam_vec[l]
    m <- fdaEstimation(XX = XX, yy = yy, optionsEst = optionsEst_ocv)
    newYpred <- fdaPrediction(m = m, newX = newX, optionsPred = list(lam_cv_type = optionsEst$lam_cv_type, intercept = optionsEst$intercept))
    if(optionsEst$model == "lm"){
      res[l] <- GET_rmsd(newY, newYpred)
    }else if(optionsEst$model == "glm"){
      res[l] <- GET_auc(y_pred = newYpred, y_true = newY, fam = optionsEst$fam)
    }
  }
  # THIS DOES TWO ROUNDS: THE 1ST WITH A LARGE AND SPARSE GRID OF LAMBDAS, THE 2ND WITH A FOCUSED GRID AROUND THE BEST LAMBDA FROM THE 1ST ROUND:
  # for(h in 1:2){
  #   if(h!=1){ lam_vec <- seq(lam_vec[max(which.min(res)-3,1)], lam_vec[min(which.min(res)+3,q)], len=q) }
  #   res <- rep(0,q)
  #   for(l in 1:length(lam_vec)){
  #     m <- fdaEstimation(XX = XX, yy = yy, optionsEst = list(model=model, penmat=penmat, lam=lam_vec[l], wgts=wgts))
  #     newYpred <- fdaPrediction(m = m, newX = newX, optionsPred = list())
  #     res[l] <- GET_rmsd(newY, newYpred)
  #   }
  #   if( which.min(res) == q ){ warning("interval for penalty parameter lambda in OCV() should be widened.") }
  # }
  if(optionsEst$model == "lm"){
    return( optionsEst_ocv$lam_vec[which.min(res)] )
  }else if(optionsEst$model == "glm"){
    return( optionsEst_ocv$lam_vec[which.max(res)] )
  }
}



#' @title Generalised Cross-Validation for penalty parameter
#'
#' @description Generalised Cross-Validation method for optimising the penalty parameter.
#'
#' @param XX Training predictor matrix.
#'
#' @param yy Training response vector.
#'
#' @param penmat Penalty matrix.
#'
#' @details
#'
#' @return The value of lambda which minimises the Generalised Cross-Validation criterion.
#'
#' @references
#'
#' @seealso
#'
#' @examples
#'
#' @export
#'
GCV <- function(XX, yy, penmat){
  XtX <- t(XX) %*% XX
  q <- 20;   lamVec <- exp(seq(log(0.001), log(20), len=q));   I <- diag(nrow(XX))
  for(h in 1:2){
    if(h != 1){ lamVec <- seq(lamVec[max(which.min(gcv)-3,1)], lamVec[min(which.min(gcv)+3,q)], len=q) }
    gcv <- rep(0,q)
    for(l in 1:length(lamVec)){
      H <- XX %*% solve( XtX + lamVec[l] * penmat ) %*% t(XX)
      ImH <- I - H
      gcv[l] <- (1/length(yy)) * sum( (ImH %*% matrix(yy, byrow=F))^2 ) / ( (1/length(yy)) * sum(diag(ImH)) )^2
    }
    if( which.min(gcv) == q ){ warning("interval for penalty parameter lambda in GCV() should be widened.") }
  }
  return( lamVec[which.min(gcv)] )
}



#' @title Compute weights
#'
#' @description Compute the weights of each observation for parameter estimation.
#'
#' @param weighted Whether to use observation weights (\code{TRUE}) or not (\code{FALSE}).
#' If \code{weighted=FALSE}, uniform weights are used.
#'
#' @param weights Externally provided weights.
#'
#' @param yy Response vector.
#'
#' @param id Randomisation IDs.
#'
#' @param uniq Unique levels in the response vector.
#'
#' @details Compute weights invesely proportional to the number of observations
#' in each level of the response variable, if a vector of weights is not provided.
#' It applies to both Regression and Classification tasks.
#'
#' @export
#'
GET_weights <- function(weighted, weights=NULL, yy, id, uniq){
  yy <- yy[id]
  ll <- length(yy)
  if(weighted){
    if(!is.null(weights)){
      wgts <- weights[id]
    }else{
      wVec <- ceiling(abs(uniq-6.5)^4)
      tab <- table(yy)
      wgts <- sapply(1:ll, function(z){ wVec[which(uniq==yy[z])] })
    }
  }else{
    wgts <- rep(1,ll)
  }
  names(wgts) <- NULL
  return(wgts/sum(wgts))
}



#' @title Compute optimal cutoff for classification
#'
#' @description Compute the optimal cutoff for Classification
#'
#' @param pred An object of class \code{prediction} generated with \code{\link[ROCR]{prediction}}.
#'
#' @param perf An object of class \code{performance} generated with \code{\link[ROCR]{performance}}.
#'
#' @param costs A named list with elements \code{cost_fp} and \code{cost_fn} giving the weights
#' associated with false positives and false negatives, respectively.
#'
#' @export
#'
GET_optROCcutoff = function(pred, perf, costs=NULL){
  # see https://hopstat.wordpress.com/2014/12/19/a-small-introduction-to-the-rocr-package
  if( is.null(costs) ){
    cut.ind = mapply(FUN=function(x, y, p){
      d = (x - 0)^2 + (y-1)^2
      ind = which(d == min(d))
      c(sensitivity = y[[ind]], specificity = 1-x[[ind]],
        cutoff = p[[ind]])
    }, perf@x.values, perf@y.values, pred@cutoffs)
  }else{
    if( !is.list(costs) | length(costs) != 2 | sum(unlist(costs)) != 1 ){
      stop("Input 'costs' in GET_optROCcutoff() must be a list of length 2 with named elements {'cost_fp','cost_fn'} which must sum up to 1.")
    }else{
      stop("Different costs for FN/FP via GET_optROCcutoff() not available yet.")
    }
  }
}



#' @title Compute roughness of coefficient function
#'
#' @description Compute roughness of a coefficient function as given by the
#' integrated squared second derivative.
#'
#' @param x Coefficient function.
#'
#' @param std Whether to standadise the coefficient function by its maximum
#' value (\code{TRUE}) or not (\code{FALSE}).
#'
#' @export
#'
GET_roughness <- function(x, std=F){
  if(std){ x <- x/max(x) }
  deriv2 <- diff(x, 1, 2)
  return(sum(deriv2^2, na.rm=T))
}



#' @title Confusion matrix plot
#'
#' @description Plots the confusion matrix.
#'
#' @param d Confusion matrix.
#'
#' @param fam Model family.
#'
#' @param avgError Average test error.
#'
#' @param binom_labels Labels for the confusion matrix in the case of binomial
#' Classification.
#'
#' @param multinom_labels Labels for the confusion matrix in the case of
#' multinomial Classification. Set to the column names of \code{d} if
#' no value is provided.
#'
#' @export
#'
confusion_plot <- function(d, fam, avgError, binom_labels=c(y0="observed class: 0", y1="observed class: 1", yhat0="predicted class: 0", yhat1="predicted class: 1"), multinom_labels=NULL){
  cex <- 0.9

  if(fam == "binomial"){

    plot(NA, xlim=c(-0.5,1.5), ylim=c(0,1)-c(-0.036,0.036), cex=cex, cex.lab=cex, bty="n", axes=F, mar=rep(0,4), mgp=c(0,0,0), xlab="", ylab="") #xlab="true class", ylab="predicted class"
    axis(side=1, at=0.5, labels=binom_labels["yhat1"], tick=F, cex=cex, cex.axis=cex, cex.lab=cex, line=-1)
    axis(side=2, at=0.5, labels=binom_labels["y0"], tick=F, cex=cex, cex.axis=cex, cex.lab=cex, line=-1)
    axis(side=3, at=0.5, labels=binom_labels["yhat0"], tick=F, cex=cex, cex.axis=cex, cex.lab=cex, line=-1)
    axis(side=4, at=0.5, labels=binom_labels["y1"], tick=F, cex=cex, cex.axis=cex, cex.lab=cex, line=-1)
    # topleft (tnr)
    rect(-0.5, d$fpr, 0.5, d$tnr+d$fpr, col="grey90")
    text(x=0, y=d$fpr+d$tnr/2, labels=paste0("tnr=",round(d$tnr,2)), cex=cex)
    # bottomright (tpr)
    rect(0.5, 0, 1.5, d$tpr, col="grey90")
    text(x=1, y=d$tpr/2, labels=paste0("tpr=",round(d$tpr,2)), cex=cex)
    # bottomleft (fpr)
    rect(-0.5, 0, 0.5, d$fpr, col="grey60")
    text(x=0, y=d$fpr/2, labels=paste0("fpr=",round(d$fpr,2)), cex=cex)
    # topright(fnr)
    rect(0.5, d$tpr, 1.5, d$tpr+d$fnr, col="grey60")
    text(x=1, y=d$tpr+d$fnr/2, labels=paste0("fnr=",round(d$fnr,2)), cex=cex)

  }else if(fam == "multinomial"){

    if(!is.null(multinom_labels)){
      colnames(d) <- multinom_labels
    }
    rescaled_d <- round(t(t(d)/colSums(d)),2) # show proportion relative to true class
    dim <- nrow(d)
    rng <- 100*range(rescaled_d[row(rescaled_d) != col(rescaled_d)])
    cols <- colorRampPalette(c("blue", "red"))(rng[2]-rng[1]+1)
    plot(0, col="white", xlim=c(0.5,dim+0.5), ylim=c(0.5,dim+1.5)+0, mar=rep(0,4), mgp=c(1,1,0), xaxs="i", xaxt="n", yaxt="n", xlab="true class", ylab="predicted class")
    axis(1, at=1:dim, labels=colnames(d), tick=F, cex=cex, cex.axis=cex, cex.lab=cex, line=-1)
    axis(2, at=1:dim, labels=rev(colnames(d)), tick=F, cex=cex, cex.axis=cex, cex.lab=cex, line=-1)
    for(d1 in 1:dim){# true class
      for(d2 in 1:dim){# predicted class
        if(d1 == d2){# diagonal
          rect(xleft=d1-0.5, ybottom=dim+0.5-d1, xright=d1+0.5, ytop=dim+1.5-d1, col="grey90", border=F)
          #text(x=d1, y=dim-d1+1, labels = d[d1,d1], cex=cex)
          text(x=d1, y=dim-d1+1, labels = sprintf("%.2f", rescaled_d[d1,d1]), cex=cex)
          text(x=d1, y=dim+0.85, labels = sprintf("%.2f", round(sum(d[-d1,d1])/sum(d[,d2]),2)), cex=cex) # class-specific error rate
        }else{# off-diagonal
          rect(xleft=d1-0.5, ybottom=dim+0.5-d2, xright=d1+0.5, ytop=dim+1.5-d2, col=alpha(cols[100*rescaled_d[d2,d1]-rng[1]+1],0.7), border=F)
          text(x=d1, y=dim-d2+1, labels=sprintf("%.2f", rescaled_d[d2,d1]), cex=cex)
        }}}
    avgError <- (sum(d) - sum(diag(d))) / sum(d)
    box(); text(x=0.35+dim/2, y=dim+1.25, labels=paste0("avg error=",round(avgError,2),"; class-specific error:"), cex=cex)
  }
}



#' @title Cross-validation levelplot
#'
#' @description Plots the cross-validation results for the number of components
#' \code{Q} and penalty parameter \code{lambda}, with the optimal and ensemble
#' models marked
#'
#' @param note_x vector of components (x-axis).
#'
#' @param note_y vector of penalty parameters (y-axis).
#'
#' @param letter A named list containing a string \code{char} to be plotted as
#' a subplot label at \code{(xloc,yloc)}.
#'
#' @details
#'
#' @export
#'
crossvalidation_levelplot <- function(note_x, note_y, letter=list(char="", xloc=1, yloc=1), ...){
  panel.levelplot(...)
  # SMOOTHEST MODEL
  panel.points(note_x[1], note_y[1], pch=20, col="black", cex=1)
  # CONTOUR FOR THE 'ACCEPTABLE MODELS'
  for(h in 1:length(note_x)){
    # bottom boundaries
    if( length(which(note_y[-h] < note_y[h])) > 0 ){
      ID_below <- which(note_y < note_y[h])
      ID_below_aligned <- ID_below[which(note_x[ID_below] == note_x[h])]
      if( length(ID_below_aligned) == 0 | ((length(ID_below_aligned) != 0) & !((note_y[h]-1) %in% note_y[ID_below_aligned])) ) # no neighbour below? draw lower line
        panel.rect(xleft = note_x[h]-0.5, xright = note_x[h]+0.5, ybottom = note_y[h]-0.5, ytop = note_y[h]-0.5, lwd=2)
    }else{
      panel.rect(xleft = note_x[h]-0.5, xright = note_x[h]+0.5, ybottom = note_y[h]-0.5, ytop = note_y[h]-0.5, lwd=2)
    }
    # top boundaries
    if( length(which(note_y[-h] > note_y[h])) > 0 ){
      ID_above <- which(note_y > note_y[h])
      ID_above_aligned <- ID_above[which(note_x[ID_above] == note_x[h])]
      if( length(ID_above_aligned) == 0 | ((length(ID_above_aligned) != 0) & !((note_y[h]+1) %in% note_y[ID_above_aligned])) ) # no neighbour above? draw lower line
        panel.rect(xleft = note_x[h]-0.5, xright = note_x[h]+0.5, ybottom = note_y[h]+0.5, ytop = note_y[h]+0.5, lwd=2)
    }else{
      panel.rect(xleft = note_x[h]-0.5, xright = note_x[h]+0.5, ybottom = note_y[h]+0.5, ytop = note_y[h]+0.5, lwd=2)
    }
    # left boundaries
    if( length(which(note_x[-h] < note_x[h])) > 0 ){
      ID_left <- which(note_x < note_x[h])
      ID_left_aligned <- ID_left[which(note_y[ID_left] == note_y[h])]
      if( length(ID_left_aligned) == 0 | ((length(ID_left_aligned) != 0) & !((note_x[h]-1) %in% note_x[ID_left_aligned])) ) # no neighbour left? draw lower line
        panel.rect(xleft = note_x[h]-0.5, xright = note_x[h]-0.5, ybottom = note_y[h]-0.5, ytop = note_y[h]+0.5, lwd=2)
    }else{
      panel.rect(xleft = note_x[h]-0.5, xright = note_x[h]-0.5, ybottom = note_y[h]-0.5, ytop = note_y[h]+0.5, lwd=2)
    }
    # right boundaries
    if( length(which(note_x[-h] > note_x[h])) > 0 ){
      ID_right <- which(note_x > note_x[h])
      ID_right_aligned <- ID_right[which(note_y[ID_right] == note_y[h])]
      if( length(ID_right_aligned) == 0 | ((length(ID_right_aligned) != 0) & !((note_x[h]+1) %in% note_x[ID_right_aligned])) ) # no neighbour right? draw lower line
        panel.rect(xleft = note_x[h]+0.5, xright = note_x[h]+0.5, ybottom = note_y[h]-0.5, ytop = note_y[h]+0.5, lwd=2)
    }else{
      panel.rect(xleft = note_x[h]+0.5, xright = note_x[h]+0.5, ybottom = note_y[h]-0.5, ytop = note_y[h]+0.5, lwd=2)
    }
  }#h
  panel.text(x=letter$xloc, y=letter$yloc, labels=bquote(bold(.(letter$char))))
}
