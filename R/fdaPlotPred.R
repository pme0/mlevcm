#' @title Prediction results plot
#'
#' @description Displays visual diagnostics information and performance measures
#' for objects of class \code{fdaModelPred}.
#'
#' @param obj An object of class \code{fdaModelPred} from \code{\link{fdaML_predict}}.
#'
#' @export
#'
fdaPlotPred <- function(obj){


  if(class(obj) != 'fdaModelPred'){
    stop('object class invalid!')
  }

  uniq <- sort(unique(obj$new_y) + obj$cenY)
  lu <- length(uniq)

  if(obj$model == "lm"){

    #  plot(new_y+cenY, rowMeans(new_y_pred)+cenY); abline(0,1)
    plot(jitter(rep(obj$new_y, obj$reps) + obj$cenY), c(obj$new_y_pred) + obj$cenY, col=alpha("black",0.2), xlab='true age', ylab='predicted age'); abline(0,1)
    legend("bottomright", c(paste0("Kopt=", obj$Kopt), paste0("Kmax=", obj$Kmax)), lty=1:2, bty="n", cex=0.9)



    ### PLOT 3 --- best cross-validation parameter for all training/testing splits
    #
    # gather all data (across 'obj$reps') in one vector for each of the unique values of 'y'
    alldata <- list()
    for(u in 1:lu){
      alldata[[u]] <- obj$new_y_pred[(obj$new_y+obj$cenY) == uniq[u]] + obj$cenY
    }
    # compute percentiles for each of the unique values of 'y'
    ofs <- (min(uniq)==0);    ofs_x <- 0.4
    percentile_values <- c(0.05, 0.15, 0.25, 0.5, 0.75, 0.85, 0.95);
    percentiles <- matrix(0, length(percentile_values), lu);   extremes <- x_vec <- matrix(0, 2, lu);   means <- rep(0,lu)
    for(u in 1:lu){
      percentiles[,u] <- quantile(alldata[[u]], probs=percentile_values)
      x_vec[,u]       <- c(uniq[u] - ofs_x, uniq[u] + ofs_x)
      means[u] <- mean(alldata[[u]])
    }
    # boxplot
    xlabels <- rep(NA,max(uniq)-min(uniq)+2); xlabels[0:(1+max(uniq)+ofs) %in% uniq] <- uniq
    range_pred <- range(percentiles) #range(c(sapply(1:obj$reps, function(r){ obj$new_y_pred[,r]+obj$cenY })) )
    ofs <- (min(uniq)==0);    ofs_x <- 0.4
    plot(0, col=NA, xlim=c(0,1+max(uniq)), ylim=c(-1,1)+range_pred, xaxt='n', xlab="true age", ylab="predicted age")
    axis(side=1, at=0:(1+max(uniq)+ofs), labels=xlabels)
    for(u in 1:lu){
      # BOXES
      rect(x_vec[1,u], percentiles[3,u], x_vec[2,u], percentiles[5,u], col="white") # fill boxes
      segments(x_vec[1,u], percentiles[4,u], x_vec[2,u], percentiles[4,u], lty=1, lwd=3, lend=1)
      # WHISKERS
      segments(uniq[u],         percentiles[5,u], uniq[u],         percentiles[7,u], lty=1, lwd=1, lend=1)
      segments(uniq[u],         percentiles[3,u], uniq[u],         percentiles[1,u], lty=1, lwd=1, lend=1)
      # WHISKER ENDS
      segments(uniq[u]-ofs_x/2, percentiles[6,u], uniq[u]+ofs_x/2, percentiles[6,u], lty=1, lwd=1, lend=1)
      segments(uniq[u]-ofs_x/2, percentiles[2,u], uniq[u]+ofs_x/2, percentiles[2,u], lty=1, lwd=1, lend=1)
      segments(uniq[u]-ofs_x/4, percentiles[7,u], uniq[u]+ofs_x/4, percentiles[7,u], lty=1, lwd=1, lend=1)
      segments(uniq[u]-ofs_x/4, percentiles[1,u], uniq[u]+ofs_x/4, percentiles[1,u], lty=1, lwd=1, lend=1)
      # MEANS
      segments(x_vec[1,u]+ofs_x/2, means[u], x_vec[2,u]-ofs_x/2, means[u], lty=1, lwd=1, col="green")
    }
    # truth
    points(x=uniq + ofs, y=uniq, pch=19, col="Gold", cex=0.8)
    #abline(0,1)
    # legend
    legend("topleft", legend=c("mean", "truth"), col=c("green","Gold"), pch=c(NA,19), lty=c(1,NA), lwd=c(1,NA), seg.len=0.75, bty="n", cex=0.8)
    legend("topright", paste0("avg pred error: ",round(mean(obj$err_test),2)), bty="n", cex=0.8, text.col="brown")
    biasVec <- colMeans(obj$bias_test)
    text(x = uniq, y = range_pred[1]-1, labels = paste0(round(biasVec,2)), cex=0.8, col="purple")
    text(uniq[lu], range_pred[1]+0.3, "bias:", cex=0.8, col="purple")

    return(list(err = round(mean(obj$err_test),2),
                largest_bias = round(biasVec[which.max(abs(biasVec))]),2)
    )

  }else if(obj$model == "glm"){

    ### PLOT 3 --- error/accuracy and optimal cutoff point (putting equal weight on sensitivity and specificity)
    #
    if(obj$family == "binomial"){

      pred <- prediction(obj$ROC_to_plot_binom$pred_test, obj$ROC_to_plot_binom$labels)
      perf <- performance(pred, "tpr", "fpr")
      measureType <- "err"
      if(measureType == "err"){ #error
        # stats
        err_perf <- performance(pred, measure = "err")
        err_ind = sapply(1:obj$reps, function(z){ which.min(err_perf@y.values[[z]]) })
        if(any(err_ind == 1)){   # get rid of infinities in 'err_perf@x.values)'
          tt <- (1:obj$reps)[err_ind == 1]; for(l in 1:length(tt)){ if( is.infinite(slot(err_perf, "x.values")[[tt[l]]][1]) ){ err_ind[tt[l]] <- 2 }}
        }
        err_val = sapply(1:obj$reps, function(z){ err_perf@y.values[[z]][err_ind[z]] })
        err_cut = sapply(1:obj$reps, function(z){ err_perf@x.values[[z]][err_ind[z]] })
        avg_err <- mean(err_val)
        avg_err_cut <- mean(err_cut)
        avg_err_avgcut <- mean((obj$new_y_pred[,1] > avg_err_cut) != rep(obj$new_y + obj$cenY, 1))
        # plot
        # err_colcode <- ifelse( all(avg_err_avgcut < c(0.5, mean(obj$y==1), mean(obj$y==0))), "green", ifelse(all(avg_err_avgcut < 0.5), "orange", "red") )
        # d_allcurves <- data.frame(x = unlist(err_perf@x.values), y = unlist(err_perf@y.values)); d_allcurves <- d_allcurves[order(d_allcurves$x),]
        # plot(err_perf, lty=1, col="grey", xlim=c(0,1), ylim=c(0,1))   # cutoff-error curve for each randomisation (1:obj$reps)
        # lines(d_allcurves$x[!is.infinite(d_allcurves$x)], predict(loess(d_allcurves$y[!is.infinite(d_allcurves$x)] ~ d_allcurves$x[!is.infinite(d_allcurves$x)])), lwd=2)
        # points(err_cut,err_val)                                       # optimal (cutoff,error) point for each randomisation (1:obj$reps)
        # abline(v = avg_err_cut); abline(h = avg_err)                  # average of optimal cutoff-error
        # #points(avg_err_cut, avg_err, pch=15, col="black")
        # points(avg_err_cut, avg_err_avgcut, pch=19, col=err_colcode)
        # legend("topleft", c(paste0("avg error (w/ specific cutoff) = ", round(avg_err,2)), paste0("avg error (w/ avg cutoff) = ", round(avg_err_avgcut,2)), paste0("Prob[y=1] = ", round(mean(obj$y==1),2))), lty=c(1,NA,NA), pch=c(NA,19,NA), col=c("black",err_colcode,NA), seg.len=0.75,  bty="n", cex=0.8)
      }
    }

    ### PLOT 4 --- ROC for best cross-validation parameter for all replications, plus average in bold
    #
    avg_auc <- mean(obj$AUC_opt[,"test"])

    if(obj$family == "binomial"){

      pred <- prediction(obj$ROC_to_plot_binom$pred_test, obj$ROC_to_plot_binom$labels)
      perf <- performance(pred, "tpr", "fpr")
      dd <- data.frame(x = unlist(perf@x.values), y = unlist(perf@y.values)); dd <- dd[order(dd[,1]),]
      xvals <- seq(0, 1, by=0.1);   lx <- length(xvals)
      intervals <- matrix(c(xvals-0.05, xvals+0.05), lx, 2);   idx <- matrix(0, lx, 2)
      for(k in 1:nrow(intervals)){ idx[k,] <- range(which(dd$x >= intervals[k,1] & dd$x < intervals[k,2])) }
      forget <- which(is.infinite(idx[,1]) * is.infinite(idx[,2]) == 1)
      idx[forget,] <- 0
      ofs_x <- 0.025
      percentile_values <- c(0.05, 0.15, 0.25, 0.5, 0.75, 0.85, 0.95)
      percentiles <- matrix(0, length(percentile_values), lx);   x_vec <- matrix(0, 2, lx);   means <- stdev <- rep(0,lx)
      for(u in 1:lx){
        percentiles[,u] <- quantile(dd$y[idx[u,1]:idx[u,2]], probs=percentile_values)
        x_vec[,u]       <- c(xvals[u] - ofs_x, xvals[u] + ofs_x)
        means[u]        <- mean(dd$y[idx[u,1]:idx[u,2]])
        stdev[u]        <- sd(dd$y[idx[u,1]:idx[u,2]])
      }
      # plot
      plot(perf, col="grey", lty=1, ylim=c(0,1.1), xlab="false positive rate", ylab="true positive rate")
      for(u in 1:lx){
        # BOXES
        rect(x_vec[1,u], percentiles[3,u], x_vec[2,u], percentiles[5,u], col="white") # fill boxes
        segments(x_vec[1,u], percentiles[4,u], x_vec[2,u], percentiles[4,u], lty=1, lwd=3, lend=1)
        # WHISKERS
        segments(xvals[u], percentiles[5,u], xvals[u], percentiles[7,u], lty=1, lwd=1, lend=1)
        segments(xvals[u], percentiles[3,u], xvals[u], percentiles[1,u], lty=1, lwd=1, lend=1)
        # WHISKER ENDS
        segments(xvals[u]-ofs_x/2, percentiles[6,u], xvals[u]+ofs_x/2, percentiles[6,u], lty=1, lwd=1, lend=1)
        segments(xvals[u]-ofs_x/2, percentiles[2,u], xvals[u]+ofs_x/2, percentiles[2,u], lty=1, lwd=1, lend=1)
        segments(xvals[u]-ofs_x/4, percentiles[7,u], xvals[u]+ofs_x/4, percentiles[7,u], lty=1, lwd=1, lend=1)
        segments(xvals[u]-ofs_x/4, percentiles[1,u], xvals[u]+ofs_x/4, percentiles[1,u], lty=1, lwd=1, lend=1)
      }
      if(length(forget) == 0){
        lines(xvals, predict(loess(means ~ xvals)), lty=1, lwd=2)
      }else{
        lines(xvals[-forget], predict(loess(means ~ xvals)), lty=1, lwd=2)
      }
      lines(seq(0,1,by=0.1), seq(0,1,by=0.1), lty=2)
      # OR just confidence intervals instead of boxplots
      #lines(xvals, predict(loess(extremes[1,] ~ xvals)), lty=2, lwd=1); lines(xvals, predict(loess(extremes[2,] ~ xvals)), lty=2, lwd=1)
      text(x = xvals, y = rep(1.1,length(xvals)), labels = round(stdev,2), col="orange", cex=0.7)
      legend("bottomright", legend = paste0("average AUC = ", round(avg_auc,2)), bty="n", cex=0.8)

    }else{
      stop("this bit not coded yet!!!!")
    }

    return(list(avg_err = round(avg_err_avgcut, 2),
                avg_auc = round(avg_auc,2)
    ))

  }

}#function
