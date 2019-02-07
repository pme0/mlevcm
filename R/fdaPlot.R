#' @title Diagnostics and performance plots
#'
#' @description Displays visual diagnostics information and performance measures
#' for objects of class \code{fdaModel} and \code{fdaModelPred}.
#'
#' @param obj An object of class \code{fdaModel} from \code{\link{fdaML_train}}
#' or class \code{fdaModelPred} from \code{\link{fdaML_predict}}.
#'
#' @param hist_range A range for the histogram plot in binomial classification.
#'
#' @param roc_smoothing A value for the parameter \code{spar} in \code{\link{smooth.spline}}.
#'
#' @param multinom_labels Labels for the multinomial classes.
#'
#' @seealso
#' \code{\link{fdaML_train}} and \code{\link{fdaML_predict}}
#'
#' @export
#'
fdaPlot <- function(obj, hist_range = NULL, roc_smoothing = NULL, multinom_labels = NULL){

  if(class(obj) == "fdaModel"){

    oldparams <- par()
    uniq <- sort(unique(obj$y))
    lu <- length(uniq)

    if(obj$model == "lm"){

      ### PLOT 1 --- cross-validation
      #
      if(!(obj$Q_len==1)){
        matplot(obj$Q_vec, t(obj$perf_cv), type="l", lty=1, col="grey", xlab="components (Q)", ylab="RMSD")
        lines(obj$Q_vec, colMeans(obj$perf_cv), lwd=2)
        points(obj$Q_vec, colMeans(obj$perf_cv), pch=1, cex=0.8)
        abline(v = c(obj$Q_max, obj$Q_opt), lty=c(2,1), lwd=1)
        legend("topright", c(paste0("Q_opt=", obj$Q_opt), paste0("Q_max=", obj$Q_max)), lty=1:2, bty="n", cex=0.9)
      }


      ### PLOT 2 --- coefficient function for all replications, plus average in bold
      #
      matplot(obj$t_range_mod, obj$beta, type="l", lty=1, col="gray", xlab="wavelength", ylab="coefficient function")
      lines(obj$t_range_mod, rowMeans(obj$beta), type="l", lwd=2 , col="black")
      legend("bottomright", c("coefficient functions", "mean coefficient function"), lty=c(1,1), lwd=c(1,2), col=c('grey','black'), bty="n", cex=0.9)


      ### PLOT 3 --- best cross-validation parameter for all training/testing splits
      #
      # gather all data (across 'obj$reps') in one vector for each of the unique values of 'y'
      alldata <- list()
      for(u in 1:lu){
        alldata[[u]] <- 0
        for(r in 1:obj$reps){
          alldata[[u]] <- c(alldata[[u]], (obj$y_testpred_opt[,r]+obj$cenY[[r]])[obj$y[obj$id_test[[r]]] == uniq[u]])
        }
        alldata[[u]] <- alldata[[u]][-1]
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
      # Boxplot
      xlabels <- rep(NA,max(uniq)-min(uniq)+2); xlabels[0:(1+max(uniq)+ofs) %in% uniq] <- uniq
      range_pred <- range(percentiles)
      ofs <- (min(uniq)==0);    ofs_x <- 0.4
      plot(0, col=NA, xlim=c(0,1+max(uniq)), ylim=c(-1,1)+range_pred, xaxt='n', xlab="true age", ylab="predicted age")
      axis(side=1, at=0:(1+max(uniq)+ofs), labels=xlabels)
      for(u in 1:lu){
        # boxes
        rect(x_vec[1,u], percentiles[3,u], x_vec[2,u], percentiles[5,u], col="white") # fill boxes
        segments(x_vec[1,u], percentiles[4,u], x_vec[2,u], percentiles[4,u], lty=1, lwd=3, lend=1)
        # whiskers
        segments(uniq[u],         percentiles[5,u], uniq[u],         percentiles[7,u], lty=1, lwd=1, lend=1)
        segments(uniq[u],         percentiles[3,u], uniq[u],         percentiles[1,u], lty=1, lwd=1, lend=1)
        # whisker ends
        segments(uniq[u]-ofs_x/2, percentiles[6,u], uniq[u]+ofs_x/2, percentiles[6,u], lty=1, lwd=1, lend=1)
        segments(uniq[u]-ofs_x/2, percentiles[2,u], uniq[u]+ofs_x/2, percentiles[2,u], lty=1, lwd=1, lend=1)
        segments(uniq[u]-ofs_x/4, percentiles[7,u], uniq[u]+ofs_x/4, percentiles[7,u], lty=1, lwd=1, lend=1)
        segments(uniq[u]-ofs_x/4, percentiles[1,u], uniq[u]+ofs_x/4, percentiles[1,u], lty=1, lwd=1, lend=1)
        # means
        segments(x_vec[1,u]+ofs_x/2, means[u], x_vec[2,u]-ofs_x/2, means[u], lty=1, lwd=1, col="green")
      }
      points(x=uniq + ofs, y=uniq, pch=19, col="Gold", cex=0.8) # truth
      legend("topleft", legend=c("mean", "truth"), col=c("green","Gold"), pch=c(NA,19), lty=c(1,NA), lwd=c(1,NA), seg.len=0.75, bty="n", cex=0.8)
      legend("topright", paste0("avg test error: ",round(mean(obj$err_test),2)), bty="n", cex=0.8, text.col="brown")
      biasVec <- colMeans(obj$bias_test)
      text(x = uniq, y = range_pred[1]-1, labels = paste0(round(biasVec,2)), cex=0.8, col="purple")
      text(uniq[lu], range_pred[1]+0.3, "bias:", cex=0.8, col="purple")

      return( list(err = round(mean(obj$err_test),2),
                   largest_bias = biasVec[which.max(abs(biasVec))]) )

    }else if(obj$model == "glm"){

      ### PLOT 1 --- cross-validation
      #
      if(!(obj$Q_len==1)){
        matplot(obj$Q_vec, t(obj$perf_cv), type="l", lty=1, col="grey", xlab="components (Q)", ylab="AUC", ylim=c(0,1))
        lines(obj$Q_vec, colMeans(obj$perf_cv), lwd=2)
        abline(v = c(obj$Q_opt, obj$Q_max), lty=c(1,2), lwd=1)
        legend("bottomright", c(paste0("Q opt = ", obj$Q_opt, " (AUC = ", round(colMeans(obj$perf_cv)[obj$Q_vec==obj$Q_opt],2),")"), paste0("Q max = ", obj$Q_max, " (AUC = ", round(colMeans(obj$perf_cv)[obj$Q_vec==obj$Q_max],2),")")), lty=1:2, bty="n", cex=0.9)
      }


      ### PLOT 2 --- coefficient function for all replications, plus average in bold
      #
      if(obj$family != "multinomial"){
        xlim_vec <- range(obj$t_range_mod)+c(-100*(obj$S+1),0)
        matplot(obj$t_range_mod, obj$beta, type="l", lty=1, col=alpha("gray",0.6), xlab="wavelength", ylab="coefficient function", xlim=xlim_vec, ylim=range(obj$beta))
        lines(obj$t_range_mod, rowMeans(obj$beta), type="l", lwd=2 , col="black")
        par(new = TRUE)
        plot(jitter(rep(xlim_vec[1], obj$reps), factor=5), obj$alpha, axes=F, col=alpha("blue",0.4), pch=5, cex=0.8, xlab="", ylab="", xlim=xlim_vec, ylim=range(rbind(obj$alpha,obj$gamma)))
        points(xlim_vec[1], mean(obj$alpha), pch="-", col="blue", cex=3, xlab="", ylab="")
        if(!is.null(obj$Z)){
          for(q in 1:obj$S){ points(jitter(rep(xlim_vec[1]+100*q, obj$reps), factor=5), obj$gamma[q,], col=alpha("blue",0.4), pch=1, cex=1.1)
            points(xlim_vec[1]+100*q, mean(obj$gamma[q,]), col="blue", pch=20, cex=1) }}
        axis(side=4, at = pretty(range(rbind(obj$alpha,obj$gamma))), ylab="?", col="blue", col.axis="blue",las=1)
      }
      legend("bottomright", c("coefficient functions", "mean coefficient function"), lty=c(1,1), lwd=c(1,2), col=c('grey','black'), bty="n", cex=0.9)


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
          error_per_rep <- sapply(1:obj$reps, function(qq){ mean((obj$y_testpred_opt[,qq] > avg_err_cut) != obj$y[obj$id_test[[qq]]]) })
          avg_err_avgcut <- mean(error_per_rep)
          stdev_err_avgcut <- sd(error_per_rep)
          err_colcode <- ifelse( all(avg_err_avgcut < c(0.5, mean(obj$y==1), mean(obj$y==0))), "green", ifelse(all(avg_err_avgcut < 0.5), "orange", "red") )
          d_allcurves <- data.frame(x = unlist(err_perf@x.values), y = unlist(err_perf@y.values)); d_allcurves <- d_allcurves[order(d_allcurves$x),]
          # plot
          plot(err_perf, lty=1, col="grey", xlim=c(0,1), ylim=c(0,1), xlab="cutoff (probability)", ylab="error rate")   # cutoff-error curve for each randomisation (1:obj$reps)
          lines(d_allcurves$x[!is.infinite(d_allcurves$x)], predict(loess(d_allcurves$y[!is.infinite(d_allcurves$x)] ~ d_allcurves$x[!is.infinite(d_allcurves$x)])), lwd=2)
          points(err_cut,err_val)                                                 # optimal (cutoff,error) point for each randomisation (1:obj$reps)
          abline(v = avg_err_cut, lty=2)                                          # average cutoff
          points(avg_err_cut, avg_err, pch=19, col="black", cex=1.3)              # avg error (w/ curve-specific cutoff)
          points(avg_err_cut, avg_err_avgcut, pch=19, col=err_colcode, cex=1.3)   # avg error (w/ avg cutoff)
          legend("topright", c("cutoff-error curves", "mean cutoff-error curve", "avg cutoff", paste0("avg error (w/ curve-specific cutoff) = ", round(avg_err,2)), paste0("avg error (w/ avg cutoff) = ", round(avg_err_avgcut,2)), paste0("freq minority class = ", round(length(obj$y[obj$y == names(table(obj$y))[which.min(table(obj$y))]]) / length(obj$y),2))), lty=c(1,1,2,NA,NA,NA), pch=c(NA,NA,NA,19,19,NA), col=c("grey","black","black","black",err_colcode,NA), seg.len=1.2,  bty="n", cex=0.9)
        }else if(measureType == "acc"){ #accuracy
          # stats
          acc_perf <- performance(pred, measure = "acc")
          acc_ind = sapply(1:obj$reps, function(z){ which.max(acc_perf@y.values[[z]]) })
          if(any(acc_ind == 1)){   # get rid of infinities in 'acc_perf@x.values'
            tt <- (1:obj$reps)[acc_ind == 1]; for(l in 1:length(tt)){ if( is.infinite(slot(acc_perf, "x.values")[[tt[l]]][1]) ){ acc_ind[tt[l]] <- 2 }}
          }
          acc_val = sapply(1:obj$reps, function(z){ acc_perf@y.values[[z]][acc_ind[z]] })
          acc_cut = sapply(1:obj$reps, function(z){ acc_perf@x.values[[z]][acc_ind[z]] })
          avg_acc <- mean(acc_val)
          avg_acc_cut <- mean(acc_cut)
          avg_acc_avgcut <- sapply(1:obj$reps, function(z){ acc_perf@y.values[[z]][which.min(abs(acc_perf@x.values[[z]] - avg_acc_cut))] }) # this is only an approximation since it seems difficult to find performance for a specific cutoff
          acc_colcode <- ifelse( all(mean(avg_acc_avgcut) > c(0.5, mean(obj$y==1), mean(obj$y==0))), "green", ifelse(all(mean(avg_acc_avgcut) > 0.5), "orange", "red") )
          # plot
          plot(acc_perf, lty=1, col="grey", xlim=c(0,1), ylim=c(0,1))   # cutoff-accuracy curve for each randomisation (1:obj$reps)
          points(acc_cut, acc_val)                                      # optimal (cutoff,accuracy) point for each randomisation (1:obj$reps)
          abline(v = avg_acc_cut); abline(h = avg_acc)                  # average of optimal cutoff-accuracy
          points(avg_acc_cut, mean(avg_acc_avgcut), pch=19, col=acc_colcode)
          legend("bottomleft", c(paste0("avg accuracy (w/ specific cutoff) = ", round(avg_acc,2)), paste0("avg accuracy (w/ avg cutoff) = ", round(mean(avg_acc_avgcut),2)), paste0("Prob[y=1] = ", round(mean(obj$y==1),2))), bty="n", cex=0.8)
        }else{
          stop("'measureType' must be either 'err' (for error rate) or 'acc' (for accuracy rate).")
        }

      }else if(obj$family == "multinomial"){
        ;
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
        plot(perf, col="grey", lty=1, ylim=c(0,1.05), xlab="false positive rate", ylab="true positive rate")
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
        text(x = xvals, y = rep(1.05,length(xvals)), labels = sapply(1:length(xvals), function(aa){ paste0(strsplit(as.character(round(stdev,2)), split="")[[aa]][2:4], collapse='') }), col="black", cex=0.7)
        # for text inside boxes use:
        #text(x = xvals, y = percentiles[5,]-0.025, labels = sapply(1:length(xvals), function(aa){ paste0(strsplit(as.character(round(stdev,2)), split="")[[aa]][2:4], collapse='') }), col="black", cex=0.7)
        labls <- c(paste0("avg AUC = ", round(avg_auc,2)), "fpr-tpr curves", "mean fpr-tpr curve", "standard deviation")
        legend("bottomright", legend=labls, lty=c(NA,1,1,NA), lwd=c(NA,1,2,NA), col=c(NA,'grey','black','black'), bty="n", cex=0.9)
        text(x=0.55, y=0.02, labels="x.xx", cex=0.9)

      }else if(obj$family == "multinomial"){

        lvls <- levels(factor(obj$y))
        d_lvl <- list(); auc_lvl <- matrix(0, obj$reps, lu);
        #pdf("~/Downloads/f1.pdf", height=6, width=6, pointsize=18)
        plot(0, col="white", xlim=c(0,1), ylim=c(0,1.04), xlab="False positive rate", ylab="True positive rate")
        for(lvls_id in 1:lu){
          pred_lvl <- prediction(obj$ROC_to_plot_multinom[[lvls[lvls_id]]]$pred_test, obj$ROC_to_plot_multinom[[lvls[lvls_id]]]$labels)
          perf_lvl <- performance(pred_lvl, "tpr", "fpr")
          auc_lvl[, lvls_id] <- unlist(slot(performance(pred_lvl, "auc"), "y.values"))
          roc_x <- roc_y <- c()
          for(kk in 1:obj$reps){
            #lines(unlist(perf_lvl@x.values[[kk]]), unlist(perf_lvl@y.values[[kk]]), lty=2, col=1+lvls_id)
            roc_x <- c(roc_x, unlist(perf_lvl@x.values[[kk]]))
            roc_y <- c(roc_y, unlist(perf_lvl@y.values[[kk]]))
          }
          d_lvl[[lvls_id]] <- cbind(roc_x, roc_y); d_lvl[[lvls_id]] <- d_lvl[[lvls_id]][order(d_lvl[[lvls_id]][,1]),]  # order points along x axis
          if(is.null(roc_smoothing)){ roc_smoothing <- 0.25 }
          smoothingSpline <- smooth.spline(d_lvl[[lvls_id]][,1], d_lvl[[lvls_id]][,2], spar=roc_smoothing)
          sSpline_x <- smoothingSpline$x
          sSpline_y <- smoothingSpline$y; sSpline_y[sSpline_y > 1] <- 1
          lines(sSpline_x, sSpline_y, col=lvls_id+1, lty=1, lwd=2)
          #points(d_lvl[[lvls_id]][,1], d_lvl[[lvls_id]][,2], pch='.', col=lvls_id+1)
          #lines(d_lvl[[lvls_id]][,1], predict(loess(d_lvl[[lvls_id]][,2] ~ d_lvl[[lvls_id]][,1], span=0.2)), col=lvls_id+1, lty=1, lwd=2)
        }
        avg_auc <- round(mean(auc_lvl),2)
        the_names <- lvls  #; the_names <- multinom_labels
        legend("bottomright", legend = paste0("y=", the_names, ": ", sprintf("%.2f",round(colMeans(auc_lvl),2))), title=expression(bold("one vs all AUC:")), lwd=2, lty=rep(1,lu), col=2+uniq, bty="n", cex=0.9)
        legend("bottomleft", legend = paste0("avg AUC: ", sprintf("%.2f", avg_auc)), bty="n", cex=0.9)
        #dev.off()
      }


      ### PLOT 5 --- positive/negative densities as a function of the linear predictor
      #
      if(obj$family == "binomial"){

        # ALTERNATIVE WAY TO FIND CUTOFF AND GIVEN DIFFERENT WAYS TO FALSE POSITIVES AND FALSE NEGATIVES
        # SEE https://www.r-bloggers.com/a-small-introduction-to-the-rocr-package
        # cost.perf = performance(pred, "cost", cost.fp = 1, cost.fn = 1)
        # pred@cutoffs[[1]][which.min(cost.perf@y.values[[1]])]

        # organise data
        id1 <- c(sapply(1:obj$reps, function(z){ obj$y[obj$id_test[[z]]]==1 }))
        d <- data.frame(label = as.factor(id1+0), val_linpred = c(obj$linear_predictor)); d$val_prob <- 1/(1+exp(-d$val_linpred))
        # plot
        nbins <- 40;   col0 <- "skyblue";   col1 <- rgb(1,0,0,0.6)
        tt <- "linpred"   # {"prob", "linpred"}
        if(tt == "linpred"){
          optcut <- -log((1-avg_err_cut)/avg_err_cut) # obtained by solving 'P = 1 / (1 + exp(-L))' for the linear predictor 'L'
          dd <- d$val_linpred
          if(is.null(hist_range)){
            n1 <- min(d$val_linpred);  n2 <- max(d$val_linpred); wd <- n2-n1;  perc <- abs(range(d$val_linpred)-optcut)/wd
          }else{
            n1 <- hist_range[1];       n2 <- hist_range[2];      wd <- n2-n1;  perc <- abs(range(hist_range)-optcut)/wd
          }
          brks <- c(seq(n1, optcut, len=perc[1]*nbins), seq(optcut, n2, len=perc[2]*nbins+1)[-1])
        }else if(tt == "prob"){
          optcut <- avg_err_cut
          dd <- d$val_prob
          brks <- seq(0, 1, len=nbins)
        }
        # confusion matrix
        # individual CM
        #confu <- confusionMatrix(data=as.factor(as.numeric(obj$y.testpred.opt[,ii] > err_cut[ii])), reference = obj$y[obj$id.test[[ii]]])$table
        # grouped CM: assumes the average cutoff for all replications
        confu <- confusionMatrix(data=factor(as.numeric(c(obj$y_testpred_opt) > avg_err_cut)), reference = factor(obj$y[unlist(obj$id_test)]))$table
        dc <- list(tnr = confu[1,1] / sum(confu[,1]),
                   tpr = confu[2,2] / sum(confu[,2]),
                   fnr = confu[1,2] / sum(confu[,2]),
                   fpr = confu[2,1] / sum(confu[,1]))
        # grouped CM: assumes specific cutoffs for each replication
        #auxM <- matrix(0, nrow(obj$y.testpred.opt), ncol(obj$y.testpred.opt))
        #for(k in 1:ncol(auxM)){ auxM[,k] <- obj$y.testpred.opt[,k] > err_cut[k] }
        #confu <- confusionMatrix(data=as.factor(as.numeric(c(auxM))), reference = obj$y[unlist(obj$id.test)])$table
        #
        #u <- par("usr"); v <- c(grconvertX(u[1:2], "user", "ndc"), grconvertY(u[3:4], "user", "ndc"))
        #v <- c( (v[1]+v[2])/2, v[2], (v[3]+v[4])/2, v[4] )
        # plot
        # fourfoldplot( round(confu / rep(colSums(confu),each=2), 2), col=c("grey50","grey80"), conf.level=0, std="all.max")
        #pdf("./nirs/typesetting/nirs_paper1_labdata/figures/qqq.pdf", height=6, width=6, pointsize=12)
        if(is.null(hist_range)){
          hist(dd[d$label==0], col=col0, freq=F, breaks=brks, border="white", xlab=ifelse(tt=="prob","probability of class 1","linear predictor"), ylab="density", main="")
          hist(dd[d$label==1], col=col1, freq=F, breaks=brks, border="white", add=T)
        }else{
          hist(dd[d$label==0 & dd > hist_range[1] & dd < hist_range[2]], col=col0, freq=F, breaks=brks, border="white", xlab=ifelse(tt=="prob","probability of class 1","linear predictor"), ylab="density", main="")
          hist(dd[d$label==1 & dd > hist_range[1] & dd < hist_range[2]], col=col1, freq=F, breaks=brks, border="white", add=T)
        }
        abline(v=optcut, lwd=2)
        #
        #l1 <- legend("topleft", "optimal cutoff", lty=1, lwd=2, col="black", seg.len=0.8, bty="n", cex=0.9)
        #l2 <- legend(x=l1$rect$left, y=l1$rect$top-l1$rect$h, c("uninfected (y=0)", "infected (y=1)"), title=expression(bold("true class:")), pch=15, col=c(col0,col1), seg.len=0.8, bty="n", cex=0.9)
        #legend("bottomleft", c(paste0("avg error = ",round(avg_err_avgcut,2)), ""), bty="n", cex=0.9)
        #
        #l1 <- legend("bottomleft", c("uninfected (y=0)", "infected (y=1)", ""), title=expression(bold("true class:")), pch=c(15,15,NA), col=c(col0,col1,NA), seg.len=0.8, bty="n", cex=0.9)
        #l2 <- legend(x=l1$rect$left, y=with(l1$rect, top+0.1), c(paste0("avg error = ",round(avg_err_avgcut,2)), "optimal cutoff"), col=c(NA,"black"), lty=c(NA,1), lwd=c(NA,2), seg.len=0.8, bty="n", cex=0.9)
        #
        legend("bottomleft", c(paste0("average error = ",round(avg_err_avgcut,2)), "optimal cutoff", expression(bold(true*" "*class*":")), "uninfected (y=0)", "infected (y=1)", ""), lty=c(NA,1,NA,1,1,NA), lwd=c(NA,2,NA,NA,NA,NA), pch=c(NA,NA,NA,15,15,NA), col=c(NA,"black",NA,col0,col1,NA), seg.len=0.8, bty="n", cex=0.9)
        par(fig = c(0.15, 0.55, 0.5, 0.9), new=TRUE, mar=2+c(0,0,0,0))
        confusion_plot(d=dc, fam=obj$family)
        options(showWarnCalls = F)
        par(oldparams)
        options(showWarnCalls = T)
        #dev.off()

      }else if(obj$family == "multinomial"){

        lvls <- levels(as.factor(obj$y)); lu <- length(lvls)
        # compute indices of test observation by true class
        lvls_id_true <- list()
        for(q in 1:lu){
          lvls_id_true[[paste0("lvl",lvls[q])]] <- unlist(sapply(1:obj$reps, function(z){ obj$id_test[[z]][obj$y[obj$id_test[[z]]]==uniq[q]] }))
        }
        # compute class probabilities by response class
        classProbs_matrix <- list()
        for(q in 1:lu){
          classProbs_matrix[[paste0("lvl",lvls[q])]] <- 0
          for(r in 1:obj$reps){
            idx <- obj$y[obj$id_test[[r]]] == uniq[q]
            classProbs_matrix[[paste0("lvl",lvls[q])]] <- rbind(classProbs_matrix[[paste0("lvl",lvls[q])]], obj$classProbs[[r]][idx,])
          }
          classProbs_matrix[[paste0("lvl",lvls[q])]] <- classProbs_matrix[[paste0("lvl",lvls[q])]][-1,]
        }
        # compute percentage of correct class as given by highest probability
        classPercents <- rep(NA,lu)
        for(q in 1:lu){
          classPercents[q] <- mean((apply(classProbs_matrix[[paste0("lvl",lvls[q])]], 1, which.max) - 1) == uniq[q])
        }
        # plot predicted class probabilities
        #flat <- unlist(lapply(classProbs_matrix, colMeans))
        #plot(uniq,rep(0,lu), col="white", xlim=range(uniq), ylim=c(0.8, 1.1)*range(flat))
        #for(t in uniq){ lines(uniq, colMeans(classProbs_matrix[[paste0("lvl",t)]]), type="l", lwd=2, col=t+1) }
        #legend("topright",paste0("y=",uniq), seg.len=0.75, bty="n", cex=0.9, col=uniq+1, lty=1, lwd=2)
        ddd <- as.table(matrix(unlist(lapply(classProbs_matrix, colMeans)), lu, lu, byrow=F)); rownames(ddd) <- paste0("pred_",lvls); colnames(ddd) <- paste0("true_",lvls)
        #pdf("~/Downloads/f2.pdf", height=6, width=6, pointsize=18)
        the_names <- lvls  ; the_names <- multinom_labels
        barplot(ddd, ylim=c(0,1.2), beside=T, col=uniq+2, names.arg=the_names, xlab="true class", ylab="average predicted class probability")
        legend("top", paste0("y=",the_names), title=expression(bold("predicted class:")), horiz=T, text.width=2.5, bty="n", cex=0.9, pch=rep(15,lu), col=uniq+2)
        #dev.off()
        # plot confusion matrix
        confu_multinom <- confusionMatrix(data=as.factor(as.numeric(c(obj$y_testpred_opt))), reference = obj$y[unlist(obj$id_test)])$table
        avg_err_avgcut <- (sum(confu_multinom) - sum(diag(confu_multinom))) / sum(confu_multinom)
        #pdf("~/Downloads/f3.pdf", height=6, width=6, pointsize=18)
        confusion_plot(d=confu_multinom, fam=obj$family, avgError=avg_err_avgcut, multinom_labels=multinom_labels)
        #dev.off()

      }

      return(list(err = avg_err_avgcut,
                  err_stdev = stdev_err_avgcut,
                  auc = avg_auc,
                  roughnessStd = mean(obj$beta_roughnessStd)) )
    }#glm

  }else if(class(obj) == "fdaModelPred"){

    if(obj$model == "lm"){

      print(paste0("avgTestErr_avgTrainCut = ", obj$avgTestErr))
      biasMat <- colMeans(obj$bias_test); names(biasMat) <- paste('age', sort(unique(obj$new_y)))
      print("Bias:")
      print(biasMat)
      plot(sort(unique(obj$new_y)), biasMat, pch=19, xlab="age", ylab="bias"); abline(0,0)
      legend('topright', paste0("avgTestErr = ", round(mean(obj$err_test), 2)), bty="n", cex=0.9)

    }else if(obj$model == "glm"){

      if(obj$family == "binomial"){
        print(paste0("avgTestErr_avgTrainCut = ", round(obj$avgTestErr, 2)))
        confusion_plot(d = obj$errorBreakdown, fam = obj$family)
      }else if(obj$family == 'multinomial'){
        print(paste0("avgTestErr = ", round(obj$avgTestErr, 2)))
        confusion_plot(d = obj$errorBreakdown, fam = obj$family) }

    }

  }else{
    stop("Class of 'obj' in fdaPlot() is invalid.")
  }

}
