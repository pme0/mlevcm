library("devtools")
library("roxygen2")

# write package documentation
setwd("~/Dropbox/research/projects/nirs/mlevcm/")
document()
?fdaML_train
?mlevcm

# install package
setwd("~/Dropbox/research/projects/nirs/")
install("mlevcm")
# or directly from GitHub:
.rs.restartR()
devtools::install_github("pmesperanca/mlevcm", force=T)
library("mlevcm")
?mlevcm

rm(list = ls())
detach("package:mlevcm", unload=TRUE)



data_zenodo <- data.frame(Scan_ID = raw_data$Scan_ID,
                    Mosquito_ID = raw_data$Mosquito_ID,
                    Position = raw_data$Position,
                    Location = raw_data$Location,
                    Generation = raw_data$Generation,
                    Species = raw_data$Species)
data_zenodo <- cbind(data_zenodo, X)


write.table(data_zenodo, file="~/Dropbox/research/projects/nirs/data/zenodo/A1T1.txt", sep=",", row.names=F)






# load data 'B_A1E1' from cleanRealData.R

# train
ii <- 1:length(y) #which(raw_data[,"Species"] == "coluzzii")
X_train <- X[ii,]
y_train <- y[ii]
set.seed(17571); obj1 <- fdaML_train(ll = list(X=X_train, y=y_train, Z=NULL, task='regr', model='lm', reduction='pls', cv='n', intercept=T, reps=30, smooth_w=rep(1,ncol(X_train)), Qlen=20, Qopt=NULL, Qvec=NULL, split_size=0.5, lam_cv_type="n", tau_Q_opt=0.05, balanced=F, estimation_w=NULL, bspline_dim=floor(0.1 * ncol(X)), t_range=wvlenghts, verbose=T) )
fdaPlot(obj1, hist_range=c(-25,15))

# spectra plot

matplot(wvlenghts, t(X_train), type="l", ylim=c(0,2.5), col=y_train, lty=1, xlab="wavelength", ylab="signal")

# predict
ii <- which(raw_data[,"Species"] == "gambiae")
X_pred <- X[ii,]
y_pred <- y[ii]
obj2 <- fdaML_predict(obj = obj1, new_x = X_pred, new_z = NULL, new_y = y_pred)
fdaPlotPred(obj2)


lb <- 1;  lc <- 2;  ll <- 3
models_list <- list()
for(nb in 1:lb){
  models_list[[nb]] <- list()
  for(nc in 1:lc){
    models_list[[nb]][[nc]] <- list()
    for(nl in 1:ll){
      models_list[[nb]][[nc]][[nl]] <- 4
}}}





