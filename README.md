# mlevcm

**mlevcm** is an R package for **ML**-based **E**pidemiological **V**ector **C**ontrol **M**onitoring using FDA techniques for NIRS data.

The package includes functions to train, predict and assess the performance of machine learning models for spectral data using techniques from Functional Data Analysis---namely functional representation of spectra, spectra smoothing and penalised parameter estimation.


## Main Functions

The main functions in this package are:
* **fdaML_train()** trains a machine learning model and assesses its performance;
* **fdaML_predict()** produces predictions for a new set of observations, given a previously trained model;
* **fdaPlot()** displays visual diagnostics information and performance measures;
