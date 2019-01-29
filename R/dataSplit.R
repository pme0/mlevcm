#' @title Data split
#'
#' @description Split dataset into training, vaidation and testing subsets.
#'
#' @param y A numeric or factor response vector of size \code{N}, where \code{N} is
#' the number of observations (spectra).
#'
#' @param train_size
#'
#' @param valid_size
#'
#' @param test_size
#'
#' @param reps
#'
#' @param balanced
#'
#' @param force_nonempty_response_class
#'
#' @return A list with the indices (running from \code{1} to \code{length(y)}) of
#' training, validation and testing observations. Each of these 3 lists is a list
#' of \code{reps} elements, each element being a vector of the length
#' \code{train_size}/\code{valid_size}/\code{test_size}, accordingly.
#'
#' @export
#'
dataSplit <- function(y, train_size, valid_size, test_size, reps, balanced, force_nonempty_response_class=TRUE){

  id_train <- id_valid <- id_test <- list()
  uniq <- sort(unique(y))
  lu <- length(uniq)
  ly <- length(y)

  if(balanced){

    nmin <- sum(y==uniq[[1]])
    ntr <- floor(train_size * nmin);   nva <- floor(valid_size * nmin);   nte <- nmin - ntr - nva  #floor(test_size  * nmin)
    positions <- matrix(0, nmin, lu)
    for(v in 1:lu){
      positions[,v] <- which(y==uniq[v])
    }
    for(r in 1:reps){
      positions <- positions[sample(nrow(positions)),]  # shuffle indices
      itr <- iva <- ite <- c()                          # define indicies for training/validation/testing sets
      for(v in 1:lu){
        itr <- c(itr, positions[1:ntr,v])
        iva <- c(iva, positions[(ntr+1):(ntr+nva),v])
        ite <- c(ite, positions[(ntr+nva+1):nmin,v])
      }
      id_train[[r]] <- itr
      id_valid[[r]] <- iva
      id_test[[r]]  <- ite
    }

  }else{

    ntr <- floor(ly * train_size)
    nva <- floor(ly * valid_size)
    nte <- floor(ly * test_size)
    while((ntr + nva + nte) != ly){
      ntr <- ntr + 1
    }
    for(r in 1:reps){
      if(force_nonempty_response_class){
        # make sure that each replication has samples from each possible response
        status <- T
        while(status){
          id_valid[[r]] <- sample(x = 1:ly, size = nva, replace = F)
          id_test[[r]]  <- sample(x = (1:ly)[!((1:ly) %in% id_valid[[r]])], size = nte, replace = F)
          id_train[[r]] <- sample(x = (1:ly)[!((1:ly) %in% c(id_valid[[r]],id_test[[r]]))], size = ntr, replace = F)
          status <- (length(unique(y[id_valid[[r]]])) != lu) || (length(unique(y[id_test[[r]]])) != lu)
        }
      }else{
        id_valid[[r]] <- sample(x = 1:ly, size = nva, replace = F)
        id_test[[r]]  <- sample(x = (1:ly)[!((1:ly) %in% id_valid[[r]])], size = nte, replace = F)
        id_train[[r]] <- sample(x = (1:ly)[!((1:ly) %in% c(id_valid[[r]],id_test[[r]]))], size = ntr, replace = F)
      }
    }
  }

  return(list(id_train = id_train,
              id_valid = id_valid,
              id_test  = id_test  ))
}


