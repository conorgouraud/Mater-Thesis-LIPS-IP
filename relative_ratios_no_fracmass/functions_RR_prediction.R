#functions
predRR21.agg.nfc <- function(x, modelOutput, conf.int = TRUE, conf.level=0.95, S=FALSE){
  
  pred <- predict(modelOutput, x, aggregate = TRUE)[[1]][,1]

  if(conf.int==TRUE){
    
    crit.val <- qnorm((1 - conf.level)/2, lower.tail=FALSE)

    if(S==FALSE){
      nObs <- 62938
      MSE <- 0.001540977
      avgMass <- 0.8993654
      predMassDiff <- (x$m0 - avgMass)^2
      sumDiff <-  20986.82
      stdErr <- sqrt(MSE*(1 + 1/nObs + predMassDiff/sumDiff))

      
    }else{
      nObs <- 2205
      MSE <- 0.001172896
      avgMass <- 1.069031
      predMassDiff <- (x$m0 - avgMass)^2
      sumDiff <- 447.149
      stdErr <- sqrt(MSE*(1 + 1/nObs + predMassDiff/sumDiff))
    }
  
    predLower <-  pred-crit.val*stdErr
    predUpper <-  pred+crit.val*stdErr
    
    output <- cbind(x, pred, predLower, predUpper)
    
  }else{
    output <- cbind(x, pred)
  }
  
  return(output)
  
}

predRR32.agg.nfc <- function(x, modelOutput, conf.int = TRUE, conf.level=0.95, S=FALSE){
  
  pred <- predict(modelOutput, x, aggregate = TRUE)[[1]][,1]

  if(conf.int==TRUE){
    
    crit.val <- qnorm((1 - conf.level)/2, lower.tail=FALSE)

    if(S==FALSE){
      nObs <- 62677
      MSE <- 0.00058191
      avgMass <- 0.9015365
      predMassDiff <- (x$m0 - avgMass)^2
      sumDiff <- 20910.98
      stdErr <- sqrt(MSE*(1 + 1/nObs + predMassDiff/sumDiff))
      
      
    }else{
      nObs <- 2190
      MSE <- 0.004651509
      avgMass <- 1.072511
      predMassDiff <- (x$m0 - avgMass)^2
      sumDiff <- 442.8918
      stdErr <- sqrt(MSE*(1 + 1/nObs + predMassDiff/sumDiff))
    }
    
    predLower <-  pred-crit.val*stdErr
    predUpper <-  pred+crit.val*stdErr
    
    output <- cbind(x, pred, predLower, predUpper)
    
  }else{
    output <- cbind(x, pred)
  }
  
  return(output)
  
}

predRR43.agg.nfc <- function(x, modelOutput, conf.int = TRUE, conf.level=0.95, S=FALSE){
  
  pred <- predict(modelOutput, x, aggregate = TRUE)[[1]][,1]

  if(conf.int==TRUE){
    
    crit.val <- qnorm((1 - conf.level)/2, lower.tail=FALSE)

    if(S==FALSE){
      nObs <- 62938
      MSE <- 0.0001490671
      avgMass <- 0.8993654
      predMassDiff <- (x$m0 - avgMass)^2
      sumDiff <- 20986.82
      stdErr <- sqrt(MSE*(1 + 1/nObs + predMassDiff/sumDiff))
      
      
    }else{
      nObs <- 2205
      MSE <- 0.0003596577
      avgMass <- 1.069031
      predMassDiff <- (x$m0 - avgMass)^2
      sumDiff <- 447.149
      stdErr <- sqrt(MSE*(1 + 1/nObs + predMassDiff/sumDiff))
    }
    
    predLower <-  pred-crit.val*stdErr
    predUpper <-  pred+crit.val*stdErr
    
    output <- cbind(x, pred, predLower, predUpper)
    
  }else{
    output <- cbind(x, pred)
  }
  
  return(output)
  
}

predRR54.agg.nfc <- function(x, modelOutput, conf.int = TRUE, conf.level=0.95, S=FALSE){
  
  pred <- predict(modelOutput, x, aggregate = TRUE)[[1]][,1]

  if(conf.int==TRUE){
    
    crit.val <- qnorm((1 - conf.level)/2, lower.tail=FALSE)

    if(S==FALSE){
      nObs <- 62677
      MSE <- 5.74326e-05
      avgMass <- 0.9015365
      predMassDiff <- (x$m0 - avgMass)^2
      sumDiff <- 20910.98
      stdErr <- sqrt(MSE*(1 + 1/nObs + predMassDiff/sumDiff))
      
    }else{
      nObs <- 2190
      MSE <- 0.0004727521
      avgMass <- 1.072511
      predMassDiff <- (x$m0 - avgMass)^2
      sumDiff <- 442.8918
      stdErr <- sqrt(MSE*(1 + 1/nObs + predMassDiff/sumDiff))
    }
    
    predLower <-  pred-crit.val*stdErr
    predUpper <-  pred+crit.val*stdErr
    
    output <- cbind(x, pred, predLower, predUpper)
    
  }else{
    output <- cbind(x, pred)
  }
  
  return(output)
  
}

