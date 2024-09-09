#functions
predRR21.agg.fc <- function(x, modelOutput, conf.int = TRUE, conf.level=0.95, S=FALSE){
  
  pred <- predict(modelOutput, x, aggregate = TRUE)[[1]][,1]

  if(conf.int==TRUE){
    
    crit.val <- qnorm((1 - conf.level)/2, lower.tail=FALSE)

    if(S==FALSE){
      nObs <- 62938
      MSE <- 0.0008415981
      avgMass <- 0.8993654
      predMassDiff <- (x$m0 - avgMass)^2
      sumDiff <-  20986.82
      stdErr <- sqrt(MSE*(1 + 1/nObs + predMassDiff/sumDiff))

      
    }else{
      nObs <- 2205
      MSE <- 0.0005122064
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

predRR32.agg.fc <- function(x, modelOutput, conf.int = TRUE, conf.level=0.95, S=FALSE){
  
  pred <- predict(modelOutput, x, aggregate = TRUE)[[1]][,1]

  if(conf.int==TRUE){
    
    crit.val <- qnorm((1 - conf.level)/2, lower.tail=FALSE)

    if(S==FALSE){
      nObs <- 62677
      MSE <- 0.0001134447
      avgMass <- 0.9015365
      predMassDiff <- (x$m0 - avgMass)^2
      sumDiff <- 20910.98
      stdErr <- sqrt(MSE*(1 + 1/nObs + predMassDiff/sumDiff))
      
      
    }else{
      nObs <- 2190
      MSE <- 0.003569434
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

predRR43.agg.fc <- function(x, modelOutput, conf.int = TRUE, conf.level=0.95, S=FALSE){
  
  pred <- predict(modelOutput, x, aggregate = TRUE)[[1]][,1]

  if(conf.int==TRUE){
    
    crit.val <- qnorm((1 - conf.level)/2, lower.tail=FALSE)

    if(S==FALSE){
      nObs <- 62938
      MSE <- 7.642264e-05
      avgMass <- 0.8993654
      predMassDiff <- (x$m0 - avgMass)^2
      sumDiff <- 20986.82
      stdErr <- sqrt(MSE*(1 + 1/nObs + predMassDiff/sumDiff))
      
      
    }else{
      nObs <- 2205
      MSE <- 0.0002348063
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

predRR54.agg.fc <- function(x, modelOutput, conf.int = TRUE, conf.level=0.95, S=FALSE){
  
  pred <- predict(modelOutput, x, aggregate = TRUE)[[1]][,1]

  if(conf.int==TRUE){
    
    crit.val <- qnorm((1 - conf.level)/2, lower.tail=FALSE)

    if(S==FALSE){
      nObs <- 62677
      MSE <- 3.327013e-05
      avgMass <- 0.9015365
      predMassDiff <- (x$m0 - avgMass)^2
      sumDiff <- 20910.98
      stdErr <- sqrt(MSE*(1 + 1/nObs + predMassDiff/sumDiff))
      
    }else{
      nObs <- 2190
      MSE <- 0.0003063309
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

