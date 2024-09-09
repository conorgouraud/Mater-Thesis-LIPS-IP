import numpy as np
import pandas as pd
from rpy2.robjects import pandas2ri, r
from rpy2.robjects.packages import importr
import rpy2.robjects as ro

base = importr('base')
utils = importr('utils')
flexmix = importr('flexmix')

# function to search df for RRs in predicted range
def calculate_element_min_max(df, lipid_and_swiss_df, element, pred_interval):
    min_col_name = f"{element}minpred_{pred_interval}"
    max_col_name = f"{element}maxpred_{pred_interval}"

    df[min_col_name] = np.nan
    df[max_col_name] = np.nan

    for i in range(len(df)):
        lower_limit = df.iloc[i]['predLower']
        upper_limit = df.iloc[i]['predUpper']

        subset_db = lipid_and_swiss_df[(lipid_and_swiss_df['RR21'] >= lower_limit) & (lipid_and_swiss_df['RR21'] <= upper_limit)]

        if not subset_db.empty:
            df.at[df.index[i], min_col_name] = subset_db[element].min()
            df.at[df.index[i], max_col_name] = subset_db[element].max()
        else:
            df.at[df.index[i], min_col_name] = np.nan
            df.at[df.index[i], max_col_name] = np.nan

    return df

# Concomitant models written in R for now
ro.r('''
#functions
predRR21_agg_fc <- function(x, modelOutput, conf.int = TRUE, conf.level=conf.level, S=FALSE){
  
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

predRR32_agg_fc <- function(x, modelOutput, conf.int = TRUE, conf.level=conf.level, S=FALSE){
  
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

predRR43_agg_fc <- function(x, modelOutput, conf.int = TRUE, conf.level=conf.level, S=FALSE){
  
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

predRR54_agg_fc <- function(x, modelOutput, conf.int = TRUE, conf.level=conf.level, S=FALSE){
  
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
     
     
# ones without frac mass:
     
predRR21_agg_nfc <- function(x, modelOutput, conf.int = TRUE, conf.level=conf.level, S=FALSE){
  
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

predRR32_agg_nfc <- function(x, modelOutput, conf.int = TRUE, conf.level=conf.level, S=FALSE){
  
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

predRR43_agg_nfc <- function(x, modelOutput, conf.int = TRUE, conf.level=conf.level, S=FALSE){
  
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

predRR54_agg_nfc <- function(x, modelOutput, conf.int = TRUE, conf.level=conf.level, S=FALSE){
  
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

''')

predRR21_agg_fc = ro.globalenv['predRR21_agg_fc']
predRR32_agg_fc = ro.globalenv['predRR32_agg_fc']
predRR43_agg_fc = ro.globalenv['predRR43_agg_fc']
predRR54_agg_fc = ro.globalenv['predRR54_agg_fc']

predRR21_agg_nfc = ro.globalenv['predRR21_agg_nfc']
predRR32_agg_nfc = ro.globalenv['predRR32_agg_nfc']
predRR43_agg_nfc = ro.globalenv['predRR43_agg_nfc']
predRR54_agg_nfc = ro.globalenv['predRR54_agg_nfc']



