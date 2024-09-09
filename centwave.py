import rpy2.robjects as robjects
from rpy2.robjects.packages import importr
import pandas as pd
import gc 
import os

# centWave implemented in R
def centwave_on_raw_data(file_path, output_file):
    xcms = importr('xcms')
    MSnbase = importr('MSnbase')

    r_script = """
    process_mass_spec_file <- function(file_path) {
      ms_data <- readMSData(file_path, mode = "onDisk")
      
      cwp <- CentWaveParam(ppm = 5, peakwidth = c(6, 60), snthresh = 1,
                           prefilter = c(3, 5000), mzCenterFun = "wMean", integrate = 1, 
                           mzdiff = -0.001, fitgauss = FALSE, noise = 0,
                           verboseColumns = FALSE, firstBaselineCheck = TRUE)
      
      res <- findChromPeaks(ms_data, param = cwp)
      
      assign("result", res, envir = .GlobalEnv)
    }

    process_mass_spec_file
    """
    robjects.r(r_script)
    process_mass_spec_file = robjects.globalenv['process_mass_spec_file']
    process_mass_spec_file(file_path)

    robjects.r('''
      peaks <- chromPeaks(result)
      peaks_df <- as.data.frame(peaks)
    ''')

    try:
        peaks_df = robjects.globalenv['peaks_df']
        peaks_data = {col: peaks_df.rx2(col) for col in peaks_df.names}
        peaks_df = pd.DataFrame(peaks_data)
        peaks_df.to_csv(output_file, index=False)
        #print(f"Centwave results saved to {output_file}")
    except Exception as e:
        print(f"Error during conversion of peaks_df: {e}")
        raise
    finally:
        # clear the R environment and gc (garbage collection)
        robjects.r('rm(list = ls(all = TRUE))')
        robjects.r('gc()')
        gc.collect()

