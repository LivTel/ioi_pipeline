import numpy as np
import statsmodels.api as sm
from utility import sf

def tukey_biweight(in_datas, scale, axis=0, max_iter=100, max_robust_unsolved_percentage=0.1, logger=None):
    # the robust estimator should be able to cope with the outliers from converting nans to num.
    location, solved = sm.robust.scale.norms.estimate_location(np.nan_to_num(in_datas), scale, norm=sm.robust.norms.TukeyBiweight(), maxiter=max_iter, axis=0, logger=logger) 
    percent_pixels_unsolved = float(float(len(np.where(solved == False)[0]))/float(solved.size))*100.
    if not np.all(solved):
        if logger is not None:
            logger.info("[statistics.tukey_biweight] " + str(sf(percent_pixels_unsolved, 2)) + "% of pixels remain unconverged of a maximum " + str(max_robust_unsolved_percentage) + "%")
        if percent_pixels_unsolved <= max_robust_unsolved_percentage:
            if logger is not None:
                logger.info("[statistics.tukey_biweight] Number of unconverged pixels within limits. Converting these pixels to NaN.") 
            location[np.where(solved==False)] = np.nan          
        else: 
            if logger is not None:
                logger.info("[statistics.tukey_biweight] Too few pixels failed to converge. Returning with Nonetype.")    
            return None
    return location
    
def huberT(in_datas, scale, axis=0, max_iter=100, max_robust_unsolved_percentage=0.1, logger=None):
    # the robust estimator should be able to cope with the outliers from converting nans to num.
    location, solved = sm.robust.scale.norms.estimate_location(np.nan_to_num(in_datas), scale, maxiter=max_iter, axis=0, logger=logger)
    percent_pixels_unsolved = float(float(len(np.where(solved == False)[0]))/float(solved.size))*100.
    if not np.all(solved):
        if logger is not None:
            logger.info("[statistics.huberT] " + str(sf(percent_pixels_unsolved)) + "% of pixels remain unconverged of a maximum " + str(max_robust_unsolved_percentage) + "%")
        if percent_pixels_unsolved <= max_robust_unsolved_percentage:
            if logger is not None:
                logger.info("[statistics.huberT] Number of unconverged pixels within limits. Converting these pixels to NaN.") 
            location[np.where(solved==False)] = np.nan          
        else: 
            if logger is not None:
                logger.info("[statistics.huberT] Too few pixels failed to converge. Returning with Nonetype.")    
            return None
    return location