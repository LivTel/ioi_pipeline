import numpy as np
import statsmodels.api as sm

def tukey_biweight(in_datas, scale, axis=0, max_iter=100, logger=None):
    # the robust estimator should be able to cope with the outliers from converting nans to num.
    location, solved = sm.robust.scale.norms.estimate_location(np.nan_to_num(in_datas), scale, norm=sm.robust.norms.TukeyBiweight(), maxiter=max_iter, axis=0, logger=logger) 
    if not np.all(solved):
        return None 
    return location
    
def huberT(in_datas, scale, axis=0, max_iter=100, logger=None):
    # the robust estimator should be able to cope with the outliers from converting nans to num.
    location, solved = sm.robust.scale.norms.estimate_location(np.nan_to_num(in_datas), scale, maxiter=max_iter, axis=0, logger=logger)
    if not np.all(solved):
        return None
    return location