import numpy as np
import scipy as sp
import pyfits 
from utility import read_FITS_file, write_FITS_file
import warnings
import scipy.stats as sps

class sky_subtraction():
    def __init__(self, logger, err):
        self.logger          = logger
        self.err             = err      
        self.data_sky        = None
        
    def interpolate_bad(self, smoothing_box_size, hard, out):
        bad_array = np.where(np.isnan(self.data_sky))
        bad_array_y = bad_array[0]
        bad_array_x = bad_array[1]
        n_bad = len(bad_array_x)
        data_interpolated = np.copy(self.data_sky)
        
        self.logger.info("[sky_subtraction.interpolate_bad] Interpolating over " + str(n_bad) + " bad sky pixel values") 
        for idx in range(n_bad):
            this_y_bad = bad_array_y[idx]
            this_x_bad = bad_array_x[idx]
         
            miny, maxy = max([this_y_bad-(smoothing_box_size-1)/2,0]), min([this_y_bad+(smoothing_box_size-1)/2,self.data_sky.shape[0]-1])
            minx, maxx = max([this_x_bad-(smoothing_box_size-1)/2,0]), min([this_x_bad+(smoothing_box_size-1)/2,self.data_sky.shape[1]-1])

            data_subset = self.data_sky[miny:maxy+1, minx:maxx+1]
            
            this_mean = sps.nanmean(data_subset.flatten())
            
            data_interpolated[this_y_bad, this_x_bad] = this_mean 

        bad_array = np.where(np.isnan(data_interpolated)) 
        bad_array_y = bad_array[0]
        bad_array_x = bad_array[1]
        n_bad = len(bad_array_x)        

        self.logger.info("[sky_subtraction.interpolate_bad] " + str(n_bad) + " bad sky pixel values remain in field") 
        if n_bad != 0:
            self.err.set_code(30, is_critical=False)      
                
        self.data_sky = data_interpolated
        
        if hard:
            self.logger.info("[sky_subtraction._interpolate_bad] Writing file out to " + out)               
            write_FITS_file(out, self.data_sky, hdr=None)   
        return self.data_sky
           
    def make_sky_frame(self, datas, sigma, hard, out):    
        data_tmp = []
        ref_sky_2d = None
        for idx, data_2d in enumerate(datas):   # for each dither (==data_2d)
            if idx == 0:                        # define "reference" sky frame as first frame
                ref_sky_2d = data_2d
            
            # make difference frame
            sky_diff = ref_sky_2d - data_2d
            
            # calculate offset as mean of clipped difference frame
            if idx == 0:        # this is the first frame (reference)
                offset = 0      # as ref_sky_2d == data_2d, offset is zero
            else:
                sky_diff_1d             = sky_diff.flatten()
                sky_diff_1d_nonan       = sky_diff_1d[np.logical_not(np.isnan(sky_diff_1d))]
                clipped_sky_values      = sp.stats.sigmaclip(sky_diff_1d_nonan, low=sigma, high=sigma)[0]        
                if len(clipped_sky_values) == 0:
                    self.err.set_code(25, is_critical=True) 
                offset = np.mean(clipped_sky_values)
                
            self.logger.info("[sky_subtraction.make_sky_frame] Offset between sky values of frame 1 and sky values of frame " + str(idx+1) + " is " + str(offset) + " ADU.")

            # apply this offset 
            data_tmp.append(data_2d+offset)
            
        self.data_sky = np.median(data_tmp, axis=0) 
        
        if hard:
            self.logger.info("[sky_subtraction.make_sky_frame] Writing file out to " + out)           
            write_FITS_file(out, self.data_sky, hdr=None)   
        return self.data_sky
        
    def sub_sky(self, data, hdr, sigma, hard, out, opt_hdr={}):
        if self.data_sky is None:
            self.err.set_code(14, is_critical=True) 
               
        # make difference frame
        sky_diff = self.data_sky - data
            
        # calculate offset as mean of clipped difference frame
        sky_diff_1d             = sky_diff.flatten()
        sky_diff_1d_nonan       = sky_diff_1d[np.logical_not(np.isnan(sky_diff_1d))]  
        clipped_sky_values      = sp.stats.sigmaclip(sky_diff_1d_nonan, low=sigma, high=sigma)[0]            
        if len(clipped_sky_values) == 0:
            self.err.set_code(25, is_critical=True) 
            
        # calculate offset as mean of clipped difference frame               
        offset = np.mean(clipped_sky_values)
        self.logger.info("[sky_subtraction.sub_sky] Offset between sky values of median sky frame and this frame is " + str(offset) + " ADU.") 
        
        # subtract off 
        data = data - (self.data_sky-offset)
        
        if hard:
            self.logger.info("[sky_subtraction.sub_sky] Writing file out to " + out)           
            write_FITS_file(out, data, hdr, opt_hdr=opt_hdr)
        return data, hdr
              