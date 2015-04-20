import copy
import pyfits 
import os
import numpy as np
from utility import read_FITS_file, write_FITS_file

class bad_pixel_mask:
    def __init__(self, logger, err):
        self.logger     = logger
        self.err        = err
        self.mask_data  = []  
        
    def _apply_bp_mask(self, data, hdr):
        if len(self.mask_data) == 0:
            self.err.set_code(28, is_critical=True)
        else:
            try: 
                data_masked = np.copy(data)
                for m in self.mask_data:
                    data_masked = data_masked * m                      
            except ValueError:
                self.err.set_code(29, is_critical=True)     
            return data_masked, hdr
          
    def combine_masks(self):
        if len(self.mask_data) < 2:
            self.err.set_code(20, is_critical=True)
        combined_mask_data = np.copy(self.mask_data[0])
        for i in range(1, len(self.mask_data)):
            if self.mask_data[0].shape != self.mask_data[i].shape:
                self.err.set_code(26, is_critical=True)
            combined_mask_data = combined_mask_data * self.mask_data[i]
        return combined_mask_data
 
    def execute(self, in_data, in_hdr, hard, out=None, opt_hdr={}):
        data = copy.deepcopy(in_data)
        hdr  = copy.deepcopy(in_hdr)
        if data is None:
            self.err.set_code(12, is_critical=True)
        data, hdr = self._apply_bp_mask(data, hdr)
        
        if hard:
            self.logger.info("[bad_pixel_mask.execute] Writing file out to " + out)
            write_FITS_file(out, data, hdr, opt_hdr=opt_hdr)
        return data, hdr           
           
    def read_mask(self, mask_path, flip=False):
        if not os.path.exists(mask_path):
            self.err.set_code(27, is_critical=True)
        data, hdr = read_FITS_file(mask_path)
        if flip:
            data = np.fliplr(data)
        self.mask_data.append(data)  
