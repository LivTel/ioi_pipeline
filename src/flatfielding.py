import copy
import pyfits 
import os
import numpy as np
from scipy import signal
from utility import read_FITS_file, write_FITS_file

class flatfielding:
    def __init__(self, logger, err):
        self.logger     = logger
        self.err        = err
        self.coeffs     = None 
        
    def execute(self, in_data, in_hdr, hard, out=None, opt_hdr={}):
        data = copy.deepcopy(in_data)
        hdr  = copy.deepcopy(in_hdr)
        if data is None:
            self.err.set_code(12, is_critical=True)
        data, hdr = self._perform_correction(data, hdr)
        
        if hard:
            self.logger.info("[flatfielding.execute] Writing file out to " + out)      
            write_FITS_file(out, data, hdr, opt_hdr=opt_hdr)
        return data, hdr    
      
    def _perform_correction(self, data, hdr):
        if self.coeffs is None:
            self.err.set_code(22, is_critical=True)
        else:
            try: 
                data_cor = data / self.coeffs
            except ValueError:
                self.err.set_code(24, is_critical=True)     
            return data_cor, hdr           
        
    def read_fcor_coeffs(self, path, flip=False):
        if not os.path.exists(path):
            self.err.set_code(21, is_critical=True)
        data, hdr = read_FITS_file(path)
        if flip:
            data = np.fliplr(data)
        self.coeffs = data        
