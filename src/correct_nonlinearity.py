import os
import numpy as np
import pyfits
from utility import read_FITS_file, write_FITS_file

class nonlinearity_correction:
    def __init__(self, logger, err):
        self.logger                     = logger
        self.err                        = err
        self.coeffs                     = []

    def _evaluate_correction(self, data, hdr, rates):
        if len(self.coeffs) == 0:
            self.err.set_code(9, is_critical=True)
        else:
            rate_corrected_data = data + rates*hdr['FRMTIME']
            data_cor = np.zeros(data.shape)
            for idx, c in enumerate(self.coeffs):
                try:
                    data_cor = data_cor + (c * pow(rate_corrected_data, idx))
                except ValueError:
                    self.err.set_code(23, is_critical=True)                
            return data_cor, hdr

    def execute(self, data, hdr, rates, hard, out=None, opt_hdr={}):
        if data is None:
            self.err.set_code(12, is_critical=True)     
        data, hdr = self._evaluate_correction(data, hdr, rates) 
        
        if hard:
            self.logger.info("[nonlinearity_correction.execute] Writing file out to " + out)          
            write_FITS_file(out, data, hdr, opt_hdr=opt_hdr)
        return data, hdr
     
    def read_lcor_coeffs(self, path, order):
        if not os.path.exists(path):
            self.err.set_code(7, is_critical=True)
        for i in range(order+1):
            self.logger.info("[nonlinearity_correction.execute] Appending order " + str(i))  
            data, hdr = read_FITS_file(path, i)
            self.coeffs.append(data)

         
      