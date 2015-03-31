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
            rate_corrected_data = data + rates*hdr['INTTIME']
            rate_linearity_corrected_data = np.zeros(data.shape)          
            for idx, c in enumerate(self.coeffs):
                try:
                    rate_linearity_corrected_data = rate_linearity_corrected_data + (c * pow(rate_corrected_data, idx))
                except ValueError:
                    self.err.set_code(23, is_critical=True)
            #
            # we need to derive correction factors (actual/observed) for rate corrected
            # data then apply this correction factor to the non-rate corrected, otherwise
            # if we were calculating it directly, we would be adding the rate to that  
            # observed.
            # 
            correction_factors = rate_linearity_corrected_data/rate_corrected_data
            data_cor = data * correction_factors
            return data_cor, hdr

    def execute(self, data, hdr, rates):
        if data is None:
            self.err.set_code(12, is_critical=True)     
        data, hdr = self._evaluate_correction(data, hdr, rates) 
        return data, hdr
     
    def read_lcor_coeffs(self, path, order, flip=False):
        if not os.path.exists(path):
            self.err.set_code(7, is_critical=True)
        for i in range(order+1):
            self.logger.info("[nonlinearity_correction.execute] Appending order " + str(i))  
            data, hdr = read_FITS_file(path, i)
            if flip:
                data = np.fliplr(data)
            self.coeffs.append(data)

         
      