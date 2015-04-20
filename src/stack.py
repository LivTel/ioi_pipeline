import copy
import numpy as np
import scipy as sp
import pyfits
from utility import read_FITS_file, write_FITS_file

class stack():
    def __init__(self, logger, err):
        self.logger          = logger
        self.err             = err     
        
    def _add_gain_keywords(self, hdr):
        if int(hdr['ASICGAIN']) == 8:
            hdr['GAIN'] = 3.05
            hdr['EPERDN'] = 3.05
        elif int(hdr['ASICGAIN']) == 10:
            hdr['ASICGAIN'] = 2.12
            hdr['EPERDN'] = 2.12
        elif int(hdr['ASICGAIN']) == 12:
            hdr['GAIN'] = 1.5
            hdr['EPERDN'] = 1.5
        elif int(hdr['ASICGAIN']) == 13:
            hdr['GAIN'] = 1.06
            hdr['EPERDN'] = 1.06            
        else:
            err.set_code(34, is_critical=False)
        return hdr

    def _coadd(self, datas, hdrs):
        data_stk = np.zeros(datas[0].shape)
        for d in datas:
            data_stk = data_stk + np.nan_to_num(d)           
        return data_stk, hdrs[0]
      
    def _correct_and_add_keywords_coadd(self, hdr, nframes):
        # EXPTIME is multiplied by number of frames
        hdr['EXPTIME'] = hdr['EXPTIME'] * nframes
        # add NCOADDS
        hdr['NCOADDS'] = (nframes, "Number of coadded frames")         
        return hdr
      
    def _correct_and_add_keywords_mean(self, hdr, nframes):
        # GAIN is multiplied by the number of dithers.
        hdr['GAIN'] = hdr['GAIN'] * nframes       
        hdr['EPERDN'] = hdr['EPERDN'] * nframes
        return hdr
      
    def _correct_and_add_keywords_mean_and_scale(self, hdr, nframes):
        # EXPTIME is multiplied by number of frames
        hdr['EXPTIME'] = hdr['EXPTIME'] * nframes
        # add NCOADDS
        hdr['NCOADDS'] = (nframes, "Number of coadded frames")
        return hdr      
        
    def execute(self, in_datas, in_hdrs, hard, method, out=None, opt_hdr={}):
        datas = copy.deepcopy(in_datas)
        hdrs  = copy.deepcopy(in_hdrs)
        if method == "MEAN":
            data, hdr = self._mean(datas, hdrs)
            self._add_gain_keywords(hdr)
            self._correct_and_add_keywords_mean(hdr, len(datas))
        elif method == "COADD":
            data, hdr = self._coadd(datas, hdrs)
            self._add_gain_keywords(hdr)  
            self._correct_and_add_keywords_coadd(hdr, len(datas))
        elif method == "MEANANDSCALE":
            data, hdr = self._mean_and_scale(datas, hdrs)
            self._add_gain_keywords(hdr)  
            self._correct_and_add_keywords_mean_and_scale(hdr, len(datas))            
        else:
            self.err.set_code(33, is_critical=True)
        
        if hard:
            self.logger.info("[stack.execute] Writing file out to " + out)
            write_FITS_file(out, data, hdr, opt_hdr=opt_hdr)
        return data, hdr
      
    def _mean_and_scale(self, datas, hdrs):
        data_stk, hdr = self._mean(datas, hdrs)     # do a mean stack
        data_stk = data_stk * len(datas)            # multiply by the number of dithers
        return data_stk, hdrs[0]
        
    def _mean(self, datas, hdrs):
        data_stk = np.nanmean(datas, axis=0)
        return data_stk, hdrs[0]  
 