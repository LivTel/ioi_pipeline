import os
import pyfits
import numpy as np
import time
from correct_nonlinearity import *
from utility import write_FITS_file

def jd_to_mjd(jd):
    return jd - 2400000.5

class combine():
    def __init__(self, logger, err):
        self.logger     = logger
        self.err        = err
       
    def execute(self, method, datas, hdrs, hard, hard_rates, out=None, out_rate=None, f_pairs=1, opt_hdr={}):   
        if method == "CDS":    
            data, hdr, rates = self._perform_CDS(datas, hdrs);  
        elif method == "FOWLER":
            data, hdr, rates = self._perform_fowler(datas, hdrs, f_pairs); 
        else:
            self.err.set_code(10, is_critical=True) 
              
        if hard:
            self.logger.info("[combine.execute] Writing file out to " + out)          
            write_FITS_file(out, data, hdr, opt_hdr=opt_hdr) 
        if hard_rates:
            self.logger.info("[combine.execute] Writing rate file out to " + out_rate)            
            write_FITS_file(out_rate, rates, hdr)             
        return data, hdr, rates           
            
    def _perform_CDS(self, datas, hdrs):
        if datas is None or len(datas) <= 1:
            self.err.set_code(11, is_critical=True)
            
        # remove first frame contribution from last
        d_s = datas[0]          
        d_e = datas[len(datas)-1]         
        data = d_e - d_s        
        
        delta = hdrs[1]['INTTIME'] - hdrs[0]['INTTIME']
        d_1 = datas[0]                   
        d_2 = datas[1]   
        rates = (d_2 - d_1)/delta       # this returns ct/s
        
        # amend headers  
        hdr = self._purge_obsolete_headers(hdrs[0])
        hdr = self._populate_DATE_keywords(hdrs[0])        
        
        return data, hdr, rates
      
    def _perform_fowler(self, datas, hdrs, f_pairs):
        if datas is None or len(datas) < 2*f_pairs:
            self.err.set_code(31, is_critical=True)
        
        # check the same number of files exist for both pedestal and read groups (1 and 2 respectively) and get indexes
        grp_1_idxs = []
        grp_2_idxs = []
        for pair in range(1, f_pairs+1):
            for idx in range(len(hdrs)):
                M = int(hdrs[idx]['SEQNUM_M'])
                N = int(hdrs[idx]['SEQNUM_N'])
                if M == 1 and N == pair:
                    grp_1_idxs.append(idx)
                if M == 2 and N == pair:
                    grp_2_idxs.append(idx)     
                    
        if len(grp_1_idxs) != len(grp_2_idxs):
            err.set_code(32, is_critical=True)                                 
           
        datas_diff = []
        for idx in range(len(grp_1_idxs)):
            g1_idx = grp_1_idxs[idx]
            g2_idx = grp_2_idxs[idx]
            datas_diff.append(datas[g2_idx] - datas[g1_idx])

        data = np.mean(datas_diff[0:f_pairs], axis=0) 
        
        delta = d_2_inttime = hdrs[1]['INTTIME'] - hdrs[0]['INTTIME']
        d_1 = datas[0]                   
        d_2 = datas[1]   
        rates = (d_2 - d_1)/delta       # this returns ct/s
        
        # amend headers  
        hdr = self._purge_obsolete_headers(hdrs[0])
        hdr = self._populate_DATE_keywords(hdrs[0])   
            
        return data, hdr, rates 
   
    def _purge_obsolete_headers(self, hdr):
        del hdr['SEQNUM_R']
        del hdr['SEQNUM_N']
        del hdr['INTTIME']
        del hdr['SEQNUM_M']
        del hdr['SEQNNAME']
        return hdr
        
    def _populate_DATE_keywords(self, hdr):    
        hdr['MJD']        = jd_to_mjd(hdr['ACQTIME'])
        hdr['DATE']       = time.strftime('%Y-%m-%d', time.strptime(hdr['ACQTIME1'], '%Y-%m-%d-%H:%M:%S.%f'))
        hdr['DATE-OBS']   = time.strftime('%Y-%m-%dT%H:%M:%S', time.strptime(hdr['ACQTIME1'], '%Y-%m-%d-%H:%M:%S.%f'))
        hdr['UTSTART']    = time.strftime('%H:%M:%S', time.strptime(hdr['ACQTIME1'], '%Y-%m-%d-%H:%M:%S.%f'))
        return hdr