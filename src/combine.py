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
        
    def _correct_nonlinearity(self, data, hdr, rates, lcor):
        rtn_data, rtn_hdr = lcor.execute(data=data, hdr=hdr, rates=rates)
        if rtn_data is None or rtn_hdr is None:
            self.err.set_code(6, is_critical=True)  
        return rtn_data, rtn_hdr
                   
    def execute(self, method, datas, hdrs, hard, lcor, out=None, f_pairs=1, opt_hdr={}):   
        if method == "CDS":    
            data, hdr, rates = self._perform_CDS(datas, hdrs, lcor);  
        elif method == "FOWLER":
            data, hdr, rates = self._perform_fowler(datas, hdrs, lcor, f_pairs); 
        else:
            self.err.set_code(10, is_critical=True) 
              
        if hard:
            self.logger.info("[combine.execute] Writing file out to " + out)          
            write_FITS_file(out, data, hdr, opt_hdr=opt_hdr)            
        return data, hdr, rates         
            
    def _perform_CDS(self, datas, hdrs, lcor):
        if datas is None or len(datas) <= 1:
            self.err.set_code(11, is_critical=True)
            
        # calculate rates
        delta = hdrs[1]['INTTIME'] - hdrs[0]['INTTIME']
        d_1 = datas[0]                   
        d_2 = datas[1]   
        rates = (d_2 - d_1)/delta       # this returns cts/s    
            
        # calculate CDS
        d_s = datas[0]          
        d_e = datas[len(datas)-1]    
        h_s = hdrs[0]
        h_e = hdrs[len(datas)-1]
        
        this_CDS  = d_e - d_s   
        hdr = hdrs[0]
        
        # add EXPTIME header
        hdr['EXPTIME'] = h_e['INTTIME'] - h_s['INTTIME']
        
        # apply nonlinearity correction if it's been requested
        if lcor is not None:
            if hdr['CLOCKING'] == 0:
                self.err.set_code(37, is_warning=True)
            this_CDS, hdr = self._correct_nonlinearity(data=this_CDS, hdr=hdr, rates=rates, lcor=lcor)  
                 
        # amend headers  
        hdr = self._purge_obsolete_headers(hdr)
        hdr = self._populate_DATE_keywords(hdr) 
        
        return this_CDS, hdr, rates
      
    def _perform_fowler(self, datas, hdrs, lcor, f_pairs):
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
            self.err.set_code(32, is_critical=True)                                 
        
        # calculate rates
        delta = d_2_inttime = hdrs[1]['INTTIME'] - hdrs[0]['INTTIME']
        d_1 = datas[0]                   
        d_2 = datas[1]   
        rates = (d_2 - d_1)/delta       # this returns cts/s
        
        # calculate CDS for each fowler pair
        datas_diff = []
        exptimes = []
        for idx in range(len(grp_1_idxs)):
            g1_idx = grp_1_idxs[idx]
            g2_idx = grp_2_idxs[idx]
            this_pair_CDS = datas[g2_idx] - datas[g1_idx]
            # apply nonlinearity correction if it's been requested
            if lcor is not None:
                if hdrs[g1_idx]['CLOCKING'] == 0:
                    self.err.set_code(37, is_warning=True)
                this_pair_CDS, hdrs[g1_idx] = self._correct_nonlinearity(data=this_pair_CDS, hdr=hdrs[g1_idx], rates=rates, lcor=lcor)  
            datas_diff.append(this_pair_CDS)
            
            exptimes.append(hdrs[g2_idx]['INTTIME'] - hdrs[g1_idx]['INTTIME']) 

        # average pairs
        data = np.mean(datas_diff[0:f_pairs], axis=0) 
        hdr = hdrs[0]
        
        # add EXPTIME key
        hdr['EXPTIME'] = np.mean(exptimes)
        
        # amend headers  
        hdr = self._purge_obsolete_headers(hdr)
        hdr = self._populate_DATE_keywords(hdr) 
     
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
        