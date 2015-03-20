import copy
import collections
import os
import numpy as np
from utility import read_FITS_file, write_FITS_file
import pyfits

class session:
    def __init__(self, logger, err):
        self.logger                     = logger
        self.err                        = err
        self.start_names                = None
        self.opt_hdr                    = collections.OrderedDict()
        self.file_data                  = None
        self.file_hdr                   = None
        self.file_data_nonss            = None
        self.file_hdr_nonss             = None
        self.file_opt_hdr_nonss         = None
        self.file_data_ss               = None
        self.file_hdr_ss                = None
        self.file_opt_hdr_ss            = None
        self.file_data_ss_stk           = None
        self.file_hdr_ss_stk            = None
        self.file_opt_hdr_ss_stk        = None
        self.file_ext                   = ''
        self.rates                      = None # cts/s
    
    def add_files(self, files):
        self.file_data   = copy.deepcopy(files)
        self.file_hdr    = copy.deepcopy(files)
        self.start_names = copy.deepcopy(files)        
        for idx_1, run in enumerate(files):          
            for idx_2, dither in enumerate(run):          
                for idx_3, f in enumerate(dither):
                    self.start_names[idx_1][idx_2][idx_3] = f
                    data, hdr = read_FITS_file(f) 
                    self.file_data[idx_1][idx_2][idx_3] = data
                    self.file_hdr[idx_1][idx_2][idx_3]  = hdr   
                    self.logger.info("[session.read_files] Added file " + f + " to this session.")
                    
    def add_amend_opt_header(self, keyword, value, comment=''):
        if keyword in self.opt_hdr:
            self.logger.info("[session.add_amend_opt_header] Amending " + keyword + " keyword. Old value was '" + str(self.file_opt_hdr[keyword]) + "', new value is '(" + str(value) + ", " + comment + ")'")
        else:
            self.logger.info("[session.add_amend_opt_header] Adding " + keyword + " keyword. New value is '(" + str(value) + ", " + comment + ")'")
        self.opt_hdr[keyword] = (value, comment)           
        
    def set_opt_header_nonss(self):
        self.file_opt_hdr_nonss = self.opt_hdr
        
    def set_opt_header_ss(self):
        self.file_opt_hdr_ss = self.opt_hdr   
        
    def set_opt_header_ss_stk(self):
        self.file_opt_hdr_ss_stk = self.opt_hdr   
   
    def set_session_vars_post_combine(self, datas, hdrs, rates):
        self.file_data_nonss = copy.deepcopy(datas)
        self.file_hdr_nonss  = copy.deepcopy(hdrs)
        self.file_data_ss    = []          
        self.file_hdr_ss     = []  
        self.rates                = rates
        for idx_1, run in enumerate(datas):
            self.file_data_ss.append([])
            self.file_hdr_ss.append([])
            for idx_2, dither in enumerate(run):
                self.file_data_ss[idx_1].append([])
                self.file_hdr_ss[idx_1].append([])
                
    def set_session_vars_post_stacking(self, datas, hdrs):
        self.file_data_ss_stk = copy.deepcopy(datas)
        self.file_hdr_ss_stk  = copy.deepcopy(hdrs)   
                                   
    def write_combined_data_as_LT(self, workingDir, extname, allow_append=False):  
        if extname == "IM_NONSS":
            data        = self.file_data_nonss
            hdr         = self.file_hdr_nonss
            opt_hdr     = self.file_opt_hdr_nonss
            self.opt_hdr['EXTNAME'] = extname
        elif extname == "IM_SS":
            data        = self.file_data_ss
            hdr         = self.file_hdr_ss
            opt_hdr     = self.file_opt_hdr_ss
            self.opt_hdr['EXTNAME'] = extname
        else:
            self.err.set_code(35, is_critical=True)
            
        for idx_1, run in enumerate(data):
            for idx_2, f in enumerate(run): 
                file_sections = os.path.basename(self.start_names[idx_1][idx_2][0]).split('_')
                this_outPath = workingDir + '_'.join(file_sections[0:5]) + "_0_1.fits"
                if data[idx_1][idx_2] is None:
                    pass   
                else:
                    write_FITS_file(data=data[idx_1][idx_2], hdr=hdr[idx_1][idx_2], out=this_outPath, opt_hdr=self.opt_hdr, allow_append=allow_append)
                                 
    def write_stacked_data_as_LT(self, workingDir, extname, allow_append=False):
        if extname == "SK_SS":  
            data        = self.file_data_ss_stk
            hdr         = self.file_hdr_ss_stk   
            opt_hdr     = self.file_opt_hdr_ss_stk
            self.opt_hdr['EXTNAME'] = extname
        else:
            self.err.set_code(36, is_critical=True)
            
        for idx_1, run in enumerate(data):
            file_sections = os.path.basename(self.start_names[idx_1][0][0]).split('_')
            this_outPath = workingDir + '_'.join(file_sections[0:4]) + "_0_0_1.fits"
            if data[idx_1] is None:
                pass
            else:
                write_FITS_file(data=data[idx_1], hdr=hdr[idx_1], out=this_outPath, opt_hdr=self.opt_hdr, allow_append=allow_append)
