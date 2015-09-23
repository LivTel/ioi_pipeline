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
        self.file_ext                   = ''
        self.opt_hdr                    = collections.OrderedDict()
        
        self.file_data                  = []
        self.file_hdr                   = []
        
        self.rates                      = []    # cts/s
        
        self.file_data_nonss            = []
        self.file_hdr_nonss             = []
        self.file_opt_hdr_nonss         = []    # per ramp
        
        self.file_data_ss               = None  # this is deepcopied
        self.file_hdr_ss                = None  # this is deepcopied
        self.file_opt_hdr_ss            = []    # per ramp
        
        self.file_data_reg              = []
        self.file_hdr_reg               = []
        self.file_opt_hdr_reg           = []    # per ramp        
        
        self.file_data_stk              = []
        self.file_hdr_stk               = []
        self.file_opt_hdr_stk           = []    # per ramp
    
    def start(self, files):
        self.start_names = copy.deepcopy(files)  
        for idx_1, run in enumerate(files):  
            self.file_data.append([])
            self.file_hdr.append([])      
            self.file_data_nonss.append([])
            self.file_hdr_nonss.append([])  
            self.rates.append([])
            self.file_data_reg.append([])
            self.file_hdr_reg.append([])  
            self.file_data_stk.append([])
            self.file_hdr_stk.append([])              
            for idx_2, dither in enumerate(run):  
                self.file_data[idx_1].append([])
                self.file_hdr[idx_1].append([])
                self.file_data_nonss[idx_1].append([])
                self.file_hdr_nonss[idx_1].append([])  
                self.rates[idx_1].append([])
                self.file_data_reg[idx_1].append([])
                self.file_hdr_reg[idx_1].append([])  
                for idx_3, f in enumerate(dither):
                    data, hdr = read_FITS_file(f)                   
                    self.file_data[idx_1][idx_2].append(data)
                    self.file_hdr[idx_1][idx_2].append(hdr)  
                    self.logger.info("[session.read_files] Added file " + f + " to this session.")
                    
    def add_amend_opt_header(self, keyword, value, comment=''):
        if keyword in self.opt_hdr:
            self.logger.info("[session.add_amend_opt_header] Amending " + keyword + " keyword. Old value was '" + str(self.file_opt_hdr[keyword]) + "', new value is '(" + str(value) + ", " + comment + ")'")
        else:
            self.logger.info("[session.add_amend_opt_header] Adding " + keyword + " keyword. New value is '(" + str(value) + ", " + comment + ")'")
        self.opt_hdr[keyword] = (value, comment)           
               
    def copy_data_hdr_nonss_to_ss(self):
        self.file_data_ss = copy.deepcopy(self.file_data_nonss)
        self.file_hdr_ss  = copy.deepcopy(self.file_hdr_nonss)

    def free_file_data_and_hdr(self):
        self.file_data = None
        self.file_hdr  = None

    def free_rates(self):
        self.rates     = None
                                   
    def write_combined_data_as_LT(self, workingDir, extname, allow_append=False):  
        if extname == "IM_NONSS":
            data        = self.file_data_nonss
            hdr         = self.file_hdr_nonss
            opt_hdr     = self.file_opt_hdr_nonss
            for idx, i in enumerate(opt_hdr):
                opt_hdr[idx]['EXTNAME'] = extname
        elif extname == "IM_SS":
            data        = self.file_data_ss
            hdr         = self.file_hdr_ss
            opt_hdr     = self.file_opt_hdr_ss
            for idx, i in enumerate(opt_hdr):
                opt_hdr[idx]['EXTNAME'] = extname
        else:
            self.err.set_code(35, is_critical=True)
            
        n_files = 0  
        if data is not None:
            for idx_1, run in enumerate(data):
                for idx_2, f in enumerate(run): 
                    file_sections = os.path.basename(self.start_names[idx_1][idx_2][0]).split('_')
                    this_outPath = workingDir + '_'.join(file_sections[0:5]) + "_0_1.fits"
                    if len(data[idx_1][idx_2]) == 0:
                        pass   
                    else:
                        write_FITS_file(data=data[idx_1][idx_2], hdr=hdr[idx_1][idx_2], out=this_outPath, opt_hdr=opt_hdr[idx_1], allow_append=allow_append)
                        n_files = n_files + 1
        return n_files
                                 
    def write_stacked_data_as_LT(self, workingDir, extname, allow_append=False):
        if extname == "SK_SS" or extname == "SK_NONSS":  
            data        = self.file_data_stk
            hdr         = self.file_hdr_stk   
            opt_hdr     = self.file_opt_hdr_stk
            for idx, i in enumerate(opt_hdr):
                opt_hdr[idx]['EXTNAME'] = extname
        else:
            self.err.set_code(36, is_critical=True)
            
        n_files = 0    
        if data is not None:
            for idx_1, run in enumerate(data):
                file_sections = os.path.basename(self.start_names[idx_1][0][0]).split('_')
                this_outPath = workingDir + '_'.join(file_sections[0:4]) + "_0_0_2.fits"
                if len(data[idx_1]) == 0:
                    pass
                else:
                    write_FITS_file(data=data[idx_1], hdr=hdr[idx_1], out=this_outPath, opt_hdr=opt_hdr[idx_1], allow_append=allow_append)
                    n_files = n_files + 1 
        return n_files
