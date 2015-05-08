import math
import numpy as np
import os
import sys
import copy
from utility import read_FITS_file, write_FITS_file
import alipy

class register():
    def __init__(self, ref_image_idx, datas, hdrs, mask_data, logger, err):
        self.ref_image_idx   = ref_image_idx
        self.datas           = datas
        self.hdrs            = hdrs
        self.mask_data       = mask_data
        self.logger          = logger
        self.err             = err
        
    def _cleanup(self, workDir):
        self.logger.info("[register._cleanup] Cleaning " + workDir + " of tmp_ files.")
        for i in os.listdir(workDir):
            if i.startswith("tmp_"):
                os.remove(workDir + i)
        
    def _set_nooverlap_region_as_nan(self, inFile, transform):
        ''' 
        required as regions that don't overlap are assigned 0 (scipy) or a v.small number (iraf)
        '''
        trans_x = int(round(transform[2]))
        trans_y = int(round(transform[3]))
        
        data, hdr = read_FITS_file(inFile)
        if trans_y >= 0:
            data[0:trans_y, :] = np.nan
        else:       
            data[data.shape[0]-1 + trans_y:data.shape[0]-1, :] = np.nan
        
        if trans_x >= 0:   
            data[:, 0:trans_x] = np.nan
        else:         
            data[:, data.shape[1]-1 + trans_x:data.shape[1]-1] = np.nan
            
        return data, hdr
      
    def _write_temporary_data_files(self, workDir):
        '''
        required as register procedure uses sExtractor (and maybe IRAF) which require hard copies of the files
        '''
        self.logger.info("[register._write_temporary_data_files] Writing out temporary files.")
        tmp_filenames = []
        for idx, d in enumerate(self.datas):
            this_tmp_outPath = workDir + "tmp_" + str(idx) + ".fits" 
            write_FITS_file(out=this_tmp_outPath, data=np.nan_to_num(self.datas[idx]), hdr=self.hdrs[idx]) # convert NaN to number or scipy breaks \
            tmp_filenames.append(this_tmp_outPath)    
        return tmp_filenames
      
    def _write_temporary_mask_files(self, workDir):
        '''
        required if bad pixel masking is set
        '''      
        self.logger.info("[register._write_temporary_data_files] Writing out temporary files.")
        tmp_mask_filenames = []
        bad_data = np.where(np.isnan(self.mask_data))
        self.mask_data.fill(0)                     # set all non-bad data to 0
        self.mask_data[bad_data] = 100             # set bad pixels to 100, that way, bp_post_reg_thresh represents the % influence that NaN has on nearby interpolated pixels
        tmp_mask_outPath = workDir + "tmp_mask.fits" 
        write_FITS_file(out=tmp_mask_outPath, data=np.nan_to_num(self.mask_data), hdr=None) # convert NaN to number or scipy breaks \
        for idx, d in enumerate(self.datas):   
            tmp_mask_filenames.append(tmp_mask_outPath)    
        return tmp_mask_filenames
      
    def execute(self, algorithm, workDir, fit_geom, bp_post_reg_thresh, hard, outs=None, mask_outs=None, keep_tmp=False):
        tmp_filenames = self._write_temporary_data_files(workDir)
        if mask_outs != None:
            tmp_mask_filenames = self._write_temporary_mask_files(workDir)
        try:
            identifications = alipy.ident.run(tmp_filenames[self.ref_image_idx], tmp_filenames, visu=False, verbose=False)
        except RuntimeError:
            self.err.set_code(18 ,is_critical=True)
            
        for id in identifications:
            if id.ok == True:
                self.logger.info("[register.execute] Found transform for " + str(id.ukn.name) + ". ")
                self.logger.info("[register.execute] " + str(id.trans))
            else:
                self.err.set_code(18, is_critical=True)

        datas = []
        hdrs  = []
        outputshape = self.datas[self.ref_image_idx].shape
        for idx, id in enumerate(identifications):
            if id.ok == True:
                tmp_outPath         = workDir + "tmp_" + str(idx) + ".trans.fits"
                tmp_mask_outPath    = workDir + "tmp_" + str(idx) + ".trans.mask.fits"
                # variant 1, using only scipy and the simple affine transform:
                if algorithm == "SCIPY":
                    # science data
                    alipy.align.affineremap(tmp_filenames[idx], id.trans, shape=outputshape, alifilepath=tmp_outPath, verbose=False, fitgeom=fit_geom)
                    self.logger.info("[register.execute] Created shifted frame " + os.path.basename(tmp_outPath) + ".")
                    # masks if set
                    if mask_outs != None:
                        alipy.align.affineremap(tmp_mask_filenames[idx], id.trans, shape=outputshape, alifilepath=tmp_mask_outPath, verbose=False, fitgeom=fit_geom)
                        self.logger.info("[register.execute] Created shifted bad pixel mask frame " + os.path.basename(tmp_mask_outPath) + ".")                    
                elif algorithm =="IRAF":
                # variant 2, using geomap/gregister:
                    import pyraf
                    alipy.align.irafalign(tmp_filenames[idx], id.uknmatchstars, id.refmatchstars, shape=outputshape, alifilepath=tmp_outPath, makepng=False, verbose=False, fitgeom=fit_geom) 
                    self.logger.info("[register.execute] Created shifted frame " + os.path.basename(tmp_outPath) + ".")
                   # masks if set
                    if mask_outs != None:
                        alipy.align.irafalign(tmp_mask_filenames[idx], id.uknmatchstars, id.refmatchstars, shape=outputshape, alifilepath=tmp_mask_outPath, makepng=False, verbose=False, fitgeom=fit_geom) 
                        self.logger.info("[register.execute] Created shifted bad pixel mask frame " + os.path.basename(tmp_mask_outPath) + ".")
                else:
                    self.err.set_code(19, is_critical=True)
  
                this_data_sci, this_hdr_sci = self._set_nooverlap_region_as_nan(tmp_outPath, id.trans.v)
                
                # reapply the transformed bad pixel mask to data
                if mask_outs != None:
                    this_data_mask, this_hdr_mask = read_FITS_file(tmp_mask_outPath)
                    new_mask = np.where(np.absolute(this_data_mask) > bp_post_reg_thresh)
                    this_data_sci[new_mask] = np.nan
                        
                datas.append(this_data_sci)
                hdrs.append(this_hdr_sci)
                if hard:
                    os.rename(tmp_outPath, outs[idx])                  
                    if len(mask_outs) != 0:
                         os.rename(tmp_mask_outPath, mask_outs[idx])    
                    
        if not keep_tmp:
            self._cleanup(workDir)             
        return datas, hdrs    
      