'''
Utility functions.
'''
import os
import logging
import numpy as np
import inspect
import time
import argparse
import ConfigParser
import collections
from FITSFile import read, write
from errors import errors
import pyfits
import math

def find_LT_files(dataPath, runNum, dithNum, minExpNum, maxExpNum, date, logger, errors):
    '''
    find FITS files in a directory matching LT nomenclature that satisfy the 
    following criteria:
     - multrun number = runNum
     - dither position = dithPos
     - date = date
     - exposure number between minExpNum and maxExpNum
    '''
    files = []
    for f in os.listdir(dataPath):
        if f.endswith(".fits"):
            file_fullpath = dataPath + f
            file_basename_noext = os.path.splitext(os.path.basename(file_fullpath))[0]
	    # assert filename consistency
            try:
                assert len(file_basename_noext.split('_')) == 7	
            except AssertionError:
                continue
            this_inst 		= str(file_basename_noext.split('_')[0])
            this_expType	= str(file_basename_noext.split('_')[1])
            this_date		= str(file_basename_noext.split('_')[2])
	    this_runNum		= int(file_basename_noext.split('_')[3])
	    this_dithNum	= int(file_basename_noext.split('_')[4])
	    this_expNum		= int(file_basename_noext.split('_')[5])
	    this_redLvl		= int(file_basename_noext.split('_')[6])
            # assert frame is from io:i and passes input criteria
            try:				
                assert this_inst == 'i'
                assert this_runNum == runNum
		assert this_dithNum == dithNum
		if date != '*':
                    assert this_date == date
		assert this_expNum >= minExpNum
                assert this_expNum <= maxExpNum
            except AssertionError:
                continue 

            files.append(os.path.abspath(file_fullpath))
            logger.info("[find_LT_files] File found: " + str(f))
    logger.info("[find_LT_files] Files found: " + str(len(files)))
    return files
  
def find_teledyne_files(dataPath, runNum, minGrpNum, maxGrpNum, minExpNum, maxExpNum, logger, errors):
    '''
    find FITS files in a directory matching Teledyne nomenclature that satisfy the 
    following criteria:
     - run number = runNum
     - group number between minGrpNum and maxGrpNum
     - exposure number between minExpNum and maxExpNum
    '''
    files = []
    for f in os.listdir(dataPath):
        if f.endswith(".fits"):
            file_fullpath = dataPath + f
            file_basename_noext = os.path.splitext(os.path.basename(file_fullpath))[0]
            # assert filename consistency
            try:
                assert len(file_basename_noext.split('_')) == 4 
            except AssertionError:
                continue
            this_inst           = str(file_basename_noext.split('_')[0])
            this_runNum         = int(file_basename_noext.split('_')[1].strip('R'))
            this_grpNum         = int(file_basename_noext.split('_')[2].strip('M'))
            this_expNum         = int(file_basename_noext.split('_')[3].strip('N'))
            # assert frame is from io:i
            try:                                
                assert this_inst == 'H2RG'
                assert this_runNum == runNum
                assert this_grpNum >= minGrpNum
                assert this_grpNum <= maxGrpNum
                assert this_expNum >= minExpNum
                assert this_expNum <= maxExpNum
            except AssertionError:
                continue 

            files.append(os.path.abspath(file_fullpath))
            logger.info("[find_teledyne_files] File found: " + str(f))
    logger.info("[find_teledyne_files] Files found: " + str(len(files)))
    return files  
  
def find_sort_files_LT(dataPath, minRunNum, maxRunNum, minDithNum, maxDithNum, minExpNum, maxExpNum, date, logger, errors):
    ''' 
    Looks for data conforming to LT nomenclature and returns a nested list with format [RUNNUM][DITHNUM].
    '''
    files_sorted = []
    
    n_files     = 0
    n_runs      = 0
    i_dither    = 0     # counter
    n_dithers   = []
    i_exp       = 0     # counter
    n_exps      = []
    for idx_1, runNum in enumerate(range(minRunNum, maxRunNum+1)):
        n_runs = idx_1
        files_sorted.append([])
        for idx_2, dithNum in enumerate(range(minDithNum, maxDithNum+1)):
            i_dithers = idx_2
            n_exps.append([])
            logger.info("[full_sort_LT_UTR] Trying runNum=" + str(runNum) + " and dithNum=" + str(dithNum))
            files = find_LT_files(dataPath, runNum, dithNum, minExpNum, maxExpNum, date, logger, errors)
            if len(files) != 0:
                files_sorted[idx_1].append(sort_files_by_FITS_key(files, "INTTIME", logger, errors))
                n_files = n_files + len(files)
                i_exp = len(files)
            n_exps[idx_2].append(i_exp)
            i_exp = 0
        n_dithers.append(i_dithers)
                
    if n_files == 0:
        errors.set_code(3, is_critical=False)
    return files_sorted, n_files, n_runs, n_dithers, n_exps
  
def find_sort_files_teledyne(dataPath, minRunNum, maxRunNum, minGrpNum, maxGrpNum, minExpNum, maxExpNum, logger, errors):
    ''' 
    Looks for data conforming to Teledyne nomenclature and returns data as a nested list with format [RUNNUM][0].
    '''
    files_sorted = []
    
    n_files     = 0
    n_runs      = 0
    i_dither    = 0          # counter
    n_dithers   = []
    i_exp       = 0          # counter
    n_exps      = []
    for idx_1, runNum in enumerate(range(minRunNum, maxRunNum+1)):
        n_runs = idx_1
        files_sorted.append([])
        logger.info("[full_sort_teledyne_UTR] Trying runNum=" + str(runNum))
        files = find_teledyne_files(dataPath, runNum, minGrpNum, maxGrpNum, minExpNum, maxExpNum, logger, errors)
        if len(files) != 0:
            files_sorted[idx_1].append(sort_files_by_FITS_key(files, "INTTIME", logger, errors))
            n_files = n_files + len(files)
            i_exp = len(files)
        n_exps.append(i_exp)  
        i_exp = 0
        
    if n_files == 0:
        errors.set_code(3, is_critical=False)   
    return files_sorted, n_files, n_runs, n_dithers, n_exps
  
def read_FITS_file(f, hdu=0):
    f_fits = read(f)
    f_fits.openFITSFile()
    f_fits.getHeaders(hdu)
    f_fits.getData(hdu)
 
    hdr     = f_fits.headers
    data    = f_fits.data 
    
    f_fits.closeFITSFile()

    return data, hdr
  
def read_ini(path):
    ini = ConfigParser.ConfigParser()
    ini.read(path)
    cfg = {}
    for section in ini.sections():
        cfg[section] = {}
        for option in ini.options(section):
            cfg[section][option] = str(ini.get(section, option))  
    return cfg
  
def sf(num, sig_figs):
    try:
        rtn = round(num, -int(math.floor(math.log10(abs(num))) - (sig_figs - 1)))
        return rtn
    except ValueError:
        return 0.  
 
def sort_files_by_FITS_key(files, key, logger, errors):
    key_values = []
    for idx, f in enumerate(files):
        logger.info("[sort_files_by_FITS_key] Sorting file " + str(idx+1) + " of " + str(len(files)))
        data, hdr = read_FITS_file(f)
        key_values.append(hdr[key])
    sorted_indexes 	= np.argsort(key_values)
    files_sorted 	= [files[i] for i in sorted_indexes]
    logger.info("[sort_files_by_FITS_key] Sorted " + str(len(files)) + " files by " + str(key) + " key.")
    return files_sorted  
  
def write_FITS_file(out, data, hdr=None, opt_hdr=None, allow_append=False): 
    if hdr is None:
        hdr = pyfits.Header(cards=[])
    if opt_hdr is not None:
        for k, v in opt_hdr.items():
            hdr[k] = v
    f_fits = write(out, data, hdr)
    f_fits.writeFITSFile(allow_append)