'''
creates a flatfield and corresponding bad pixel map from sky frames.
'''
import sys
import argparse
import logging
import os
import numpy as np

sys.path.append("../")
from utility import read_FITS_file, write_FITS_file

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--f', action='append', dest='paths', help='paths to sky files')
    parser.add_argument('--bt', action='store', dest='badThresh', default=0.35, help='coefficient threshold to be flagged as bad pixel')
    parser.add_argument('--log', action='store', default='INFO', dest='logLevel', type=str, help='level (DEBUG|INFO|WARNING|ERROR|CRITICAL)') 
    parser.add_argument('--fits', action='store_true', dest='fits', help='make flat and bad pixel map FITS files')
    
    args = parser.parse_args()
    params = {
        'paths' : args.paths,
        'badThresh' : float(args.badThresh),
        'logLevel' : str(args.logLevel),
        'fits' : str(args.fits)
    }    
    
    # console logging
    logger = logging.getLogger('measure_flatfield')
    logger.setLevel(getattr(logging, params['logLevel']))
    
    ## console handler
    ch = logging.StreamHandler()
    ch.setLevel(getattr(logging, params['logLevel']))

    ## set logging format
    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s", "%Y-%m-%d %H:%M:%S")
    ch.setFormatter(formatter)

    ## add handlers to logging object
    logger.addHandler(ch)    
    
    if params['paths'] == None:
        logger.critical("[paths] is empty. quitting.")
        exit(0)
    
    datas_cor = []
    for i in params['paths']:
        if os.path.exists(i):
            logger.info("reading file " + i)
            data, hdr = read_FITS_file(i)
            print i
            data_mean = np.nanmean(data)
            logger.info("mean of data pre-normalisation is " + str(data_mean))
            logger.info("normalising data") 
            data_cor = data/data_mean
            logger.info("mean of data post-normalisation is " + str(np.nanmean(data_cor))) 
            datas_cor.append(data_cor)
        else:
            logger.warning("path " + i + " does not exist. skipping.")        
    
    # make flat and set "bad pixels" to 1
    flat = np.nanmean(datas_cor, axis=0) 
    bad_array   = np.where(flat < params['badThresh'])   # establish bad pixels
    flat[bad_array] = 1 
    if params['fits']:
        logger.info("writing flat")    
        write_FITS_file(out="flat.fits", data=flat, hdr=None)
    
    # make bad pixel map
    bad         = np.ones(flat.shape)   
    bad[bad_array] = np.nan
    if params['fits']:   
        logger.info("writing bad pixel map")    
        write_FITS_file(out="flat_bad.fits", data=bad, hdr=None)    
        
