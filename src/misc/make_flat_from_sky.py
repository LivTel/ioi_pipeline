'''
creates a flatfield and corresponding bad pixel map from sky frames.
'''
import sys
import argparse
import logging
import os
import numpy as np
import scipy.stats as sps

sys.path.append("/usr/local/bin/ioi_pipeline/src")
sys.path.append("../")
from utility import read_FITS_file, write_FITS_file

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--f', action='append', dest='paths', help='paths to sky files')
    parser.add_argument('--btl', action='store', dest='badThreshLower', default=0.35, help='lower coefficient threshold to be flagged as bad pixel')
    parser.add_argument('--btu', action='store', dest='badThreshUpper', default=2.0, help='upper coefficient threshold to be flagged as bad pixel')
    parser.add_argument('--log', action='store', default='INFO', dest='logLevel', type=str, help='level (DEBUG|INFO|WARNING|ERROR|CRITICAL)') 
    parser.add_argument('--fits', action='store_true', dest='fits', help='make flat and bad pixel map FITS files')
    parser.add_argument('--fl', action='store_true', dest='flip', help='flip array on write?')
    parser.add_argument('--hdu', action='store', dest='hdu', default=0, type=int, help='hdu index to use')
    parser.add_argument('--o', action='store_true', dest='clobber', help='clobber existing files on write?')
    parser.add_argument('--sig', action='store', dest='sigclip', default=10, type=int, help='sigma to use for clipping dataset prior to mean calculation')
    
    args = parser.parse_args()
    params = {
        'paths' : args.paths,
        'badThreshLower' : float(args.badThreshLower),
        'badThreshUpper' : float(args.badThreshUpper),
        'logLevel' : str(args.logLevel),
        'fits' : bool(args.fits),
        'flip' : bool(args.flip),
        'hdu' : int(args.hdu),
        'clobber' : bool(args.clobber),
        'sigclip' : int(args.sigclip)
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
            logger.info("reading file " + i + " at HDU " + str(params['hdu']))
            data, hdr = read_FITS_file(i, hdu=params['hdu'])
            logger.info("converting NaNs to numbers")
            data_nonan = np.nan_to_num(data)
            logger.info("clipping data with sigma of " + str(params['sigclip']))
            data_nonan_sigclip = sps.sigmaclip(data_nonan, low=params['sigclip'], high=params['sigclip'])[0]
            data_nonan_sigclip_mean = np.mean(data_nonan_sigclip)
            logger.info("mean of data pre-normalisation is " + str(data_nonan_sigclip_mean))
            logger.info("normalising data") 
            data_cor = data_nonan/data_nonan_sigclip_mean
            logger.info("mean of data post-normalisation is " + str(np.mean(data_cor))) 
            datas_cor.append(data_cor)
        else:
            logger.warning("path " + i + " does not exist. skipping.")        

    # create master flat if more than one file specified
    flat = np.mean(datas_cor, axis=0) 
        
    # make flat and set "bad pixels" to 1        
    bad_array   = np.where(np.logical_or(flat < params['badThreshLower'], flat > params['badThreshUpper']))   # establish bad pixels
    flat[bad_array] = 1 
    if params['fits']:
        if params['clobber']:
            if os.path.exists("flat.fits"):
                os.remove("flat.fits")
                logger.info("clobber set. removed flat.fits file")
        logger.info("writing flat")   
        if params['flip']:
            flat = np.fliplr(flat)
        try:
            write_FITS_file(out="flat.fits", data=flat, hdr=None)
        except IOError:
            logger.info("flat.fits file already exists, and hasn't been clobbered")
           
    
    # make bad pixel map
    bad         = np.ones(flat.shape)   
    bad[bad_array] = np.nan
    if params['fits']:   
        if params['clobber']:
            if os.path.exists("flat_bad.fits"):
                os.remove("flat_bad.fits")
                logger.info("clobber set. removed flat_bad.fits file") 
        logger.info("writing bad pixel map")    
        if params['flip']:
            bad = np.fliplr(bad)
        try:
            write_FITS_file(out="flat_bad.fits", data=bad, hdr=None)
        except IOError:
            logger.info("flat_bad.fits file already exists, and hasn't been clobbered") 