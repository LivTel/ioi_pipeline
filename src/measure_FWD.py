'''
calculate FWD a single UTR sequence.
'''
import sys
import optparse
import os
import logging
import math

import scipy
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

from utility import read_FITS_file
from run_pipe import run_pipe
from errors import errors

# util functions
def sf(num, sig_figs):
    try:
        rtn = round(num, -int(math.floor(math.log10(abs(num))) - (sig_figs - 1)))
        return rtn
    except ValueError:
        return 0.

if __name__ == "__main__":
    parser = optparse.OptionParser()
    group1 = optparse.OptionGroup(parser, "General")
    group1.add_option('--p', action='store', default='/mnt/NAS/devel/IOI/images_and_analysis/remote_5/images/9/', dest='dataPath', type=str, help='path to data')
    group1.add_option('--wd', action='store', default='test', dest='workingDir', type=str, help='path to working directory')
    group1.add_option('--o', action='store_true', dest='clobber', help='clobber working directory?')
    group1.add_option('--pa', action='store', default='../config/paths_rmb.ini', type=str, dest='pathsCfgPath', help='path to paths config file')
    group1.add_option('--pi', action='store', default='../config/pipeline_fwd_rmb.ini', type=str, dest='pipeCfgPath', help='path to pipeline config file')     
    group1.add_option('--log', action='store', default='DEBUG', dest='logLevel', type=str, help='log level (DEBUG|INFO|WARNING|ERROR|CRITICAL)')      
    group1.add_option('--glo', action='store', default=0, type=int, dest='minGrpNum', help='lowest group number to use')
    group1.add_option('--ghi', action='store', default=20, type=int, dest='maxGrpNum', help='highest group number to use')
    parser.add_option_group(group1)
    
    options, args = parser.parse_args()
    params = {
        'dataPath' : str(options.dataPath).rstrip("/") + "/",
        'workingDir' : str(options.workingDir).rstrip("/") + "/",
        'clobber' : bool(options.clobber),       
        'pathsCfgPath' : str(options.pathsCfgPath),
        'pipeCfgPath' : str(options.pipeCfgPath),
        'minRunNum' : 1,
        'maxRunNum' : 1,
        'minDithNum' : 0,
        'maxDithNum' : 0,
        'minGrpNum' : int(options.minGrpNum),
        'maxGrpNum' : int(options.maxGrpNum),
        'minExpNum' : 1,
        'maxExpNum' : 1,        
        'logLevel' : str(options.logLevel.upper())
    }
    
    # console logging
    logger = logging.getLogger('run_pipe')
    logger.setLevel(getattr(logging, params['logLevel']))
    
    ## console handler
    ch = logging.StreamHandler()
    ch.setLevel(getattr(logging, params['logLevel']))

    ## set logging format
    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s", "%Y-%m-%d %H:%M:%S")
    ch.setFormatter(formatter)

    ## add handlers to logging object
    logger.addHandler(ch)
    
    # error handler

    err = errors(logger)   

    # create res directory to store metadata
    if os.path.exists(params['workingDir']) is True:
        if params['clobber'] is True:
            for i in os.listdir(params['workingDir']):    
                os.remove(params['workingDir'] + i)
            os.rmdir(params['workingDir'])
        else:
            err.set_code(1, is_critical=True)
    os.mkdir(params['workingDir'])

    # file loggingReduction
    fh = logging.FileHandler(params['workingDir'] + "res.log")
    fh.setLevel(logging.DEBUG)

    ## set logging format
    fh.setFormatter(formatter)

    ## add handlers to logging object
    logger.addHandler(fh)
    
    f_r1_data  = []
    f_r1_hdr   = []
        
    pipe = run_pipe(params, logger, err)
    pipe.go()   
        
    data        = pipe.session.file_data[0][0]
    hdr         = pipe.session.file_hdr[0][0]
    
    ## get ASIC gain and frmtime
    gain        = float(hdr[0].comments['ASICGAIN'].split('(')[1].split('dB')[0])
    frmtime     = float(hdr[0]['FRMTIME']) 
    
    for d in data:
        print np.mean(d)