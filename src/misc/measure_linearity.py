'''
    Analyse UTR sequence with and without nonlinearity correction.
    
    n.b. the counts measured here are after CDS subtraction.
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

sys.path.append("../")
from utility import read_FITS_file
from run_pipe import run_pipe
from errors import errors

if __name__ == "__main__":
    np.seterr(all='raise')
    
    parser = optparse.OptionParser()
    group1 = optparse.OptionGroup(parser, "General")
    group1.add_option('--p', action='store', default='/mnt/NAS/devel/IOI/images_and_analysis/remote_5/images/15/', dest='dataPath', type=str, help='path to data')
    group1.add_option('--wd', action='store', default='test', dest='workingDir', type=str, help='path to working directory')
    group1.add_option('--o', action='store_true', dest='clobber', help='clobber working directory?')
    group1.add_option('--pa', action='store', default='../../config/paths_rmb.ini', type=str, dest='pathsCfgPath', help='path to paths config file')    
    group1.add_option('--log', action='store', default='DEBUG', dest='logLevel', type=str, help='log level (DEBUG|INFO|WARNING|ERROR|CRITICAL)')      
    group1.add_option('--glo', action='store', default=0, type=int, dest='minGrpNum', help='lowest group number to use')
    group1.add_option('--ghi', action='store', default=50, type=int, dest='maxGrpNum', help='highest group number to use ')
    group1.add_option('--pl', action='store_true', dest='plt', help='make hard plot?') 
    parser.add_option_group(group1)
    
    options, args = parser.parse_args()
    params = {
        'dataPath' : str(options.dataPath).rstrip("/") + "/",
        'workingDir' : str(options.workingDir).rstrip("/") + "/",
        'clobber' : bool(options.clobber),       
        'pathsCfgPath' : str(options.pathsCfgPath),
        'minRunNum' : 1,
        'maxRunNum' : 1,
        'minDithNum' : 0,
        'maxDithNum' : 0,
        'minGrpNum' : int(options.minGrpNum),
        'maxGrpNum' : int(options.maxGrpNum),
        'minExpNum' : 1,
        'maxExpNum' : 16,        
        'logLevel' : str(options.logLevel.upper()),
        'plt' : bool(options.plt),
        'flip' : False
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
    
    # UNCORRECTED
    params['pipeCfgPath'] = "config/pipeline_lin_1.ini"
    data_mean_nolin     = []
    exptime_nolin       = []
    for i in range(params['minGrpNum'], params['maxGrpNum']-1):
        params['maxGrpNum'] = i+2
            
        pipe = run_pipe(params, logger, err)
        pipe.go()   
            
        data = pipe.session.file_data_nonss[0][0]
        rates = pipe.session.rates[0][0]
        exptime_nolin.append(pipe.session.file_hdr_nonss[0][0]['EXPTIME'])
        frmtime = pipe.session.file_hdr_nonss[0][0]['FRMTIME']
        missing_flux = rates*frmtime
        data_mean_nolin.append(np.nanmean(data+missing_flux))

    # CORRECTED
    params['pipeCfgPath'] = "config/pipeline_lin_2.ini"        
    data_mean_lin       = []
    exptime_lin         = []
    for i in range(params['minGrpNum'], params['maxGrpNum']-1):
        params['maxGrpNum'] = i+2
           
        pipe = run_pipe(params, logger, err)
        pipe.go()   
            
        data = pipe.session.file_data_nonss[0][0]
        rates = pipe.session.rates[0][0]
        exptime_lin.append(pipe.session.file_hdr_nonss[0][0]['EXPTIME'])
        frmtime = pipe.session.file_hdr_nonss[0][0]['FRMTIME']
        missing_flux = rates*frmtime
        data_mean_lin.append(np.nanmean(data+missing_flux))
     
    fig = plt.figure()
    plt.subplot(2,1,1)
    plt.plot(exptime_nolin, data_mean_nolin, 'rx-', label="no linearity correction")
    plt.plot(exptime_lin, data_mean_lin, 'bx-', label="linearity correction")

    # linearity line
    fitted_coeffs = np.polyfit([exptime_lin[0], exptime_lin[1]], [data_mean_nolin[0], data_mean_nolin[1]], 1)
    plt.plot([0, exptime_lin[-1]], np.polyval(fitted_coeffs, [0, exptime_nolin[-1]]), 'k--', label='linear')
    
    plt.ylim([0,66000])
    plt.title("Flux v EXPTIME")
    plt.xlabel("EXPTIME (s)")
    plt.ylabel("Mean counts after CDS (ADU)")
    plt.legend(loc='upper left', fontsize=10)
    
    plt.subplot(2,1,2)   
    plt.plot(data_mean_lin, (abs(np.polyval(fitted_coeffs, exptime_nolin)-data_mean_nolin)/np.polyval(fitted_coeffs, exptime_nolin))*100, 'rx-', label='residual nonlinearity % (uncorrected)') 
    plt.plot(data_mean_lin, (abs(np.polyval(fitted_coeffs, exptime_lin)-data_mean_lin)/np.polyval(fitted_coeffs, exptime_lin))*100, 'bx-', label='residual nonlinearity % (corrected)')   
    
    # 1% RESIDUAL nonlinearity line
    plt.plot(data_mean_lin, [1 for i in range(len(data_mean_lin))] , 'k--', label="1% nonlinearity")   
     
    plt.title("Residual nonlinearity")
    plt.xlabel("Mean counts after CDS (ADU)")
    plt.ylabel("Nonlinearity %")
    plt.legend(loc='upper left', fontsize=10)
    plt.yscale('log')
    plt.ylim([0.1,100])
     
    if params['plt']: 
        fig.tight_layout()
        plt.savefig("linearity.png")
    else:
        plt.show()
        