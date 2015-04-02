'''
    Analyse an UTR sequence for FWD.
'''
import sys
import optparse
import os
import logging
import math

import scipy
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

sys.path.append("../")
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

class measure_FWD:
    def __init__(self, params):
        self.params = params
        
    def go(self):
        # console logging
        logger = logging.getLogger('run_pipe')
        logger.setLevel(getattr(logging, self.params['logLevel']))
        
        ## console handler
        ch = logging.StreamHandler()
        ch.setLevel(getattr(logging, self.params['logLevel']))

        ## set logging format
        formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s", "%Y-%m-%d %H:%M:%S")
        ch.setFormatter(formatter)

        ## add handlers to logging object
        logger.addHandler(ch)
        
        # error handler
        err = errors(logger)   

        # create res directory to store metadata
        if os.path.exists(self.params['workingDir']) is True:
            if self.params['clobber'] is True:
                for i in os.listdir(self.params['workingDir']):    
                    os.remove(self.params['workingDir'] + i)
                os.rmdir(self.params['workingDir'])
            else:
                err.set_code(1, is_critical=True)
        os.mkdir(self.params['workingDir'])

        # file logging
        fh = logging.FileHandler(self.params['workingDir'] + "res.log")
        fh.setLevel(logging.DEBUG)

        ## set logging format
        fh.setFormatter(formatter)

        ## add handlers to logging object
        logger.addHandler(fh)
        
        # run pipeline
        self.params['pipeCfgPath'] = "config/pipeline_fwd.ini"
        pipe = run_pipe(self.params, logger, err)
        pipe.go()   
        
        # retrieve ref subtracted data/hdrs
        data        = pipe.session.file_data[0][0]
        hdr         = pipe.session.file_hdr[0][0]
        
        # get ASICGAIN and FRMTIME
        gain        = float(hdr[0].comments['ASICGAIN'].split('(')[1].split('dB')[0])
        frmtime     = float(hdr[0]['FRMTIME']) 
        
        # establish means of frames for each INTTIME
        data_mean   = []
        inttime     = []
        for idx, d in enumerate(data):
            data_mean.append(np.mean(data[idx]))
            inttime.append(hdr[idx]['INTTIME'])
    
        # plot
        font = {'family' : 'normal',
                'weight' : 'normal',
                'size'   : 12}
        matplotlib.rc('font', **font)
        
        plt.plot(inttime, data_mean, 'k-')
        plt.plot([0, max(inttime)], [np.min(data_mean) for i in [0, max(inttime)]], 'k--', label='Mean Bias = ' + str(np.min(data_mean)) + 'ADU')
        plt.plot([0, max(inttime)], [np.max(data_mean) for i in [0, max(inttime)]], 'k:', label='Mean FWD = ' + str(np.max(data_mean)) + 'ADU') 
            
        plt.title("FWD (" + str(gain) + "dB)")
        plt.xlabel("INTTIME (s)")
        plt.ylabel("Mean counts")
        plt.xlim([0, max(inttime)])    
        plt.ylim([0,66000])
        plt.legend(loc='upper left', fontsize=11)
        
        if self.params['fits']:     
            plt.savefig("fwd.png")
        
        logger.info("Calculated FWD is " + str(np.max(data_mean) - np.min(data_mean)) + "ADU.")

if __name__ == "__main__":
    parser = optparse.OptionParser()
    group1 = optparse.OptionGroup(parser, "General")
    group1.add_option('--p', action='store', default='/mnt/NAS/devel/IOI/images_and_analysis/remote_5/images/15/', dest='dataPath', type=str, help='path to data')
    group1.add_option('--wd', action='store', default='test', dest='workingDir', type=str, help='path to working directory')
    group1.add_option('--o', action='store_true', dest='clobber', help='clobber working directory?')
    group1.add_option('--pa', action='store', default='../../config/paths_rmb.ini', type=str, dest='pathsCfgPath', help='path to paths config file')
    group1.add_option('--log', action='store', default='DEBUG', dest='logLevel', type=str, help='log level (DEBUG|INFO|WARNING|ERROR|CRITICAL)')      
    group1.add_option('--glo', action='store', default=0, type=int, dest='minGrpNum', help='lowest group number to use')
    group1.add_option('--ghi', action='store', default=20, type=int, dest='maxGrpNum', help='highest group number to use')
    group1.add_option('--hard', action='store_true', dest='fits', help='make plot?') 
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
        'maxExpNum' : 1,        
        'logLevel' : str(options.logLevel.upper()),
        'fits' : bool(options.fits),
        'flip' : False
    }
    
    m = measure_FWD(params)
    m.go()
        