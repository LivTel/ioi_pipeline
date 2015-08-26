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
import pylab as plt
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
        
        ## process section argument
        x1 = int(self.params['section'].split(',')[0])
        y1 = int(self.params['section'].split(',')[1])
        x2 = int(self.params['section'].split(',')[2])
        y2 = int(self.params['section'].split(',')[3])        
        
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
            data_mean.append(np.mean(data[idx][y1:y2,x1:x2]))
            inttime.append(hdr[idx]['INTTIME'])
    
        # plot
        font = {'family' : 'normal',
                'weight' : 'normal',
                'size'   : 12}
        matplotlib.rc('font', **font)
        
        plt.plot(inttime, data_mean, 'kx-')
        plt.plot([0, max(inttime)], [np.min(data_mean) for i in [0, max(inttime)]], 'k--', linewidth=3, label='Mean Bias = ' + str(np.min(data_mean)) + 'ADU')
        plt.plot([0, max(inttime)], [np.max(data_mean) for i in [0, max(inttime)]], 'k--', label='Mean Sat. = ' + str(np.max(data_mean)) + 'ADU') 
        
        if self.params['linLine']:
            # linearity line
            fitted_coeffs = np.polyfit([inttime[0], inttime[1]], [data_mean[0], data_mean[1]], 1)
            plt.plot([0, inttime[-1]], np.polyval(fitted_coeffs, [0, inttime[-1]]), 'k:', label='Linearity')
            
        plt.title("FWD (" + str(gain) + "dB)")
        plt.xlabel("EXPTIME (s)")
        plt.ylabel("Mean counts (ADU)")
        plt.xlim([0, max(inttime)])    
        plt.ylim([0,60000])
        plt.legend(loc='lower right', fontsize=11)
        
        plt.ticklabel_format(axis='y',style='sci',scilimits=(1,3))   
        
        if self.params['hard']: 
            plt.savefig("fwd.eps")        
        
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
    group1.add_option('--hard', action='store_true', dest='hard', help='make hard plot?') 
    group1.add_option('--s', action='store', dest='section', default='0,0,2047,2047', help='csv section of image to use [x1,y1,x2,y2]')
    group1.add_option('--l', action='store_true', dest='linLine', help='add linearity line?')
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
        'hard' : bool(options.hard),
        'flip' : False,
        'section' : str(options.section),
        'linLine' : bool(options.linLine)
    }
    
    m = measure_FWD(params)
    m.go()
        