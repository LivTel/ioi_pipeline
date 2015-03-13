'''
calculate gain and read noise from two UTR sequences by making a variance v signal plot.
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

def calc_theoretical_gain(gain):
    A = pow(10, (gain/20))
    g = 3 * (40*pow(10, -15))/(pow(2, 16) * A * 1.6 * pow(10, -19))
    return g

if __name__ == "__main__":
    parser = optparse.OptionParser()
    group1 = optparse.OptionGroup(parser, "General")
    group1.add_option('--p', action='store', default='/mnt/NAS/devel/IOI/images_and_analysis/remote_5/images/7/', dest='dataPath', type=str, help='path to data')
    group1.add_option('--wd', action='store', default='test', dest='workingDir', type=str, help='path to working directory')
    group1.add_option('--o', action='store_true', dest='clobber', help='clobber working directory?')
    group1.add_option('--pa', action='store', default='../config/paths_rmb.ini', type=str, dest='pathsCfgPath', help='path to paths config file')
    group1.add_option('--pi', action='store', default='../config/pipeline_gain.ini', type=str, dest='pipeCfgPath', help='path to pipeline config file')     
    group1.add_option('--log', action='store', default='DEBUG', dest='logLevel', type=str, help='log level (DEBUG|INFO|WARNING|ERROR|CRITICAL)')  
    group1.add_option('--s', action='store', default='900,900,1000,1000', type=str, dest='sect', help='section (x1,y1,x2,y2)')    
    group1.add_option('--glo', action='store', default=0, type=int, dest='minGrpNum', help='lowest group number to use')
    group1.add_option('--ghi', action='store', default=25, type=int, dest='maxGrpNum', help='highest group number to use')
    parser.add_option_group(group1)
    
    options, args = parser.parse_args()
    params = {
        'dataPath' : str(options.dataPath).rstrip("/") + "/",
        'workingDir' : str(options.workingDir).rstrip("/") + "/",
        'clobber' : bool(options.clobber),       
        'pathsCfgPath' : str(options.pathsCfgPath),
        'pipeCfgPath' : str(options.pipeCfgPath),
        'minRunNum' : 1,
        'maxRunNum' : 2,
        'minDithNum' : 0,
        'maxDithNum' : 0,
        'minGrpNum' : int(options.minGrpNum),
        'maxGrpNum' : int(options.maxGrpNum),
        'minExpNum' : 1,
        'maxExpNum' : 1,        
        'logLevel' : str(options.logLevel.upper()),
        'sect' : str(options.sect),
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
    f_r1_rates = []
    f_r2_data  = []
    f_r2_hdr   = []
    f_r2_rates = []
    for i in range(params['minGrpNum'], params['maxGrpNum']-1):
        params['minGrpNum'] = params['minGrpNum']
        params['maxGrpNum'] = i+2
        
        pipe = run_pipe(params, logger, err)
        pipe.go()   
        
        f_r1_data.append(pipe.session.file_data[0][0])
        f_r1_hdr.append(pipe.session.file_hdr[0][0])
        f_r1_rates.append(pipe.session.rates[0][0])
        f_r2_data.append(pipe.session.file_data[1][0])
        f_r2_hdr.append(pipe.session.file_hdr[1][0]) 
        f_r2_rates.append(pipe.session.rates[1][0])
     
    ## get ASIC gain and frmtime
    gain        = float(pipe.session.file_hdr[0][0].comments['ASICGAIN'].split('(')[1].split('dB')[0])
    frmtime     = float(pipe.session.file_hdr[0][0]['FRMTIME']) 
        
    # section of image to use
    sect_x = (int(params['sect'].split(',')[0].strip()), int(params['sect'].split(',')[2].strip()))
    sect_y = (int(params['sect'].split(',')[1].strip()), int(params['sect'].split(',')[3].strip()))

    if len(f_r1_data) != len(f_r2_data):
        print "Ramps must have equal number of frames."
        exit(0)      
   
    # calculate signal mean and variance of difference frame 
    means_average       = []
    shot_and_read_noise = []
    for i in range(0, len(f_r1_data), 1):
        print "Processing file " + str(i+1) + " of " + str(len(f_r1_data))

        ## section files
        f_sect_r1 = f_r1_data[i][sect_y[0]:sect_y[1]+1, sect_x[0]:sect_x[1]+1] + f_r1_rates[i][sect_y[0]:sect_y[1]+1, sect_x[0]:sect_x[1]+1]*frmtime
        f_sect_r2 = f_r2_data[i][sect_y[0]:sect_y[1]+1, sect_x[0]:sect_x[1]+1] + f_r2_rates[i][sect_y[0]:sect_y[1]+1, sect_x[0]:sect_x[1]+1]*frmtime   
        
        ## get means/stds of sections for each file
        this_mean_r1 = np.mean(f_sect_r1)
        this_mean_r2 = np.mean(f_sect_r2)
        this_std_r1  = np.std(f_sect_r1)
        this_std_r2  = np.std(f_sect_r2)

        ## calculate average signal and average (total) noise
        this_mean_average = (this_mean_r1 + this_mean_r2)/2.    # this needs to be rate adjusted. CDS loses one frame time of signal!
        this_std_average = (this_std_r1 + this_std_r2)/2.
        means_average.append(this_mean_average)
 
        ## calculate std of diff frame to remove FPN component
        diff = (f_sect_r1 + 1000) - f_sect_r2   # add 1000 ADU to avoid negative offsets
        this_std_diff = np.std(diff)

        ## calculate shot_and_read_noise component for single frame
        this_shot_and_read_noise = this_std_diff/pow(2, 0.5)    # taking the difference (CDS) increases noise by sqrt(2)
        shot_and_read_noise.append(this_shot_and_read_noise)
        
        print '\t' + 'mean:\t\t' + str(this_mean_average)
        print '\t' + 'sd(shot+rd):\t' + str(this_shot_and_read_noise)
    
    # get variance 
    vars_diff = [pow(i, 2) for i in shot_and_read_noise]

    # plot the result
    plt.figure()

    fitted_coeffs = np.polyfit(means_average, vars_diff, 1)
    fitted_total_noise = np.polyval(fitted_coeffs, means_average)

    calc_gain           = 1./fitted_coeffs[0]
    y_intercept         = fitted_coeffs[1]
    y_intercept_e       = y_intercept * calc_gain
    rn                  = pow(y_intercept, 0.5)*calc_gain
    rn_cds              = rn * pow(2, 0.5)
    
    plt.plot(means_average, vars_diff, 'rx')
    plt.plot(means_average, fitted_total_noise, 'r-', label='gain: ' + str(gain) + 'dB')

    plt.title("gain plot")
    plt.legend(loc='upper left')
    plt.xlabel("Signal (ADU)")
    plt.ylabel("Variance (ADU)")
    plt.annotate('calc_gain = ' + str(sf(calc_gain, 2)) + ' e-/ADU', fontsize=10, xy=(0.5, 0.3), xycoords='axes fraction')
    plt.annotate('theory_gain = ' + str(sf(calc_theoretical_gain(gain), 3)) + ' e-/ADU', fontsize=10, xy=(0.5, 0.25), xycoords='axes fraction')
    plt.annotate('y_intercept = ' + str(sf(y_intercept, 2)) + ' ADU^2 = ' + str(sf(y_intercept_e, 2)) + 'e-^2', fontsize=10, xy=(0.5, 0.2), xycoords='axes fraction')
    plt.annotate('RN = ' + str(sf(rn, 2)) + 'e-', fontsize=10, xy=(0.5, 0.15), xycoords='axes fraction')
    plt.annotate('CDS RN = ' + str(sf(rn_cds, 2)) + 'e-', fontsize=10, xy=(0.5, 0.1), xycoords='axes fraction')

    plt.show()

