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
      
def gauss(x,a,x0,sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))

def calc_theoretical_gain(gain):
    A = pow(10, (gain/20))
    g = 3 * (40*pow(10, -15))/(pow(2, 16) * A * 1.6 * pow(10, -19))
    return g

if __name__ == "__main__":
  
    np.seterr(all='raise')
    
    parser = optparse.OptionParser()
    group1 = optparse.OptionGroup(parser, "General")
    group1.add_option('--p', action='store', default='/mnt/NAS/devel/IOI/images_and_analysis/remote_5/images/9/', dest='dataPath', type=str, help='path to data')
    group1.add_option('--wd', action='store', default='test', dest='workingDir', type=str, help='path to working directory')
    group1.add_option('--o', action='store_true', dest='clobber', help='clobber working directory?')
    group1.add_option('--pa', action='store', default='../../config/paths_rmb.ini', type=str, dest='pathsCfgPath', help='path to paths config file')
    group1.add_option('--pi', action='store', default='../../config/pipeline_gain_rmb.ini', type=str, dest='pipeCfgPath', help='path to pipeline config file')     
    group1.add_option('--log', action='store', default='DEBUG', dest='logLevel', type=str, help='log level (DEBUG|INFO|WARNING|ERROR|CRITICAL)')      
    group1.add_option('--glo', action='store', default=0, type=int, dest='minGrpNum', help='lowest group number to use')
    group1.add_option('--ghi', action='store', default=20, type=int, dest='maxGrpNum', help='highest group number to use')
    group1.add_option('--l', action='store', dest='maxLinADU', default=10000, type=float, help='maximum ADU to consider for linear fit')
    group1.add_option('--win', action='store', dest='windowSize', default=50, type=int, help='size of window to use (px)')  
    group1.add_option('--pl', action='store_true', dest='plt', help='make plot?') 
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
        'maxLinADU' : float(options.maxLinADU),
        'windowSize' : int(options.windowSize),
        'plt' : bool(options.plt)
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
        
        f_r1_data.append(pipe.session.file_data_nonss[0][0])
        f_r1_hdr.append(pipe.session.file_hdr_nonss[0][0])
        f_r1_rates.append(pipe.session.rates[0][0])
        f_r2_data.append(pipe.session.file_data_nonss[1][0])
        f_r2_hdr.append(pipe.session.file_hdr_nonss[1][0]) 
        f_r2_rates.append(pipe.session.rates[1][0])

    ## get ASIC gain and frmtime
    gain        = float(pipe.session.file_hdr_nonss[0][0].comments['ASICGAIN'].split('(')[1].split('dB')[0])
    frmtime     = float(pipe.session.file_hdr_nonss[0][0]['FRMTIME']) 
    
    if len(f_r1_data) != len(f_r2_data):
        print "Ramps must have equal number of frames."
        exit(0)       
   
    RN          = []
    RN_CDS      = []
    CALC_GAIN   = []
    
    # calculate signal mean and variance of difference frame for 50x100 pixel different windows
    for s1 in range(100, 1900, params['windowSize']):
        for s2 in range(100, 1900, params['windowSize']):
            means_average       = []
            shot_and_read_noise = []
            ## section files
            s1_lo = int(s1)
            s1_hi = int(s1+100+1)
            s2_lo = int(s2)
            s2_hi = int(s2+100+1)
            for i in range(0, len(f_r1_data), 1):
                f_sect_r1 = f_r1_data[i][s1_lo:s1_hi, s2_lo:s2_hi] + np.median(f_r1_rates[i][s1_lo:s1_hi, s2_lo:s2_hi])   # use median otherwise we add another component of readout noise!
                f_sect_r2 = f_r2_data[i][s1_lo:s1_hi, s2_lo:s2_hi] + np.median(f_r2_rates[i][s1_lo:s1_hi, s2_lo:s2_hi])   # use median otherwise we add another component of readout noise!

                ## get means/stds of sections for each file
                this_mean_r1 = np.mean(f_sect_r1)
                this_mean_r2 = np.mean(f_sect_r2)
                this_std_r1  = np.std(f_sect_r1)
                this_std_r2  = np.std(f_sect_r2)

                ## calculate average signal and average (total) noise
                this_mean_average = (this_mean_r1 + this_mean_r2)/2.
                means_average.append(this_mean_average)
        
                ## calculate std of diff frame to remove FPN component
                diff = (f_sect_r1 + 1000) - f_sect_r2   # add 1000 ADU to avoid negative offsets
                this_std_diff = np.std(diff)

                ## calculate shot_and_read_noise component for single frame
                this_shot_and_read_noise = this_std_diff/pow(2, 0.5)    # taking the difference (CDS) increases noise by sqrt(2)
                shot_and_read_noise.append(this_shot_and_read_noise)
                
            # get variance 
            vars_diff = [pow(i, 2) for i in shot_and_read_noise]
            
            ## calculate coeffs of linear fit to early integration times (where accumulation of signal ~ linear)
            idx_lo = np.where(np.array(means_average)<params['maxLinADU'])[0][0]
            idx_hi = np.where(np.array(means_average)<params['maxLinADU'])[0][-1]
            if idx_lo == 0 and idx_hi == 0:
                continue 
            fitted_coeffs = np.polyfit(means_average[idx_lo:idx_hi], vars_diff[idx_lo:idx_hi], 1)
            fitted_total_noise = np.polyval(fitted_coeffs, means_average)
            
            try: 
                calc_gain           = 1./fitted_coeffs[-2]
                y_intercept         = np.polyval(fitted_coeffs, 0)
                y_intercept_e       = y_intercept * calc_gain
                rn                  = pow(y_intercept, 0.5)*calc_gain       # e-
                rn_cds              = rn * pow(2, 0.5)                      # e-
            except FloatingPointError:
                continue
            
            RN.append(rn)
            RN_CDS.append(rn_cds)
            CALC_GAIN.append(calc_gain) 
            
    # CDS RN    
    if gain == 21:
        fit = [15,40]
    elif gain == 18:
       fit = [15,45]
    elif gain == 15:
        fit = [20,45]
    elif gain == 12:
        fit = [20,50]

    RN_CDS_TRIM = []
    for i in np.where(np.logical_and(np.array(RN_CDS)>=fit[0], np.array(RN_CDS)<=fit[1]))[0].tolist():
        RN_CDS_TRIM.append(RN_CDS[i])        
    n, bins, patches = plt.hist(RN_CDS_TRIM, bins=np.arange(fit[0],fit[1],1.0), label='data')
    popt,pcov = curve_fit(gauss,[bins[i]+((bins[i+1]-bins[i])/2) for i in range(len(bins)-1)],n,p0=[1,np.nanmean(RN_CDS_TRIM),np.nanstd(RN_CDS_TRIM)])
    
    logger.info("Calculated CDS RN is " + str(sf(popt[1], 3)) + "e-.")
        
    if params['plt']:
        plt.figure()
        
        plt.title("Histogram CDS RN Plot (" + str(gain) + "dB)")
        plt.xlabel("CDS RN (e-)")
        plt.ylabel("Number")
        
        plt.plot(np.arange(fit[0],fit[1],1.0),gauss(np.arange(fit[0],fit[1],1.0),*popt),'r-', linewidth=2, label='fit')
        plt.plot([popt[1] for i in range(0, int(np.ceil(max(n))))], range(0, int(np.ceil(max(n)))), 'r--', linewidth=2, label='CDS RN = ' + str(sf(popt[1], 3)) + 'e-')
        plt.xlim(fit)
        plt.legend(loc='upper left', fontsize=10)
        plt.savefig("cdsrn.png")    
     
    # GAIN   
    if gain == 21:
        fit = [0.9,1.2]
    elif gain == 18:
        fit = [1.25,1.7] 
    elif gain == 15:
        fit = [1.85,2.25]
    elif gain == 12:
        fit = [2.5,3.2]   
        
    CALC_GAIN_TRIM = []
    for i in np.where(np.logical_and(np.array(CALC_GAIN)>=fit[0], np.array(CALC_GAIN)<=fit[1]))[0].tolist():
        CALC_GAIN_TRIM.append(CALC_GAIN[i])      
    n, bins, patches = plt.hist(CALC_GAIN_TRIM, bins=np.arange(fit[0], fit[1], 0.02), label='data')
    popt,pcov = curve_fit(gauss,[bins[i]+((bins[i+1]-bins[i])/2) for i in range(len(bins)-1)],n,p0=[1,np.nanmean(CALC_GAIN_TRIM),np.nanstd(CALC_GAIN_TRIM)])
    
    logger.info("Calculated gain is " + str(sf(popt[1], 3)) + "e-/ADU.")    
        
    if params['plt']:        
        plt.figure() 
        
        plt.title("Histogram Gain Plot (" + str(gain) + "dB)")
        plt.xlabel("Gain (e-)")
        plt.ylabel("Number")
        
        plt.plot(np.arange(fit[0], fit[1], 0.02),gauss(np.arange(fit[0], fit[1], 0.02),*popt),'r-', linewidth=2, label='fit')
        plt.plot([popt[1] for i in range(0, int(np.ceil(max(n))))], range(0, int(np.ceil(max(n)))), 'r--', linewidth=2, label='Calc gain = ' + str(sf(popt[1], 3)) + 'e-')
        plt.xlim(fit)
        plt.legend(loc='upper left', fontsize=10)
        plt.savefig("gain.png")
