'''
calculate gain and read noise by making a variance v signal plot from either:
a) two UTR sequences (R01 and R02) in the same directory
b) a series of fowler exposures, each with two ramps.

b) requires a directory structure like

  PARENT/
         SUBDIR_1
         SUBDIR_2
         SUBDIR_3
         .
         .
         SUBDIR_N
   where each SUBDIR_? directory (different FS exptimes) contain two ramps.
   
n.b. for FS, note that you will need to set the .ini file to reflect the number of
pairs!
n.b.2 sometimes you'll need to change the number of bins around if the fit fails.
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
import matplotlib

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
    
class measure_gain_RN:
    def __init__(self, params):
        self.params = params        
        
    def go(self):
        np.seterr(all='raise')
        
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

        # file loggingReduction
        fh = logging.FileHandler(self.params['workingDir'] + "res.log")
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
        # run pipeline with either UTR sequence
        self.params['pipeCfgPath'] = "config/pipeline_gain_UTR.ini"
        if self.params['method'] == 'UTR':
            for i in range(self.params['minGrpNum'], self.params['maxGrpNum']-1):
                self.params['maxGrpNum'] = i+2
                
                pipe = run_pipe(self.params, logger, err)
                pipe.go()   
                
                f_r1_data.append(pipe.session.file_data_nonss[0][0])
                f_r1_hdr.append(pipe.session.file_hdr_nonss[0][0])
                f_r1_rates.append(pipe.session.rates[0][0])
                f_r2_data.append(pipe.session.file_data_nonss[1][0])
                f_r2_hdr.append(pipe.session.file_hdr_nonss[1][0]) 
                f_r2_rates.append(pipe.session.rates[1][0])
                
        # run pipeline with FS sequence (see note at top for expected directory format)   
        self.params['pipeCfgPath'] = "config/pipeline_gain_FS.ini"
        if self.params['method'] == 'FOWLER':             
            self.params['minGrpNum'] = 1
            self.params['maxGrpNum'] = 2    
            root_params = self.params['dataPath']
            for d in os.listdir(self.params['dataPath']):
                self.params['dataPath'] = root_params + d + "/"
                
                pipe = run_pipe(self.params, logger, err)
                pipe.go()   
                
                f_r1_data.append(pipe.session.file_data_nonss[0][0])
                f_r1_hdr.append(pipe.session.file_hdr_nonss[0][0])
                f_r1_rates.append(pipe.session.rates[0][0])
                f_r2_data.append(pipe.session.file_data_nonss[1][0])
                f_r2_hdr.append(pipe.session.file_hdr_nonss[1][0]) 
                f_r2_rates.append(pipe.session.rates[1][0])

        # get ASIC gain and frmtime
        gain        = float(pipe.session.file_hdr_nonss[0][0].comments['ASICGAIN'].split('(')[1].split('dB')[0])
        frmtime     = float(pipe.session.file_hdr_nonss[0][0]['FRMTIME']) 
        
        if len(f_r1_data) != len(f_r2_data):
            print "Ramps must have equal number of frames."
            exit(0)       
      
        RN          = []
        CALC_GAIN   = []
        # calculate signal mean and variance of difference frame for different windows
        for s1 in range(100, 1900, self.params['windowSize']):
            for s2 in range(100, 1900, self.params['windowSize']):
                means_average       = []
                shot_and_read_noise = []
                ## section files
                s1_lo = int(s1)
                s1_hi = int(s1+self.params['windowSize']+1)
                s2_lo = int(s2)
                s2_hi = int(s2+self.params['windowSize']+1)
                for i in range(0, len(f_r1_data), 1):
                    f_sect_r1 = f_r1_data[i][s1_lo:s1_hi, s2_lo:s2_hi] + np.median(f_r1_rates[i][s1_lo:s1_hi, s2_lo:s2_hi])*frmtime   # use median otherwise we add another component of readout noise!
                    f_sect_r2 = f_r2_data[i][s1_lo:s1_hi, s2_lo:s2_hi] + np.median(f_r2_rates[i][s1_lo:s1_hi, s2_lo:s2_hi])*frmtime   # use median otherwise we add another component of readout noise!

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
                    this_shot_and_read_noise = this_std_diff/pow(2, 0.5)   # taking the difference increases noise by sqrt(2). THIS GIVES US BACK THE CDS NOISE.
                    shot_and_read_noise.append(this_shot_and_read_noise)
                    
                ## get variance 
                vars_diff = [pow(i, 2) for i in shot_and_read_noise]
                
                ## calculate coeffs of fit to integration times
                ADUrange = np.where(np.logical_and(np.array(means_average)>=float(self.params['fit'][0]), np.array(means_average)<=float(self.params['fit'][1])))
                idx_lo = ADUrange[0][0]
                idx_hi = ADUrange[0][-1]
                if idx_lo == idx_hi:
                    continue
                fitted_coeffs = np.polyfit(means_average[idx_lo:idx_hi], vars_diff[idx_lo:idx_hi], self.params['fitCoeff'])
                fitted_total_noise = np.polyval(fitted_coeffs, means_average) 

                try: 
                    calc_gain           = 1./fitted_coeffs[-2]
                    y_intercept         = np.polyval(fitted_coeffs, 0)
                    rn                  = pow(y_intercept, 0.5)
                    rn_e                = rn*calc_gain
                except FloatingPointError:
                    continue
                RN.append(rn_e)
                CALC_GAIN.append(calc_gain)  

        font = {'family' : 'normal',
                'weight' : 'normal',
                'size'   : 12}
        matplotlib.rc('font', **font)              
        # plot RN        
        ## set up some RN histogram attributes
        if self.params['calcRead']:
            rn_median = np.median(RN)
            lim = (rn_median-40, rn_median+40)
            n_bins = 25
            y_pad = 100
            bin_width = (lim[1]-lim[0])/n_bins  
              
            n, bins, patches = plt.hist(RN, bins=np.arange(lim[0],lim[1],bin_width))                                                                     # histogram
            plt.setp(patches, 'facecolor', 'w', 'alpha', 1)                                   
            popt, pcov = curve_fit(gauss,[bins[i]+((bins[i+1]-bins[i])/2) for i in range(len(bins)-1)],n,p0=[1,rn_median,np.nanstd(RN)])                 # gauss fit
                
            plt.title("Read Noise (" + str(gain) + "dB)")
            plt.xlabel("CDS RN (e-)")
            plt.ylabel("Number")
                
            #plt.plot(np.arange(lim[0],lim[1],bin_width),gauss(np.arange(lim[0],lim[1],bin_width),*popt),'r-', linewidth=2, label='fit')
            plt.plot([popt[1] for i in range(0, int(np.ceil(max(n)))+y_pad)], range(0, int(np.ceil(max(n)))+y_pad), 'k--', linewidth=1, label='CDS RN = ' + str(sf(popt[1], 3)) + 'e-')
            plt.xlim(lim)
            plt.ylim([0, np.ceil(max(n) + y_pad)])
            plt.legend(loc='upper left', fontsize=10)
            if self.params['fits']:
                plt.savefig("rn.eps")   

            logger.info("Calculated CDS RN is " + str(sf(popt[1], 3)) + "e-.")            
              
        # plot GAIN    
        ## set up some GAIN histogram attributes
        if self.params['calcGain']:    
            gain_median = np.median(CALC_GAIN)
            lim = [gain_median-1, gain_median+1]
            n_bins = 30
            y_pad = 100
            bin_width = (lim[1]-lim[0])/n_bins   

            n, bins, patches = plt.hist(CALC_GAIN, bins=np.arange(lim[0], lim[1], bin_width))                                                                    # histogram
            plt.setp(patches, 'facecolor', 'w', 'alpha', 1)
            popt, pcov = curve_fit(gauss,[bins[i]+((bins[i+1]-bins[i])/2) for i in range(len(bins)-1)],n,p0=[1,gain_median,np.nanstd(CALC_GAIN)])                # gauss fit
                
            plt.title("Gain (" + str(gain) + "dB)")
            plt.xlabel("Gain (e-/ADU)")
            plt.ylabel("Number")
                
            #plt.plot(np.arange(lim[0], lim[1], bin_width),gauss(np.arange(lim[0], lim[1], bin_width),*popt),'r-', linewidth=2, label='fit')
            plt.plot([popt[1] for i in range(0, int(np.ceil(max(n)))+y_pad)], range(0, int(np.ceil(max(n)))+y_pad), 'k--', linewidth=1, label='Calc gain = ' + str(sf(popt[1], 3)) + 'e-/ADU')
            plt.xlim(lim)
            plt.ylim([0, np.ceil(max(n) + y_pad)])
            plt.legend(loc='upper left', fontsize=10)
            
            if self.params['hard']:
                plt.savefig("gain.eps")   
            
            logger.info("Calculated gain is " + str(sf(popt[1], 3)) + "e-/ADU.")    
            
if __name__ == "__main__":
    parser = optparse.OptionParser()
    group1 = optparse.OptionGroup(parser, "General")
    group1.add_option('--p', action='store', default='/mnt/NAS/devel/IOI/images_and_analysis/remote_5/images/12/', dest='dataPath', type=str, help='path to data')
    group1.add_option('--wd', action='store', default='test', dest='workingDir', type=str, help='path to working directory')
    group1.add_option('--o', action='store_true', dest='clobber', help='clobber working directory?')
    group1.add_option('--pa', action='store', default='../../config/paths_rmb.ini', type=str, dest='pathsCfgPath', help='path to paths config file')    
    group1.add_option('--log', action='store', default='DEBUG', dest='logLevel', type=str, help='log level (DEBUG|INFO|WARNING|ERROR|CRITICAL)')      
    group1.add_option('--glo', action='store', default=0, type=int, dest='minGrpNum', help='lowest group number to use (method==UTR only)')
    group1.add_option('--ghi', action='store', default=20, type=int, dest='maxGrpNum', help='highest group number to use (method==UTR only)')
    group1.add_option('--r', action='store', dest='fit', default="0,10000", help='ADU range to consider when fitting')
    group1.add_option('--c', action='store', default=1, type=float, dest='fitCoeff', help='Coeffient of fit to signal')
    group1.add_option('--win', action='store', dest='windowSize', default=50, type=int, help='size of window to use (px)')  
    group1.add_option('--hard', action='store_true', dest='hard', help='make hard plot?')    
    group1.add_option('--m', action='store', dest='method', default='UTR', type=str, help='method (UTR||FOWLER)')
    group1.add_option('--gain', action='store_true', dest='calc_gain', help='calculate gain?')
    group1.add_option('--read', action='store_true', dest='calc_read', help='calculate read?')
    parser.add_option_group(group1)
    
    options, args = parser.parse_args()
    params = {
        'dataPath' : str(options.dataPath).rstrip("/") + "/",
        'workingDir' : str(options.workingDir).rstrip("/") + "/",
        'clobber' : bool(options.clobber),       
        'pathsCfgPath' : str(options.pathsCfgPath),
        'minRunNum' : 1,
        'maxRunNum' : 2,
        'minDithNum' : 0,
        'maxDithNum' : 0,
        'minGrpNum' : int(options.minGrpNum),
        'maxGrpNum' : int(options.maxGrpNum),
        'minExpNum' : 1,
        'maxExpNum' : 16,        
        'logLevel' : str(options.logLevel.upper()),
        'fit' : str(options.fit).split(','),
        'fitCoeff' : float(options.fitCoeff),
        'windowSize' : int(options.windowSize),
        'hard' : bool(options.hard),
        'method' : str(options.method),
        'calcGain': bool(options.calc_gain),
        'calcRead': bool(options.calc_read),
        'flip' : False
    }
    
    ## execute
    m = measure_gain_RN(params)
    m.go()
