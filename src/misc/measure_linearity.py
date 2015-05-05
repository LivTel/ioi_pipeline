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
import pylab as plt
from scipy.optimize import curve_fit

sys.path.append("../")
from utility import read_FITS_file
from run_pipe import run_pipe
from errors import errors

class measure_linearity:
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

        # file loggingReduction
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
        
        # UNCORRECTED
        self.params['pipeCfgPath'] = "config/pipeline_lin_1.ini"
        data_mean_nolin     = []
        exptime_nolin       = []
        for i in range(self.params['minGrpNum'], self.params['maxGrpNum']-1):
            self.params['maxGrpNum'] = i+2
                
            pipe = run_pipe(self.params, logger, err)
            pipe.go()   
                
            data = pipe.session.file_data_nonss[0][0]
            rates = pipe.session.rates[0][0]
            exptime_nolin.append(pipe.session.file_hdr_nonss[0][0]['EXPTIME'])
            frmtime = pipe.session.file_hdr_nonss[0][0]['FRMTIME']
            missing_flux = rates*frmtime
            data_mean_nolin.append(np.nanmean(data[y1:y2,x1:x2]+missing_flux[y1:y2,x1:x2]))

        if self.params['doLinearity']:
            # CORRECTED
            self.params['pipeCfgPath'] = "config/pipeline_lin_2.ini"        
            data_mean_lin       = []
            exptime_lin         = []
            for i in range(self.params['minGrpNum'], self.params['maxGrpNum']-1):
                self.params['maxGrpNum'] = i+2
                  
                pipe = run_pipe(self.params, logger, err)
                pipe.go()   
                    
                data = pipe.session.file_data_nonss[0][0]
                rates = pipe.session.rates[0][0]
                exptime_lin.append(pipe.session.file_hdr_nonss[0][0]['EXPTIME'])
                frmtime = pipe.session.file_hdr_nonss[0][0]['FRMTIME']
                missing_flux = rates*frmtime
                data_mean_lin.append(np.nanmean(data[y1:y2,x1:x2]+missing_flux[y1:y2,x1:x2]))
                
        # linearity line
        fitted_coeffs = np.polyfit([exptime_nolin[0], exptime_nolin[1]], [data_mean_nolin[0], data_mean_nolin[1]], 1)         
        
        if self.params['p1']:
            #plt.clf()
            plt.plot([0, exptime_nolin[-1]], np.polyval(fitted_coeffs, [0, exptime_nolin[-1]]), 'k:', label='Linear fit')       
            plt.plot(exptime_nolin, data_mean_nolin, 'kx-', label="No linearity correction")
            if self.params['doLinearity']:
                plt.plot(exptime_lin, data_mean_lin, 'ko--', label="With linearity correction")
            
            plt.ylim([0,66000])
            plt.title("Flux v EXPTIME")
            plt.xlabel("EXPTIME (s)")
            plt.ylabel("CDS mean counts (ADU)")
            plt.legend(loc='lower right', fontsize=10)
            plt.ticklabel_format(style='sci', axis='y', scilimits=(0,66000))
              
            plt.ticklabel_format(axis='y',style='sci',scilimits=(1,3))           
            
            if self.params['hard']: 
                plt.savefig("linearity1.eps")
                
        if self.params['p2']:    
            #plt.clf()
            plt.plot(data_mean_nolin, (abs(np.polyval(fitted_coeffs, exptime_nolin)-data_mean_nolin)/np.polyval(fitted_coeffs, exptime_nolin))*100, 'kx-', label='Nonlinearity % (uncorrected)') 
            if self.params['doLinearity']:
                plt.plot(data_mean_lin, (abs(np.polyval(fitted_coeffs, exptime_lin)-data_mean_lin)/np.polyval(fitted_coeffs, exptime_lin))*100, 'ko--', label='Residual nonlinearity % (corrected)')   
            
            # 5% RESIDUAL nonlinearity line
            plt.plot([0, 70000], [5, 5] , 'k--', label="5% nonlinearity")   
            
            plt.title("Residual Nonlinearity")
            plt.xlabel("CDS mean counts (ADU)")
            plt.ylabel("Nonlinearity %")
            plt.legend(loc='lower right', fontsize=10)
            plt.yscale('log')
            plt.ylim([0.1,10])
            
            if self.params['hard']: 
                plt.savefig("linearity2.eps")        
        
if __name__ == "__main__":
    np.seterr(all='raise')
    
    parser = optparse.OptionParser()
    group1 = optparse.OptionGroup(parser, "General")
    group1.add_option('--p', action='store', default='/mnt/NAS/devel/IOI/images_and_analysis/remote_5/images/14/', dest='dataPath', type=str, help='path to data')
    group1.add_option('--wd', action='store', default='test', dest='workingDir', type=str, help='path to working directory')
    group1.add_option('--o', action='store_true', dest='clobber', help='clobber working directory?')
    group1.add_option('--pa', action='store', default='../../config/paths_rmb.ini', type=str, dest='pathsCfgPath', help='path to paths config file')    
    group1.add_option('--log', action='store', default='DEBUG', dest='logLevel', type=str, help='log level (DEBUG|INFO|WARNING|ERROR|CRITICAL)')      
    group1.add_option('--glo', action='store', default=0, type=int, dest='minGrpNum', help='lowest group number to use')
    group1.add_option('--ghi', action='store', default=32, type=int, dest='maxGrpNum', help='highest group number to use ')
    group1.add_option('--cor', action='store_true', dest='doLinearity', help='compute corrections?')
    group1.add_option('--s', action='store', dest='section', default='287,1600,703,1803', help='csv section of image to use [x1,y1,x2,y2]')    
    group1.add_option('--p1', action='store_true', dest='p1', help='make plot of exptime v counts?')    
    group1.add_option('--p2', action='store_true', dest='p2', help='make plot of counts v linearity?') 
    group1.add_option('--hard', action='store_true', dest='hard', help='make hard plot?')
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
        'doLinearity' : bool(options.doLinearity),
        'section' : str(options.section),
        'p1' : bool(options.p1),
        'p2' : bool(options.p2),
        'hard' : bool(options.hard),
        'flip' : False
    }
    
    m = measure_linearity(params)
    m.go()