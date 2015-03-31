'''
    measure non-linearity using UTR sequences.
    
    when setting the yFit parameter, you should be aware that
    this is defined after the frame has been CDS'ed, so the 
    upper limit is restricted by bias level of the first frame.
    e.g. consider:
    
    MAX BIAS ~ 25k, fills to 65K where the ADC saturates.
    YFIT upper limit for these pixels will never pass 65k-25k 
    = 40k.
    
    it should also be noted that close to FWD, the polyfit will not 
    adequately match the real functional shape and the nonlinearity %
    derived will be large. as it happens, we don't really need
    to match all the way up to FWD, because anything past 
    intrinsic 5% nonlinearity is probably garbage anyway.
'''
import sys
import optparse
import os
import numpy as np
import matplotlib.pyplot as plt
import pyfits
import warnings
import logging
import math

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

if __name__ == "__main__":
    parser = optparse.OptionParser()
    group1 = optparse.OptionGroup(parser, "General")
    group1.add_option('--p', action='store', default='/mnt/NAS/devel/IOI/images_and_analysis/remote_5/images/19/', dest='dataPath', type=str, help='path to data')
    group1.add_option('--wd', action='store', default='test', dest='workingDir', type=str, help='path to working directory')
    group1.add_option('--o', action='store_true', dest='clobber', help='clobber working directory?')
    group1.add_option('--pa', action='store', default='../../config/paths_rmb.ini', type=str, dest='pathsCfgPath', help='path to paths config file')
    group1.add_option('--log', action='store', default='DEBUG', dest='logLevel', type=str, help='log level (DEBUG|INFO|WARNING|ERROR|CRITICAL)')   
    group1.add_option('--rlo', action='store', default=1, type=int, dest='minRunNum', help='lowest ramp number to use')
    group1.add_option('--rhi', action='store', default=1, type=int, dest='maxRunNum', help='highest ramp number to use')    
    group1.add_option('--glo', action='store', default=0, type=int, dest='minGrpNum', help='lowest group number to use')
    group1.add_option('--ghi', action='store', default=25, type=int, dest='maxGrpNum', help='highest group number to use')
    group1.add_option('--l', action='store', default='20,100', type=str, dest='linFit', help='Time range over which to establish linear fit to slope')
    group1.add_option('--f', action='store', default='2500,40000', type=str, dest='yFit', help='ADU range over which to establish polyfit.')
    group1.add_option('--c', action='store', default=3, type=float, dest='fitCoeff', help='Coeffient of fit to signal')
    group1.add_option('--b', action='store', default=0.8, type=float, dest='badThresh', help='Maximum percentage of non-linearity post-correction to flag as bad pixel')
    group1.add_option('--sig', action='store', default=3, type=float, dest='sig', help='Number of sigma which mean rate must lie above 0 to be classed as good')    
    group1.add_option('--plb', action='store_true', dest='plotB', help='plot bad pixels')
    group1.add_option('--pla', action='store_true', dest='plotA', help='plot all pixels')
    group1.add_option('--fits', action='store_true', dest='fits', help='make flat and bad pixel map FITS files')
    group1.add_option('--fl', action='store_true', dest='flip', help='flip array on write?')
    parser.add_option_group(group1)
    
    options, args = parser.parse_args()
    params = {
        'dataPath' : str(options.dataPath).rstrip("/") + "/",
        'workingDir' : str(options.workingDir).rstrip("/") + "/",
        'clobber' : bool(options.clobber),       
        'pathsCfgPath' : str(options.pathsCfgPath),
        'minRunNum' : int(options.minRunNum),
        'maxRunNum' : int(options.maxRunNum),
        'minDithNum' : 0,
        'maxDithNum' : 0,
        'minGrpNum' : int(options.minGrpNum),
        'maxGrpNum' : int(options.maxGrpNum),
        'minExpNum' : 1,
        'maxExpNum' : 1,        
        'logLevel' : str(options.logLevel.upper()),
        'linFit' : str(options.linFit),
        'yFit' : str(options.yFit),
        'fitCoeff' : int(options.fitCoeff),
        'badThresh' : float(options.badThresh),
        'sig' : float(options.sig),        
        'plotA' : bool(options.plotA),
        'plotB' : bool(options.plotB),
        'fits' : bool(options.fits),
        'flip' : bool(options.flip)
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

    # file logging
    fh = logging.FileHandler(params['workingDir'] + "res.log")
    fh.setLevel(logging.DEBUG)

    ## set logging format
    fh.setFormatter(formatter)

    ## add handlers to logging object
    logger.addHandler(fh)
    
    # run pipeline
    data        = []
    inttime     = []
    params['pipeCfgPath'] = "config/pipeline_make_lin_cube.ini"
    for i in range(params['minGrpNum'], params['maxGrpNum']-1):
        params['maxGrpNum'] = i+2
        data_this_grp = []
        for j in range(params['minRunNum'], params['maxRunNum']+1, 1):
            pipe = run_pipe(params, logger, err)
            pipe.go()   
            frmtime = pipe.session.file_hdr_nonss[j-1][0]['FRMTIME']
            data_this_grp.append(pipe.session.file_data_nonss[j-1][0] + pipe.session.rates[j-1][0]*frmtime) 
            if j == 1:
                inttime.append(pipe.session.file_hdr_nonss[j-1][0]['EXPTIME'])
        data.append(np.median(data_this_grp, axis=0))

    # create coeffs, nonlinearities and badpixmap arrays
    coeffs = []
    for i in range(params['fitCoeff']+1):
        coeffs.append(np.zeros(data[0].shape))
    nonlinearities = np.zeros(data[0].shape)    # create nonlinearities array
    badpixmap = np.ones(data[0].shape)          # create bad pixel map array 

    bad_pix_count = 0
    for jj in range(0, 1):
        for ii in range(0, 5):
            this_pixel_values = []
            flagged_pixel = False
            ## obtain all values for this pixel
            for t in data:
                this_pixel_values.append(t[jj, ii])
                
            ## get linear fit to slope within [linFit] limits
            this_pixel_vals_in_fit_lin	 	= [i for i,j in zip(this_pixel_values, inttime) if j > float(params['linFit'].split(',')[0]) and j < float(params['linFit'].split(',')[1])]
            this_pixel_inttimes_in_fit_lin	= [j for i,j in zip(this_pixel_values, inttime) if j > float(params['linFit'].split(',')[0]) and j < float(params['linFit'].split(',')[1])]    
            linfit_coeffs 			= np.polyfit(this_pixel_inttimes_in_fit_lin, this_pixel_vals_in_fit_lin, 1)

            ## fit slope to uncorrected data within [yfit] limit
            this_pixel_vals_in_fit 	= [i for i,j in zip(this_pixel_values, inttime) if i > float(params['yFit'].split(',')[0]) and i < float(params['yFit'].split(',')[1])]
            this_pixel_inttimes_in_fit	= [j for i,j in zip(this_pixel_values, inttime) if i > float(params['yFit'].split(',')[0]) and i < float(params['yFit'].split(',')[1])]

            if len(this_pixel_vals_in_fit) == 0:        # TEST: if count rate is so low that array is empty..
                flagged_pixel = True
                
                ideal_signal_vals_in_fit 	= []
                this_pixel_correction_coeffs 	= []
                this_pixel_values_cor 		= []
                this_pixel_nonlinearity_uncor 	= []
                this_pixel_nonlinearity_cor	= []
                
                ## store coeffs as nan
                for idx in range(len(coeffs)):  
                    coeffs[idx][jj, ii] = 1
                ## store nonlinearity as nan
                nonlinearities[jj, ii] = 0
                ## store bad pixel as nan
                badpixmap[jj, ii] = np.nan
                
                print "FAILED FLUX: Bad pixel found @ " + str(ii) + ", " + str(jj) + " with no increase of flux with time."
                bad_pix_count = bad_pix_count + 1
                print "NUMBER OF BAD PIXELS: " + str(bad_pix_count)
            else:
                ## find ideal signal for all inttimes available
                ideal_signal_vals_in_fit = np.polyval(linfit_coeffs, this_pixel_inttimes_in_fit)

                ## find correction (mapping of observed to "ideal")
                this_pixel_correction_coeffs = np.polyfit(this_pixel_vals_in_fit, ideal_signal_vals_in_fit, params['fitCoeff'])
          
                ## store coeffs
                for idx in range(len(coeffs)):  
                    coeffs[idx][jj, ii] = this_pixel_correction_coeffs[idx]
 
                ## linearise and store this pixel's data
                this_pixel_values_cor = []
                for t in data:
                    this_corrected_value = np.polyval(this_pixel_correction_coeffs, t[jj,ii])
                    this_pixel_values_cor.append(this_corrected_value)

                ## find nonlinearity as percentage pre/post and store

                this_pixel_nonlinearity_uncor = abs((this_pixel_values - np.polyval(linfit_coeffs, inttime))/np.polyval(linfit_coeffs, inttime))*100
                this_pixel_nonlinearity_cor = abs((this_pixel_values_cor - np.polyval(linfit_coeffs, inttime))/np.polyval(linfit_coeffs, inttime))*100
                avg_nonlinearity = np.median([this_pixel_nonlinearity_cor[i] for i in range(len(this_pixel_vals_in_fit))])
                nonlinearities[jj, ii] = avg_nonlinearity
                
                ### TEST: if the residual nonlinearity is below the "bad" threshold:
                if avg_nonlinearity > params['badThresh']:
                    flagged_pixel = True
                    print "FAILED LINEARITY: Bad pixel found @ " + str(ii) + ", " + str(jj) + " with a nonlinearity of " + str(sf(np.median(avg_nonlinearity), 2)) + "%"
                    badpixmap[jj, ii] = np.nan
                    bad_pix_count = bad_pix_count + 1
                    print "NUMBER OF BAD PIXELS: " + str(bad_pix_count)
                    
                ### TEST: if the rate is effectively zero:
                rates = [(this_pixel_values[i]-this_pixel_values[i-1])/(inttime[i]-inttime[i-1]) for i in range(len(this_pixel_vals_in_fit))]
                if np.mean(rates) - (params['sig'] * np.std(rates)) < 0:
                    flagged_pixel = True
                    print "FAILED RATE: Bad pixel found @ " + str(ii) + ", " + str(jj) + " with a rate of " + str(sf(np.mean(rates), 2)) + "ADU/s"
                    badpixmap[jj, ii] = np.nan
                    bad_pix_count = bad_pix_count + 1
                    print "NUMBER OF BAD PIXELS: " + str(bad_pix_count)
                    
       	    ## plot graphs if requested
            if params['plotA'] or (params['plotB'] and flagged_pixel):
                fig = plt.figure()

                fig.add_subplot(321)
                plt.plot(inttime, this_pixel_values, 'rx-', label="pixel " + str(ii) + ", " + str(jj) + " uncorrected")
                try:
                    plt.plot(inttime, this_pixel_values_cor, 'gx-', label="pixel " + str(ii) + ", " + str(jj) + " corrected")
                except ValueError:
                    pass
                plt.plot(inttime, np.polyval(linfit_coeffs, inttime), 'b-', label="Linear fit")

                plt.plot(inttime, [65536 for i in inttime], 'k-', label="16bit ADC limit")
                plt.annotate('ADC Limit = 65536 ADU', xy=(200, 65536+1500), xytext=(200, 65536+1500))
                plt.plot(inttime, [float(params['yFit'].split(',')[1]) for i in inttime], 'k--', label="Fit limit")
                plt.plot(inttime, [linfit_coeffs[1] for i in inttime], 'k-', label="Bias")
                plt.annotate('Bias = ' + str(round(linfit_coeffs[1])) + ' ADU', xy=(200, min(this_pixel_values)+1500), xytext=(200, min(this_pixel_values)+1500))

                plt.ylabel("Signal (ADU)")
                plt.xlabel("Time (s)")
                plt.ylim([0,100000])
                leg = plt.legend(loc='lower right', prop={'size':11})
                leg.get_frame().set_alpha(0.75)

                fig.add_subplot(322)
                try:
                    plt.plot(this_pixel_vals_in_fit, ideal_signal_vals_in_fit, 'rx-')
                except ValueError:
                    pass
                try:
                    plt.plot(this_pixel_vals_in_fit, np.polyval(this_pixel_correction_coeffs, this_pixel_vals_in_fit), 'b-', label="Polyfit(" + str(params['fitCoeff'])+")")
                except ValueError:
                    pass
                plt.ylabel("Detected signal (ADU)")
                plt.xlabel("\"Linear   \" signal (ADU)")
                leg = plt.legend(loc='lower right', prop={'size':11})
                leg.get_frame().set_alpha(0.75)

                fig.add_subplot(323)
                try:
                    plt.plot(this_pixel_values, this_pixel_nonlinearity_uncor, 'rx-', label="uncorrected")
                except ValueError:
                    pass
                try:
                    plt.plot(this_pixel_values_cor, this_pixel_nonlinearity_cor, 'gx-', label="corrected")
                except ValueError:
                    pass
                plt.ylabel("% nonlinearity")
                plt.xlabel("Signal (ADU)")
                plt.xlim([float(params['yFit'].split(',')[0]), float(params['yFit'].split(',')[1])])
                plt.ylim([0, 5])
                leg = plt.legend(loc='upper right', prop={'size':11})
                leg.get_frame().set_alpha(0.75)

                fig.add_subplot(324)
                plt.plot([this_pixel_values[i] for i in range(1, len(this_pixel_values)-2)], [(this_pixel_values[i]-this_pixel_values[i-1])/(inttime[i]-inttime[i-1]) for i in range(1, len(this_pixel_values)-2)], 'rx-', label="uncorrected")
                try:
                    plt.plot([this_pixel_values_cor[i] for i in range(1, len(this_pixel_values_cor)-2)], [(this_pixel_values_cor[i]-this_pixel_values_cor[i-1])/(inttime[i]-inttime[i-1]) for i in range(1, len(this_pixel_values_cor)-2)], 'gx-', label="corrected")
                except ValueError:
                    pass
                plt.ylabel("Rate (ADU/s)")
                plt.xlabel("Signal (ADU)")
                leg = plt.legend(loc='lower left', prop={'size':11})
                leg.get_frame().set_alpha(0.75)

                ## plot image postage stamp (end integration)
                ax = fig.add_subplot(325)
                win_y_s = 0 if jj-5 < 0 else jj-5
                win_y_e = data[0].shape[1] if jj+5 > data[0].shape[1] else jj+5
                win_x_s = 0 if ii-5 < 0 else ii-5
                win_x_e = data[0].shape[0] if ii+5 > data[0].shape[1] else ii+5
                cax = ax.imshow(t[win_y_s:win_y_e, win_x_s:win_x_e], interpolation="nearest", aspect='auto')
                ax.scatter(ii-win_x_s,jj-win_y_s, marker='x', color='white')
                cb_ax = fig.colorbar(cax)

                plt.show()
        print "row " + str(jj+1) + " of " + str((data[0].shape[1]))

    if params['fits']:
        if params['flip']:
            nonlinearities = np.fliplr(nonlinearities)
            badpixmap = np.fliplr(badpixmap)
        pyfits.writeto("lin.fits", data=nonlinearities) 
        pyfits.writeto("lin_bad.fits", data=badpixmap)         
        for i in range(params['fitCoeff']+1):  
            this_coeffs_32 = coeffs[params['fitCoeff']-i].astype(np.float32, copy=False)    # params['fitCoeff']-i to reverse so that coeff_0 is the constant..  
            if params['flip']:
                this_coeffs_32 = np.fliplr(this_coeffs_32)
                
            if os.path.exists("lin_coeffs.fits"):
                pyfits.append("lin_coeffs.fits", this_coeffs_32)    
            else:
                pyfits.writeto("lin_coeffs.fits", this_coeffs_32)               
