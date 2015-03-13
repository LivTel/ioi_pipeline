'''
measure non-linearity using a single ramp of data.

naming assumes Teledyne's convention (H2RG_RXX_MYY_NZZ).
'''
import sys
import argparse
import os
import numpy as np
import matplotlib.pyplot as plt
import pyfits
import warnings
import logging
from utility import read_FITS_file

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', action='store', default='/mnt/NAS/devel/IOI/images_and_analysis/remote_4/analysis/12/noref/', dest='dataPath', help='path to data directory')
    parser.add_argument('-rn', action='store', default=1, type=int, dest='rampNum', help='Ramp number to use')
    parser.add_argument('-elo', action='store', default=0, type=int, dest='expNumLo', help='Lowest exposure number to use')
    parser.add_argument('-ehi', action='store', default=20, type=int, dest='expNumHi', help='Highest exposure number to use')
    parser.add_argument('-glo', action='store', default=0, type=int, dest='grpNumLo', help='lowest group number to use')
    parser.add_argument('-ghi', action='store', default=250, type=int, dest='grpNumHi', help='highest group number to use')
    parser.add_argument('-s', action='store', default='0,0,2039,2039', type=str, dest='sect', help='section (x1,y1,x2,y2)')
    parser.add_argument('-linf', action='store', default='20,50', type=str, dest='linFit', help='Time range over which to establish linear fit to slope')
    parser.add_argument('-yf', action='store', default='2500,50000', type=str, dest='yFit', help='ADU range over which to establish polyfit. be aware that any pixel which has a well depth short of the upper value will most likely be flagged as bad')
    parser.add_argument('-fc', action='store', default=3, type=float, dest='fitCoeff', help='Coeffient of fit to signal')
    parser.add_argument('-bt', action='store', default=0.8, type=float, dest='badThresh', help='Maximum percentage of non-linearity post-correction to flag as bad pixel')
    parser.add_argument('-sig', action='store', default=3, type=float, dest='sig', help='Number of sigma which mean rate must lie above 0 to be classed as good')    
    parser.add_argument('-plb', action='store_true', dest='plotB', help='plot bad pixels')
    parser.add_argument('-pla', action='store_true', dest='plotA', help='plot all pixels')
    parser.add_argument('-fits', action='store_true', dest='fits', help='make coeff, linearity and bad pixel map FITS files')
    parser.add_argument('-log', action='store', default='INFO', dest='logLevel', type=str, help='level (DEBUG|INFO|WARNING|ERROR|CRITICAL)')  
    parser.add_argument('-suf', action='store', default='', dest='suffix', type=str, help='file suffix')

    args = parser.parse_args()
    params = {
        'dataPath' : str(args.dataPath),
	'rampNum' : int(args.rampNum),
	'expNumLo' : int(args.expNumLo),
	'expNumHi' : int(args.expNumHi),
	'grpNumLo' : int(args.grpNumLo),
	'grpNumHi' : int(args.grpNumHi),
        'sect' : str(args.sect),
        'linFit' : str(args.linFit),
        'yFit' : str(args.yFit),
	'fitCoeff' : int(args.fitCoeff),
        'sig' : float(args.sig),	
	'badThresh' : float(args.badThresh),
        'plotA' : bool(args.plotA),
        'plotB' : bool(args.plotB),
        'fits' : bool(args.fits),
        'logLevel' : str(args.logLevel),
        'suffix' : str(args.suffix)
    }
    
    # console logging
    logger = logging.getLogger('measure_linearity')
    logger.setLevel(getattr(logging, params['logLevel']))
    
    ## console handler
    ch = logging.StreamHandler()
    ch.setLevel(getattr(logging, params['logLevel']))

    ## set logging format
    formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s", "%Y-%m-%d %H:%M:%S")
    ch.setFormatter(formatter)

    ## add handlers to logging object
    logger.addHandler(ch)    

    # some filename filters viz. ramp number, exposure number
    files_to_process = []
    for i in os.listdir(params['dataPath']):
        rnum = int(i.split('_')[1].strip('R'))
        if rnum == params['rampNum']:
            enum = int(i.split('_')[3].strip('.noref.fits').strip('N'))
            gnum = int(i.split('_')[2].strip('M'))
            if enum >= params['expNumLo'] and enum <= params['expNumHi'] and gnum >= params['grpNumLo'] and gnum <= params['grpNumHi']:
                files_to_process.append(params['dataPath'] + i)

    # section of the image to use
    sect_1_x = (int(params['sect'].split(',')[0].strip()), int(params['sect'].split(',')[2].strip()))
    sect_1_y = (int(params['sect'].split(',')[1].strip()), int(params['sect'].split(',')[3].strip()))

    # process
    inttime 	= [] 
    sect_1_data = []
    for idx, r_file in enumerate(files_to_process):
        print "file " + str(idx) + " of " + str(len(files_to_process)-1)
        data, hdr = read_FITS_file(r_file)

        r_1_section = data[sect_1_y[0]:sect_1_y[1]+1, sect_1_x[0]:sect_1_x[1]+1]
        sect_1_data.append(r_1_section)
        inttime.append(hdr['INTTIME'])

    # sort both sect_1_data and inttime by inttime ascending
    sort_indexes 	= np.argsort(inttime)
    inttime_sorted 	= [inttime[i] for i in sort_indexes]
    sect_1_data_sorted 	= [sect_1_data[i] for i in sort_indexes]

    # determine correction
    ## create coeffs array
    coeffs = []
    for i in range(params['fitCoeff']+1):
        coeffs.append(np.zeros(((sect_1_y[1]-sect_1_y[0]) + 1, (sect_1_x[1]-sect_1_x[0]) + 1)))
   
    ## create nonlinearities array
    nonlinearities = np.zeros(((sect_1_y[1]-sect_1_y[0]) + 1, (sect_1_x[1]-sect_1_x[0]) + 1))
    
    ## create bad pixel map array
    badpixmap = np.ones(((sect_1_y[1]-sect_1_y[0]) + 1, (sect_1_x[1]-sect_1_x[0]) + 1))    
    
    '''print badpixmap.shape
    print nonlinearities.shape
    for i in range(params['fitCoeff']+1):
      print coeffs[i].shape'''
    
    ## CDS files. this means that the flux from 0->1.4s is "lost", i.e. a 0 recorded signal at the start should actually be non-zero, see *** 
    first_file_data = sect_1_data_sorted[0]
    sect_1_data_sorted = sect_1_data_sorted - sect_1_data_sorted[0]

    bad_pix_count = 0
    for jj in range(0, (sect_1_y[1]-sect_1_y[0]) + 1):
        for ii in range(0, (sect_1_x[1]-sect_1_x[0]) + 1):
            this_pixel_values = []
            flagged_pixel = False
            ## obtain all values for this pixel
            for t in sect_1_data_sorted:
                this_pixel_values.append(t[jj, ii])
                
            ## get linear fit to slope
            this_pixel_vals_in_fit_lin	 	= [i for i,j in zip(this_pixel_values, inttime_sorted) if j > float(params['linFit'].split(',')[0]) and j < float(params['linFit'].split(',')[1])]
            this_pixel_inttimes_in_fit_lin	= [j for i,j in zip(this_pixel_values, inttime_sorted) if j > float(params['linFit'].split(',')[0]) and j < float(params['linFit'].split(',')[1])] 
            
            ## *** get the average flux accumulation rate for this pixel (per readout time) for times which have been linearly fitted
            this_pixel_rate_in_fit_lin = [(this_pixel_vals_in_fit_lin[i]-this_pixel_vals_in_fit_lin[i-1])/(this_pixel_inttimes_in_fit_lin[i]-this_pixel_inttimes_in_fit_lin[i-1]) for i in range(1, len(this_pixel_vals_in_fit_lin))]
            average_rate = np.mean(this_pixel_rate_in_fit_lin)
            
            ## *** add this rate on to account for the single readout missed (CDS) 
            this_pixel_vals_in_fit_lin + average_rate
            
            ## work out linear fit
            linfit_coeffs 			= np.polyfit(this_pixel_inttimes_in_fit_lin, this_pixel_vals_in_fit_lin, 1)

            ## fit slope to uncorrected data within fit limit
            this_pixel_vals_in_fit 	= [i for i,j in zip(this_pixel_values, inttime_sorted) if i > float(params['yFit'].split(',')[0]) and i < float(params['yFit'].split(',')[1])]
            this_pixel_inttimes_in_fit	= [j for i,j in zip(this_pixel_values, inttime_sorted) if i > float(params['yFit'].split(',')[0]) and i < float(params['yFit'].split(',')[1])]

            if len(this_pixel_vals_in_fit) == 0: # if count rate is too low that array is empty..
                flagged_pixel = True
                ideal_signal_vals_in_fit 	= []
                this_pixel_correction_coeffs 	= []
                this_pixel_values_cor 		= []
                this_pixel_nonlinearity_uncor 	= []
                this_pixel_nonlinearity_cor	= []
                ## store coeffs as nan
                for idx in range(len(coeffs)):  
                    coeffs[idx][jj, ii] = np.nan
                ## store nonlinearity as nan
                nonlinearities[jj, ii] = np.nan
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
                for t in sect_1_data_sorted:
                    this_corrected_value = np.polyval(this_pixel_correction_coeffs, t[jj,ii])
                    this_pixel_values_cor.append(this_corrected_value)

                ## find nonlinearity as percentage pre/post and store
                this_pixel_nonlinearity_uncor = abs((this_pixel_values - np.polyval(linfit_coeffs, inttime_sorted))/np.polyval(linfit_coeffs, inttime_sorted))*100
                this_pixel_nonlinearity_cor = abs((this_pixel_values_cor - np.polyval(linfit_coeffs, inttime_sorted))/np.polyval(linfit_coeffs, inttime_sorted))*100
                nonlinearities[jj, ii] = np.median(this_pixel_nonlinearity_cor)        
                
                ### test if the corrected linearity is below the "bad" threshold:
                if np.median([this_pixel_nonlinearity_cor[i] for i in range(len(this_pixel_vals_in_fit))]) > params['badThresh']:
                    flagged_pixel = True
                    print "FAILED LINEARITY: Bad pixel found @ " + str(ii) + ", " + str(jj) + " with a linearity of " + str(np.median(this_pixel_nonlinearity_cor))
                    badpixmap[jj, ii] = np.nan
                    bad_pix_count = bad_pix_count + 1
                    print "NUMBER OF BAD PIXELS: " + str(bad_pix_count)
                    
                ### test if the rate is effectively zero
                rates = [(this_pixel_values[i]-this_pixel_values[i-1])/(inttime_sorted[i]-inttime_sorted[i-1]) for i in range(len(this_pixel_vals_in_fit))]
                if np.mean(rates) - (params['sig'] * np.std(rates)) < 0:
                    flagged_pixel = True
                    print "FAILED RATE: Bad pixel found @ " + str(ii) + ", " + str(jj) + " with a rate of " + str(np.mean(rates))
                    badpixmap[jj, ii] = np.nan
                    bad_pix_count = bad_pix_count + 1
                    print "NUMBER OF BAD PIXELS: " + str(bad_pix_count)
       	    ## plot graphs if requested
            if params['plotA'] or (params['plotB'] and flagged_pixel):
                fig = plt.figure()

                fig.add_subplot(321)
                plt.plot(inttime_sorted, this_pixel_values, 'rx-', label="pixel " + str(ii + sect_1_x[0]) + ", " + str(jj + sect_1_y[0]) + " uncorrected")
                try:
                    plt.plot(inttime_sorted, this_pixel_values_cor, 'gx-', label="pixel " + str(ii + sect_1_x[0]) + ", " + str(jj + sect_1_y[0]) + " corrected")
                except ValueError:
                    pass
                plt.plot(inttime_sorted, np.polyval(linfit_coeffs, inttime_sorted), 'b-', label="Linear fit")

                plt.plot(inttime_sorted, [65536 for i in inttime_sorted], 'k-', label="16bit ADC limit")
                plt.annotate('ADC Limit = 65536 ADU', xy=(200, 65536+1500), xytext=(200, 65536+1500))
                plt.plot(inttime_sorted, [float(params['yFit'].split(',')[1]) for i in inttime_sorted], 'k--', label="Fit limit")
                plt.plot(inttime_sorted, [linfit_coeffs[1] for i in inttime_sorted], 'k-', label="Bias")
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
                plt.xlabel("\"Linear\" signal (ADU)")
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
                plt.plot([this_pixel_values[i] for i in range(1, len(this_pixel_values)-2)], [(this_pixel_values[i]-this_pixel_values[i-1])/(inttime_sorted[i]-inttime_sorted[i-1]) for i in range(1, len(this_pixel_values)-2)], 'rx-', label="uncorrected")
                try:
                    plt.plot([this_pixel_values_cor[i] for i in range(1, len(this_pixel_values_cor)-2)], [(this_pixel_values_cor[i]-this_pixel_values_cor[i-1])/(inttime_sorted[i]-inttime_sorted[i-1]) for i in range(1, len(this_pixel_values_cor)-2)], 'gx-', label="corrected")
                except ValueError:
                    pass
                plt.ylabel("Rate (ADU/s)")
                plt.xlabel("Signal (ADU)")
                leg = plt.legend(loc='lower left', prop={'size':11})
                leg.get_frame().set_alpha(0.75)

                ## plot image postage stamp (end integration)
                ax = fig.add_subplot(325)
                win_y_s = 0 if jj-5 < 0 else jj-5
                win_y_e = sect_1_y[1]-sect_1_y[0] if jj+5 > sect_1_y[1]-sect_1_y[0] else jj+5
                win_x_s = 0 if ii-5 < 0 else ii-5
                win_x_e = sect_1_x[1]-sect_1_x[0] if ii+5 > sect_1_x[1]-sect_1_x[0] else ii+5
                cax = ax.imshow(t[win_y_s:win_y_e, win_x_s:win_x_e], interpolation="nearest", aspect='auto')
                ax.scatter(ii-win_x_s,jj-win_y_s, marker='x', color='white')
                cb_ax = fig.colorbar(cax)

                plt.show()

        print "row " + str(jj+1) + " of " + str((sect_1_y[1]-sect_1_y[0]) + 1)

    if params['fits']:
        pyfits.writeto(params['suffix'] + "lin.fits", data=nonlinearities) 
        pyfits.writeto(params['suffix'] + "lin_bad.fits", data=badpixmap)         
        for i in range(params['fitCoeff']+1):   
            pyfits.writeto(params['suffix'] + "lin_coeffs_" + str(i) + ".fits", coeffs[params['fitCoeff']-i])          # params['fitCoeff']-i to reverse so that coeff_0 is the constant..        
        #pyfits.writeto(params['suffix'] + "data.fits", data=sect_1_data_sorted[len(sect_1_data_sorted)-1])   
    #n, bins, patches = plt.hist(nonlinearities.flatten(), 1000)        # now has NaN.
    #print np.median(nonlinearities.flatten())
    #plt.xlabel("Nonlinearity (%)")
    #plt.ylabel("Number")
    #plt.show()            

