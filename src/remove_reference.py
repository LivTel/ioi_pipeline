import pyfits 
import os
import numpy as np
from scipy import signal
from utility import read_FITS_file, write_FITS_file  

def gauss_kern(size, sizey=None):
    """ Returns a normalized 2D gauss kernel array for convolutions """
    size = int(size)
    if not sizey:
        sizey = size
    else:
        sizey = int(sizey)
    x, y = np.mgrid[-size:size+1, -sizey:sizey+1]
    g = np.exp(-(x**2/float(size)+y**2/float(sizey)))
    return g / g.sum()

def smooth_image(im, n, ny=None) :
    """ blurs the image by convolving with a gaussian kernel of typical
        size n. The optional keyword argument ny allows for a different
        size in the y direction.
    """
    g = gauss_kern(n, sizey=ny)
    improc = signal.convolve2d(im,g, mode='same', boundary="symm", fillvalue=0)		# 'same' = same size as input array, boundary is mirrored.
    return(improc)

class ref_subtraction:
    def __init__(self, logger, err):
        self.logger     = logger
        self.err        = err  
        
    def execute(self, method, data, hdr, hard, strip, out=None, opt_hdr={}):
        if data is None:
            self.err.set_code(12, is_critical=True)
            
        # section of the image to use
        self.sections = []
        self.sections.append(0)
        try:
            if hdr['NOUTPUTS'] == 1:
                sections_interval = 2048
            elif hdr['NOUTPUTS'] == 4:
                sections_interval = 512
            elif hdr['NOUTPUTS'] == 32:
                sections_interval = 64
        except KeyError:
            self.err.set_code(4, is_warning=True)
            return None
            
        for i in range(sections_interval, 2048+1, sections_interval):
            self.sections.append(i)   
            
        if method=="TRIM":
            pass
        elif method=="CONSTANT":
            self._remove_constant(data, hdr)
        elif method=="RAMP":
            self._remove_ramp(data, hdr)
        elif method=="COLUMN":
            self._remove_smoothed_column(data, hdr)
        elif method=="FULL":
            self._remove_ramp(data, hdr)
            self._remove_smoothed_column(data, hdr)
        else:
            self.err.set_code(5, is_critical=True)
          
        if strip:
            data = data[4:2044, 4:2044]
        if hard:
            self.logger.info("[ref_subtraction.execute] Writing file out to " + out)  
            write_FITS_file(out, data, hdr, opt_hdr=opt_hdr)          
        return data, hdr  
       
    def _remove_constant(self, data, hdr):     
        # REMOVE CONSTANT PER OUTPUT OFFSET
        # --------------------------------
        # Remove constant from the upper ref to lower ref pixels
        # on a per output, odd-and-even basis
        for idx in range(len(self.sections)-1):
            x_lo = self.sections[idx]
            x_hi = self.sections[idx+1]

            this_out_ref_pix_odd = []
            this_out_ref_pix_even = []

            for i in range(x_lo,x_hi):
                if i % 2 == 0:
                    this_out_ref_pix_even.append(data[1:4, i])   	# miss off border reference pixels, these could be wrong
                    this_out_ref_pix_even.append(data[2044:2047, i])	# miss off border reference pixels, these could be wrong
                elif i % 2 == 1:
                    this_out_ref_pix_odd.append(data[1:4, i])		# miss off border reference pixels, these could be wrong
                    this_out_ref_pix_odd.append(data[2044:2047, i])	# miss off border reference pixels, these could be wrong     

            this_out_ref_pix_even_mean 	= np.mean(this_out_ref_pix_even)
            this_out_ref_pix_odd_mean 	= np.mean(this_out_ref_pix_odd) 

            for i in range(x_lo,x_hi):
                if i % 2 == 0:
                    data[0:2048, i] -= this_out_ref_pix_even_mean
                elif i % 2 == 1:
                    data[0:2048, i] -= this_out_ref_pix_odd_mean
        return data, hdr

    def _remove_ramp(self, data, hdr):
        # REMOVE PER OUTPUT OFFSET
        # ------------------------
        # Remove a slope from the upper ref to lower ref pixels
        # on a per output, odd-and-even basis.
        for idx in range(len(self.sections)-1):
            x_lo = self.sections[idx]
            x_hi = self.sections[idx+1]

            this_out_lower_ref_pix_odd = []
            this_out_lower_ref_pix_even = []
            this_out_upper_ref_pix_odd = []
            this_out_upper_ref_pix_even = []

            for i in range(x_lo,x_hi):
                if i % 2 == 0:
                    this_out_lower_ref_pix_even.append(data[1:4, i])   	# miss off border reference pixels, these could be wrong
                    this_out_upper_ref_pix_even.append(data[2044:2047, i])	# miss off border reference pixels, these could be wrong
                elif i % 2 == 1:
                    this_out_lower_ref_pix_odd.append(data[1:4, i])	# miss off border reference pixels, these could be wrong
                    this_out_upper_ref_pix_odd.append(data[2044:2047, i])	# miss off border reference pixels, these could be wrong                        
 
            this_out_lower_ref_pix_even_mean = np.mean(this_out_lower_ref_pix_even)
            this_out_upper_ref_pix_even_mean = np.mean(this_out_upper_ref_pix_even)
            this_out_lower_ref_pix_odd_mean = np.mean(this_out_lower_ref_pix_odd)
            this_out_upper_ref_pix_odd_mean = np.mean(this_out_upper_ref_pix_odd)

            coeffs_even = np.polyfit([1.5, 2055.5], [this_out_lower_ref_pix_even_mean, this_out_upper_ref_pix_even_mean], 1)
            coeffs_odd = np.polyfit([1.5, 2055.5], [this_out_lower_ref_pix_odd_mean, this_out_upper_ref_pix_odd_mean], 1)
        
            thisOffsetEven = np.polyval(coeffs_even, range(0,2048))
            thisOffsetOdd = np.polyval(coeffs_odd, range(0,2048))
            for i in range(x_lo,x_hi):
                if i % 2 == 0:
                    data[0:2048, i] -= thisOffsetEven
                elif i % 2 == 1:
                    data[0:2048, i] -= thisOffsetOdd
        return data, hdr                    
                        
    def _remove_smoothed_column(self, data, hdr):
	# REMOVE COLUMN SMOOTHED OFFSET
	# -----------------------------
	# Take a number of full reference columns, smooth with a 2D gaussian smoothing kernel
	# and subtract this contribution off.
        left_ref_pix = data[0:2048, 1:4]					# we want 5 columns as the kernel can only be an odd size. don't want 7 as that'd use border reference pixels.
        right_ref_pix = data[0:2048, 2044:2046]	
        col_ref_pix = np.append(left_ref_pix, right_ref_pix, axis=1)

 	# produce same size array as input because of global "mirrored" boundary condition which is necessary to retain same size y. therefore need to take CENTRAL column of output array.
        smoothed_image = smooth_image(col_ref_pix, 100, ny=2)			# note reversed n and ny.
        smoothed_column = smoothed_image[0:2048, 2]

        for i in range(0,2048, 1):
            data[0:2048, i] -= smoothed_column
        return data, hdr
      