import pylab as plt
import numpy as np
import sys
import pyfits
from scipy.optimize import curve_fit
import math
import matplotlib

fname = sys.argv[1]

fits = pyfits.open(fname)

# util functions
def sf(num, sig_figs):
    try:
        rtn = round(num, -int(math.floor(math.log10(abs(num))) - (sig_figs - 1)))
        return rtn
    except ValueError:
        return 0.

def gauss(x,a,x0,sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))
  
lim_0 = (-1000, 1000)
lim_1 = (0.7, 1.3)
lim_2 = (-0.000006, 0.000009)
lim_3 = (-0.1*pow(10, -9), 0.1*pow(10, -9))
nbins = 50

lim = [lim_0, lim_1, lim_2, lim_3]

font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 12}
matplotlib.rc('font', **font)

fig = plt.figure()
for idx, ext in enumerate(fits): 
    data          = ext.data[0].flatten()
    data_filtered = data[np.where(np.logical_and(data > lim[idx][0], data < lim[idx][1]))]
    
    n, bins = np.histogram(data_filtered, bins=nbins)
    popt, pcov = curve_fit(gauss,[bins[i]+((bins[i+1]-bins[i])/2) for i in range(len(bins)-1)],n,p0=[1,np.median(data),np.median(data)]) 
    
    plt.subplot(2,2,idx+1)
    plt.plot(bins,gauss(bins,*popt),'k--', linewidth=2, label='fit (centre=' + str(sf(popt[1],3)) +')')
    n, bins, patches = plt.hist(data_filtered, bins=nbins, label='data')
    plt.setp(patches, 'facecolor', 'w', 'alpha', 1)
    plt.ticklabel_format(axis='x', style='sci', scilimits=lim[idx])
    plt.xlim(lim[idx])
    plt.title("Order " + str(idx))
    plt.xlabel("Value")
    plt.ylabel("Number")
    plt.legend(loc="lower right", fontsize=10)
    
plt.tight_layout()
plt.show()