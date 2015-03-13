ioi_pipeline
============

# Installation on lt-qc

Inside /usr/local/bin/:

$ git clone https://github.com/LivTel/ioi_pipeline

and recursively chown the resulting directory to data.

Requires:

* pyfits > wget https://pypi.python.org/packages/source/p/pyfits/pyfits-3.3.tar.gz, gunzip, untar and "python setup.py install" into src/ directory. Installing into the default dir breaks astrometry.net (reliant on an older version of pyfits) on lt-qc
* asciidata > wget http://www.stecf.org/software/PYTHONtools/astroasciidata/source/asciidata-1.1.1.tar.gz, gunzip, untar, "python setup.py install" into default dir.
* sExtractor > yum install sextractor
* Any additional python modules will be indicated on execution of the pipeline. Use "pip install [package]" to install them into the default dir

# Running from lt-qc

[data@lt-qc ~/rmb]$ set IOI\_PIPE\_BASE = /usr/local/bin/ioi\_pipeline

[data@lt-qc ~/rmb]$ python $IOI\_PIPE\_BASE/src/run\_pipe.py --pa $IOI\_PIPE\_BASE/config/paths.ini --pi $IOI\_PIPE\_BASE/config/pipeline.ini --rlo 1 --rhi 1 --d 20150227 --wd test --p $IOI\_PIPE\_BASE/test/

A breakdown of the parameters:

[data@lt-qc ~/rmb]$ python $IOI\_PIPE\_BASE/src/run\_pipe.py --h

Usage: run\_pipe.py [options]

Options:
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;-h, --help           show this help message and exit

  General:<br/>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;--p=DATAPATH       path to data<br/>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;--wd=WORKINGDIR    path to working directory<br/>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;--o                clobber working directory?<br/>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;--pa=PATHSCFGPATH  path to paths config file<br/>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;--pi=PIPECFGPATH   path to pipeline config file<br/>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;--rlo=MINRUNNUM    minimum multrun number<br/>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;--rhi=MAXRUNNUM    maximum multrun number<br/>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;--elo=MINEXPNUM    lowest exposure number to use<br/>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;--ehi=MAXEXPNUM    highest exposure number to use<br/>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;--log=LOGLEVEL     log level (DEBUG|INFO|WARNING|ERROR|CRITICAL)<br/>

  LT specific general parameters:<br/>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;--d=DATE           date (YYYYMMDD)<br/>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;--dlo=MINDITHNUM   minimum dither number<br/>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;--dhi=MAXDITHNUM   maximum dither number

  Teledyne specific general parameters:<br/>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;--glo=MINGRPNUM    lowest group number to use<br/>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;--ghi=MAXGRPNUM    highest group number to use<br/>

Some notes:

* The --pa and --pi flags are mandatory, and shouldn't change from above.
* The --*lo and --*hi flags are inclusive and can be used to selected single multrun/dithers, i.e. --rlo 1 --rhi 1 (as above)
* The date flag (--d) must be set, otherwise it'll use all dates in a given directory

# An overview of the pipeline process

i) **REFERENCE PIXEL SUBTRACTION**.

The array contains 4 rows (and columns, but I don't use these in the subtraction scheme i'm using at the moment) adjacent to each border of the array that are connected to readout, but are not active (reference pixels). We are operating in 32 output mode, so the 2048 pixel array "width" is split into 64 bit sections, each read out through one output.

For each section, I fit a slope to the mean of the top and bottom reference pixels (in actuality, I create two ramps per section, one for odd numbered pixels and one for even). For each pixel in the section, I evaluate the "bias" contribution using this slope, and subtract it off. I call this mode "RAMP" in the config/pipeline.ini file.

There is also a "CONSTANT" mode, which as the name suggests, takes the mean for each section of BOTH the top and bottom reference pixels, and subtracts.

"COLUMN" mode was an attempt at trying to remove low frequency banding along the vertical axis of the array. In short, it uses the left and right border reference pixels, smooths the data with a gaussian kernel, and subtracts. It worked to an extent but only had a marginal improvement data quality, so in the interests of simplicity I dropped it.

"FULL" does both "RAMP" and "COLUMN" modes.

"TRIM" just subtracts off the reference pixels without subtraction.


ii) **FRAME COMBINATION**.

Either CDS or Fowler. Not much to note here, except that in addition to performing the combination, I calculate the rate of charge accumulation. This is done (at the moment) by subtracting the second file in the sequence by the first, and dividing by the difference between their INTTIME values. This yields a value for counts/s, and is necessary for accurate nonlinearity correction (see next).


iii) **NONLINEARITY CORRECTION**.

I have generated a series of "correction coefficient files" by using an up-the-ramp (UTR) sequence. For each pixel, the integration time is plotted against the counts observed. A linear fit is made to the early part of this sequence (when the charge accumulation is assumed to be linear), and another fit made between the values derived from the linear fit (the "actual" counts) and the observed counts. At the moment, the order of this fit is 3, hence there are three filenames for the linearity correction (paths\_*.ini > nonlinearity\_corrections > coeff\_filenames: lin\_coeff\_0\_011214.fits,lin\_coeffs\_1\_011214.fits,lin\_coeffs\_2\_011214.fits,lin\_coeffs\_3\_011214.fits). This makes for easy application of the correction in the pipeline (it's just array multiplication..).

A caveat to this is that because we are combining the frames, we "lose" 1.4s (one readout time in 32 output mode) of signal. If the rate of accumulation is the same between the UTR sequence used to find the nonlinearity coefficients and the target you're observing, then this isn't a problem. However, if they're different, which they obviously will be in any real scenario, then you could possibly end up on the wrong part of the linearity slope, e.g. consider:

* nonlinearity correction derived with count rate of 50ADU/s, so 50*1.4=70ADU of flux "lost" from frame combination.
* target count rate of 250ADU/s, so 250*1.4=350ADU of flux "lost" from frame combination.

This means that in your combined frames, a "0" value actually correspond to a difference in signal of 280ADU.

In the process above I correct for this by adding 1.4s worth of flux for each pixel. This 1.4s of flux is derived in the same way the rate was derived for FRAME COMBINATION above.

When generating the correction files, I also output a bad pixel mask. Bad pixels are defined as both a) no accumulation of charge with time, and b) having a percentage nonlinearity ((observed&corrected-actual)/actual)*100 > some threshold value (currently 0.8% - arbitrary but derived by looking at the distribution of linearity values). A .FITS file with the linearity values is outputted and currently resides as "lin\_[DATE]" in the config/nonlinearity\_corrections dir. OK pixels are stored with a value of "1", bad pixels are stored with a value of "NaN". This assignment makes it very easy (and is also necessary) for the pipeline to combine/use masks. This bad pixel mask is named "lin\_bad\_[DATE]" and is also currently in the config/nonlinearity\_corrections dir.

The script I used to measure the linearity is found in src/measure\_linearity.py. I am currently in the process of refactoring it.

iv) **FLATFIELDING**.

Don't think I need to write much about this, each combined frame is divided through by the flat. An extra point to note is that when I generated flats previously, I also generated a bad pixel mask to find coefficients < some threshold value (0.35). OK pixels are stored with a value of "1", bad pixels are stored with a value of "NaN". This bad pixel mask is named "flat\_bad\_[DATE]" and is currently in the config/flatfielding dir.

The script I used to make a flat is found in src/measure\_flatfield\_from\_sky.py.

v) **BAD PIXEL MASKING**.

A cursory glance over raw array data and you'll see that there are a few notable blotches. There is an option to mask this data by feeding in any number of bad pixel masks (I used two, one from flatfielding, one from nonlinearity corrections). As mentioned before, OK pixels are stored with a value of "1", bad pixels are stored with a value of "NaN".


vi) **SKY SUBTRACTION**.

A median sky is constructed from ALL the combined frames available (corresponding to one per dither). When making the median sky image, each combined frame is first scaled. Scaling is done by iteratively sigma clipping the science images (sigma = 3 = pipeline.ini > sky\_subtraction > bg\_sigma\_clip) to find the background "sky" value, taking the mean of this clipped dataset (should this be a median?), and dividing it by the mean (median?) of reference sky frame values to find a scaling factor which is then applied to the science image. The reference sky frame values are taken as the sky values of the first combined frame, i.e. the scaling factor for the first frame == 1.

When subtracting the sky contribution from each combined frame, the median sky is scaled again: sigma clip the science frame to find the sky values, take the mean, then divide this mean by the median of the constructed median sky. The scaled median sky is then subtracted off.

Any cold/hot pixels in the array will propagate through to the final median sky. If bad pixel masking has been selected, then each NaN in the array is locally median smoothed (box size = 5 = pipeline.ini > sky\_subtraction > smoothing\_box\_size).


vii) **IMAGE REGISTRATION/STACKING**.

Uses alipy (http://obswww.unige.ch/~tewes/alipy/) which in turn uses sExtractor to make quads (similar to astrometry.net) of sources and find per frame transformations. It can either use IRAF imalign to perform alignment given each transformation, or python's affine\_transform (pipeline.ini > registration > registration\_algorithm). It's currently set to use the latter.


**OTHER BITS**:

In pipeline.ini, you'll find:

* each routine has a "do" parameter, allowing toggling of that particular stage. Note that if you set frame\_combination.do to 0, it will not be able to proceed past reference subtraction.
* "quit" parameters for most routines, telling the pipeline to stop after this stage. This was a requirement when we were considering on-site processing, you should never need to set quit for any routine.
* a "write\_LT\_file" parameter, which, if quit is NOT set, will make \_1.fits images. If it IS set, it will make intermediate images (\_0.fits) with the data processed up to the point at which it was told to quit.
* "hard" parameters to produce temporary files showing the output of the processed data up to the point.


**IF YOU WANT TO CHANGE THE ORDER OF PROCESSES**:

I suspect that I should do this as it involves fiddling with the code.

For reference, the code blocks containing nonlinearity correction, flatfielding, bad pixel masking and sky subtraction (delineated by "# -------------------------") *should* all be switchable between each other. They take in the same format data and output the same format data. This is NOT true of reference subtraction, frame combination and registration/stacking, e.g. reference subtraction expects all the fowler pair files for a given dither, whereas nonlinearity correction expects only the combined frame.

