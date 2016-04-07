ioi_pipeline
============

# Overview

This package provides a pipeline to reduce data taken with IO:I.

# Installation

1. Clone the repository

2. Install dependencies. Rather than list them all, it is probably more efficient to install 
them as and when they are called for at runtime (there aren't that many!). Use `sudo pip install` 
or `sudo easy_install`, although see the notes about pyfits in **Known Issues**.

3. Move the altered src/sm/robust/norms.py to your local copy of statsmodels (`robust/`).

# Configuration Files

The configuration files are kept in `config/`. There are two files: `paths.ini` and `pipeline.ini`. 

The former defines paths to the base and various subdirectories. Typically only [**path\_base**] 
should need to be altered as the other paths are defined relative to this. The two other subheadings, 
[**flatfielding**] and [**nonlinearity_corrections**] define paths and files to be used in each of 
the corresponding processes. The generation and use of these files is further discussed in the **An Overview of the Pipeline Process** 
section.

The latter, `pipeline.ini`, defines the various pipeline specific parameters. Most, if not all, routines have:

* a [**do**] parameter allowing toggling of that particular stage. Note that if you set [**frame\_combination.do**] to 0, it will not be able to proceed
past reference subtraction
* a [**quit**] parameter telling the pipeline to stop after this stage. You should never 
need to set quit for any routine, except maybe during testing
* a [**hard**] parameter to produce temporary files showing the output of the processed data up to that point

The routine-specific parameters are commented, and should generally be left alone unless you know what they do.

# Invoking the Pipeline

You will be required to specify both configuration files using the `--pa` (paths.ini) and `--pi` (pipeline.ini) 
flags at the command prompt. The help blurb will tell you more about the various parameters available, but generally 
for LT specific files, you should be setting the ramp low/hi (`--rlo`, `--rhi` respectively) inclusive range, the date of the 
run (`--d`), the working directory (`--wd`) and the data path (`--p`). The date flag *must* be set, otherwise it'll use all dates 
in a given directory. An example command is:

`python run_pipe.py --pa paths.ini --pi pipeline.ini --rlo 1 --rhi 1 --d 20150227 --wd test --p test/`

If you're working with raw Teledyne files, the group low/hi (`--glo`, `--ghi`) inclusive range needs to be set instead of the ramp 
range and date. 

The help command can be invoked by `python run_pipe.py --h`.

# Known Issues

Installing pyfits into the default directory will break older versions of astrometry.net. Instead, local install into the 
`src/` directory. 

# An Overview of the Pipeline Process

The pipeline process is different for Teledyne and LT specific reduction modes. For example, the concept of dithering 
is not understood when processing raw Teledyne files as there is no position information held in the 
filename or header. This makes some of the stages (e.g. sky subtraction) somewhat ambiguous.

The process defined here relates more to the LT specific reduction mode. A more in depth description of the methodology 
is given in Barnsley, 2016.

## i) **REFERENCE PIXEL SUBTRACTION**.

The array contains four rows and four columns adjacent to each border of the array that are connected to readout, 
but are not active (reference pixels). Operating in 32 output mode, the 2048 pixel array "width" is split into 64 sections, 
each read out through one output.

The [**reference_subtraction.method**] parameter in `pipeline.ini` defines how the 
reference pixels are used. These modes are discussed in the paper.

## ii) **FRAME COMBINATION**.

The [**frame_combination.method**] in `pipeline.ini` defines how the resulting frame 
is generated. This can be set to either **CDS** or **FOWLER**. At this stage, the rate of charge accumulation is also calculated. 
This is done by subtracting the second file in the sequence by the first, and dividing by the difference between their INTTIME values. 
This yields a value for counts/s, and is necessary for accurate nonlinearity correction.

Checks are done to ensure that there are equal numbers of pedestal/readout frames when performing Fowler sampling. 

Some headers are also appended/updated in this section, removing Teledyne run (**SEQNUM\_R**), group (**SEQNUM\_M**) and exposure 
(**SEQNUM\_N**) numbers as well as **SEQNNAME** and **INTTIME** keywords, all of which have no meaning after frame combination. **MJD**, 
**DATE**, **DATE-OBS** and **UTSTART** are all also amended, in addition to setting the correct **EXPTIME**.

## iii) **NONLINEARITY CORRECTION**.

A series of "correction coefficient files" have been generated using an up-the-ramp (UTR) sequence. For each pixel, the integration time is plotted against the counts observed. A linear fit is made to the early part of this sequence (when the charge accumulation is assumed to be linear), and another fit made between the values derived from the linear fit (the "actual" counts) and the observed counts. At the moment, the order of this fit is 3, producing a single MEF file with each order as an extension ([**nonlinearity_correction.coeff\_filename**] in `paths.ini`). This makes for easy application of the correction in the pipeline (it's just array multiplication..).

A caveat to this is that because we are combining the frames, we "lose" 1.4s (one readout time in 32 output mode) of signal. If the rate of accumulation is the same between the UTR sequence used to find the nonlinearity coefficients and the target you're observing, then this isn't a problem. However, if they're different, which they obviously will be in any real scenario, then you could possibly end up on the wrong part of the linearity slope, e.g. consider:

* nonlinearity correction derived with count rate of 50ADU/s, so 50*1.4=70ADU of flux "lost" from frame combination.
* bright target with count rate of 250ADU/s, so 250*1.4=350ADU of flux "lost" from frame combination.

This would mean that in your combined frames, a "0" value actually correspond to a difference in signal of 350-70 = 280ADU.

In the process above I correct for this by adding 1.4s worth of flux for each pixel (this has been done to the calibration frames also). This 1.4s of flux is derived in the same way the rate was derived for FRAME COMBINATION above.

When generating the correction files, a bad pixel mask is also outputted. Bad pixels are defined as either a) no accumulation of charge with time, or b) having a percentage nonlinearity ((observed&corrected-actual)/actual)\*100 is greater than some threshold value (currently 0.8% - arbitrary but derived by looking at the distribution of linearity values). This file is defined in ([**nonlinearity\_correction.bad\_pix\_filename**] in `paths.ini`.

The script used to make the current nonlinearity cube is `src/make_nonlinearity_cube.py`. 

## iv) **FLATFIELDING**.

Each frame can be divided through by a flat, the path to which is defined by the [**flatfielding.coeff\_filename**] in `pipeline.ini`. There is also a bad pixel mask for this stage. Bad pixels are defined as having a coefficient less than some value (0.35 currently). The path to this bad pixel mask is [**flatfielding.bad\_pix\_filename**] in `pipeline.ini`.

## v) **BAD PIXEL MASKING**.

There is an option to mask this data by feeding in any number of bad pixel masks (currently two - one from flatfielding, one from nonlinearity corrections). In these files, OK pixels are stored with a value of "1", bad pixels are stored with a value of "NaN", making it easy to combine masks. Which bad pixel masks are used is defined by the [**bad\_pix\_masking.which**] parameter in `pipeline.ini`.

## vi) **SKY SUBTRACTION**.

The sky is constructed using either a median or more robust estimators (Tukey, HuberT) to combine either all frames or only a frame's peers. Sky images are constructed on a per dither basis by offsetting each background to a reference average value.

## vii) **IMAGE REGISTRATION/STACKING**.

This routine uses alipy, which in turn uses sExtractor to make quads (similar to astrometry.net) of sources and find per frame transformations. It can either use IRAF imalign to perform alignment given each transformation, or python's affine\_transform ([**registration.registration\_algorithm**] in `pipeline.ini`).

Stacking can either be done by taking the mean (**MEAN**), coadding (**COADD**) or taking the mean and multiplying by the number of dithers (**MEANANDSCALE**). This can be specified in [**stacking.method**] in `pipeline.ini`.


