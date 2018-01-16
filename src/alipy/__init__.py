"""
This is alipy, a python package to quickly, automatically, and robustly identify geometrical transforms between optical astronomical images.
Originally written by Malte Tewes, (v 2.0, 2012) and modifed to support IO:I by Robert Smith.
All the hard work is from Malte Tewes, with IO:I additions are primarily data handling, formatting and wrappers.

:Authors: Malte Tewes, Robert Smith
:License: GPLv3


"""

__author__ = "Robert Smith"
__copyright__ = "2012, Malte Tewes"
__version__ = "2.0ioi"

import imgcat, pysex, star, quad, ident, align


