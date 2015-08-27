'''
name:		FITSFile.py
author:		rmb

description: 	A class to handle FITS files
'''

import pyfits
import os
import numpy as np

class read:
    '''
    a class for reading FITS files
    '''
    def __init__(self, filepath):
        self.filePath	= filepath
	self.hduList	= None
	self.headers	= None
	self.data	= None
	self.fileOpen	= False
	
    def closeFITSFile(self):
        '''
        close FITS file cleanly
        '''
        if self.hduList is not None:
            self.hduList.close()
            self.fileOpen = False
            return True
        else:
            return False	

    def getData(self, HDU):
        '''
        get data from a specified HDU
        '''
        if self.hduList is not None:
            self.data = self.hduList[HDU].data.astype(np.float16)
        else:
            return False   
          
    def getHeaders(self, HDU):
        '''
        retrieve headers from a specified HDU
        '''
        if self.hduList is not None:
            self.headers = self.hduList[HDU].header
        else:
            return False       
          
    def openFITSFile(self):
        '''
        check a FITS file exists and try to open it
        '''
        if os.path.exists(self.filePath):
            try:
                self.hduList = pyfits.open(self.filePath)
            except IOError:
                return False  
            self.fileOpen = True
            return True
        else:
            return False
          
class write:
    '''
    a class for writing FITS files
    '''         
    def __init__(self, filepath, data, headers=None):
        self.filePath   = filepath
        self.data       = data
        self.headers    = headers
        
    def setData(self, data):
        self.data = data
        
    def setHeaders(self, hdr):
        self.headers = hdr        
        
    def writeFITSFile(self, allow_append):
        self.data = self.data.astype(np.float32, copy=False)    # 64->32 bit (float)
        #hdu.scale('int16', bzero=32768)                        # 64->16 bit (int)
        if allow_append and os.path.exists(self.filePath):
            pyfits.append(self.filePath, self.data, self.headers)        
        else:
            hdu = pyfits.PrimaryHDU(self.data, self.headers)
            hdu.writeto(self.filePath)    
