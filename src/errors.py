import inspect
import time

class errors():
    def __init__(self, logger):
        self.codes = []
        self.logger = logger
        self.current_code = 0
        self.err_dict = {
            0   :       "Success.",
            1   :       "Working directory not empty and clobber not set.",
	    2	:	"File doesn't conform to LT nomenclature. Skipping.",
	    3	:	"Unable to find any files to process.",
	    4   :       "Couldn't find NOUTPUTS key in header.",
	    5   :       "Unknown reference subtraction method.",
	    6   :       "Return data or header is Nonetype.",
	    7   :       "The nonlinearity correction coefficient file doesn't exist.",
            8   :       "Couldn't find key in .ini files. Either the path to the file or the key doesn't exist.",
            9   :       "A nonlinearity correction coefficient file has not been defined.",
            10  :       "Unknown frame combination method.",
            11  :       "Too few files to perform CDS combination, or [files] is NoneType.",
            12  :       "[data] is NoneType.",
            13  :       "Failed to make stack",
            14  :       "[data_sky] has not been set.",
            15  :       "Insufficient files to perform sky subtraction.",
            16  :       "Insufficient files to perform stacking.",
            17  :       "No suitable files found in path.",
            18  :       "Unable to obtain id for registering.",
            19  :       "Unknown registration method.",
            20  :       "Insufficient bad pixel masks set to combine.",
            21  :       "Flatfielding correction coefficient file doesn't exist",
            22  :       "The flatfielding correction coefficient file has not been defined.",
            23  :       "Error occurred while performing linear correction operation, it is possible that the coefficient or bad pixel mask files are not consistent in size with the files being processed.",
            24  :       "Error occurred while performing flatfielding operation, it is possible that the coefficient file is not consistent in size with the files being processed.",
            25  :       "Sky clipping was unsuccessful, no values returned.",
            26  :       "Bad pixel mask size mismatch.",
            27  :       "Bad pixel mask file does not exist.", 
            28  :       "Bad pixel mask file(s) not defined.",  
            29  :       "Error occurred while performing bad pixel masking operation, it is possible that a mask file is not consistent in size with the files being processed.",
            30  :       "Bad pixels remain in sky field.",
            31  :       "Too few files to perform fowler combination, or [files] is NoneType.",
            32  :       "Unequal number of pedestal/readout files or exposure number mismatch.",
            33  :       "Unknown stacking method.",
            34  :       "Value of ASICGAIN gives an unmapped GAIN/EPERDN, ignoring.",
            35  :       "EXTNAME is not valid, expecting either IM_NONSS or IM_SS.",
            36  :       "EXTNAME is not valid, expecting SK_SS.",
            37  :       "Frames not taken with enhanced horizontal clocking. Linearity correction may not be optimal."
        }
        
    def set_code(self, code, is_warning=False, is_critical=False):
        self.current_code = code
        self.codes.append(code)
        parent_calling_func = inspect.stack()[1][3]
        if not is_warning and not is_critical:
            self.logger.info("[" + parent_calling_func + "] " + "Code: " + str(self.current_code))
            self.logger.info("[" + parent_calling_func + "] " + "Message: \"" + self.err_dict[self.current_code] + "\"")
        elif is_warning:
            self.logger.warning("[" + parent_calling_func + "] " + "Code:" + str(self.current_code))
            self.logger.warning("[" + parent_calling_func + "] " + "Message: \"" + self.err_dict[self.current_code] + "\"")
        elif is_critical:
            self.logger.critical("[" + parent_calling_func + "] " + "Code: " + str(self.current_code))
            self.logger.critical("[" + parent_calling_func + "] " + "Message: \"" + self.err_dict[self.current_code] + "\"")
            exit(self.current_code)
            
    def get_code(self):
        return self.current_code
        
