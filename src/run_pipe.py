from session import *
from remove_reference import *
from correct_nonlinearity import *
from flatfielding import *
from utility import *
from combine import *
from bad_pixel_masking import *
from sky_subtraction import *
from register import *
from stack import *
import logging
import optparse
    
class run_pipe():
    def __init__(self, params, logger, err):
        self.params = params
        self.logger = logger
        self.err    = err
        self.session = None
        
    def go(self):
        params = self.params
        logger = self.logger
        err = self.err
        
        self.session = session(logger, err)    # create session object   

        # ----------------------
        # ---- read configs ----
        # ---------------------- 
        paths_cfg       = read_ini(params['pathsCfgPath'])
        pipe_cfg        = read_ini(params['pipeCfgPath'])
            
        ## paths
        try:
            path_base                           = str(paths_cfg['general']['path_base'].rstrip("/") + "/")
            full_path_config                    = str(path_base + paths_cfg['general']['path_config'].rstrip("/") + "/")
            full_path_man                       = str(path_base + paths_cfg['general']['path_man'].rstrip("/") + "/")
            full_path_src                       = str(path_base + paths_cfg['general']['path_src'].rstrip("/") + "/")
            full_path_nonlinearity_corrections  = str(paths_cfg['nonlinearity_corrections']['path_nonlinearity_corrections'].rstrip("/") + "/")    
            full_path_flatfielding_corrections  = str(paths_cfg['flatfielding']['path_flatfielding_corrections'].rstrip("/") + "/")      
            file_naming                         = str(pipe_cfg['general']['file_naming'])
            write_lt_file                       = bool(int(pipe_cfg['general']['write_lt_file']))  
            ## reference_subtraction
            do_refsub                           = bool(int(pipe_cfg['reference_subtraction']['do']))   
            refsub_method                       = str(pipe_cfg['reference_subtraction']['method'].upper()) 
            refsub_hard                         = bool(int(pipe_cfg['reference_subtraction']['hard']))  
            refsub_quit                         = bool(int(pipe_cfg['reference_subtraction']['quit']))  
            ## frame_combination
            do_fc                               = bool(int(pipe_cfg['frame_combination']['do']))          
            combination_method                  = str(pipe_cfg['frame_combination']['method'].upper())
            combination_f_pairs                 = int(pipe_cfg['frame_combination']['fowler_pairs'])
            combination_hard                    = bool(int(pipe_cfg['frame_combination']['hard'])) 
            combination_hard_rates              = bool(int(pipe_cfg['frame_combination']['hard_rates']))
            combination_quit                    = bool(int(pipe_cfg['frame_combination']['quit']))         
            ## nonlinearity_corrections
            do_nonlincor                        = bool(int(pipe_cfg['nonlinearity_corrections']['do']))       
            nonlincor_order                     = int(pipe_cfg['nonlinearity_corrections']['order'])
            nonlincor_coeff_path                = str(full_path_nonlinearity_corrections + paths_cfg['nonlinearity_corrections']['coeff_filename'])
            nonlincor_badpix_path               = str(full_path_nonlinearity_corrections + paths_cfg['nonlinearity_corrections']['bad_pix_filename'])        
            nonlincor_hard                      = bool(int(pipe_cfg['nonlinearity_corrections']['hard'])) 
            nonlincor_quit                      = bool(int(pipe_cfg['nonlinearity_corrections']['quit'])) 
            ## flatfielding
            do_ff                               = bool(int(pipe_cfg['flatfielding']['do']))   
            ff_coeff_path                       = str(full_path_flatfielding_corrections + paths_cfg['flatfielding']['coeff_filename'])
            ff_badpix_path                      = str(full_path_flatfielding_corrections + paths_cfg['flatfielding']['bad_pix_filename'])
            ff_hard                             = bool(int(pipe_cfg['flatfielding']['hard'])) 
            ff_quit                             = bool(int(pipe_cfg['flatfielding']['quit'])) 
            ## bad_pix_correction
            do_bp_masking                       = bool(int(pipe_cfg['bad_pix_masking']['do']))
            bp_which                            = str(pipe_cfg['bad_pix_masking']['which']).split(',')
            bp_post_reg_thresh                  = float(pipe_cfg['bad_pix_masking']['bad_post_register_threshold']) 
            bp_hard                             = bool(int(pipe_cfg['bad_pix_masking']['hard'])) 
            bp_quit                             = bool(int(pipe_cfg['bad_pix_masking']['quit'])) 
            ## sky_subtraction
            do_SS                               = bool(int(pipe_cfg['sky_subtraction']['do']))
            ss_smoothing_box_size               = int(pipe_cfg['sky_subtraction']['smoothing_box_size'])
            ss_bg_sigma_clip                    = float(pipe_cfg['sky_subtraction']['bg_sigma_clip'])
            ss_hard                             = bool(int(pipe_cfg['sky_subtraction']['hard']))
            ss_quit                             = bool(int(pipe_cfg['sky_subtraction']['quit'])) 
            ## registration
            do_registration                     = bool(int(pipe_cfg['registration']['do'])) 
            registration_algorithm              = str(pipe_cfg['registration']['registration_algorithm'].upper())
            registration_hard                   = bool(int(pipe_cfg['registration']['hard']))
            registration_fit_geometry           = str(pipe_cfg['registration']['fit_geometry'].lower())
            registration_quit                   = bool(int(pipe_cfg['registration']['quit'])) 
            ## stacking
            do_stacking                         = bool(int(pipe_cfg['stacking']['do']))
            stacking_hard                       = bool(int(pipe_cfg['stacking']['hard']))   
            stacking_method                     = str(pipe_cfg['stacking']['method'])
        except KeyError:
            err.set_code(8, is_critical=True)
            
        # ----------------------
        # ---- print header ----
        # ----------------------   
        if params['logLevel'] == "DEBUG" or params['logLevel'] == "INFO":
            with open(full_path_man + "HEADER") as f:
                for line in f:
                     print line.strip('\r\n')
                
        logger.info("[run_pipe.go] Using paths config: " + os.path.abspath(params['pathsCfgPath']))
        logger.info("[run_pipe.go] Using pipeline config: " + os.path.abspath(params['pipeCfgPath']))
   
        # ----------------------------
        # ---- find files in path ----
        # ---------------------------- 
        logger.info("[run_pipe.go] Searching for files in " + params['dataPath']) 
        if file_naming == "TELEDYNE":
            files = find_sort_files_teledyne(params['dataPath'], params['minRunNum'], params['maxRunNum'], params['minGrpNum'], params['maxGrpNum'], params['minExpNum'], params['maxExpNum'], logger, err)  
        elif file_naming == "LT":
            files = find_sort_files_LT(params['dataPath'], params['minRunNum'], params['maxRunNum'], params['minDithNum'], params['maxDithNum'], params['minExpNum'], params['maxExpNum'], params['date'], logger, err)
        else:
            err.set_code(17, is_critical=True)
        self.session.add_files(files)
    
        # -------------------------------   
        # ---- reference subtraction ----   
        # -------------------------------      
        if do_refsub:
            logger.info("[run_pipe.go] Beginning reference subtraction with method: " + refsub_method)
            self.session.add_amend_opt_header('L1REFSUB', 1, 'reference subtracted') 
            self.session.file_ext = ".noref"
            sub = ref_subtraction(logger, err)
            for idx_1, run in enumerate(self.session.start_names):
                for idx_2, dither in enumerate(run):
                    for idx_3, f in enumerate(dither):
                        logger.info("[run_pipe.go] Removing reference from run:" + str(idx_1+params['minRunNum']) + ", dither:" + str(idx_2+params['minDithNum']) + ", exp:" + str(idx_3+params['minExpNum']))
                        this_outPath = params['workingDir'] + str(idx_1+params['minRunNum']) + "_" + str(idx_2+params['minDithNum']) + "_" + str(idx_3+params['minExpNum']) + self.session.file_ext + ".fits"
                        rtn_data, rtn_hdr = sub.execute(method=refsub_method, data=self.session.file_data[idx_1][idx_2][idx_3], hdr=self.session.file_hdr[idx_1][idx_2][idx_3], \
                                                        hard=refsub_hard, strip=True, out=this_outPath, opt_hdr=self.session.opt_hdr) 
                        if rtn_data is not None and rtn_hdr is not None:
                            self.session.file_data[idx_1][idx_2][idx_3] = rtn_data
                            self.session.file_hdr[idx_1][idx_2][idx_3]  = rtn_hdr 
                        else:
                            err.set_code(6, is_critical=True)              
            if registration_quit:
                logger.info("[run_pipe.go] Returning with code: " + str(err.current_code))
                return err.current_code
        else:
              logger.info("[run_pipe.go] No reference subtraction requested.")   
              self.session.add_amend_opt_header('L1REFSUB', 0, 'reference subtracted')
        
        # ---------------------------   
        # ---- frame combination ----   
        # ---------------------------       
        if do_fc:
            logger.info("[run_pipe.go] Beginning frame combination with method: " + combination_method) 
            self.session.add_amend_opt_header('L1FMCO', 1, 'frame combined')
            self.session.add_amend_opt_header('L1FMCOME', combination_method, 'frame combination method')
            self.session.add_amend_opt_header('L1FMCOPA', combination_f_pairs, 'num pairs used (fowler only)')
            comb = combine(logger, err)
            file_data_post_combine     = []
            file_hdr_post_combine      = []
            rates                      = []
            for idx_1, run in enumerate(self.session.file_data):
                file_data_post_combine.append([])
                file_hdr_post_combine.append([])
                rates.append([])
                for idx_2, dither in enumerate(run): 
                    logger.info("[run_pipe.go] Combining files for run:" + str(idx_1+params['minRunNum']) + ", dither:" + str(idx_2+params['minDithNum']))  
                    this_outPath_comb = params['workingDir'] + "comb_" + str(idx_1+params['minRunNum']) + "_" + str(idx_2+params['minDithNum']) + self.session.file_ext + ".fits"   
                    this_outPath_rate = params['workingDir'] + "comb_" + str(idx_1+params['minRunNum']) + "_" + str(idx_2+params['minDithNum']) + self.session.file_ext + ".rate.fits"  
                    rtn_data, rtn_hdr, rtn_rates = comb.execute(method=combination_method, datas=self.session.file_data[idx_1][idx_2], hdrs=self.session.file_hdr[idx_1][idx_2], \
                                                                hard=combination_hard, hard_rates=combination_hard_rates, out=this_outPath_comb, out_rate=this_outPath_rate, \
                                                                f_pairs=combination_f_pairs, opt_hdr=self.session.opt_hdr) 
                    if rtn_data is not None and rtn_hdr is not None and rtn_rates is not None: 
                        file_data_post_combine[idx_1].append(rtn_data)
                        file_hdr_post_combine[idx_1].append(rtn_hdr)
                        rates[idx_1].append(rtn_rates)
                    else:
                        err.set_code(6, is_critical=True)
            if combination_quit:
                logger.info("[run_pipe.go] Returning with code: " + str(err.current_code))
                return err.current_code
              
            self.session.set_session_vars_post_combine(file_data_post_combine, file_hdr_post_combine, rates)          # this sets file_[data||hdr]_[ss||nonss] and rates vars.   

            # ---------------------------------   
            # ---- nonlinearity correction ----   
            # ---------------------------------            
            if do_nonlincor:    
                logger.info("[run_pipe.go] Beginning nonlinearity correction.")
                self.session.add_amend_opt_header('L1LCOR', 1, 'linearity corrected')
                self.session.add_amend_opt_header('L1LCORF', os.path.basename(nonlincor_coeff_path), 'nonlinearity coeff file used')
                self.session.file_ext = self.session.file_ext + ".lcor"
                cor = nonlinearity_correction(logger, err) 
                cor.read_lcor_coeffs(nonlincor_coeff_path, nonlincor_order)
                for idx_1, run in enumerate(self.session.file_data_nonss):
                  for idx_2, f in enumerate(run):
                      logger.info("[run_pipe.go] Applying nonlinearity correction to run:" + str(idx_1+params['minRunNum']) + ", dither:" + str(idx_2+params['minDithNum']))  
                      this_outPath = params['workingDir'] + "comb_" + str(idx_1+params['minRunNum']) + "_" + str(idx_2+params['minDithNum']) + self.session.file_ext + ".fits"
                      rtn_data, rtn_hdr = cor.execute(data=self.session.file_data_nonss[idx_1][idx_2], hdr=self.session.file_hdr_nonss[idx_1][idx_2], \
                                                      rates=self.session.rates[idx_1][idx_2], hard=nonlincor_hard, out=this_outPath, opt_hdr=self.session.opt_hdr)
                      if rtn_data is not None and rtn_hdr is not None:
                          self.session.file_data_nonss[idx_1][idx_2] = rtn_data
                          self.session.file_hdr_nonss[idx_1][idx_2]  = rtn_hdr
                      else:
                          err.set_code(6, is_critical=True)       
                if nonlincor_quit:
                    logger.info("[run_pipe.go] Returning with code: " + str(err.current_code))
                    return err.current_code                       
            else:
                logger.info("[run_pipe.go] No linearity correction requested.") 
                self.session.add_amend_opt_header('L1LCOR', 0, 'linearity corrected')          
                      
            # -------------------------------   
            # -------- flatfielding ---------   
            # -------------------------------      
            if do_ff:
                logger.info("[run_pipe.go] Beginning flatfielding.")
                self.session.add_amend_opt_header('L1FLAT', 1, 'flatfielded')
                self.session.add_amend_opt_header('L1FLATF', os.path.basename(ff_coeff_path), 'flatfield file used')
                self.session.file_ext = self.session.file_ext + ".fcor"
                ff = flatfielding(logger, err)
                ff.read_fcor_coeffs(ff_coeff_path)
                for idx_1, run in enumerate(self.session.file_data_nonss):
                  for idx_2, f in enumerate(run):
                      logger.info("[run_pipe.go] Applying flatfield correction to run:" + str(idx_1+params['minRunNum']) + ", dither:" + str(idx_2+params['minDithNum']))  
                      this_outPath = params['workingDir'] + "comb_" + str(idx_1+params['minRunNum']) + "_" + str(idx_2+params['minDithNum']) + self.session.file_ext + ".fits"
                      rtn_data, rtn_hdr = ff.execute(data=self.session.file_data_nonss[idx_1][idx_2], hdr=self.session.file_hdr_nonss[idx_1][idx_2], hard=ff_hard, out=this_outPath, \
                                                        opt_hdr=self.session.opt_hdr) 
                      if rtn_data is not None and rtn_hdr is not None:
                          self.session.file_data_nonss[idx_1][idx_2] = rtn_data
                          self.session.file_hdr_nonss[idx_1][idx_2]  = rtn_hdr       
                      else:
                          err.set_code(6, is_critical=True)
                if ff_quit:
                    logger.info("[run_pipe.go] Returning with code: " + str(err.current_code))
                    return err.current_code    
            else:
                logger.info("[run_pipe.go] No flatfielding requested.") 
                self.session.add_amend_opt_header('L1FLAT', 0, 'flatfielded')           
            
            # -------------------------------------  
            # -------- bad pixel masking  ---------   
            # -------------------------------------
            if do_bp_masking:
                logger.info("[run_pipe.go] Masking science files using bad pixel masks.")
                self.session.add_amend_opt_header('L1BADPX', 1, 'bad pixel masks applied')
                self.session.add_amend_opt_header('L1BADPXF', os.path.basename(ff_badpix_path), 'bad pixel mask used (flat)')
                self.session.add_amend_opt_header('L1BADPXL', os.path.basename(nonlincor_badpix_path), 'bad pixel mask used (linearity)')            
                self.session.file_ext = self.session.file_ext + ".bp"
                bp = bad_pixel_mask(logger, err)
                if "FLAT" in bp_which:
                    bp.read_mask(nonlincor_badpix_path) # nonlinearity bad pixel mask
                    logger.info("[run_pipe.go] Using bad pixel mask: " + os.path.basename(nonlincor_badpix_path))
                if "LIN" in bp_which:
                    bp.read_mask(ff_badpix_path)        # flatfield bad pixel mask
                    logger.info("[run_pipe.go] Using bad pixel mask: " + os.path.basename(ff_badpix_path))                
                for idx_1, run in enumerate(self.session.file_data_nonss):
                    for idx_2, f in enumerate(run):         
                      logger.info("[run_pipe.go] Applying bad pixel mask to run:" + str(idx_1+params['minRunNum']) + ", dither:" + str(idx_2+params['minDithNum']))  
                      this_outPath = params['workingDir'] + "comb_" + str(idx_1+params['minRunNum']) + "_" + str(idx_2+params['minDithNum']) + self.session.file_ext + ".fits"
                      rtn_data, rtn_hdr = bp.execute(data=self.session.file_data_nonss[idx_1][idx_2], hdr=self.session.file_hdr_nonss[idx_1][idx_2], \
                                                      hard=bp_hard, out=this_outPath, opt_hdr=self.session.opt_hdr) 
                      if rtn_data is not None and rtn_hdr is not None:
                          self.session.file_data_nonss[idx_1][idx_2] = rtn_data    
                          self.session.file_hdr_nonss[idx_1][idx_2]  = rtn_hdr
                      else:
                          err.set_code(6, is_critical=True)
                if bp_quit:
                    logger.info("[run_pipe.go] Returning with code: " + str(err.current_code))
                    return err.current_code   
            else:
                logger.info("[run_pipe.go] No bad pixel masking requested.") 
                self.session.add_amend_opt_header('L1BADPX', 0, 'bad pixel mask applied')    
            self.session.set_opt_header_nonss()  


            if file_naming == "LT":     # this only makes sense in the context of LT files
                if write_lt_file and file_naming:
                    logger.info("[run_pipe.go] Writing IM_NONSS session data out...")
                    self.session.write_combined_data_as_LT(params['workingDir'], extname="IM_NONSS")              

                # -------------------------   
                # ---- sky subtraction ----   
                # -------------------------  
                if do_SS:
                    logger.info("[run_pipe.go] Beginning sky subtraction process.") 
                    self.session.add_amend_opt_header('L1SKYSUB', 1, 'sky subtracted')
                    for idx_1, run in enumerate(self.session.file_data_nonss): 
                        ## subtract sky from combined frames if more than one dither position exists
                        if len(run) == 1:
                            logger.info("[run_pipe.go] Only one dither position found. Data will not be sky subtracted.")
                            self.session.file_data_ss[idx_1][idx_2] = None
                            self.session.file_hdr_ss[idx_1][idx_2]  = None
                        else:
                            logger.info("[run_pipe.go] Run:" + str(idx_1+params['minRunNum']))
                            self.session.file_ext = self.session.file_ext + ".ss"
                            ss = sky_subtraction(logger, err)
                        
                            ## make sky
                            logger.info("[run_pipe.go] Making sky frame.")
                            this_outPath = params['workingDir'] + "comb_sky_" + str(idx_1+params['minRunNum']) + ".fits"
                            data_sky = ss.make_sky_frame(datas=self.session.file_data_nonss[idx_1], sigma=ss_bg_sigma_clip, hard=ss_hard, out=this_outPath)
                        
                            ## interpolate over bad pixels if bad pixel masking has been requested
                            if do_bp_masking:
                                ### interpolate over bad pixels
                                this_outPath = params['workingDir'] + "comb_sky_" + str(idx_1+params['minRunNum']) + ".interp.fits" 
                                ss.interpolate_bad(smoothing_box_size=ss_smoothing_box_size, hard=all([bp_hard, ss_hard]), out=this_outPath)
                                
                            for idx_2, f in enumerate(run):
                                logger.info("[run_pipe.go] Subtracting sky from run:" + str(idx_1+params['minRunNum']) + ", dither:" + str(idx_2+params['minDithNum']))  
                                this_outPath = params['workingDir'] + "comb_" + str(idx_1+params['minRunNum']) + "_" + str(idx_2+params['minDithNum']) + self.session.file_ext + ".fits"                   
                                rtn_data, rtn_hdr = ss.sub_sky(data=self.session.file_data_nonss[idx_1][idx_2], hdr=self.session.file_hdr_nonss[idx_1][idx_2], \
                                                              sigma=ss_bg_sigma_clip, hard=ss_hard, out=this_outPath, opt_hdr=self.session.opt_hdr)
                                if rtn_data is not None and rtn_hdr is not None:
                                    self.session.file_data_ss[idx_1][idx_2] = rtn_data
                                    self.session.file_hdr_ss[idx_1][idx_2]  = rtn_hdr 
                                else:
                                    err.set_code(6, is_critical=True)    
                        self.session.set_opt_header_ss()  
                        if ss_quit: 
                            logger.info("[run_pipe.go] Returning with code: " + str(err.current_code))
                            return err.current_code
                          
                    if write_lt_file:
                        logger.info("[run_pipe.go] Writing IM_SS session data out...")
                        self.session.write_combined_data_as_LT(params['workingDir'], extname="IM_SS", allow_append=True)

                    # -----------------------------------------  
                    # ---- SS image registration (LT only) ----   
                    # -----------------------------------------
                    if file_naming == "LT":
                        if do_registration:   
                            logger.info("[run_pipe.go] Beginning registration process.")
                            self.session.add_amend_opt_header('L1REG', 1, 'registered')
                            for idx_1, run in enumerate(self.session.file_data_ss): 
                                ## subtract sky from combined frames if more than one dither position exists
                                if len(run) == 1:
                                    logger.info("[run_pipe.go] Only one dither position found. Data will not be registered.")
                                    pass
                                else:              
                                    logger.info("[run_pipe.go] Run:" + str(idx_1+params['minRunNum']))
                                    self.session.file_ext = self.session.file_ext + ".reg"
                                    bp_combined_mask_data = None
                                    if do_bp_masking:
                                        bp_combined_mask_data = bp.combine_masks()
                                    reg = register(0, self.session.file_data_ss[idx_1], self.session.file_hdr_ss[idx_1], bp_combined_mask_data, logger, err) 
                                    this_outPaths = []
                                    this_mask_outPaths = []
                                    for idx_2, d in enumerate(run):
                                        this_outPaths.append(params['workingDir'] + "comb_" + str(idx_1+params['minRunNum']) + "_" + str(idx_2+params['minDithNum']) + self.session.file_ext + ".fits")
                                        this_mask_outPaths.append(params['workingDir'] + "comb_" + str(idx_1+params['minRunNum']) + "_" + str(idx_2+params['minDithNum']) + self.session.file_ext + ".mask.fits")
                                    rtn_datas, rtn_hdrs = reg.execute(registration_algorithm, params['workingDir'], fit_geom=registration_fit_geometry, bp_post_reg_thresh=bp_post_reg_thresh, \
                                                                      hard=registration_hard, outs=this_outPaths, mask_outs=this_mask_outPaths)                   
                                    if rtn_datas is not None and rtn_hdrs is not None:                    
                                        self.session.file_data_ss[idx_1] = rtn_datas                
                                        self.session.file_hdr_ss[idx_1]  = rtn_hdrs   
                                    else:
                                        err.set_code(6, is_critical=True)                 
                                if registration_quit: 
                                    logger.info("[run_pipe.go] Returning with code: " + str(err.current_code))
                                    return err.current_code

                            # ---------------------------  
                            # ---- SS image stacking ----   
                            # --------------------------- 
                            if do_stacking:
                                logger.info("[run_pipe.go] Beginning stacking process.")
                                self.session.add_amend_opt_header('L1STK', 1, 'stacked')
                                try:
                                    file_data_post_stacking = []
                                    file_hdr_post_stacking  = []
                                    for idx_1, run in enumerate(self.session.file_data_ss): 
                                        if len(run) == 1:
                                            logger.info("[run_pipe.go] Only one dither position found. Data will not be stacked.")
                                            pass
                                        else: 
                                            stk = stack(logger, err)
                                            this_outPath = params['workingDir'] + "comb_stk_" + str(idx_1+params['minRunNum']) + ".fits"                    
                                            rtn_data, rtn_hdr = stk.execute(datas=self.session.file_data_ss[idx_1], hdrs=self.session.file_hdr_ss[idx_1], \
                                                                            hard=stacking_hard, method=stacking_method, out=this_outPath, opt_hdr=self.session.opt_hdr)  
                                            if rtn_data is not None and rtn_hdr is not None:
                                                file_data_post_stacking.append(rtn_data)
                                                file_hdr_post_stacking.append(rtn_hdr)
                                            else:
                                              err.set_code(6, is_critical=True)             
                                    self.session.set_opt_header_ss_stk()    
                                except RuntimeError:
                                    err.set_code(13, is_critical=True)
                                self.session.set_session_vars_post_stacking(file_data_post_stacking, file_hdr_post_stacking)          # this sets file_[data||hdr]_[nonss]_stk vars.  
                                
                                # --------------------------  
                                # ---- write STK output ----   
                                # -------------------------- 
                                if write_lt_file:
                                    logger.info("[run_pipe.go] Writing SK_SS session data out...")
                                    self.session.write_stacked_data_as_LT(params['workingDir'], extname="SK_SS")
                        else:
                            logger.info("[run_pipe.go] No stacking requested or the absence of a previous operation mean it cannot be performed.")  
                            self.session.add_amend_opt_header('L1STK', 0, 'stacked')  
                    else:
                        logger.info("[run_pipe.go] No registration requested.")  
                        self.session.add_amend_opt_header('L1REG', 0, 'registered')         
                else:
                    logger.info("[run_pipe.go] No sky subtraction requested.")     
                    self.session.add_amend_opt_header('L1SKYSUB', 0, 'sky subtracted') 
        else:
            logger.info("[run_pipe.go] No frame combination requested.") 
            self.session.add_amend_opt_header('L1FMCO', 0, 'frame combined')
             
        logger.info("[run_pipe.go] Returning with code: " + str(err.current_code))
        return err.current_code
        
if __name__ == "__main__":
    parser = optparse.OptionParser()
    group1 = optparse.OptionGroup(parser, "General")
    group1.add_option('--p', action='store', default='../test/', dest='dataPath', type=str, help='path to data')
    group1.add_option('--wd', action='store', default='test', dest='workingDir', type=str, help='path to working directory')
    group1.add_option('--o', action='store_true', dest='clobber', help='clobber working directory?')
    group1.add_option('--pa', action='store', default='../config/paths_rmb.ini', type=str, dest='pathsCfgPath', help='path to paths config file')
    group1.add_option('--pi', action='store', default='../config/pipeline_rmb.ini', type=str, dest='pipeCfgPath', help='path to pipeline config file')
    group1.add_option('--rlo', action='store', dest='minRunNum', default=1, type=int, help='minimum multrun number')
    group1.add_option('--rhi', action='store', dest='maxRunNum', default=1, type=int, help='maximum multrun number')          
    group1.add_option('--elo', action='store', default=1, type=int, dest='minExpNum', help='lowest exposure number to use')
    group1.add_option('--ehi', action='store', default=4, type=int, dest='maxExpNum', help='highest exposure number to use')    
    group1.add_option('--log', action='store', default='INFO', dest='logLevel', type=str, help='log level (DEBUG|INFO|WARNING|ERROR|CRITICAL)')
    parser.add_option_group(group1)
    
    group2 = optparse.OptionGroup(parser, "LT specific general parameters")    
    group2.add_option('--d', action='store', dest='date', default='*', type=str, help='date (YYYYMMDD)')
    group2.add_option('--dlo', action='store', dest='minDithNum', default=1, type=int, help='minimum dither number')
    group2.add_option('--dhi', action='store', dest='maxDithNum', default=9, type=int, help='maximum dither number')  
    parser.add_option_group(group2)

    group3 = optparse.OptionGroup(parser, "Teledyne specific general parameters")    
    group3.add_option('--glo', action='store', default=0, type=int, dest='minGrpNum', help='lowest group number to use')
    group3.add_option('--ghi', action='store', default=1, type=int, dest='maxGrpNum', help='highest group number to use')
    parser.add_option_group(group3)
    
    options, args = parser.parse_args()
    params = {
        'dataPath' : str(options.dataPath).rstrip("/") + "/",
        'workingDir' : str(options.workingDir).rstrip("/") + "/",
        'clobber' : bool(options.clobber),       
        'pathsCfgPath' : str(options.pathsCfgPath),
        'pipeCfgPath' : str(options.pipeCfgPath),
        'minRunNum' : int(options.minRunNum),
        'maxRunNum' : int(options.maxRunNum),
        'date' : str(options.date),
        'minDithNum' : int(options.minDithNum),         
        'maxDithNum' : int(options.maxDithNum), 
        'minExpNum' : int(options.minExpNum),
        'maxExpNum' : int(options.maxExpNum),
        'minGrpNum' : int(options.minGrpNum),
        'maxGrpNum' : int(options.maxGrpNum),
        'logLevel' : str(options.logLevel.upper()),
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

    # file loggingReduction
    fh = logging.FileHandler(params['workingDir'] + "res.log")
    fh.setLevel(logging.DEBUG)

    ## set logging format
    fh.setFormatter(formatter)

    ## add handlers to logging object
    logger.addHandler(fh)
    
    pipe = run_pipe(params, logger, err)
    rtn = pipe.go()        
    exit(rtn)      
