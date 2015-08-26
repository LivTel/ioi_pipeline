from session import *
from remove_reference import *
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
            write_output                        = bool(int(pipe_cfg['general']['write_output']))  
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
            combination_quit                    = bool(int(pipe_cfg['frame_combination']['quit'])) 
            combination_fowler_bodge            = float(pipe_cfg['frame_combination']['fowler_bodge'])
            ## nonlinearity_corrections
            do_nonlincor                        = bool(int(pipe_cfg['nonlinearity_corrections']['do']))       
            nonlincor_order                     = int(pipe_cfg['nonlinearity_corrections']['order'])
            nonlincor_coeff_path                = str(full_path_nonlinearity_corrections + paths_cfg['nonlinearity_corrections']['coeff_filename'])
            nonlincor_badpix_path               = str(full_path_nonlinearity_corrections + paths_cfg['nonlinearity_corrections']['bad_pix_filename'])        
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
            ss_peers_only                       = bool(int(pipe_cfg['sky_subtraction']['peers_only']))
            ss_make_algorithm                   = str(pipe_cfg['sky_subtraction']['make_algorithm'])
            ss_max_robust_iterations            = int(pipe_cfg['sky_subtraction']['max_robust_iterations'])
            ss_max_robust_unsolved_percentage   = float(pipe_cfg['sky_subtraction']['max_robust_unsolved_percentage'])
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
        except KeyError, e:
            logger.info("[run_pipe.go] Key/section " + str(e) + " appears to be missing.")
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
        logger.info("[run_pipe.go] Searching for LT files in " + params['dataPath']) 
        try:
            files, n, n_runs, n_dithers, n_exp = find_sort_files_LT(params['dataPath'], params['minRunNum'], params['maxRunNum'], params['minDithNum'], params['maxDithNum'], params['minExpNum'], params['maxExpNum'], params['date'], logger, err)
        except KeyError:        # this may happen if relevant LT input params haven't been set when instantiating, e.g. date.
            n = 0
        if n == 0:
            logger.info("[run_pipe.go] No LT files found in " + params['dataPath'] + "or relevant LT parameters haven't been set.") 
            logger.info("[run_pipe.go] Searching for Teledyne files in " + params['dataPath'])
            try:
                files, n, n_runs, n_dithers, n_exp = find_sort_files_teledyne(params['dataPath'], params['minRunNum'], params['maxRunNum'], params['minGrpNum'], params['maxGrpNum'], params['minExpNum'], params['maxExpNum'], logger, err)
            except KeyError:    # this may happen if relevant Teledyne input params haven't been set when instantiating.
                n = 0
            if n == 0:
                logger.info("[run_pipe.go] No Teledyne files found in " + params['dataPath'] + "or relevant Teledyne parameters haven't been set.") 
                err.set_code(17, is_critical=True)
            else:
                logger.info("[run_pipe.go] Assuming Teledyne nomenclature.")
                are_LT_files = False
        else:
            logger.info("[run_pipe.go] Assuming LT nomenclature.")  
            are_LT_files = True
            
        self.session.start(files)

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
                        rtn_data, rtn_hdr = sub.execute(method=refsub_method, in_data=self.session.file_data[idx_1][idx_2][idx_3], in_hdr=self.session.file_hdr[idx_1][idx_2][idx_3], \
                                                        hard=refsub_hard, strip=True, out=this_outPath, opt_hdr=self.session.opt_hdr) 
                        if rtn_data is not None and rtn_hdr is not None:
                            self.session.file_data[idx_1][idx_2][idx_3]      = rtn_data
                            self.session.file_hdr[idx_1][idx_2][idx_3]       = rtn_hdr
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
            
            # ----------------------------------------   
            # ---- set up nonlinearity correction ----   
            # ----------------------------------------  
            #
            # this needs to be done on a per CDS frame basis 
            # if fowler mode is requested as the linearity
            # correction needs to be separately applied to
            # the different pairs (each pair has a different
            # "zero" value corresponding to a multiple of the
            # readout time). It therefore has to be done after
            # each individual pair is made, i.e. in the 
            # combine routine.
            #
            # painful, I know.
            #
            if do_nonlincor:
                logger.info("[run_pipe.go] Nonlinearity correction requested.")
                self.session.add_amend_opt_header('L1LCOR', 1, 'linearity corrected')
                self.session.add_amend_opt_header('L1LCORF', os.path.basename(nonlincor_coeff_path), 'nonlinearity coeff file used')
                self.session.file_ext = self.session.file_ext + ".lcor"
                lcor = nonlinearity_correction(logger, err) 
                lcor.read_lcor_coeffs(nonlincor_coeff_path, nonlincor_order, flip=params['flip'])
            else:
                logger.info("[run_pipe.go] No linearity correction requested.") 
                self.session.add_amend_opt_header('L1LCOR', 0, 'linearity corrected')  
                lcor = None
                
            for idx_1, run in enumerate(self.session.file_data):
                for idx_2, dither in enumerate(run): 
                    logger.info("[run_pipe.go] Combining files for run:" + str(idx_1+params['minRunNum']) + ", dither:" + str(idx_2+params['minDithNum']))  
                    this_outPath_comb = params['workingDir'] + "comb_" + str(idx_1+params['minRunNum']) + "_" + str(idx_2+params['minDithNum']) + self.session.file_ext + ".fits"   
                    rtn_data, rtn_hdr, rtn_rates = comb.execute(method=combination_method, in_datas=self.session.file_data[idx_1][idx_2], in_hdrs=self.session.file_hdr[idx_1][idx_2], \
                                                                hard=combination_hard, lcor=lcor, out=this_outPath_comb, f_pairs=combination_f_pairs, opt_hdr=self.session.opt_hdr, \
                                                                fowler_bodge = combination_fowler_bodge) 
                    if rtn_data is not None and rtn_hdr is not None and rtn_rates is not None: 
                        self.session.file_data_nonss[idx_1][idx_2]      = rtn_data
                        self.session.file_hdr_nonss[idx_1][idx_2]       = rtn_hdr
                        self.session.rates[idx_1][idx_2]                = rtn_rates
                    else:
                        err.set_code(6, is_critical=True)        
            
            if combination_quit:
                logger.info("[run_pipe.go] Returning with code: " + str(err.current_code))
                return err.current_code   
              
            # -------------------------------   
            # -------- flatfielding ---------   
            # -------------------------------      
            if do_ff:
                logger.info("[run_pipe.go] Beginning flatfielding.")
                self.session.add_amend_opt_header('L1FLAT', 1, 'flatfielded')
                self.session.add_amend_opt_header('L1FLATF', os.path.basename(ff_coeff_path), 'flatfield file used')
                self.session.file_ext = self.session.file_ext + ".fcor"
                ff = flatfielding(logger, err)
                ff.read_fcor_coeffs(ff_coeff_path, flip=params['flip'])
                for idx_1, run in enumerate(self.session.file_data_nonss):
                  for idx_2, f in enumerate(run):
                      logger.info("[run_pipe.go] Applying flatfield correction to run:" + str(idx_1+params['minRunNum']) + ", dither:" + str(idx_2+params['minDithNum']))  
                      this_outPath = params['workingDir'] + "comb_" + str(idx_1+params['minRunNum']) + "_" + str(idx_2+params['minDithNum']) + self.session.file_ext + ".fits"
                      rtn_data, rtn_hdr = ff.execute(in_data=self.session.file_data_nonss[idx_1][idx_2], in_hdr=self.session.file_hdr_nonss[idx_1][idx_2], hard=ff_hard, out=this_outPath, \
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
                if "LIN" in bp_which:
                    bp.read_mask(nonlincor_badpix_path, flip=params['flip']) # nonlinearity bad pixel mask
                    logger.info("[run_pipe.go] Using bad pixel mask: " + os.path.basename(nonlincor_badpix_path))
                if "FLAT" in bp_which:
                    bp.read_mask(ff_badpix_path, flip=params['flip'])        # flatfield bad pixel mask
                    logger.info("[run_pipe.go] Using bad pixel mask: " + os.path.basename(ff_badpix_path))                
                for idx_1, run in enumerate(self.session.file_data_nonss):
                    for idx_2, f in enumerate(run):         
                      logger.info("[run_pipe.go] Applying bad pixel mask to run:" + str(idx_1+params['minRunNum']) + ", dither:" + str(idx_2+params['minDithNum']))  
                      this_outPath = params['workingDir'] + "comb_" + str(idx_1+params['minRunNum']) + "_" + str(idx_2+params['minDithNum']) + self.session.file_ext + ".fits"
                      rtn_data, rtn_hdr = bp.execute(in_data=self.session.file_data_nonss[idx_1][idx_2], in_hdr=self.session.file_hdr_nonss[idx_1][idx_2], \
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
            self.session.file_opt_hdr_nonss.append(collections.OrderedDict(self.session.opt_hdr))

            if are_LT_files:     # these processes only make sense in the context of LT files. 
              
                self.session.copy_data_hdr_nonss_to_ss()        # make a copy of nonss session vars to ss.
                
                data    = self.session.file_data_nonss          # pointer to allow for option of ss or nonss data in registration/stacking.
                hdr     = self.session.file_hdr_nonss           # pointer to allow for option of ss or nonss data in registration/stacking. 
                
                # -----------------------------------   
                # ---- sky subtraction (LT only) ----   
                # -----------------------------------
                if do_SS:
                    logger.info("[run_pipe.go] Beginning sky subtraction process.") 
                    for idx_1, run in enumerate(self.session.file_data_nonss): 
                        ## omit sky subtraction if only one dither position exists.
                        if len(run) == 1:
                            logger.info("[run_pipe.go] Only one dither position found. Data will not be sky subtracted.")
                            self.session.add_amend_opt_header('L1SKYSUB', 0, 'sky subtracted')
                            self.session.file_data_ss[idx_1][idx_2] = None
                            self.session.file_hdr_ss[idx_1][idx_2]  = None
                        else:  
                            logger.info("[run_pipe.go] Run:" + str(idx_1+params['minRunNum']))
                            self.session.add_amend_opt_header('L1SKYSUB', 1, 'sky subtracted')
                            self.session.file_ext = self.session.file_ext + ".ss"  
                            for idx_2, f in enumerate(run):
                                logger.info("[run_pipe.go] Subtracting sky from run:" + str(idx_1+params['minRunNum']) + ", dither:" + str(idx_2+params['minDithNum']))   
                                
                                ss = sky_subtraction(logger, err)
                                logger.info("[run_pipe.go] Constructing sky frame.")
                                this_outPath = params['workingDir'] + "comb_sky_" + str(idx_1+params['minRunNum']) + '_' + str(idx_2+params['minDithNum']) + ".fits"
                                data_sky = ss.make_sky_frame(in_datas=self.session.file_data_nonss[idx_1], idx_of_current_frame=idx_2, sigma=ss_bg_sigma_clip, interpolate_bad=do_bp_masking, \
                                                             smoothing_box_size=ss_smoothing_box_size, make_algorithm=ss_make_algorithm, max_iter=ss_max_robust_iterations, peers_only=ss_peers_only, \
                                                             max_robust_unsolved_percentage=ss_max_robust_unsolved_percentage, hard=ss_hard, out=this_outPath)                         
                                logger.info("[run_pipe.go] Subtracting sky frame.")      
                                ### subtract sky
                                this_outPath = params['workingDir'] + "comb_sky_sub_" + str(idx_1+params['minRunNum']) + '_' + str(idx_2+params['minDithNum']) + ".fits"
                                rtn_data, rtn_hdr = ss.sub_sky(in_data=self.session.file_data_nonss[idx_1][idx_2], in_hdr=self.session.file_hdr_nonss[idx_1][idx_2], \
                                                               sigma=ss_bg_sigma_clip, hard=ss_hard, out=this_outPath, opt_hdr=self.session.opt_hdr)                     
                                if rtn_data is not None and rtn_hdr is not None:
                                    self.session.file_data_ss[idx_1][idx_2] = rtn_data
                                    self.session.file_hdr_ss[idx_1][idx_2]  = rtn_hdr 
                                else:
                                    err.set_code(6, is_critical=True)    
                        self.session.file_opt_hdr_ss.append(collections.OrderedDict(self.session.opt_hdr))  
                        if ss_quit: 
                            logger.info("[run_pipe.go] Returning with code: " + str(err.current_code))
                            return err.current_code        
                    data = self.session.file_data_ss            # ss has been applied, redirect pointer  
                    hdr  = self.session.file_hdr_ss             # ss has been applied, redirect pointer
                else:
                    logger.info("[run_pipe.go] No sky subtraction requested.")     
                    self.session.add_amend_opt_header('L1SKYSUB', 0, 'sky subtracted') 
                    
                # -----------------------------------------  
                # ---- SS image registration (LT only) ----   
                # -----------------------------------------
                if do_registration:   
                    logger.info("[run_pipe.go] Beginning registration process.")
                    for idx_1, run in enumerate(data): 
                        ## subtract sky from combined frames if more than one dither position exists
                        if len(run) == 1:
                            logger.info("[run_pipe.go] Only one dither position found. Data will not be registered.")
                            self.session.add_amend_opt_header('L1REG', 0, 'registered')
                        else:              
                            logger.info("[run_pipe.go] Run:" + str(idx_1+params['minRunNum']))
                            self.session.add_amend_opt_header('L1REG', 1, 'registered')
                            self.session.file_ext = self.session.file_ext + ".reg"
                            bp_combined_mask_data = None
                            if do_bp_masking:
                                bp_combined_mask_data = bp.combine_masks()
                            reg = register(0, data[idx_1], hdr[idx_1], bp_combined_mask_data, logger, err) 
                            this_outPaths = []
                            this_mask_outPaths = []
                            for idx_2, d in enumerate(run):
                                this_outPaths.append(params['workingDir'] + "comb_reg_" + str(idx_1+params['minRunNum']) + "_" + str(idx_2+params['minDithNum']) + self.session.file_ext + ".fits")
                                this_mask_outPaths.append(params['workingDir'] + "comb_reg_" + str(idx_1+params['minRunNum']) + "_" + str(idx_2+params['minDithNum']) + self.session.file_ext + ".mask.fits")
                            rtn_datas, rtn_hdrs = reg.execute(registration_algorithm, params['workingDir'], fit_geom=registration_fit_geometry, bp_post_reg_thresh=bp_post_reg_thresh, \
                                                              hard=registration_hard, outs=this_outPaths, mask_outs=this_mask_outPaths)                   
                            if rtn_datas is not None and rtn_hdrs is not None:
                                    self.session.file_data_reg[idx_1] = rtn_datas            
                                    self.session.file_hdr_reg[idx_1]  = rtn_hdrs                          
                            else:
                                err.set_code(6, is_critical=True)                 
                        if registration_quit: 
                            logger.info("[run_pipe.go] Returning with code: " + str(err.current_code))
                            return err.current_code    
                          
                    # -------------------------------------  
                    # ---- SS image stacking (LT only) ----   
                    # ------------------------------------- 
                    if do_stacking:
                        logger.info("[run_pipe.go] Beginning stacking process.")
                        try:
                            for idx_1, run in enumerate(data): 
                                if len(run) == 1:
                                    logger.info("[run_pipe.go] Only one dither position found. Data will not be stacked.")
                                    self.session.add_amend_opt_header('L1STK', 0, 'stacked')
                                else: 
                                    self.session.add_amend_opt_header('L1STK', 1, 'stacked')
                                    stk = stack(logger, err)
                                    this_outPath = params['workingDir'] + "comb_stk_" + str(idx_1+params['minRunNum']) + ".fits"                    
                                    rtn_data, rtn_hdr = stk.execute(in_datas=self.session.file_data_reg[idx_1], in_hdrs=self.session.file_hdr_reg[idx_1], \
                                                                    hard=stacking_hard, method=stacking_method, out=this_outPath, opt_hdr=self.session.opt_hdr)                            
                                    if rtn_data is not None and rtn_hdr is not None:
                                        self.session.file_data_stk[idx_1] = rtn_data
                                        self.session.file_hdr_stk[idx_1]  = rtn_hdr
                                    else:
                                        err.set_code(6, is_critical=True)             
                            self.session.file_opt_hdr_stk.append(collections.OrderedDict(self.session.opt_hdr))  
                        except RuntimeError:
                            err.set_code(13, is_critical=True)   

                    else:
                        logger.info("[run_pipe.go] No stacking requested.")  
                        self.session.add_amend_opt_header('L1STK', 0, 'stacked')  
                else:
                    logger.info("[run_pipe.go] No registration requested.")  
                    self.session.add_amend_opt_header('L1REG', 0, 'registered') 
            else:
                logger.info("[run_pipe.go] File naming does not conform to LT nomenclature, processing will not proceed any further.")  
        else:
            logger.info("[run_pipe.go] No frame combination requested.") 
            self.session.add_amend_opt_header('L1FMCO', 0, 'frame combined')
           
        # ----------------------  
        # ---- write output ----   
        # ----------------------  
        # 
        # Each dither position:
        # a) do_SS set   = IM_SS, IM_NONSS
        # b) do_SS unset = IM_NONSS
        #
        # If stacking successful:
        # a) do_SS set   = SK_SS
        # b) do_SS unset = SK_NONSS
        #
        if write_output and are_LT_files:
            if do_SS:
                logger.info("[run_pipe.go] Writing IM_SS extension session data out...")
                n = self.session.write_combined_data_as_LT(params['workingDir'], extname="IM_SS")
                logger.info("[run_pipe.go] " + str(n) + " files written or appended to.")
            logger.info("[run_pipe.go] Writing IM_NONSS extension session data out...")
            n = self.session.write_combined_data_as_LT(params['workingDir'], extname="IM_NONSS", allow_append=True)
            logger.info("[run_pipe.go] " + str(n) + " files written.")            
            if do_SS:
                logger.info("[run_pipe.go] Writing SK_SS extension session data out...")
                n = self.session.write_stacked_data_as_LT(params['workingDir'], extname="SK_SS")
                logger.info("[run_pipe.go] " + str(n) + " files written or appended to.")
            else:
                logger.info("[run_pipe.go] Writing SK_NONSS extension session data out...")
                n = self.session.write_stacked_data_as_LT(params['workingDir'], extname="SK_NONSS")
                logger.info("[run_pipe.go] " + str(n) + " files written.")
                
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
    group1.add_option('--f', action='store_true', dest='flip', help='flip calibration data about x axis?')
    parser.add_option_group(group1)
    
    group2 = optparse.OptionGroup(parser, "LT specific general parameters")    
    group2.add_option('--d', action='store', dest='date', default='*', type=str, help='date (YYYYMMDD)')
    group2.add_option('--dlo', action='store', dest='minDithNum', default=1, type=int, help='minimum dither number')
    group2.add_option('--dhi', action='store', dest='maxDithNum', default=999, type=int, help='maximum dither number')  
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
