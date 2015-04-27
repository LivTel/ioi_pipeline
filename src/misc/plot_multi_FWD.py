import os
import pylab as plt
from measure_FWD import measure_FWD as mfwd
from measure_linearity import measure_linearity as mlin

paths = ["/mnt/NAS/devel/IOI/images_and_analysis/remote_5/images/22/",
         "/mnt/NAS/devel/IOI/images_and_analysis/remote_5/images/14/",
         "/mnt/NAS/devel/IOI/images_and_analysis/remote_5/images/15/",
         "/mnt/NAS/devel/IOI/images_and_analysis/remote_5/images/21/"
        ]

ghis = [31,40,32,29]
#ghis = [5]
for idx, p in enumerate(paths):
    # parameters are shared between measure_FWD and measure_linearity.
    params = {                  
        'dataPath' : p,
        'workingDir' : "tmp",
        'clobber' : True,       
        'pathsCfgPath' : "../../config/paths_rmb.ini",
        'minRunNum' : 1,
        'maxRunNum' : 1,
        'minDithNum' : 0,
        'maxDithNum' : 0,
        'minGrpNum' : 0,
        'maxGrpNum' : ghis[idx],
        'minExpNum' : 1,
        'maxExpNum' : 1,        
        'logLevel' : "INFO",
        'doLinearity' : False,
        'p1' : True,
        'p2' : False,
        'hard' : False,
        'flip' : False,
        'section' : '287,1600,703,1803',
        'linLine' : False
    }
    
    # make exptime v counts plot
    plt.subplot(4,3,(idx*3)+1)
    m = mfwd(params)
    m.go()  
    plt.title('')
    if idx < len(paths)-1:
        plt.xlabel('')
    if idx == 0 or idx < len(paths)-1:
        plt.xlabel('')
        
    # make exptime v counts plot (CDS)
    plt.subplot(4,3,(idx*3)+2)
    m = mlin(params)
    m.go()
    plt.title('')  
    if idx < len(paths)-1:
        plt.xlabel('')    
    if idx == 0 or idx < len(paths)-1:
        plt.xlabel('')
        
    # make exptime v linearity plot   
    params['p1'] = False
    params['p2'] = True
    plt.subplot(4,3,(idx*3)+3)
    m = mlin(params)
    m.go()    
    plt.title('')  
    if idx == 0 or idx < len(paths)-1:
        plt.xlabel('')   
            
plt.show()