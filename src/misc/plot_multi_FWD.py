import os
import pylab as plt
from measure_FWD import measure_FWD as mfwd

paths = ["/mnt/NAS/devel/IOI/images_and_analysis/remote_5/images/13/",
         "/mnt/NAS/devel/IOI/images_and_analysis/remote_5/images/14/",
         "/mnt/NAS/devel/IOI/images_and_analysis/remote_5/images/15/"
        ]

for idx, p in enumerate(paths):
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
        'maxGrpNum' : 50,
        'minExpNum' : 1,
        'maxExpNum' : 1,        
        'logLevel' : "INFO",
        'fits' : False,
        'flip' : False
    }
    plt.subplot(3,1,idx+1)
    m = mfwd(params)
    m.go()
    
    # amend plots as desired
    ## remove x labels for all apart from last plot
    if idx < len(paths)-1:
        plt.xlabel("")
    
plt.show()