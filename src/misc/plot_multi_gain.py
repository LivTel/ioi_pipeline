import os
import pylab as plt
from measure_gain_RN import measure_gain_RN as mgrn

paths = ["/mnt/NAS/devel/IOI/images_and_analysis/remote_5/images/9/",
         "/mnt/NAS/devel/IOI/images_and_analysis/remote_5/images/10/",
         "/mnt/NAS/devel/IOI/images_and_analysis/remote_5/images/11/",
         "/mnt/NAS/devel/IOI/images_and_analysis/remote_5/images/12/",
        ]

for idx, p in enumerate(paths):
    params = {
        'dataPath' : p,
        'workingDir' : "tmp",
        'clobber' : True,       
        'pathsCfgPath' : "../../config/paths_rmb.ini",
        'minRunNum' : 1,
        'maxRunNum' : 2,
        'minDithNum' : 0,
        'maxDithNum' : 0,
        'minGrpNum' : 0,
        'maxGrpNum' : 20,
        'minExpNum' : 1,
        'maxExpNum' : 16,        
        'logLevel' : "INFO",
        'fit' : [0,10000],
        'fitCoeff' : 1,
        'windowSize' : 50,
        'fits' : False,
        'method' : "UTR",
        'calcGain': True,
        'calcRead': False,
        'flip' : False
    }
    plt.subplot(2,2,idx+1)
    m = mgrn(params)
    m.go()
    
plt.show()