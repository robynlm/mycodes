"""
Author: Robyn Munoz
Date: 2nd July 2020

This calculates the Richardson expansion errors as a function of time of the files:
 - h5_data.csv -> h5_data_error.csv
 - asc_data.csv -> asc_data_error.csv

From the param.py file the following variables are needed:
 sim_name, dx_ref, dtfac
 
To use:
python RError.py simulation1 simulation2 simulation3 location

The 3 simulations are the same but with different grid size.
The grid spacing MUST be doubled every simulation.
(The order in which they are listed does not matter)
The error is calculated for the simulation of the smallest grid points 
but with the iterations of the simulation with the biggest grid points.
"""
import re
import sys
import numpy as np
import pandas as pd
from tools import ChooseParam
from tools import ReadingTools as RRead

def order(f, i):
    return f[i[0]], f[i[1]], f[i[2]]

if __name__ == "__main__":
    
    print("\n#################################")
    print("\n        Richardson Error         \n")
    print("#################################\n")
    
    simA, pathA, OGpathA = ChooseParam.ChooseP(sys.argv[2], sys.argv[1])
    simB, pathB, OGpathB = ChooseParam.ChooseP(sys.argv[3], sys.argv[1])
    simC, pathC, OGpathC = ChooseParam.ChooseP(sys.argv[4], sys.argv[1])
    sims = [simA, simB, simC]
    paths = [pathA, pathB, pathC]
    
    # a lot of change is needed here
    dts = [sim.dx*sim.dtfac for sim in sims]
    idxs = np.argsort(dts)[::-1]

    dt1, dt2, dt4 = order(dts, idxs)
    if dt4*2!=dt2 or dt4*4!=dt1:
        print('ERROR: dts are not x2 from each other')
        sys.exit()
    
    sim1, sim2, sim4 = order(sims, idxs)
    path1, path2, path4 = order(paths, idxs)
    
    # =====================================================
    Files = re.split('\n|/|.csv', RRead.BASH('ls '+path4+'*.csv'))
    FileName = [f for f in Files if f in ['h5_data', 'asc_data', 'TA_radius']]
    print(FileName)
    for filename in FileName:
        print('\n ===== Calc Errors for file:',filename+'.csv')
        dummy = RRead.BASH('rm '+path4+filename+'_error.csv')
        dummy = RRead.BASH('rm '+path4+filename+'_conv.csv')

        file1 = pd.read_csv(path1+filename+'.csv', delimiter=',')
        cols = list(file1.columns)
        file2 = pd.read_csv(path2+filename+'.csv', delimiter=',')
        file4 = pd.read_csv(path4+filename+'.csv', delimiter=',')
        t1, t2, t4 = np.array(file1['t']), np.array(file2['t']), np.array(file4['t'])

        if 'h5' in filename or 'TA' in filename:
            filedatconv = [0.0, 0.0]
            filedaterr = [file2['it'], file2['t']]
            colnames = ['it','t']
            ini_iteration = 2
        elif 'asc' in filename:
            filedatconv = [0.0]
            filedaterr = [file2['t']]
            colnames = ['t']
            ini_iteration = 1
        
        for i in range(ini_iteration, len(cols)):
            f1, f2, f4 = np.array(file1[cols[i]]), np.array(file2[cols[i]]), np.array(file4[cols[i]])

            len1 = np.min([len(f1), int(len(f2)/2), int(len(f4)/4)])
            f1_len1 = f1[:len1]
            f2_len1 = f2[:2*len1]
            f4_len1 = f4[:4*len1]

            len2 = np.min([len(f2), int(len(f4)/2)])
            f2_len2 = f2[:len2]
            f4_len2 = f4[:2*len2]
            
            # Convergence
            c_len1 = RRead.safe_division((f1_len1-f2_len1[0::2]), (f2_len1[0::2]-f4_len1[0::4]))
            cmean = abs(np.mean(c_len1))

            c_len2 = []
            i1 = 0
            i2 = 0
            while i2<len2:
                c_len2 += [c_len1[i1]]
                i2 += 1
                if i1<len1-1:
                    if t2[i2]>t1[i1+1]:
                    	i1 += 1

            print(len(c_len1), len(c_len2))                    
            # Calc Error
            print(len(f2_len2), len(f4_len2[0::2]), len(c_len2))
            RError = abs(f2_len2-f4_len2[0::2])/(np.array(c_len2)-1)
            print(cols[i], ': Convergence = ', cmean, ' Error = ', np.mean(RError))

            # Record results
            filedatconv += [cmean]
            filedaterr += [RError]
            colnames += [cols[i]]
           
        colnameserr = [col+'_error' for col in colnames]
        lenmin = np.min([len(i) for i in filedaterr])
        data = np.array([i[:lenmin] for i in filedaterr]).T
        pd.DataFrame(data, columns=colnameserr).to_csv(path4+filename+'_error.csv')
        
        filedat = [colnames, filedatconv]
        pd.DataFrame(np.array(filedat).T, columns=['Var', 'Convergence']).to_csv(path4+filename+'_conv.csv')
    print('\nDone!')
                                
                                
