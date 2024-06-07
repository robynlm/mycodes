import sys
import numpy as np
import pandas as pd
from tools import LinData
from tools import ReadingTools as RRead

if __name__ == "__main__":
    
    simname = sys.argv[1]
    param = RRead.read_parameters(simname)
    Lin = LinData.LinData_Class(param)
    
    all_it, all_t = RRead.collect_iteration_and_time(param['datapath'])
    all_a = Lin.evo.a(all_t)
    all_an = Lin.an_initial(all_t)
    all_H = Lin.evo.Hprop(all_t)
    all_z = Lin.evo.z(all_t)
    
    temporal_step_variation = [abs(np.sum(np.diff(np.diff(temp_var)))) 
                               for temp_var in [all_t[:-2], all_a[:-2], all_H[:-2], all_z[:-2]]]
    temporal_step = ['dt', 'da', 'dH', 'dz'][np.argmin(temporal_step_variation)]
    print(temporal_step, 'constant in this simulation', flush=True)
    
    filedat = np.array([all_it, all_t, all_a, all_an, all_H, all_z]).T
    pd.DataFrame(filedat).to_csv(param['datapath']+'Time_'+temporal_step+'.csv', 
                                 header=['it', 't', 'a', 'an', 'H', 'z'], index=False)
