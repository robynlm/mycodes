"""
Author: Robyn Munoz
Date: 9th June 2020

This Class:
 - links the program to the appropriate parameter class
 - links the path to the relevant location
 - collects all iteration steps
 
 From the param.py file the following variables are needed:
 sim_name, path_H and/or path_S
"""

import inspect
import re
from . import ReadingTools as RRead
import param
import sys
import numpy as np
from . import LinData

def ChooseP(simname, HorS):
        
    # identify the simulation
    method_to_call = getattr(param, simname)
    sim = method_to_call()
    print('simulation: '+sim.sim_name)
    
    # identify the path 
    if HorS=='home':
        HorSloc = '~/simulations/'
    elif HorS=='sciama':
        HorSloc = '/mnt/lustre2/ET_sims/'
    else:
        print("HorS should be 'home' or 'sciama'")
        sys.exit()
    path = HorSloc+sim.sim_name+'/output-0000/'+sim.sim_name+'/'
    
    return sim, path

def ChooseIt(var_wanted_str, var_wanted_float_arg, Lin, sim, OGpath):
    var_wanted_float = [float(var) for var in var_wanted_float_arg.split(',')]
    if var_wanted_str=='it':
        it_wanted = [int(it) for it in var_wanted_float]
    elif var_wanted_str in ['a', 'an', 'H', 't', 'z']:
        it_wanted = [int(Lin.temp_from_temp('it', var_wanted_str, var)) for var in var_wanted_float]
        
    try:    
        all_it_file_names = re.split('\n', RRead.BASH('ls '+OGpath+'all_iterations/'+sim.sim_name+'*'))
        all_it = [int(re.split('_|.h', file_i)[-2]) for file_i in all_it_file_names]
        all_hdf5it = all_it
        all_hdf5it = Lin.temporal_file['it'][0::Lin.param['IOHDF5::out_every']]
        it_plotted = [int(all_hdf5it[np.argmin(abs(np.array(all_hdf5it)-it))]) for it in it_wanted]
        it_plotted = np.sort(list(set(it_plotted)))
    except:
        it_plotted = np.sort(list(set([int(it) for it in it_wanted])))
    print('Iterations wanted:', it_wanted)
    print('Iterations that will be plotted:', it_plotted)
    return it_plotted


