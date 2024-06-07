"""
Author: Robyn Munoz
Date: 11th June 2020

After a simulation is run and the data is split per iteration 
this code is used to make plots of the different variables available in GetVars_Plot2d.py
the images will be along the desired plane: xd, xy, yz
and will be recorded in the directory: path+'plots/var'.

To run this type in a terminal:
>> python MakeVid.py sim_name HorS var plane iterations
sim_name: check param.py
HorS: home, sciama (from ChooseParam)
var: check GetVars_Plot2d
plane: xd, xy, yz
iterations: all, 1,2,3

From the param.py file the following variables are needed:
dx_ref, L, ti

"""

import sys
import numpy as np
from tools import Plot2d
from tools import ChooseParam
from tools import ReadingTools as RRead
from tools import GetVars_Plot2d as GVar
from tools import LinData
        
if __name__ == "__main__":
    
    print("\n#################################")
    print("\n           Make Plots            \n")
    print("#################################\n")
    
    ###############################################
    # Read user input parameters
    ###############################################
    
    sim, path, OGpath = ChooseParam.ChooseP(sys.argv[2], sys.argv[1])
    Lin = LinData.LinData_Class(sim, path)
    print('path: ', path)
    # Get the function that will collect var
    var = sys.argv[3]
    print('variable: '+var)
    get_var = GVar.Get_var(Lin, OGpath)
    try:
        f_var = getattr(get_var, 'f_'+var)
    except:
        print("\n Error: desired variable is not an option, those available are: "+', '.join(get_var.options)+'\n')
        
    # Make the folders where the plots will be recorded
    folder = path+'plots'
    foldervar = folder+'/'+var
    RRead.MakeDir(folder)
    RRead.MakeDir(foldervar)
    
    # Axis of plane to be plotted
    plane = sys.argv[4]
    print('plane :', plane)
    beginningfilename = foldervar+'/'+var
    makeplot = Plot2d.MyPlot2dclass(sim, Lin, title=get_var.labels[var], filename=beginningfilename)
    try:
        f_plot = getattr(makeplot, 'f_'+plane)
    except:
        print("\n Error: desired plane is not an option, those available are: "+', '.join(makeplot.options)+'\n')
    
    
    # Iterations to be plotted
    itplotted = ChooseParam.ChooseIt(sys.argv[5], sys.argv[6], Lin, sim, OGpath)
    
        
    ###############################################
    #               Make the plots
    ###############################################

    print(' -- nbr of iterations : ', len(itplotted))
    print(' i = ', sep=' ', end=' ', flush=True)
    for i, it in enumerate(itplotted, start=1):
        print(i, ', ', sep=' ', end=' ', flush=True)
        data = f_var(it)
        f_plot(data, it)        
    print('\n Done!')

