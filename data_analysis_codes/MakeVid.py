"""
Author: Robyn Munoz
Date: 11th June 2020

After a simulation is run and the data is split per iteration 
this code is used to make a video of the different variables available in GetVars_Plot2d.py
the images will be along the x-d (y=z) plane
and will be recorded in the directory: path+'videos'.

To run this type in a terminal:
>> python MakeVid.py sim_name var fixcolorbar

The 'fixcolorbar' keyword will make the colorbar fixed during the simulation.

From the param.py file the following variables are needed:
dx_ref, L

"""

import sys
import h5py
import numpy as np
from tools import Plot2d
from tools import ChooseParam
from tools import ReadingTools as RRead
from tools import GetVars_Plot2d as GVar
from tools import LinData
        
if __name__ == "__main__":
    
    print("\n#################################")
    print("\n           Make Video            \n")
    print("#################################\n")
    
    # Read user input parameters
    simname = sys.argv[1]
    param = RRead.read_parameters(simname)
    Lin = LinData.LinData_Class(param)
        
    # Get the function that will collect var
    var = sys.argv[2] # Variable that will be calculated
    print('variable: '+var)
    get_var = GVar.Get_var(param, Lin)
    try:
        f_var = getattr(get_var, 'f_'+var)
    except:
        print("\n Error: desired variable is not an option, those available are: "+', '.join(get_var.options)+'\n')
        
    # Make the folder where video will be created
    folder = param['datapath']+'videos'
    RRead.MakeDir(folder)
    # In case var.mp4 already exists delete it
    try:
        RRead.BASH("rm "+folder+"/"+var+".mp4")
        print("Remove existing "+folder+"/"+var+".mp4")
    except:
        print("No existing "+folder+"/"+var+".mp4")
        
    #==================================
    print('\n Collecting data')
    recfile = h5py.File(folder + "/data.hdf5", 'w')
    all_it, all_t = RRead.collect_iteration_and_time(param['datapath'])
    all_hdf5it = all_it[0::param['IOHDF5::out_every']]
    print(' -- nbr of iterations : ', len(all_hdf5it))
    print(' i = ', sep=' ', end=' ', flush=True)
    print('no chunk distinction', flush=True)
    for i, it in enumerate(all_hdf5it, start=1):
        print(i, ', ', sep=' ', end=' ', flush=True)
        data = f_var(it)
            
        #Record
        recfile.create_dataset('data_{}'.format(it), data=data)
    recfile.close()
    

    #==================================
    print('\n Making images')
        
    # Would you like the ticks of the colorbar to change during the simulation
    # or to be fixed? This is controled with the following logical
    if len(sys.argv)>5:
        CBarlogical = sys.argv[5] == 'fixcolorbar'
    else:
        CBarlogical = False
    print('fixcolorbar :', CBarlogical)
    Dmin = np.min(data)
    Dmax = np.max(data)
    print('min/max :', Dmin, Dmax)

    # Axis of plane to be plotted
    plane = sys.argv[3]
    beginningfilename = folder+'/'+var
    makeplot = Plot2d.MyPlot2dclass(param, Lin, 
                                    title = get_var.labels[var], 
                                    filename = beginningfilename, 
                                    fixcolorbar = CBarlogical, 
                                    CBarLimits = [Dmin, Dmax])
    try:
        f_plot = getattr(makeplot, 'f_'+plane)
    except:
        print("\n Error: desired plane is not an option, those available are: "
              + ', '.join(makeplot.options)+'\n')
    
    # Now make plots
    recfile = h5py.File(folder + "/data.hdf5", 'r')
    lab = get_var.labels[var]
    print(' -- nbr of iterations : ', len(all_hdf5it))
    print(' i = ', sep=' ', end=' ', flush=True)
    for i, it in enumerate(all_hdf5it, start=1):
        print(i, ', ', sep=' ', end=' ', flush=True)
        pdata = recfile['data_{}'.format(it)]
        f_plot(pdata, it)
    recfile.close()
    
    
    #==================================
    print('\n All images formed, now make video')
    RRead.BASH("ffmpeg -framerate 5 -pattern_type glob -i '"+folder+"/*.png'   -c:v libx264 -r 30 "+folder+"/"+var+".mp4")
    #RRead.BASH("ffmpeg -framerate 1 -i '"+folder+"/file%03d.png'   -c:v libx264 -r 30 "+folder+"/"+var+".mp4")
    
    
    #==================================
    print('\n Video made, now remove images')
    RRead.BASH("rm "+folder+"/data.hdf5")
    RRead.BASH("rm "+folder+"/*.png")
    print('\n Done!')

