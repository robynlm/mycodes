import h5py
import sys
import numpy as np
import pandas as pd
#from multiprocessing import Pool
from tools import ChooseParam
from tools import LinData
from tools import TAradius
from tools import MassCalc
from tools import GetVars_Plot2d as GVar
from tools import ODUDLoc
from tools import ReadingTools as RRead

class MyClass:
    def __init__(self, param, Lin):
        self.param = param
        self.cell_vol = param['dx']**3
        self.Lin = Lin
        locfinder = ODUDLoc.ODUDLocClass(param)
        self.defined_locs = locfinder.findlocations()
        self.TAr = TAradius.TA_Class(param, Lin, self.defined_locs[0])
        self.get_var = GVar.Get_var(param, Lin)
        self.it_file_name = (param['HorSpath'] + param['simname']
                             + '/output-0000/' + param['simname'] 
                             + '/all_iterations/' + param['simname'])
        self.MC = MassCalc.MassCalcClass(param, Lin, self.defined_locs[0])
        
    def average(self, phi, weight):
        return np.sum(phi * weight) / np.sum(weight)
        
    def getvals(self, Radius, all_hdf5it):
        data = []
        for r in Radius:
            weight = self.MC.getWeightsOfWholeBox(r)
            
            for i, it in enumerate(all_hdf5it):
                f = h5py.File('{}_it_{:06d}.hdf5'.format(self.it_file_name, 
                                                         int(it)), 'r')
                t = self.Lin.temp_from_temp('t', 'it', it)
            
                # Get the variable
                Kdict = self.get_var.get_the_curv(f, int(it))
                var = Kdict['A2']
                gdet = Kdict['metric_dic']['gdet']
                del Kdict
                
                # Compute averages
                var_average = self.average(var, weight)
                domain_weight = weight * np.sqrt(gdet) * self.cell_vol
                var_domain_average = self.average(var, domain_weight)
                
                if r==Radius[0]:
                    data += [[it, t, var_domain_average, var_average]]
                else:
                    data[i] += [var_domain_average, var_average]
            
            print('evo = {:.2f}'.format(100 * 
                                        np.argmin(abs(np.array(Radius)-r)) 
                                        / len(Radius)), flush=True)
        return data
        
if __name__ == "__main__":
    
    print("\n#################################")
    print("\n       Turn Around Radius        \n")
    print("#################################\n")
    
    # simulation_name
    
    simname = sys.argv[1]
    param = RRead.read_parameters(simname)
    Lin = LinData.LinData_Class(param)
    funcs = MyClass(param, Lin)

    all_it = Lin.temporal_file['it']
    all_hdf5it = np.array(all_it[0::param['IOHDF5::out_every']])
    print(' -- number of iterations : ', len(all_hdf5it), flush=True)
    
    radius_denominator = np.append(np.arange(3, 11, 1), np.arange(20, 128, 10))
    Radius = [param['Lx']/i for i in radius_denominator]
    data_values = funcs.getvals(Radius, all_hdf5it)
    
    data_header = ['it', 't']+list(np.ravel([[var+str(i) 
                                              for var in ['A2_dav_', 'A2_av_']] 
                                             for i in radius_denominator]))
    path = (param['HorSpath'] + '/' + param['simname'] 
            + '/output-0000/' + param['simname'] + '/')
    pd.DataFrame(data_values).to_csv(path+'A2_in_sphere.csv', 
                                     header=data_header, index=False)
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
            
