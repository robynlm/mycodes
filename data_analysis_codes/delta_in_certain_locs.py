import h5py
import sys
import numpy as np
import pandas as pd
from multiprocessing import Pool
from tools import ChooseParam
from tools import LinData
from tools import TAradius
from tools import GetVars_Plot2d as GVar
from tools import ODUDLoc

class MyClass:
    def __init__(self, p, OGpath, Lin):
    
        self.P4 = np.pi/4
        self.P2 = np.pi/2
        self.p = p
        self.Lin = Lin
        locfinder = ODUDLoc.ODUDLocClass(p)
        self.defined_locs = locfinder.findlocations()
        self.TAr = TAradius.TA_Class(p, Lin, self.defined_locs[0])
        self.get_var = GVar.Get_var(Lin, OGpath)
        self.it_file_name = OGpath+'all_iterations/'+p.sim_name
    
    def getvals(self, it):
        f = h5py.File('{}_it_{:06d}.hdf5'.format(self.it_file_name, int(it)), 'r')
        t = self.Lin.temp_from_temp('t', 'it', it)
        a_lin = self.Lin.temp_from_temp('a', 'it', it)
        
        delta = self.get_var.get_the_rho(f, int(it))['drho']
        
        delta_defined_locs = [delta[i][0] for i in self.defined_locs]
        
        T0_P0   = [[0.0, self.P2, self.P2], [0.0, 0.0, self.P2], 1.0]
        T45_P0  = [[self.P4, self.P4, self.P4, self.P4, self.P2, self.P2], [0, self.P2, np.pi, 3*self.P2, self.P4, 3*self.P4]]
        T45_P45 = [[self.P4, 3*self.P4, self.P4, 3*self.P4], [self.P4, self.P4, 3*self.P4, 3*self.P4]]
        
        delta_t0p0_locs = []
        for i in [2.0, 3.0, 4.0, 6.0, 8.0]:
            wanted_length = self.p.L/i
            delta_loc = []
            for theta, phi in zip(T0_P0[0], T0_P0[1]):
                delta_loc += self.TAr.get_delta_in_loc(theta, phi, 1.0, delta, wanted_length)
            delta_t0p0_locs += [np.average(delta_loc)]
            
        delta_t45p0_locs = []
        for i in [2.0, 3.0, 4.0, 6.0, 8.0, 4.0/np.sqrt(2)]:
            wanted_length = self.p.L/i
            delta_loc = []
            for theta, phi in zip(T45_P0[0], T45_P0[1]):
                delta_loc += self.TAr.get_delta_in_loc(theta, phi, np.sqrt(2), delta, wanted_length)
            delta_t45p0_locs += [np.average(delta_loc)]
               
        delta_t45p45_locs = []
        for i in [2.0, 3.0, 4.0, 6.0, 8.0, 4.0/np.sqrt(3)]:
            wanted_length = self.p.L/i
            delta_loc = []
            for theta, phi in zip(T45_P45[0], T45_P45[1]):
                delta_loc += self.TAr.get_delta_in_loc(theta, phi, np.sqrt(3), delta, wanted_length)
            delta_t45p45_locs += [np.average(delta_loc)]
        
        
        print('It = {}'.format(int(it/self.p.h5_every)), flush=True)
        
        return [it, t]+delta_defined_locs+delta_t0p0_locs+delta_t45p0_locs+delta_t45p45_locs

if __name__ == "__main__":
    
    print("\n#################################")
    print("\n       Turn Around Radius        \n")
    print("#################################\n")
    
    sim, path, OGpath = ChooseParam.ChooseP(sys.argv[3], sys.argv[1])
    pool = Pool(processes=int(sys.argv[2]))
    Lin = LinData.LinData_Class(sim, path)
    funcs = MyClass(sim, OGpath, Lin)
    all_it = Lin.temporal_file['it']
    all_hdf5it = np.array(all_it[0::sim.h5_every])
    print(' -- number of iterations : ', len(all_hdf5it), flush=True)
    data_values = np.array(pool.map(funcs.getvals,  all_hdf5it))
    data_values = data_values[data_values[:,0].argsort()]
    pd.DataFrame(data_values).to_csv(path+'delta_at_loc.csv', header=['it', 't',
                'delta_OD', 'delta_midOD', 'delta_center', 'delta_midUD', 'delta_UD',
        'delta_T0_P0_L2', 'delta_T0_P0_L3', 'delta_T0_P0_L4', 'delta_T0_P0_L6', 'delta_T0_P0_L8', 
        'delta_T45_P0_L2', 'delta_T45_P0_L3', 'delta_T45_P0_L4', 'delta_T45_P0_L6', 'delta_T45_P0_L8', 'delta_T45_P0_dxfacL4', 
   'delta_T45_P45_L2', 'delta_T45_P45_L3', 'delta_T45_P45_L4', 'delta_T45_P45_L6', 'delta_T45_P45_L8', 'delta_T45_P45_dxfacL4'], index=False)
    
    
    
    