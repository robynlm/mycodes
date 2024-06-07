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
        self.MC = MassCalc.MassCalcClass(p.dx, Lin, self.defined_locs[0])
        
    def getvals(self, Radius, all_hdf5it):
        data = []
        for r in Radius:
            weight = self.MC.getWeightsOfWholeBox(r)
            
            for i, it in enumerate(all_hdf5it):
                f = h5py.File('{}_it_{:06d}.hdf5'.format(self.it_file_name, int(it)), 'r')
                t = self.Lin.temp_from_temp('t', 'it', it)
            
                rho_dict = self.get_var.get_the_rho(f, int(it))
                gdet = self.get_var.get_the_metric(f, int(it))['gdet']
        
                Mass_prop = np.sum(rho_dict['rho']*weight*np.sqrt(gdet))*(self.p.dx**3)
                Mass = np.sum(weight*rho_dict['rho'])*self.p.dx**3
                V = np.sum(weight*np.sqrt(gdet))*self.p.dx**3
                delta_domain_average = np.sum(rho_dict['drho']*weight*np.sqrt(gdet))*(self.p.dx**3)/V
                delta_average = np.sum(rho_dict['drho']*weight)/np.sum(weight)
            
                iX, iY, iZ = np.where(self.TAr.radius_in_grid<r)
                region = [(ix, iy, iz) for ix, iy, iz in zip(iX, iY, iZ)]
                rAc,rAp=self.TAr.get_region_Radius(r,gdet,region,[0.0,self.P2,self.P2], [0.0,0.0,self.P2], 1.0)
                rBc,rBp=self.TAr.get_region_Radius(r,gdet,region,[self.P4,self.P4,self.P4,self.P4,self.P2,self.P2],
                                                                  [0.0,self.P2,np.pi,3*self.P2,self.P4,3*self.P4], np.sqrt(2))
                rCc,rCp=self.TAr.get_region_Radius(r,gdet,region,[self.P4,3*self.P4,self.P4,3*self.P4],
                                                                  [self.P4,self.P4,3*self.P4,3*self.P4], np.sqrt(3))
                if r==Radius[0]:
                    data += [[it, t, Mass_prop, Mass, delta_domain_average, delta_average, rAc, rAp, rBc, rBp, rCc, rCp]]
                else:
                    data[i] += [Mass_prop, Mass, delta_domain_average, delta_average, rAc, rAp, rBc, rBp, rCc, rCp]
            print('evo = {:.2f}'.format(100*np.argmin(abs(np.array(Radius)-r))/len(Radius)), flush=True)
        return data
        
if __name__ == "__main__":
    
    print("\n#################################")
    print("\n       Turn Around Radius        \n")
    print("#################################\n")
    
    # arg: 1 2
    # home_or_sciama simulation_name
    
    sim, path, OGpath = ChooseParam.ChooseP(sys.argv[2], sys.argv[1])
    Lin = LinData.LinData_Class(sim, path)
    
    funcs = MyClass(sim, OGpath, Lin)

    all_it = Lin.temporal_file['it']
    all_hdf5it = np.array(all_it[0::sim.h5_every])
    print(' -- number of iterations : ', len(all_hdf5it), flush=True)
    
    #data_values = np.array(pool.map(funcs.getvals,  all_hdf5it))
    #data_values = data_values[data_values[:,0].argsort()]
    radius_denominator = np.append(np.arange(3, 11, 1), np.arange(20, 128, 10))
    Radius = [sim.L/i for i in radius_denominator]
    data_values = funcs.getvals(Radius, all_hdf5it)
    
    data_header = ['it', 't']+list(np.ravel([[var+str(i) for var in ['Mass_prop_', 'Mass_', 'delta_dav_', 'delta_av_', 'r_vertice_com_', 'r_vertice_prop_', 'r_edge_com_', 'r_edge_prop_', 'r_face_com_', 'r_face_prop_']] for i in radius_denominator]))
    pd.DataFrame(data_values).to_csv(path+'mass_radius_evo.csv', header=data_header, index=False)
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
            
