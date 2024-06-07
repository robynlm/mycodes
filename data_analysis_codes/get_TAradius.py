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
        self.OD_loc = locfinder.findlocations()[0]
        self.TAr = TAradius.TA_Class(p, Lin, self.OD_loc)
        self.get_var = GVar.Get_var(Lin, OGpath)
        self.it_file_name = OGpath+'all_iterations/'+p.sim_name
    
    def getvals(self, it):
        f = h5py.File('{}_it_{:06d}.hdf5'.format(self.it_file_name, int(it)), 'r')
        t = self.Lin.temp_from_temp('t', 'it', it)
        K, dK, metric_dic = [self.get_var.get_the_curv(f, int(it))[var] for var in ['K', 'dK', 'metric_dic']]
        gdet = metric_dic['gdet']
        del metric_dic

        if dK[self.OD_loc]+1<0:
            del dK
            rt0p0c,  rt0p0p  =self.TAr.get_TA_Radius(gdet,K,[0.0,self.P2,self.P2], 
                                                     [0.0,0.0,self.P2], 1.0)
            rt45p0c, rt45p0p =self.TAr.get_TA_Radius(gdet,K,
                                                     [self.P4,self.P4,self.P4,self.P4,self.P2,self.P2],
                                                     [0.0,self.P2,np.pi,3*self.P2,self.P4,3*self.P4], 
                                                     np.sqrt(2))
            rt45p45c,rt45p45p=self.TAr.get_TA_Radius(gdet,K,
                                                     [self.P4,3*self.P4,self.P4,3*self.P4],
                                                     [self.P4,self.P4,3*self.P4,3*self.P4], 
                                                     np.sqrt(3))
            print('It = {}, TA radius = {:.2f}, {:.2f}, {:.2f}, {:.2f}, {:.2f}, {:.2f}'.format(int(it/self.p.h5_every), rt0p0c, rt0p0p, rt45p0c, rt45p0p, rt45p45c, rt45p45p), flush=True)
            return [it, t, rt0p0c, rt0p0p, rt45p0c, rt45p0p, rt45p45c, rt45p45p]
        else:
            print('It = {}, Before TA'.format(int(it/self.p.h5_every)), flush=True)
            return [it, t, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
            
if __name__ == "__main__":
    
    print("\n#################################")
    print("\n       Turn Around Radius        \n")
    print("#################################\n")
    
    sim, path, OGpath = ChooseParam.ChooseP(sys.argv[3], sys.argv[1])
    pool = Pool(processes=int(sys.argv[2]))
    Lin = LinData.LinData_Class(sim, path)
    funcs = MyClass(sim, OGpath, Lin)
    all_it = Lin.temporal_file['it']
    
    it_start = Lin.temp_from_temp('it', 'an', Lin.an_TA*0.9)
    
    all_hdf5it = np.array(all_it[0::sim.h5_every])
    it_used = all_hdf5it[np.argmin(abs(all_hdf5it-it_start)):]
    print(' -- number of iterations : ', len(all_hdf5it), flush=True)
    data_values = np.array(pool.map(funcs.getvals,  it_used))
    data_values = data_values[data_values[:,0].argsort()]
    pd.DataFrame(data_values).to_csv(path+'TA_radius.csv', header=['it', 't', 'R_T0_P0_COM', 'R_T0_P0_PROP', 'R_T45_P0_COM', 'R_T45_P0_PROP', 'R_T45_P45_COM', 'R_T45_P45_PROP'], index=False)
