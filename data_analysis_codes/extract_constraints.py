"""
Author: Robyn Munoz
Date: 26th Oct 2022
"""

import h5py
import sys
import re
import numpy as np
import pandas as pd
import os, psutil
from multiprocessing import Pool, Value
from tools import LinData
from tools import ReadingTools as RRead
from tools import Ricci_CoGrad_Weyl as RCW_file
from tools import FD as FD_File
from tools import ODUDLoc

def init(args):
    global it_counter
    it_counter = args

class MyClass:
    def __init__(self, param, Lin, nbrit):
    
        self.param = param
        self.it_file_name = param['h5datapath'] + param['simname']
        self.nbrit = nbrit
        
        self.Lin = Lin
        self.FD = FD_File.FD_Class(param['dx'], param['dy'], param['dz'])
        self.RCW = RCW_file.Ricci_CoGrad_Weyl_Class(self.FD)
        
        locfinder = ODUDLoc.ODUDLocClass(param)
        self.Locations, self.loclab = locfinder.findlocations()
        print('Locations recorded:', self.Locations)
        
        self.header = []
        
    def get_key(self, name, key, it):
        return name + '::{} it={} tl=0 rl=0'.format(key, it)
            
    def get_1ddata(self, f, name, key, it):
        data = RRead.fixij(f[self.get_key(name, key, it)])
        data = RRead.cut0(data, self.param['ghost_size'], self.param['Nx'])
        return data

    # Averaging scheme
    def average(self, V, gdet, phi):
        try:
            cell_vol = param['dx'] * param['dy'] * param['dz']
            return np.sum(phi*np.sqrt(gdet))*cell_vol/V
        except:
            return 0.0
        
    def recval(self, V, gdet, F):
        if isinstance(F, (int, float, np.ndarray, np.generic)):
            F = [F]
        f_rec = []
        for f in F:
            f_av = self.average(V, gdet, f)
            f_L1 = np.average(f)
            f_rec  += [f_av, f_L1]
            var = abs(f)
            f_rec += [np.nanmin(var), 
                      np.nanpercentile(var, 25), 
                      np.nanpercentile(var, 50), 
                      np.nanpercentile(var, 75), 
                      np.nanmax(var)]
            for loc in self.Locations:
                f_rec += [np.average(f[loc])]
        return f_rec

    def update_header(self, varnames, withrec=False):
        if withrec:
            suffix = ['_av', '_L1', '_minabs', '_lowQabs', 
                      '_medianabs', '_highQabs', '_maxabs'] + self.loclab
            for ivar in varnames:
                self.header += [ivar + measure for measure in suffix]
        else:
            self.header += varnames

    """
    The __main__ loops over this function for each iteration
    """
    def getvals(self, it):
        it = int(it)
        filename = '{}_it_{:06d}.hdf5'.format(self.it_file_name, it)
        f = h5py.File(filename, 'r')
        NEWROW = [it]
        if it==0: self.update_header(['it'])
            
        #---------------------------------------------------
        # Hamiltonian Constraint
        #---------------------------------------------------
                    
        # Get gdown, gdet, Volume, gup
        gij_keys = [self.get_key('ADMBASE', gij, it) 
                    for gij in ['gxx', 'gxy', 'gxz', 'gyy', 'gyz', 'gzz']]
        metric, gdown_full = RRead.read_xyz(f, gij_keys)
        gdown = RRead.cut2(gdown_full, self.param['ghost_size'], self.param['Nx'])
        gdet = RRead.det3(gdown)
        Volume = np.sum(np.sqrt(gdet))*self.param['dx']**3
        gup = RRead.inv3(gdown)
        del gij_keys, metric, gdown_full
                    
        # Get Gudd, RicciS
        Gudd, Gddd = self.RCW.Christoffel_symbol(gdown, gup)
        RicciTdown, RicciS = self.RCW.Ricci_TandS(gup, Gudd)
        del RicciTdown, Gddd
                    
        # Get K, A2
        kij_keys = [self.get_key('ADMBASE', kij, it) 
                    for kij in ['kxx', 'kxy', 'kxz', 'kyy', 'kyz', 'kzz']]
        curv, Kdown_full = RRead.read_xyz(f, kij_keys)
        Kdown = RRead.cut2(Kdown_full, self.param['ghost_size'], self.param['Nx'])
        K = np.einsum('ij...,ij... -> ...', gup, Kdown)
        Adown = Kdown - (K/3)*gdown
        Aup = np.einsum('ib...,ja...,ab... -> ij...', gup, gup, Adown)
        A2 = np.einsum('ij...,ij... -> ...', Adown, Aup)/2
        del kij_keys, curv, Kdown_full, gdown, Adown, Aup
                    
        # Get rho
        if ((self.get_key('ADMBASE', 'alp', it) in f.keys()) 
            and (self.get_key('TMUNUBASE', 'eTtt', it) in f.keys())):
            # rho in the frame normal to the hypersurface
            alpha = self.get_1ddata(f, 'ADMBASE', 'alp', it)
            T00 = self.get_1ddata(f, 'TMUNUBASE', 'eTtt', it) # down
            if ((self.get_key('ADMBASE', 'betax', it) in f.keys()) 
                and (self.get_key('TMUNUBASE', 'eTtt', it) in f.keys())):
                betaup = np.array([self.get_1ddata(f, 'ADMBASE', 'betax', it),
                                   self.get_1ddata(f, 'ADMBASE', 'betay', it),
                                   self.get_1ddata(f, 'ADMBASE', 'betaz', it)])
                
                T0i = np.array([self.get_1ddata(f, 'TMUNUBASE', 'eTtx', it),
                                self.get_1ddata(f, 'TMUNUBASE', 'eTty', it),
                                self.get_1ddata(f, 'TMUNUBASE', 'eTtz', it)])
                
                Txy = self.get_1ddata(f, 'TMUNUBASE', 'eTxy', it)
                Txz = self.get_1ddata(f, 'TMUNUBASE', 'eTxz', it)
                Tyz = self.get_1ddata(f, 'TMUNUBASE', 'eTyz', it)
                Tij = np.array([[self.get_1ddata(f, 'TMUNUBASE', 'eTxx', it), Txy, Txz],
                               [Txy, self.get_1ddata(f, 'TMUNUBASE', 'eTyy', it), Tyz],
                               [Txz, Tyz, self.get_1ddata(f, 'TMUNUBASE', 'eTzz', it)]])
                rho = (T00 - 2 * np.einsum('i...,i...->...', T0i, betaup) 
                       + np.einsum('ij...,i...,j...->...', Tij, betaup, betaup)) / (alpha**2)
                Eflux = (np.einsum('ij...,j...->i...', Tij, betaup) - T0i) / alpha
                del betaup, T0i, Tij, Txy, Txz, Tyz
            else:
                rho = T00 / (alpha**2)
                Eflux = np.zeros((3, self.param['Nx'], self.param['Ny'], self.param['Nz']))
            del alpha, T00
        else:
            # rho in fluid frame equivalent to the prior in synchronous comoving gauge
            if self.get_key('CT_DUST', 'rho', it) in f.keys():
                rho = self.get_1ddata(f, 'CT_DUST', 'rho', it)
                Eflux = np.zeros((3, self.param['Nx'], self.param['Ny'], self.param['Nz']))
            else:
                rho = self.get_1ddata(f, 'HYDROBASE', 'rho', it)
                Eflux = np.zeros((3, self.param['Nx'], self.param['Ny'], self.param['Nz']))
                    
        # Calc Ham
        Ham = (RicciS + (2/3)*K**2 - 2*A2 
               - 2*self.Lin.evo.kappa*rho - 2*self.Lin.evo.Lambda)
        
        # Calc Ham energy scale
        Ham_energy_scale = np.sqrt(abs((RicciS)**2 + ((2/3)*K**2)**2 + (2*A2)**2 
                                       + (2*self.Lin.evo.kappa*rho)**2 
                                       + (2*self.Lin.evo.Lambda)**2))
        del RicciS, K, A2, rho
                    
        # Record
        NEWROW += self.recval(Volume, gdet, 
                              [Ham, Ham_energy_scale, 
                               RRead.safe_division(Ham, Ham_energy_scale)])
        if it==0: self.update_header(['Ham', 'HamEScale', 
                                      'Ham/HamEScale'], withrec=True)
        del Ham, Ham_energy_scale
            
        #---------------------------------------------------
        # Momentum Constraint
        #---------------------------------------------------
        
        # Calc Mom
        # TO DO: include extra term!
        DdKdd = self.RCW.CovD3_tensor2down(Gudd, Kdown)
        DdKud = np.einsum('bc..., bca... -> a...', gup, DdKdd)
        DdK = np.einsum('bc..., acb... -> a...', gup, DdKdd)
        Mom = DdKud - DdK - self.Lin.evo.kappa*Eflux
        del DdKdd, Gudd, Kdown
        
        # Calc Mom energy scale
        DdKud_DdKuu = np.einsum('a..., ad..., d... -> ...', DdKud, gup, DdKud)
        DdK_DuK = np.einsum('a..., ad..., d... -> ...', DdK, gup, DdK)
        Eflux2 = ((self.Lin.evo.kappa**2)
                  * np.einsum('a..., ad..., d... -> ...', Eflux, gup, Eflux))
        Mom_energy_scale = np.sqrt(abs(DdKud_DdKuu + DdK_DuK + Eflux2))
        del DdKud, DdK, gup, DdKud_DdKuu, DdK_DuK, Eflux, Eflux2
        
        # Record
        NEWROW += self.recval(Volume, gdet, 
                              [Mom[0], Mom[1], Mom[2], Mom_energy_scale])
        if it==0: self.update_header(['Mom1', 'Mom2', 'Mom3', 
                                      'MomEScale'], withrec=True)
        NEWROW += self.recval(Volume, gdet, 
                              [RRead.safe_division(Mom[0], Mom_energy_scale), 
                               RRead.safe_division(Mom[1], Mom_energy_scale), 
                               RRead.safe_division(Mom[2], Mom_energy_scale)])
        if it==0: self.update_header(['Mom1/MomEScale', 'Mom2/MomEScale', 
                                      'Mom3/MomEScale'], withrec=True)
        del Mom, Mom_energy_scale, Volume, gdet
            
        #---------------------------------------------------
        # Close iteration file and return results
        #---------------------------------------------------
        
        f.close()
        iteration = int(it / self.param['IOHDF5::out_every'])
        process = psutil.Process(os.getpid())
        Mem = process.memory_info()[0] / (1024**2)
        global it_counter
        with it_counter.get_lock():
            it_counter.value += 1
        percent_done = (it_counter.value * 100 
                        / self.nbrit)
        print(self.param['simname'] 
              + ', Iteration = {:5d},'.format(iteration) 
              + ' EndMemory = {:.2f} MB,'.format(Mem)
              + ' Progress = {:.2f}%'.format(percent_done), flush=True)
        return NEWROW
    
if __name__ == "__main__":
    
    print("\n#################################", flush=True)
    print("\n     Extract Constraint Data   \n", flush=True)
    print("#################################\n", flush=True)
    
    simname = sys.argv[1]
    nbr_processes = int(sys.argv[2])
    
    param = RRead.read_parameters(simname)
    Lin = LinData.LinData_Class(param)
    all_it, all_t = RRead.collect_iteration_and_time(Lin.param)
    all_hdf5it = all_it[0::param['IOHDF5::out_every']]
    Lin = LinData.LinData_Class(param)
    nbrit = len(all_hdf5it)
    print(' -- number of iterations : ', nbrit, flush=True)
    
    # Set up multiprocessing
    it_counter = Value('i', 0)
    pool = Pool(processes = nbr_processes, 
                initializer = init, 
                initargs = (it_counter, ))
    
    # Define classes
    funcs = MyClass(param, Lin, nbrit)
    
    # Make and record the header
    NEWROW_it0 = funcs.getvals(0)
    saveheader = funcs.header
    print(saveheader, flush=True)
    
    # Compute the data
    data_values = np.array(pool.map(funcs.getvals,  all_hdf5it[1:]))
    data_values = np.append(data_values, np.array([NEWROW_it0]), axis=0)
    data_values = data_values[data_values[:, 0].argsort()]
    
    # Record the data
    pd.DataFrame(data_values).to_csv(param['datapath']+'constraints.csv', 
                                     header=saveheader, index=False)
    

    
    
    
    
    
    
    
    
    
    
    
    
    
    
