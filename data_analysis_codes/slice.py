"""
Author: Robyn Munoz
Date: 

sys.arg : HorS simname it

"""

import re
import h5py
import sys
import numpy as np
import pandas as pd
from multiprocessing import Pool
from tools import ChooseParam
from tools import LinData
from tools import ReadingTools as RRead
from tools import GetVars_Plot2d as GVar
from tools import FD as FD_file
from tools import Ricci_CoGrad_Weyl as RCW_file

class MyClass:
    def __init__(self, p, HorS, OGpath, save_path, str_var_wanted):
        self.p = p
        self.Lin = Lin
        self.get_var = GVar.Get_var(param, Lin)
        self.save_path = save_path
        self.str_var_wanted = str_var_wanted
        self.it_file_name = OGpath+'all_iterations/'+p.sim_name
        self.FD = FD_file.FD_Class(p.dx, p.dx, p.dx)
        self.RCW = RCW_file.Ricci_CoGrad_Weyl_Class(self.FD)

    def diag(self, data):
        dim = np.shape(data)[0]
        dataxyz = np.zeros(dim)
        for j in range(dim):
            dataxyz[j] = data[j, j, j]
        return dataxyz

    def getcontrib(self, it):
        f = h5py.File('{}_it_{:06d}.hdf5'.format(self.it_file_name, it), 'r')   #Open iteration file

        t = self.Lin.temp_from_temp('t', 'it', it)
        var_wanted = self.Lin.temp_from_temp(self.str_var_wanted, 'it', it)
        K2_lin = (self.Lin.K(t)**2)

        metr_dic = self.get_var.get_the_metric(f, it)
        Christoffeludd = self.RCW.Christoffel_symbol(metr_dic['gdown'], metr_dic['gup'])
        RicciTdown, RicciS = self.RCW.Ricci_TandS(metr_dic['gup'], Christoffeludd)
        rho_dic = self.get_var.get_the_rho(f, it)
        curv_dic = self.get_var.get_the_curv(f, it)
        
        gdet_diag = self.diag(metr_dic['gdet'])
        dgdet_diag = self.diag(metr_dic['dgdet'])
        RicciS_diag = self.diag(RicciS)
        K_diag = self.diag(curv_dic['K'])
        dK_diag = self.diag(curv_dic['dK'])
        A2_diag = self.diag(curv_dic['A2'])
        rho_diag = self.diag(rho_dic['rho'])
        drho_diag = self.diag(rho_dic['drho'])
        
        # Electric and Magnetic parts of the Weyl tensor
        #Eup, E2   = self.RCW.Weyl_E(metr_dic['gdown'], metr_dic['gup'], RicciTdown, curv_dic['K'], curv_dic['Kdown'], curv_dic['Kmixed'], self.Lin.evo.kappa, rho_dic['rho'], self.Lin.evo.Lambda)
        #Bdown, B2 = self.RCW.Weyl_B(self.Lin.LeviCivita(), metr_dic['gup'], metr_dic['gmixed'], Christoffeludd, curv_dic['Kdown'])
        del RicciS, metr_dic, RicciTdown, curv_dic, rho_dic, Christoffeludd
        #E2_diag   = self.diag(E2)
        #B2_diag   = self.diag(B2)
        #M_diag    = self.diag(self.RCW.Weyl_M(Eup, Bdown))
        #L_diag    = self.diag(self.RCW.Weyl_L(E2, B2))
        #del Eup, E2, Bdown, B2
        
        # Hamiltonian Constraint
        Ham_diag = RicciS_diag + (2/3)*(K_diag**2) - 2*A2_diag - 2*self.Lin.evo.kappa*rho_diag - 2*self.Lin.evo.Lambda
        #HamPert_diag = RicciS_diag + (2/3)*((K_diag**2)-K2_lin) - 2*A2_diag - 2*self.Lin.kappa*(rho_diag-self.Lin.rho(t))

        # Raychaudhuri equation
        dtK_diag = (K_diag**2)/3 + 2*A2_diag + self.Lin.evo.kappa*rho_diag/2 - self.Lin.evo.Lambda
        #dtKPert_diag  = dtK_diag - (K2_lin/3 + self.Lin.kappa*self.Lin.rho(t)/2)

        colnames = ['gdet', 'dgdet', 'RicciS', 'K', 'dK', 'A2', 'rho', 'drho', 'Ham', 'dtK']
        data = np.array([gdet_diag, dgdet_diag, RicciS_diag, K_diag, dK_diag, A2_diag, rho_diag, drho_diag, Ham_diag, dtK_diag]).T
        pd.DataFrame(data, columns=colnames).to_csv(self.save_path+'Slice_'+str(self.str_var_wanted)+'={:.2f}.csv'.format(var_wanted))
        f.close()
    
if __name__ == "__main__":
    
    print("\n#################################")
    print("\n     Get Slice along x=y=z       \n")
    print("#################################\n")
    
    #horsloc, simname, strwanted, timewanted
    
    str_var_wanted = sys.argv[3]
    if str_var_wanted=='it':
        print('WARNING: Correspondence to lower res simulations will not be good')    
    
    sim_name = sys.argv[2]
    if 'N128' in sim_name:
        sim_name_cutup = re.split('N128', sim_name)
        sim_names = [sim_name_cutup[0]+N+sim_name_cutup[1] for N in ['N32', 'N64', 'N128']]
    elif 'N64' in sim_name:
        sim_name_cutup = re.split('N64', sim_name)
        sim_names = [sim_name_cutup[0]+N+sim_name_cutup[1] for N in ['N32', 'N64']]
    elif 'N32' in sim_name:
        sim_name_cutup = re.split('N32', sim_name)
        sim_names = [sim_name_cutup[0]+N+sim_name_cutup[1] for N in ['N32']]
    else:
        sim_names = [sim_name]
        
        
    for sim_name in sim_names:
        sim, path = ChooseParam.ChooseP(sim_name, sys.argv[1])
        param = RRead.read_parameters(sim_name)
        Lin = LinData.LinData_Class(param)
        itplotted = ChooseParam.ChooseIt(str_var_wanted, sys.argv[4], Lin, sim, path)

        save_path = path+'DataSlice/'
        RRead.MakeDir(save_path)

        ###############################################
        #         Get and rec data along slice
        ###############################################

        print(' -- nbr of iterations : ', len(itplotted))
        print(' i = ', sep=' ', end=' ', flush=True)
        funcs = MyClass(sim, sys.argv[1], path, save_path, str_var_wanted)
        for i, it in enumerate(itplotted, start=1):
            print(i, ', ', sep=' ', end=' ', flush=True)
            funcs.getcontrib(it)
        print('\n Done!')


