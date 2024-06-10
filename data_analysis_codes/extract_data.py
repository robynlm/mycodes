"""
Author: Robyn Munoz
Date: 30th June 2020

This code collects parameters from:

 - hdf5 files: the data needs to be split per iteration before use 
 then this code will
             - get the location of Over/Under density
             - determine which variables are available to collect
             - run each iteration file in MyClass, this will extract 
             the available parameters listed in data_names
             - record the results in path+'h5_data.csv'
             
 - ascii files: if possible this code will collect the varibables: 
 t, Ham_av, Mom_av, a_av, gdet_av, K_av, rho_av, press_av
                and record them in the file path+'asc_data.csv'

From the param.py file the following variables are needed:
 sim_name, refinement_level, nbr_ghost, Amp_pert, ti

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
from tools import NumMethods

def init(args):
    global it_counter
    it_counter = args

class MyClass:
    def __init__(self, param, Lin, nbr_iterations, verbose=False):
        self.verbose = verbose
        self.param = param
        self.cell_vol = param['dx'] * param['dy'] * param['dz']
        self.it_file_name = param['h5datapath'] + param['simname']
        self.nbr_iterations = nbr_iterations
        
        self.Lin = Lin
        self.FD = FD_File.FD_Class(param['dx'], param['dy'], param['dz'])
        self.RCW = RCW_file.Ricci_CoGrad_Weyl_Class(self.FD)

        locfinder = ODUDLoc.ODUDLocClass(param)        
        self.Locations, loclab = locfinder.findlocations()
        self.locsuffix = ['_av', '_L1', '_var', '_max', '_min'] + loclab
        self.print_header_val = True
        self.header = []

    # Averaging scheme
    def average(self, V, gdet, phi):
        try:
            return np.sum(phi*np.sqrt(gdet))*self.cell_vol/V
        except:
            return 0.0

    # Variance scheme
    def variance(self, V, gdet, phi, phi_av):
        return self.average(V, gdet, (phi-phi_av)**2)
    
    #L1 error
    def L1_average(self, var):
        return np.average(var)
    def L1_error(self, var, lin=0):
        if lin != 0:
            return self.L1_average(abs(var/lin-1))
        else:
            return self.L1_average(abs(var))
    def L1_constraint(self, T, B):
        return RRead.safe_division(self.L1_error(T), self.L1_error(B))
    
    def pert(self, var, var_str, t):
        var_th = getattr(self.Lin, var_str)(t)
        return RRead.safe_division(var, var_th)-1    
    def pert_withloc(self, var, var_str, t, loc):
        var_th = getattr(self.Lin, var_str)(t, loc)
        dvar = RRead.safe_division(var, var_th)-1 
        return dvar
    
    # Provide the average, variance, 
    # and the values at the under and over density
    def recval(self, V, gdet, F):
        if isinstance(F, (int, float, np.ndarray, np.generic)):
            F = [F]
        f_rec = []
        for f in F:
            f_av = self.average(V, gdet, f)
            f_L1 = self.L1_average(f)
            f_rec  += [f_av, f_L1, self.variance(V, gdet, f, f_av)]
            f_rec += [np.nanmax(abs(f)), np.nanmin(abs(f))]
            for loc in self.Locations:
                f_rec += [np.average(f[loc])]
        return f_rec

    def update_header(self, varnames, withrec=False):
        if withrec:
            for ivar in varnames:
                self.header += [ivar+measure for measure in self.locsuffix]
            if self.print_header_val:
                print(varnames, end=" ", flush=True)
        else:
            self.header += varnames
            if self.print_header_val:
                print(varnames, ' global value', end=" ", flush=True)
            
    def collect_var(self, key, f, NEWROW, Volume, gdet):
        if '::' in key:
            var_name = key.split('::')[1].split(' ')[0]
            it = int(key.split('it=')[1].split(' ')[0])
            var = self.get_1ddata(f, key)
            NEWROW += self.recval(Volume, gdet, [var])
            if it==0: self.update_header([var_name], withrec=True)
            del var
            
    def get_key(self, thorn_name, var_name, it):
        return thorn_name + '::{} it={} tl=0 rl=0'.format(var_name, it)
            
    def get_1ddata(self, f, key):
        if type(key)==list:
            key = self.get_key(key[0], key[1], key[2])
        data = RRead.fixij(f[key])
        data = RRead.cut0(data, self.param['ghost_size'], self.param['Nx'])
        return data
    
    def CalcMass(self, rho_u, gdet, region):
        return [np.sum(rho_u[region]*np.sqrt(gdet[region]))*self.cell_vol]

    """
    This function calculates the different values that will be recorded in 
    path+'data.asc'
    All values listed in return are recorded as is
    The __main__ loops over this function for each iteration
    """
    def getvals(self, it):       
        # Open iteration file
        it = int(it)
        f = h5py.File('{}_it_{:06d}.hdf5'.format(self.it_file_name, it), 'r')
        
        # working keys
        working_keys = [k for k in list(f.keys()) if ' rl=0' in k]
        
        NEWROW = [it]
        if self.verbose: print('it = ', it, flush=True)
        if it==0: self.update_header(['it'])

        t = f[working_keys[0]].attrs['time']
        NEWROW += [t]
        if self.verbose: print('t = ', t, flush=True)
        if it==0: self.update_header(['t'])

        #---------------------------------------------------
        # Spatial Metric
        #---------------------------------------------------
        # Component values
        if self.verbose: print('get metric', flush=True)
        gij_keys = [self.get_key('ADMBASE', gij, it) 
                    for gij in ['gxx', 'gxy', 'gxz', 'gyy', 'gyz', 'gzz']]
        metric, gdown_full = RRead.read_xyz(f, gij_keys)
        gdown = RRead.cut2(gdown_full, self.param['ghost_size'], self.param['Nx'])
        gup = RRead.inv3(gdown)
        del gij_keys, gdown_full, metric

        # Determinant and Volume
        gdet = RRead.det3(gdown)
        aloc = gdet**(1/6)
        Volume = np.sum(np.sqrt(gdet))*self.cell_vol
        NEWROW += self.recval(Volume, gdet, [gdet, aloc])
        if it==0: self.update_header(['gdet', 'aloc'], withrec=True)
        del aloc

        NEWROW += self.recval(Volume, gdet, 
                              [gup[0, 0], gup[0, 1], gup[0, 2], 
                               gup[1, 1], gup[1, 2], gup[2, 2]])
        if it==0: self.update_header(['guxx', 'guxy', 'guxz', 
                                      'guyy', 'guyz', 'guzz'], 
                                     withrec=True)

        #---------------------------------------------------
        # All variables
        #---------------------------------------------------
        for key in working_keys:
            self.collect_var(key, f, NEWROW, Volume, gdet)

        #---------------------------------------------------
        #  Curvature
        #---------------------------------------------------
        Christoffeludd, Christoffelddd = self.RCW.Christoffel_symbol(
            gdown, gup)
        RicciTdown, RicciS = self.RCW.Ricci_TandS(gup, Christoffeludd)
        NEWROW += self.recval(Volume, gdet, [RicciS])
        NEWROW += [self.L1_error(RicciS)]
        if it==0: self.update_header(['RicciS'], withrec=True)
        if it==0: self.update_header(['RicciS_L1error'])

        #---------------------------------------------------
        #  Lapse
        #---------------------------------------------------
        if self.get_key('ADMBASE', 'alp', it) in working_keys:
            if self.verbose: print('get lapse', flush=True)
            alpha = self.get_1ddata(f, ['ADMBASE', 'alp', it])
        else:
            alpha = np.ones(np.shape(gdown[0,0]))
        alpha2 = alpha**2

        if self.get_key('ADMBASE', 'dtalp', it) in working_keys:
            dtalpha = self.get_1ddata(f, ['ADMBASE', 'dtalp', it])
        else:
            dtalpha = np.zeros(np.shape(gdown[0,0]))
        alpha2 = alpha**2

        #---------------------------------------------------
        #  Shift
        #---------------------------------------------------
        if self.get_key('ADMBASE', 'betax', it) in working_keys:
            if self.verbose: print('get shift', flush=True)
            betax = self.get_1ddata(f, ['ADMBASE', 'betax', it])
            betay = self.get_1ddata(f, ['ADMBASE', 'betay', it])
            betaz = self.get_1ddata(f, ['ADMBASE', 'betaz', it])
            betaup = np.array([betax, betay, betaz])
            betadown = np.einsum('ij...,i...->j...', gdown, betaup)
            betasquare = np.einsum('i...,i...->...', betadown, betaup)
            NEWROW += self.recval(Volume, gdet, [betasquare])
            if it==0: self.update_header(['betasquare'], withrec=True)
            del betax, betay, betaz
        else:
            betaup = np.zeros(np.shape(gdown[0]))
            betadown = np.zeros(np.shape(gdown[0]))
            betasquare = np.zeros(np.shape(gdown[0,0]))

        #---------------------------------------------------
        #  Spacetime metric
        #---------------------------------------------------
        Box1 = np.ones(
            (self.param['Nx'], self.param['Ny'], self.param['Nz']))
        Box0 = np.zeros(
            (self.param['Nx'], self.param['Ny'], self.param['Nz']))
        g4down4 = np.array(
            [[-alpha2 + betasquare, 
              betadown[0], betadown[1], betadown[2]],
            [betadown[0], gdown[0,0], gdown[0,1], gdown[0,2]],
            [betadown[1], gdown[1,0], gdown[1,1], gdown[1,2]],
            [betadown[2], gdown[2,0], gdown[2,1], gdown[2,2]]])
        g4up4 = RRead.inv4(g4down4)

        #---------------------------------------------------
        #  Density
        #---------------------------------------------------
        if self.verbose: print('get rho_u', flush=True)
        # rho along fluidgauge
        if self.get_key('CT_DUST', 'rho', it) in working_keys:
            rho_u = self.get_1ddata(f, ['CT_DUST', 'rho', it])
        else:
            # rest mass energy density
            rho0 = self.get_1ddata(f, ['HYDROBASE', 'rho', it]) 
            # specific internal energy
            eps = self.get_1ddata(f, ['HYDROBASE', 'eps', it]) 
            rho_u = rho0 * (1 + eps)

        #---------------------------------------------------
        #  Fluid velocity
        #---------------------------------------------------
        if not param['synchronous']:
            if self.verbose: print('get fluid velocity', flush=True)
            if self.get_key('CT_DUST', 'u1', it) in working_keys:
                u1down = self.get_1ddata(f, ['CT_DUST', 'u1', it])
                u2down = self.get_1ddata(f, ['CT_DUST', 'u2', it])
                u3down = self.get_1ddata(f, ['CT_DUST', 'u3', it])
                W = self.get_1ddata(f, ['CT_DUST', 'W', it])
                W = W / np.sqrt(gdet)
                u0up = W / alpha
                udown3 = np.array([u1down, u2down, u3down])
                u0down = (-alpha2 * u0up 
                          + np.einsum('i...,i...->...', udown3, betaup))
                udown4 = np.array([u0down, u1down, u2down, u3down])
                del u0down, u1down, u2down, u3down, udown3
            else:
                W = self.get_1ddata(f, ['HYDROBASE', 'w_lorentz', it])
                u0up = W / alpha
                u1up = W * (self.get_1ddata(f, ['HYDROBASE', 'vel[0]', it]) 
                            - (betaup[0]/alpha))
                u2up = W * (self.get_1ddata(f, ['HYDROBASE', 'vel[1]', it]) 
                            - (betaup[1]/alpha))
                u3up = W * (self.get_1ddata(f, ['HYDROBASE', 'vel[2]', it]) 
                            - (betaup[2]/alpha))
                uup4 = np.array([u0up, u1up, u2up, u3up])
                udown4 = np.einsum('ab...,b...->a...', g4down4, uup4)
                u0down = udown4[0]
                del u0down, u1up, u2up, u3up, uup4

            hdown4 = g4down4 + np.einsum('a...,b...->ab...', 
                                         udown4, udown4)
            hdet4 = RRead.det4(hdown4)
            hdet3 = RRead.det3(hdown4[1:,1:])
            VolumeF = np.sum(np.sqrt(hdet3))*self.cell_vol
            NEWROW += self.recval(Volume, gdet, 
                                  [W, u0up, 
                                   udown4[1], udown4[2], udown4[3], 
                                   hdet4, hdet3])
            if it==0: self.update_header(['W', 'u0up', 
                                          'u1down', 'u2down', 'u3down', 
                                          'hdet4', 'hdet3'], 
                                         withrec=True)
            del W, u0up, hdet4
        else:
            udown4 = np.array([-Box1, Box0, Box0, Box0])
            hdown4 = np.array([[Box0, Box0, Box0, Box0],
                               [Box0, gdown[0,0], gdown[0,1], gdown[0,2]],
                               [Box0, gdown[1,0], gdown[1,1], gdown[1,2]],
                               [Box0, gdown[2,0], gdown[2,1], gdown[2,2]]])
            hdet3 = gdet
            VolumeF = Volume

        #---------------------------------------------------
        #  T00
        #---------------------------------------------------
        if self.get_key('TMUNUBASE', 'eTtt', it) in working_keys:
            if self.verbose: print('get T_00', flush=True)
            T00 = self.get_1ddata(f, ['TMUNUBASE', 'eTtt', it])
            rho_n = T00 / alpha2
            NEWROW += self.recval(Volume, gdet, [rho_n])
            if it==0: self.update_header(['rho_n'], withrec=True)
            del rho_n, T00, rho_u

        #---------------------------------------------------
        # Extrinsic Curvature
        #---------------------------------------------------
        # Component values
        if self.verbose: print('get the curv', flush=True)
        kij_keys = [self.get_key('ADMBASE', kij, it) 
                    for kij in ['kxx', 'kxy', 'kxz', 'kyy', 'kyz', 'kzz']]
        curv, Kdown_full = RRead.read_xyz(f, kij_keys)
        Kdown = RRead.cut2(Kdown_full, self.param['ghost_size'], 
                           self.param['Nx'])
        del kij_keys, curv, Kdown_full

        # Kup
        Kup = np.einsum('ia...,jb...,ab... -> ij...', gup, gup, Kdown)
        NEWROW += self.recval(Volume, gdet, 
                              [Kup[0, 0], Kup[0, 1], Kup[0, 2], 
                               Kup[1, 1], Kup[1, 2], Kup[2, 2]])
        if it==0: self.update_header(['Kuxx', 'Kuxy', 'Kuxz', 
                                      'Kuyy', 'Kuyz', 'Kuzz'], 
                                     withrec=True)
        del Kup

        # Trace 
        K = np.einsum('ij...,ij... -> ...', gup, Kdown)
        NEWROW += self.recval(Volume, gdet, [K])
        if it==0: self.update_header(['K'], withrec=True)

        # Traceless part of extrinsic curvature
        Adown = Kdown - (K/3)*gdown
        Aup = np.einsum('ib...,ja...,ab... -> ij...', gup, gup, Adown)
        A2 = np.einsum('ij...,ij... -> ...', Adown, Aup)/2
        NEWROW += self.recval(Volume, gdet, A2)
        if it==0: self.update_header(['A2'], withrec=True)
        del Adown, Aup, K, A2, gdet, Volume
        
        # Close iteration file and return results
        f.close()
        iteration = int(it / self.param['IOHDF5::out_every'])
        process = psutil.Process(os.getpid())
        Mem = process.memory_info()[0]/1024**2
        global it_counter
        with it_counter.get_lock():
            it_counter.value += 1
        percent_done = (it_counter.value * 100 / self.nbr_iterations)
        print(self.param['simname'] 
              + ', Iteration = {:4d},'.format(iteration)
              + ' EndMemory = {:.2f} MB,'.format(Mem)
              + ' Progress = {:.2f}%'.format(percent_done), flush=True)
        return NEWROW
    
if __name__ == "__main__":
    
    print("\n#################################", flush=True)
    print("\n          Extract Data           \n", flush=True)
    print("#################################\n", flush=True)
    
    # Input: simname nbr_processes
    
    # Simulation to analyse
    simname = sys.argv[1]
    param = RRead.read_parameters(simname)
    print(param['simname'], flush=True)
    Lin = LinData.LinData_Class(param)
    
    ###############################################
    print('\n===== Collecting iterations', flush=True)
    ###############################################
    if param['nbr_restarts']==1:
        all_it, all_t = RRead.collect_iteration_and_time(Lin.param)
        all_hdf5it = all_it[0::param['IOHDF5::out_every']
                            * (param['max_refinement_levels'] - 1)]
    else:
        all_hdf5it = RRead.collect_h5iteration(Lin.param)
        all_it = all_hdf5it
        all_it_wrl0 = []
        all_t = []
        for it in all_hdf5it:
            f = h5py.File('{}_it_{:06d}.hdf5'.format(
                param['h5datapath'] + param['simname'], it), 'r')
            all_t += [f[list(f.keys())[0]].attrs['time']]
            working_keys = [k for k in list(f.keys()) if ' rl=0' in k]
            if len(working_keys)!=0:
                all_it_wrl0 += [it]
        all_hdf5it = all_it_wrl0
        
    # Create time file
    filedat = np.array([all_it, all_t]).T
    pd.DataFrame(filedat).to_csv(param['datapath']+'Time_dt.csv', 
                                 header=['it', 't'], index=False)
    Lin = LinData.LinData_Class(param)
    
    nbr_iterations = len(all_hdf5it)
    print('number of iterations : ', int(nbr_iterations), flush=True)
    
    ###############################################
    print('\n===== Set up code', flush=True)
    ###############################################
    
    # Set up multiprocessing
    nbr_processes = int(sys.argv[2])
    print('with ', nbr_processes, ' processe(s)', flush=True)
    it_counter = Value('i',0)
    pool = Pool(processes = int(sys.argv[2]), 
                initializer = init, 
                initargs = (it_counter, ))
    
    # Define computing class
    funcs = MyClass(param, Lin, nbr_iterations, verbose=False)
    
    ###############################################
    print('\n===== Creating data header', flush=True)
    ###############################################
    
    print('the estimators and locations recorded are: ', 
          funcs.locsuffix, flush=True)
    print('unless the variable is already a global value', flush=True)
    NEWROW_it0 = funcs.getvals(0)
    save_header = funcs.header.copy()
    
    ###############################################
    print('\n===== Extract values', flush=True)
    ###############################################
    
    data_values = np.array(pool.map(funcs.getvals,  all_hdf5it[1:]))
    data_values = np.append(data_values, np.array([NEWROW_it0]), axis=0)
    data_values = data_values[data_values[:,0].argsort()]
    
    ###############################################
    print('\n===== Extract data function of time', flush=True)
    ###############################################
    
    # Already have proper time?
    if not param['synchronous']:
        if 'tau_av' in funcs.header:
            print('Got proper time from simulation output', flush=True)
            CalcPropTime = False
        else:
            print('Calculate proper time', flush=True)
            CalcPropTime = True
    
    # Update Header
    funcs.print_header_val = False
    if param['synchronous']:
        funcs.update_header(['tau', 'a', 'an', 'H', 'z'], 
                            withrec=False)
        keys = []
    else:
        if CalcPropTime:
            keys = ['tau', 'a', 'an', 'H', 'z']
        else:
            keys = ['a', 'an', 'H', 'z']
            
    keys += ['drho_u', 'ddrho_u']
    if 'theta' in funcs.header:
        keys += ['dTheta', 'ddTheta']
    
    for loc in funcs.locsuffix:
        for key in keys:
            funcs.update_header([key+loc], withrec=False)
    print(funcs.header[len(save_header):], flush=True)
    
    # Functions to call previously extracted variables
    def get_prev_vali(i, key):
        idx = np.where(np.array(funcs.header)==key)[0][0]
        return data_values[i, idx]
    def get_prev_val(key):
        idx = np.where(np.array(funcs.header)==key)[0][0]
        return data_values.T[idx]
    def func_of_tau(tau):
        return [Lin.evo.a(tau), Lin.an_initial(tau), 
                Lin.evo.Hprop(tau), Lin.evo.z(tau)]
    
    # Calculate proper time
    tau_dict = {}
    if not param['synchronous']:
        if CalcPropTime:
            for loc in funcs.locsuffix:
                tau_dict[loc] = RK4_tau(get_prev_val('t'), 
                                        get_prev_val('alpha'+loc))
        else:
            for loc in funcs.locsuffix:
                tau_dict[loc] = get_prev_val('tau'+loc)
    
    # Collect data
    new_data_values = []
    for i in range(len(all_hdf5it)):
        NEWROW = list(data_values[i]).copy()
        
        # Homogeneous
        if param['synchronous']:
            tau = get_prev_vali(i, 't')
            NEWROW += [tau] + func_of_tau(tau)
            
        # Inhomogeneous
        for loc in funcs.locsuffix:
            # Functions of proper time
            if not param['synchronous']:
                tau = tau_dict[loc][i]
                if CalcPropTime:
                    NEWROW += [tau] + func_of_tau(tau)
                else:
                    NEWROW += func_of_tau(tau)
                
            # Contrasts
            drho_u_f = funcs.pert(get_prev_vali(i, 'rho'+loc), 'rho', tau)
            drho_u = get_prev_vali(i, 'rho'+loc) / Lin.rho_u(tau) - 1
            if i==0: Lin.delta_initial[loc] = drho_u
            NEWROW += [drho_u]
            NEWROW += [funcs.pert_withloc(drho_u, 'drho', tau, loc)]
            
            if 'theta' in funcs.header:
                dTheta = funcs.pert(get_prev_vali(i, 'theta'+loc), 'Theta', tau)
                NEWROW += [dTheta]
                NEWROW += [funcs.pert_withloc(dTheta, 'dTheta', tau, loc)]
        new_data_values += [NEWROW]
            
    
    ###############################################
    print('\n===== Record data', flush=True)
    ###############################################
    pd.DataFrame(new_data_values).to_csv(param['datapath']+'h5_data.csv', 
                                     header=funcs.header, index=False)

    print(param['simname'], flush=True)
    print('Done :D', flush=True)
