"""
Author: Robyn Munoz
Date: 23rd April 2024

I didn't finish writing this, the only part mising is:
partial_t u^mu (or) partial_t u_mu

"""

import h5py
import sys
import re
import numpy as np
import pandas as pd
import os, psutil
from tools import LinData
from multiprocessing import Pool, Value
from tools import ReadingTools as RRead
from tools import Ricci_CoGrad_Weyl as RCW_file
from tools import FD as FD_File
from tools import NumMethods

def init(args):
    global it_counter
    it_counter = args

class MyClass:
    def __init__(self, param, nbr_iterations, verbose=False):
        self.verbose = verbose
        self.param = param
        self.cell_vol = param['dx'] * param['dy'] * param['dz']
        self.it_file_name = param['h5datapath'] + param['simname']
        self.nbr_iterations = nbr_iterations
        
        self.FD = FD_File.FD_Class(param['dx'], param['dy'], param['dz'])
        self.RCW = RCW_file.Ricci_CoGrad_Weyl_Class(self.FD)
        self.recalculate = False
            
    def get_key(self, thorn_name, var_name, it, rl):
        return thorn_name + '::{} it={} tl=0 rl={}'.format(var_name, it, rl)
            
    def get_1ddata(self, f, key):
        if type(key)==list:
            key = self.get_key(key[0], key[1], key[2], key[3])
        data = RRead.fixij(f[key])
        return data
    
    def getTheta(self, it):
        
        # Open iteration file
        it = int(it)
        if self.verbose: print('it = ', it, flush=True)
        f = h5py.File('{}_it_{:06d}.hdf5'.format(self.it_file_name, it), 'r+') 
        
        first_key = list(f.keys())[0]
        varname = first_key.split(' ')[0].split('::')[1]
        all_rls = set([int(k.split('rl=')[1]) 
                       for k in list(f.keys()) 
                       if varname in k])
        for rl in all_rls:
            working_keys = [k for k in list(f.keys()) if 'rl='+str(rl) in k]
            Nx, Ny, Nz = np.shape(np.array(f[working_keys[0]]))
            Box1 = np.ones((Nx, Ny, Nz))
            Box0 = np.zeros((Nx, Ny, Nz))
            
            self.FD.dx = self.param['dx']/(2**rl)
            self.FD.dy = self.param['dy']/(2**rl)
            self.FD.dz = self.param['dz']/(2**rl)
            
            theta_key = 'COFLUID::theta it={:06d} rl={}'.format(it, rl)
            
            if (theta_key not in working_keys) or self.recalculate:

                #---------------------------------------------------
                # Spatial Metric
                #---------------------------------------------------
                # Component values
                if self.verbose: print('get metric', flush=True)
                gij_keys = [self.get_key('ADMBASE', gij, it, rl) 
                            for gij in ['gxx', 'gxy', 'gxz', 
                                        'gyy', 'gyz', 'gzz']]
                metric, gdown = RRead.read_xyz(f, gij_keys)
                gup = RRead.inv3(gdown)
                gdet = RRead.det3(gdown)
                sgdet = np.sqrt(gdet)
                del gij_keys, metric, gdet

                #---------------------------------------------------
                #  Curvature
                #---------------------------------------------------
                Gudd3, Gddd3 = self.RCW.Christoffel_symbol(gdown, gup)

                #---------------------------------------------------
                #  Lapse
                #---------------------------------------------------
                if self.get_key('ADMBASE', 'alp', it, rl) in working_keys:
                    if self.verbose: print('get lapse', flush=True)
                    alpha = self.get_1ddata(f, ['ADMBASE', 'alp', it, rl])
                else:
                    alpha = Box1
                alpha2 = alpha**2

                if self.get_key('ADMBASE', 'dtalp', it, rl) in working_keys:
                    dtalpha = self.get_1ddata(f, ['ADMBASE', 'dtalp', it, rl])
                else:
                    dtalpha = Box0

                #---------------------------------------------------
                #  Shift
                #---------------------------------------------------
                if self.get_key('ADMBASE', 'betax', it, rl) in working_keys:
                    if self.verbose: print('get shift', flush=True)
                    betax = self.get_1ddata(f, ['ADMBASE', 'betax', it, rl])
                    betay = self.get_1ddata(f, ['ADMBASE', 'betay', it, rl])
                    betaz = self.get_1ddata(f, ['ADMBASE', 'betaz', it, rl])
                    betaup = np.array([betax, betay, betaz])
                    betadown = np.einsum('ij...,i...->j...', gdown, betaup)
                    betasquare = np.einsum('i...,i...->...', betadown, betaup)
                    del betax, betay, betaz
                else:
                    betaup = np.array([Box0, Box0, Box0])
                    betadown = np.array([Box0, Box0, Box0])
                    betasquare = Box0

                if self.get_key('ADMBASE', 'dtbetax', it, rl) in working_keys:
                    if self.verbose: print('get dtshift', flush=True)
                    dtbetax = self.get_1ddata(f, ['ADMBASE', 'dtbetax', it, rl])
                    dtbetay = self.get_1ddata(f, ['ADMBASE', 'dtbetay', it, rl])
                    dtbetaz = self.get_1ddata(f, ['ADMBASE', 'dtbetaz', it, rl])
                    dtbetaup = np.array([dtbetax, dtbetay, dtbetaz])
                    del dtbetax, dtbetay, dtbetaz
                else:
                    dtbetaup = np.array([Box0, Box0, Box0])

                #---------------------------------------------------
                #  Spacetime metric
                #---------------------------------------------------
                g4down400 = -alpha2 + betasquare
                g4down4 = np.array([[g4down400,
                                    betadown[0], betadown[1], betadown[2]],
                                   [betadown[0], gdown[0,0], 
                                    gdown[0,1], gdown[0,2]],
                                   [betadown[1], gdown[1,0], 
                                    gdown[1,1], gdown[1,2]],
                                   [betadown[2], gdown[2,0], 
                                    gdown[2,1], gdown[2,2]]])

                g4up4 = RRead.inv4(g4down4)
                del betasquare

                #---------------------------------------------------
                #  Fluid velocity
                #---------------------------------------------------
                if not param['synchronous']:
                    if self.verbose: print('get fluid velocity', flush=True)
                    if self.get_key('CT_DUST', 'u1', it, rl) in working_keys:
                        u1down = self.get_1ddata(f, ['CT_DUST', 'u1', it, rl])
                        u2down = self.get_1ddata(f, ['CT_DUST', 'u2', it, rl])
                        u3down = self.get_1ddata(f, ['CT_DUST', 'u3', it, rl])
                        W = self.get_1ddata(f, ['CT_DUST', 'W', it, rl])
                        W = W / sgdet
                        u0up = W / alpha
                        udown3 = np.array([u1down, u2down, u3down])
                        u0down = (-alpha2 * u0up 
                                  + np.einsum('i...,i...->...', udown3, betaup))
                        udown4 = np.array([u0down, u1down, u2down, u3down])
                        uup4 = np.einsum('ab...,b...->a...', g4up4, udown4)
                        del u0up, u0down, u1down, u2down, u3down, udown3
                    else:
                        W = self.get_1ddata(f, ['HYDROBASE', 'w_lorentz', it, rl])
                        u0up = W / alpha
                        u1up = W * (self.get_1ddata(f, ['HYDROBASE', 'vel[0]', it, rl]) 
                                    - (betaup[0]/alpha))
                        u2up = W * (self.get_1ddata(f, ['HYDROBASE', 'vel[1]', it, rl]) 
                                    - (betaup[1]/alpha))
                        u3up = W * (self.get_1ddata(f, ['HYDROBASE', 'vel[2]', it, rl]) 
                                    - (betaup[2]/alpha))
                        uup4 = np.array([u0up, u1up, u2up, u3up])
                        udown4 = np.einsum('ab...,b...->a...', g4down4, uup4)
                        del u0up, u1up, u2up, u3up
                else:
                    W = Box1
                    udown4 = np.array([-Box1, Box0, Box0, Box0])
                    uup4 = np.einsum('ab...,b...->a...', g4up4, udown4)
                del Box1, alpha2

                #---------------------------------------------------
                #  Density
                #---------------------------------------------------

                # rest mass energy density
                if self.verbose: print('get rho_u', flush=True)
                if self.get_key('CT_DUST', 'rho', it, rl) in working_keys:
                    rho0 = self.get_1ddata(f, ['CT_DUST', 'rho', it, rl])
                else:
                    rho0 = self.get_1ddata(f, ['HYDROBASE', 'rho', it, rl])

                # specific internal energy
                if self.verbose: print('get eps', flush=True)
                if self.get_key('CT_DUST', 'eps', it, rl) in working_keys:
                    eps = self.get_1ddata(f, ['CT_DUST', 'eps', it, rl])
                elif self.get_key('HYDROBASE', 'eps', it, rl) in working_keys:
                    eps = self.get_1ddata(f, ['HYDROBASE', 'eps', it, rl])
                else:
                    eps = Box0

                eosw = 0
                enthalpy = (1 + eps) * (1 + eosw)
                press = eosw * rho0
                del Box0, eosw

                #---------------------------------------------------
                # Extrinsic Curvature
                #---------------------------------------------------
                # Component values
                if self.verbose: print('get the curv', flush=True)
                kij_keys = [self.get_key('ADMBASE', kij, it, rl) 
                            for kij in ['kxx', 'kxy', 'kxz', 'kyy', 'kyz', 'kzz']]
                curv, Kdown = RRead.read_xyz(f, kij_keys)
                del kij_keys, curv

                #---------------------------------------------------
                # Metric time derivative
                #---------------------------------------------------
                # dt of spatial metric indices up
                Kup = np.einsum('im...,jn...,ij...->mn...', gup, gup, Kdown)
                dbetaup = self.FD.D3_tensor1(betaup)
                Lbgup = (np.einsum('s...,sij...->ij...', 
                                   betaup, self.FD.D3_tensor2(gup))
                         - np.einsum('si...,sj...->ij...', dbetaup, gup)
                         - np.einsum('sj...,is...->ij...', dbetaup, gup))
                dtgup = Lbgup - 2*alpha*Kup
                del Kup, dbetaup, Lbgup

                # dt of square(determinant(spatial metric))
                Ktrace = np.einsum('ij...,ij...->...', gup, Kdown)
                CovDbeta = self.RCW.CovD3_tensor1up(Gudd3, betaup)
                dtsgdet = sgdet * ( - alpha * Ktrace
                                    + np.einsum('ii...->...', CovDbeta))
                del Gudd3, Ktrace, CovDbeta

                #---------------------------------------------------
                # Fluid time derivative
                #---------------------------------------------------

                # define conserved quantities in Wilson formalism
                D = rho0 * W * sgdet
                E = rho0 * eps * W * sgdet
                S4down = rho0 * enthalpy * udown4 * W * sgdet
                del rho0, eps, enthalpy

                S4up = np.einsum('im...,mn...->i...', g4up4, S4down)
                V = uup4[1:] / uup4[0]

                # define time derivatives of conserved quantities
                dtD = - np.einsum('ii...->...', self.FD.D3_tensor1(D * V))
                dtS = (- np.einsum('jij...->...',
                                   self.FD.D3_tensor2(np.einsum('i...,j...->ij...', 
                                                                S4down[1:], V))) 
                       + (S4down[0]*self.FD.D3_scalar(g4down400)/2)
                       + np.einsum('j...,ij...->i...',
                                   S4up[1:], self.FD.D3_tensor1(betadown))
                       + (np.einsum('j...,k...,ijk...->i...', 
                                    S4up[1:], S4up[1:], 
                                    self.FD.D3_tensor2(gdown))/(2*S4down[0]))
                       - alpha*sgdet*self.FD.D3_scalar(press))
                del betadown, g4down400, S4up

                # dtW
                W2m1 = W**2 - 1
                W2m1oDpE = W2m1 / (D + E)
                fac1 = RRead.safe_division(W2m1,
                                           2*W*np.einsum('ij...,i...,j...->...',
                                                         gup, S4down[1:], S4down[1:]))
                par1 = (np.einsum('ij...,i...,j...->...',
                                  dtgup, S4down[1:], S4down[1:])
                        + 2 * np.einsum('ij...,i...,j...->...',
                                        gup, S4down[1:], dtS))
                fac2 = - W2m1oDpE / W
                par2 = (dtD
                        - np.einsum('ii...->...', self.FD.D3_tensor1(E * V))
                        - press * np.einsum('ii...->...',
                                            self.FD.D3_tensor1(W * sgdet * V)))
                fac3 = press * W2m1oDpE * dtsgdet
                dtW = ((fac1 * par1 + fac2 * par2 + fac3)
                       / ( 1 - (press * sgdet * W2m1oDpE / W )))
                del W2m1, W2m1oDpE, fac1, par1, fac2, par2, fac3

                # dtE
                dtsgdetW = dtsgdet * W + sgdet * dtW
                dtE = (- np.einsum('ii...->...', self.FD.D3_tensor1(E * V))
                       - press*(dtsgdetW + np.einsum('ii...->...',
                                                     self.FD.D3_tensor1(sgdet*W*V))))
                del dtW, V, dtsgdet, dtsgdetW, sgdet, press

                # define time derivative of fluid velocity, \partial_t(u_\mu)
                dtudown3 = udown4[1:] * (RRead.safe_division(dtS, S4down[1:])
                                         - ((dtD + dtE) / (D + E)))
                dtudown0 = (((np.einsum('ij...,i...,j...->...',
                                        dtgup, udown4[1:], udown4[1:])
                              + 2*np.einsum('i...,i...->...', 
                                            uup4[1:], dtudown3))/(2*alpha*W)) 
                           - (uup4[0] * dtalpha / alpha))
                dtudown4 = np.array([dtudown0, dtudown3[0], dtudown3[1], dtudown3[2]])
                # delete variables unused in the rest of the code for memory
                del W, uup4, dtgup
                del D, E, S4down, dtD, dtS, dtE, dtudown3, dtudown0

                #---------------------------------------------------
                # 4D christoffel symbols
                #---------------------------------------------------
                Gudd4 = self.RCW.Christoffel_symbol4beta(alpha, dtalpha, 
                                                         betaup, dtbetaup,
                                                         Kdown, gup, gdown,
                                                         Gddd3)
                del alpha, dtalpha, betaup, dtbetaup, Kdown, gup, gdown, Gddd3

                #---------------------------------------------------
                # Fluid projection tensor
                #---------------------------------------------------
                hdown4 = g4down4 + np.einsum('a...,b...->ab...', udown4, udown4)
                hmixed4 = np.einsum('ac...,cb...->ab...', g4up4, hdown4)
                hup4 = np.einsum('ac...,bd...,cd...->ab...', g4up4, g4up4, hdown4)
                del g4down4

                #---------------------------------------------------
                # Fluid expansion
                #---------------------------------------------------
                # covariant derivative of fluid velocity
                CovDu = self.RCW.CovD4_tensor1down(Gudd4, udown4, dtudown4)
                CovariantCovDu = np.einsum('ab...,ac...->bc...', hmixed4, CovDu)
                del Gudd4, udown4, dtudown4, hmixed4, CovDu

                # expansion tensor
                Thetadown = self.RCW.symmetric_tensor(CovariantCovDu)

                # trace and traceless parts of expansion tensor 
                # -> expansion scalar and shear
                Theta = np.einsum('ab...,ab...->...', hup4, Thetadown)
                sheardown = Thetadown - (1/3) * Theta * hdown4
                del Thetadown, hdown4, hup4

                #---------------------------------------------------
                # Vorticity
                #---------------------------------------------------

                omegadown = self.RCW.antisymmetric_tensor(CovariantCovDu)
                omegaup = np.einsum('im...,jn...,mn...->ij...',
                                    g4up4, g4up4, omegadown)
                omegamag = np.einsum('ij...,ij...->...', omegaup, omegadown) / 2
                del CovariantCovDu, omegadown, omegaup, g4up4

                #---------------------------------------------------
                # Save data
                #---------------------------------------------------

                added_keys = ['Theta', 'sigmaxx', 'sigmaxy', 'sigmaxz', 'sigmayy', 
                              'sigmayz', 'sigmazz', 'omegamag']
                for k in added_keys:
                    if k in working_keys:
                        del f[k]

                f.create_dataset('COFLUID::theta it={:06d} rl={}'.format(it, rl), 
                                 data=Theta)
                f.create_dataset('COFLUID::sigmaxx it={:06d} rl={}'.format(it, rl), 
                                 data=sheardown[0,0])
                f.create_dataset('COFLUID::sigmaxy it={:06d} rl={}'.format(it, rl), 
                                 data=sheardown[0,1])
                f.create_dataset('COFLUID::sigmaxz it={:06d} rl={}'.format(it, rl), 
                                 data=sheardown[0,2])
                f.create_dataset('COFLUID::sigmayy it={:06d} rl={}'.format(it, rl), 
                                 data=sheardown[1,1])
                f.create_dataset('COFLUID::sigmayz it={:06d} rl={}'.format(it, rl), 
                                 data=sheardown[1,2])
                f.create_dataset('COFLUID::sigmazz it={:06d} rl={}'.format(it, rl), 
                                 data=sheardown[2,2])
                f.create_dataset('COFLUID::omegamag it={:06d} rl={}'.format(it, rl), 
                                 data=omegamag)
                del Theta, sheardown, omegamag
        
        #---------------------------------------------------
        # Close things up
        #---------------------------------------------------    

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
    
if __name__ == "__main__":
    
    print("\n#################################", flush=True)
    print("\n        Calculate Theta         \n", flush=True)
    print("#################################\n", flush=True)
    
    # Input: simname nbr_processes
    
    # Simulation to analyse
    simname = sys.argv[1]
    param = RRead.read_parameters(simname)
    print(param['simname'], flush=True)
    Lin = LinData.LinData_Class(param)
    
    ###############################################
    print('\n===== Iterations', flush=True)
    ###############################################
    all_hdf5it = RRead.collect_h5iteration(Lin.param)
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
    funcs = MyClass(param, nbr_iterations, verbose=False)
    
    ###############################################
    print('\n===== Calculating Theta', flush=True)
    ###############################################
    
    for it in all_hdf5it:
        funcs.getTheta(it)
    #pool.map(funcs.getTheta,  all_hdf5it)

    print(param['simname'], flush=True)
    print('Done :D', flush=True)
