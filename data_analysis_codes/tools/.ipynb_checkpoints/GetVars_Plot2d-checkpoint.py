"""
Author: Robyn Munoz
Date: 9th June 2020

This class is to be called by MakeVid.py and MakePlot.py
For a given iteration it will calculate the 3d following variables:
 - f_gdet : normalised determinant of spatial metric
 - f_K    : normalised trace of extrinsic curvature
 - f_rho  : normalised density
 - f_RicciS   : Ricci Scalar
 
  From the param.py file the following variables are needed:
  sim_name, nbr_ghost

"""

import h5py
import numpy as np
from . import ReadingTools as RRead
from . import Ricci_CoGrad_Weyl as RCW_file
from . import FD as FD_file

class Get_var():
    def __init__(self, param, Lin):
        self.param = param
        self.FD = FD_file.FD_Class(param['dx'], param['dy'], param['dz'])
        self.RCW = RCW_file.Ricci_CoGrad_Weyl_Class(self.FD)
        self.Lin = Lin
        
        
        self.it_file_name = (param['h5datapath'] + param['simname'])
        
        # What is currently available with this class
        self.options = ['alpha', 'gdet', 'dgdet', 'K', 'absK', 'shear', 'rho_u', 'drho_u', 
                        'RicciS', 'Ca', 'Da', 'Za', 
                        'EWeyl', 'Etrace', 'BWeyl', 'Btrace', 
                        'Ham', 'RdiffBE', 'divE', 'divB', 'curlE', 'curlB', 
                        'L_invar', 'M_invar', 'ReI_invar', 'ImI_invar', 
                        'ReJ_invar', 'ImJ_invar', 'ReS_invar', 'ImS_invar']
        self.labels = {'alpha':r"$\alpha$", 'gdet':r"$\gamma$", 'dgdet':r"$\delta\gamma$", 
                       'K':r"$\delta K$", 
                       'absK':r"$|K|$", 'shear':r"$|\sigma^2|$", 
                       'rho_u':r"$\rho$", 
                       'drho_u':r"$\delta$", 'RicciS':r"${}^{(3)}R$", 
                       'Ca':r"$|\mathcal{C}_a|$", 'Da':r"$|\mathcal{D}_a|$", 
                       'Za':r"$|\mathcal{Z}_a|$", 'EWeyl':r"$E^2$", 
                       'Etrace':r"$E^T$", 'BWeyl':r"$B^2$", 'Btrace':r"$B^T$", 
                       'Ham':'Ham', 'RdiffBE':r'$B^2/E^2$', 
                       'divE':r'$|\nabla \cdot E|$', 'divB':r'$|\nabla \cdot B|$', 
                       'curlE':r'$|\nabla \times E|$', 
                       'curlB':r'$|\nabla \times B|$', 'L_invar':r'$L$', 
                       'M_invar':r'$M$', 'ReI_invar':r'$Re(I)$', 
                       'ImI_invar':r'$Im(I)$', 'ReJ_invar':r'$Re(J)$', 
                       'ImJ_invar':r'$Im(J)$', 'ReS_invar':r'$Re(S)$', 
                       'ImS_invar':r'$Im(S)$'}
        
    def record(self, varname, var, it):
        #print('{}_it_{:06d}.hdf5'.format(self.it_file_name, it))
        with h5py.File('{}_it_{:06d}.hdf5'.format(self.it_file_name, it),
                       'a') as fnew:
            if varname not in fnew.keys():
                fnew[varname]=var
                
    def get_the_metric(self, f, it):
        #t = self.Lin.temp_from_temp('t', 'it', it)
        # Read spatial metric components
        gij_keys = ['ADMBASE::{} it={} tl=0 rl=0'.format(
            gij, it) for gij in 
                      ['gxx', 'gxy', 'gxz', 'gyy', 'gyz', 'gzz']]
        metric, gdown_full = RRead.read_xyz(f, gij_keys)
        gdown = RRead.cut2(gdown_full, self.param['ghost_size'], self.param['Nx'])
        
        # Calculate the determinant
        gdet_full = RRead.det3(metric)
        gdet = RRead.cut0(gdet_full, self.param['ghost_size'], self.param['Nx']) 
        # Cut ghost off
        #dgdet = RRead.safe_division(gdet, self.Lin.gdet(t))-1    
        
        # Calculate gup
        gup_full = RRead.inv3(metric)
        gup = RRead.cut2(gup_full, self.param['ghost_size'], self.param['Nx'])
        
        # Calculate gmixed
        gmixed = np.einsum('ij...,jk...->ik...', gup, gdown)
        
        #return {"gdown":gdown, "gdet":gdet, "dgdet":dgdet, 
        #        "gup":gup, "gmixed":gmixed}
        return {"gdown":gdown, "gdet":gdet, 
                "gup":gup, "gmixed":gmixed}
    
    def get_the_curv(self, f, it):
        #t = self.Lin.temp_from_temp('t', 'it', it)
        metric_dic = self.get_the_metric(f, it)
        gdown = metric_dic['gdown']
        gup = metric_dic['gup']
    
        # Get extrinsic curvature
        kij_keys = ['ADMBASE::{} it={} tl=0 rl=0'.format(
            kij, it) for kij in 
                    ['kxx', 'kxy', 'kxz', 'kyy', 'kyz', 'kzz']]
        curv, Kdown_full = RRead.read_xyz(f, kij_keys)
        Kdown = RRead.cut2(Kdown_full, self.param['ghost_size'], self.param['Nx'])
        
        # Calculate the trace
        Kmixed = np.einsum('ik...,kj... -> ij...', gup, Kdown)
        K = np.einsum('ij...,ij... -> ...', gup, Kdown)
        #dK = RRead.safe_division(K, self.Lin.K(t))-1    
        KijKji = np.einsum('ij...,ji...->...', Kmixed, Kmixed)

        # Calculate traceless part
        Adown = Kdown - (K/3)*gdown
        Aup   = np.einsum('ib...,ja...,ab... -> ij...', gup, gup, Adown)
        A2    = np.einsum('ij...,ij... -> ...', Adown, Aup)/2
        
        #return {"Kdown":Kdown, "Kmixed":Kmixed, "K":K, "Adown":Adown, "A2":A2, 
        #        "dK":dK, "KijKji":KijKji, "metric_dic":metric_dic}
        return {"Kdown":Kdown, "Kmixed":Kmixed, "K":K, "Adown":Adown, "A2":A2, 
                "KijKji":KijKji, "metric_dic":metric_dic}
    
    def get_the_rho(self, f, it):
        t = self.Lin.temp_from_temp('t', 'it', it)
        
        rho_key = 'CT_DUST::rho it={} tl=0 rl=0'.format(it)
        rho_full = RRead.fixij(np.array(f[rho_key]))
        rho = RRead.cut0(rho_full, self.param['ghost_size'], self.param['Nx'])
        drho = RRead.safe_division(rho, self.Lin.rho_u(t))-1
        return {"rho":rho, "drho":drho}
    
    def get_the_ricci(self, f, it):
        metric_dic = self.get_the_metric(f, it)
        Christoffeludd = self.RCW.Christoffel_symbol(metric_dic['gdown'], 
                                                     metric_dic['gup'])
        RicciTdown, RicciS = self.RCW.Ricci_TandS(metric_dic['gup'], 
                                                  Christoffeludd)
        return {"RicciS":RicciS, "RicciTdown":RicciTdown, 
                "Christoffeludd":Christoffeludd}
    
    def get_the_EWeyl(self, f, it):
        metric_dic = self.get_the_metric(f, it)
        gdown = metric_dic['gdown']
        gup = metric_dic['gup']
        del metric_dic
        
        Box_0 = np.zeros(np.shape(gdown[0,0]))
        gdown4 = np.array([[-np.ones(np.shape(gdown[0,0])), Box_0, Box_0, Box_0],
                           [Box_0, gdown[0,0], gdown[0,1], gdown[0,2]],
                           [Box_0, gdown[1,0], gdown[1,1], gdown[1,2]],
                           [Box_0, gdown[2,0], gdown[2,1], gdown[2,2]]])
        LCuud3 = np.einsum('ae...,bf...,efc...->abc...', gup, 
                           gup, self.RCW.LeviCivita4(gdown4)[0,1:,1:,1:])
        del gdown4, Box_0
        
        Christoffeludd     = self.RCW.Christoffel_symbol(gdown, gup)
        RicciTdown, RicciS = self.RCW.Ricci_TandS(gup, Christoffeludd)
        Kdown = self.get_the_curv(f, it)['Kdown']
        rho = self.get_the_rho(f, it)['rho']
        Tdown3 = np.zeros(np.shape(gdown)) 
        # not really but valid, for what is computed
        Edict = self.RCW.Weyl_E(gdown, gup, LCuud3, Christoffeludd, RicciS, 
                                RicciTdown, Kdown, Tdown3)
        return Edict
    
    def get_the_BWeyl(self, f, it):
        metric_dic = self.get_the_metric(f, it)
        gdown = metric_dic['gdown']
        gup = metric_dic['gup']
        del metric_dic
        
        Box_0 = np.zeros(np.shape(gdown[0,0]))
        gdown4 = np.array([[-np.ones(np.shape(gdown[0,0])), Box_0, Box_0, Box_0],
                           [Box_0, gdown[0,0], gdown[0,1], gdown[0,2]],
                           [Box_0, gdown[1,0], gdown[1,1], gdown[1,2]],
                           [Box_0, gdown[2,0], gdown[2,1], gdown[2,2]]])
        LCuud3 = np.einsum('ae...,bf...,efc...->abc...', gup, gup, 
                           self.RCW.LeviCivita4(gdown4)[0,1:,1:,1:])
        nup = np.array([-np.ones(np.shape(gdown[0,0])), Box_0, Box_0, Box_0])
        del gdown4, Box_0
        
        
        Christoffeludd = self.RCW.Christoffel_symbol(gdown, gup)
        Kdown = self.get_the_curv(f, it)['Kdown']
        
        Bdict = self.RCW.Weyl_B(gdown, gup, nup, LCuud3, Christoffeludd, Kdown)
        return Bdict
    
    def get_the_invar(self, f, it):
        metric_dic = self.get_the_metric(f, it)
        gdown = metric_dic['gdown']
        gup = metric_dic['gup']
        del metric_dic

        Box_0 = np.zeros(np.shape(gdown[0,0]))
        gdown4 = np.array([[-np.ones(np.shape(gdown[0,0])), Box_0, Box_0, Box_0],
                           [Box_0, gdown[0,0], gdown[0,1], gdown[0,2]],
                           [Box_0, gdown[1,0], gdown[1,1], gdown[1,2]],
                           [Box_0, gdown[2,0], gdown[2,1], gdown[2,2]]])
        LCuud3 = np.einsum('ae...,bf...,efc...->abc...', gup, gup, 
                           self.RCW.LeviCivita4(gdown4)[0,1:,1:,1:])
        nup = np.array([-np.ones(np.shape(gdown[0,0])), Box_0, Box_0, Box_0])
        del gdown4, Box_0

        Christoffeludd     = self.RCW.Christoffel_symbol(gdown, gup)
        RicciTdown, RicciS = self.RCW.Ricci_TandS(gup, Christoffeludd)
        Kdown = self.get_the_curv(f, it)['Kdown']
        rho = self.get_the_rho(f, it)['rho']
        Tdown3 = np.zeros(np.shape(gdown)) 
        # not really but valid, for what is computed
        Edict = self.RCW.Weyl_E(gdown, gup, LCuud3, Christoffeludd, RicciS, 
                                RicciTdown, Kdown, self.Lin.evo.kappa, Tdown3)
        Bdict = self.RCW.Weyl_B(gdown, gup, nup, LCuud3, Christoffeludd, Kdown)
        
        return self.RCW.Weyl_Invar(gdown, Edict, Bdict)

    """ 
    ================================================================
    Function to calculate the normalised determinant of the 
    spatial metric at a given iteration 
    """
    def f_alpha(self, it):
        f = h5py.File('{}_it_{:06d}.hdf5'.format(self.it_file_name, it), 'r')
        alpha_key = 'ADMBASE::alp it={} tl=0 rl=0'.format(it)
        alpha_full = RRead.fixij(np.array(f[alpha_key]))
        alpha = RRead.cut0(alpha_full, self.param['ghost_size'], self.param['Nx'])
        f.close()
        return alpha
    
    def f_gdet(self, it):
        f = h5py.File('{}_it_{:06d}.hdf5'.format(self.it_file_name, it), 'r')
        gdet = self.get_the_metric(f, it)['gdet']
        f.close()
        return gdet
    
    def f_dgdet(self, it):
        gdet = self.f_gdet(it)
        t = self.Lin.temp_from_temp('t', 'it', it)
        gdet_flrw = self.Lin.gdet(t)
        return gdet/gdet_flrw-1
    
    """ 
    ================================================================
    Function to calculate the normalised trace of the 
    extrinsic curvature at a given iteration 
    """
    def f_K(self, it):
        f = h5py.File('{}_it_{:06d}.hdf5'.format(self.it_file_name, it), 'r')
        K = self.get_the_curv(f, it)['K']
        t = self.Lin.temp_from_temp('t', 'it', it)
        K_flrw = self.Lin.K(t)
        f.close()
        return K/K_flrw-1
    
    """ 
    ================================================================
    Function to calculate the normalised trace of the 
    extrinsic curvature at a given iteration 
    """
    def f_absK(self, it):
        f = h5py.File('{}_it_{:06d}.hdf5'.format(self.it_file_name, it), 'r')
        K = self.get_the_curv(f, it)['K']
        f.close()
        return abs(K)
    
    """ 
    ================================================================
    Function to calculate the normalised trace of the 
    extrinsic curvature at a given iteration 
    """
    def f_shear(self, it):
        f = h5py.File('{}_it_{:06d}.hdf5'.format(self.it_file_name, it), 'r')
        A2 = self.get_the_curv(f, it)['A2']
        f.close()
        return abs(A2)
    
    """ 
    ================================================================
    Function to calculate the normalised density at a given iteration 
    """
    def f_rho_u(self, it):
        f = h5py.File('{}_it_{:06d}.hdf5'.format(self.it_file_name, it), 'r')
        rho = self.get_the_rho(f, it)['rho']
        f.close()
        return rho
    
    def f_drho_u(self, it):
        f = h5py.File('{}_it_{:06d}.hdf5'.format(self.it_file_name, it), 'r')
        drho = self.get_the_rho(f, it)['drho']
        f.close()
        return drho
    
    """ 
    ================================================================
    Function to calculate the Ricci scalar at a given iteration 
    """
    def f_RicciS(self, it):
        f = h5py.File('{}_it_{:06d}.hdf5'.format(self.it_file_name, it), 'r')
        RicciS = self.get_the_ricci(f, it)['RicciS']
        f.close()
        #self.record('RicciS', RicciS, it)
        return RicciS
    
    """ 
    ================================================================
    Function to calculate Ca at a given iteration 
    """
    def f_Ca(self, it):
        # Get spatial components
        f = h5py.File('{}_it_{:06d}.hdf5'.format(self.it_file_name, it), 'r')
        metric_dic = self.get_the_metric(f, it)
        gdown = metric_dic['gdown']
        gup = metric_dic['gup']
        
        t = self.Lin.temp_from_temp('t', 'it', it)
        a = self.Lin.evo.a(t)
        
        # Pass to specialised class to get RicciS
        Christoffeludd     = self.RCW.Christoffel_symbol(gdown, gup)
        RicciTdown, RicciS = self.RCW.Ricci_TandS(gup, Christoffeludd)  
        
        Ca = self.RCW.CoGrad_Ca(gmixed, a_lin, RicciS)
        Ca_norm = RRead.norm(Ca)
        f.close()
        return Ca_norm

    """ 
    ================================================================
    Function to calculate Da at a given iteration 
    """
    def f_Da(self, it):
        # Get spatial components
        f = h5py.File('{}_it_{:06d}.hdf5'.format(self.it_file_name, it), 'r')
        
        gmixed = self.get_the_metric(f, it)['gmixed']
        t = self.Lin.temp_from_temp('t', 'it', it)
        a = self.Lin.evo.a(t)
        
        rho = self.get_the_rho(f, it)
        
        Da = self.RCW.CoGrad_Da(gmixed, a_lin, rho)
        Da_norm = RRead.norm(Da)
        f.close()
        return Da_norm
    
    """ 
    ================================================================
    Function to calculate Za at a given iteration 
    """
    def f_Za(self, it):
        # Get spatial components
        f = h5py.File('{}_it_{:06d}.hdf5'.format(self.it_file_name, it), 'r')
        
        gmixed = self.get_the_metric(f, it)['gmixed']
        K = self.get_the_curv(f, it)['K']
        
        theta = -K  
        
        t = self.Lin.temp_from_temp('t', 'it', it)
        a = self.Lin.evo.a(t)
        
        Za = self.RCW.CoGrad_Za(gmixed, a_lin, K)
        Za_norm = RRead.norm(Za)
        f.close()
        return Za_norm
    
    """ 
    ================================================================
    Function to calculate the electric and magnetic part of the 
    Weyl tensor at a given iteration 
    """
    def f_EWeyl(self, it):
        # Get spatial components
        f = h5py.File('{}_it_{:06d}.hdf5'.format(self.it_file_name, it), 'r')
        data = self.get_the_EWeyl(f, it)['E2']
        f.close()
        return data
    
    def f_Etrace(self, it):
        # Get spatial components
        f = h5py.File('{}_it_{:06d}.hdf5'.format(self.it_file_name, it), 'r')
        data = self.get_the_EWeyl(f, it)['Etrace']
        f.close()
        return data
    
    def f_divE(self, it):
        # Get spatial components
        f = h5py.File('{}_it_{:06d}.hdf5'.format(self.it_file_name, it), 'r')
        data = self.get_the_EWeyl(f, it)['divE_norm']
        f.close()
        return data
    
    def f_curlE(self, it):
        # Get spatial components
        f = h5py.File('{}_it_{:06d}.hdf5'.format(self.it_file_name, it), 'r')
        data = self.get_the_EWeyl(f, it)['curlE_norm']
        f.close()
        return data
    
    def f_BWeyl(self, it):
        # Get spatial components
        f = h5py.File('{}_it_{:06d}.hdf5'.format(self.it_file_name, it), 'r')
        data = self.get_the_BWeyl(f, it)['B2']
        f.close()
        return data
    
    def f_Btrace(self, it):
        # Get spatial components
        f = h5py.File('{}_it_{:06d}.hdf5'.format(self.it_file_name, it), 'r')
        data = self.get_the_BWeyl(f, it)['Btrace']
        f.close()
        return data
    
    def f_divB(self, it):
        # Get spatial components
        f = h5py.File('{}_it_{:06d}.hdf5'.format(self.it_file_name, it), 'r')
        data = self.get_the_BWeyl(f, it)['divB_norm']
        f.close()
        return data
    
    def f_curlB(self, it):
        # Get spatial components
        f = h5py.File('{}_it_{:06d}.hdf5'.format(self.it_file_name, it), 'r')
        data = self.get_the_BWeyl(f, it)['curlB_norm']
        f.close()
        return data
    
    def f_Ham(self, it):
        f = h5py.File('{}_it_{:06d}.hdf5'.format(self.it_file_name, it), 'r')
        curv_dic = self.get_the_curv(f, it)
        ricci_dic = self.get_the_ricci(f, it)
        rho = self.get_the_rho(f, it)['rho']
        Ham = (ricci_dic['RicciS'] + (2/3)*curv_dic['K']**2 
               - 2*curv_dic['A2'] - 2*self.Lin.evo.kappa*rho 
               - 2*self.Lin.evo.Lambda)
        return Ham
        
        