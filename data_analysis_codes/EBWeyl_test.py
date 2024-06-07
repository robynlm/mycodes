"""
Author: Robyn Munoz
Date: 12th June 2021

This class calculates the electric and magnetic parts of the Weyl tensor purely geometrically (never calling the Einstein equations) by assuming the comoving synchronous gauge.

"""

import numpy as np
import h5py
from tools import FD as FD_file
from tools import ReadingTools as RRead
import matplotlib.pyplot as plt
from multiprocessing import Pool, Value
import EBWeylparam

def init(args):
    global it_counter
    it_counter = args

class MyClass():
    def __init__(self, OGpath, all_it):
        self.N = EBWeylparam.N
        self.dx = EBWeylparam.L/self.N
        self.dt = self.dx*EBWeylparam.dtfac
        self.nbr_ghost = 0
        self.FD = FD_file.FD_Class(self.dx, periodic_boundary=EBWeylparam.boundary)#, order6=True)
        self.Rsym = False
        self.box1 = np.ones((self.N, self.N, self.N))
        self.box0 = np.zeros((self.N, self.N, self.N))
        self.Abox0 = np.array([self.box0]*(len(all_it)-self.FD.mask*2))
        self.OGpath = OGpath
        self.it_file_name = OGpath+'all_iterations/'+EBWeylparam.simname
        self.all_it = all_it
        
    def read(self, var, it):
        if it=='all':
            dstr = ['dt', 'dx', 'dy', 'dz']
            gstr = ['gxx', 'gxy', 'gxz', 'gyy', 'gyz', 'gzz']
            dgstr = np.ravel([[d+g for g in gstr] for d in dstr])
            if var in dgstr: # all it
                f = h5py.File(self.OGpath+var+'.hdf5', 'r')
                return np.array([f['it={:06d}'.format(it)] for it in self.all_it])
        else:
            if var=='Kdown':
                f = h5py.File('{}_it_{:06d}.hdf5'.format(self.it_file_name, it), 'r')
                kij_keys = ['{} it={}'.format(kij, it) for kij in 
                        ['kxx', 'kxy', 'kxz', 'kyy', 'kyz', 'kzz']]
                curv, Kdown_full = RRead.read_xyz(f, kij_keys, fixindexes=False)
                Kdown = RRead.cut2(Kdown_full, self.nbr_ghost, self.N)
                return Kdown
            if var=='gdown' or var=='gup':
                f = h5py.File('{}_it_{:06d}.hdf5'.format(self.it_file_name, it), 'r')
                gij_keys = ['{} it={}'.format(gij, it) for gij in 
                        ['gxx', 'gxy', 'gxz', 'gyy', 'gyz', 'gzz']]
                metric, gdown_full = RRead.read_xyz(f, gij_keys, fixindexes=False)
                gdown = RRead.cut2(gdown_full, self.nbr_ghost, self.N)
                if var=='gdown':
                    return gdown
                else:
                     return RRead.cut2(RRead.inv(metric, RRead.det(metric)), self.nbr_ghost, self.N), gdown
            if var in ['dgdown4', 'dtgdown4', 'dxgdown4', 'dygdown4', 'dzgdown4', 
                       'dtdtgdown4', 'dtdxgdown4', 'dtdygdown4', 'dtdzgdown4']:
                f = h5py.File(self.OGpath+'EBWeyl/it={:06d}.hdf5'.format(it), 'r')
                return f[var]
        
        
    def record(self, varname, var, it, separate_iterations=True):
        if separate_iterations:
            if it=='all':
                for i, it in enumerate(self.all_it[self.FD.mask*2:]):         #!!! dt FD here 2 or 4
                    with h5py.File(self.OGpath+'EBWeyl/it={:06d}.hdf5'.format(it), 'a') as fnew:
                        if varname not in fnew.keys():
                            fnew[varname]=var[:,:,i,:,:,:]
            else:
                with h5py.File(self.OGpath+'EBWeyl/it={:06d}.hdf5'.format(it), 'a') as fnew:
                        if varname not in fnew.keys():
                            fnew[varname]=var
        else:
            with h5py.File(self.OGpath+varname+'.hdf5', 'a') as fnew:
                if 'it={:06d}'.format(it) not in fnew.keys():
                    fnew['it={:06d}'.format(it)]=var
        
    def calc_dtgdown(self, it):
        Kdown = self.read('Kdown', it)
        self.record('dtgxx', -2*Kdown[0,0], it, separate_iterations=False)
        self.record('dtgxy', -2*Kdown[0,1], it, separate_iterations=False)
        self.record('dtgxz', -2*Kdown[0,2], it, separate_iterations=False)
        self.record('dtgyy', -2*Kdown[1,1], it, separate_iterations=False)
        self.record('dtgyz', -2*Kdown[1,2], it, separate_iterations=False)
        self.record('dtgzz', -2*Kdown[2,2], it, separate_iterations=False)
        self.record('dtgdown4',  np.array([[self.box0, self.box0, self.box0, self.box0],
                                           [self.box0, -2*Kdown[0,0], -2*Kdown[0,1], -2*Kdown[0,2]],
                                           [self.box0, -2*Kdown[1,0], -2*Kdown[1,1], -2*Kdown[1,2]],
                                           [self.box0, -2*Kdown[2,0], -2*Kdown[2,1], -2*Kdown[2,2]]]), it)
    def calc_dxgdown(self, it):
        gdown = self.read('gdown', it)
        dxgxx = self.FD.D3x(gdown[0,0])
        dxgxy = self.FD.D3x(gdown[0,1])
        dxgxz = self.FD.D3x(gdown[0,2])
        dxgyy = self.FD.D3x(gdown[1,1])
        dxgyz = self.FD.D3x(gdown[1,2])
        dxgzz = self.FD.D3x(gdown[2,2])
        self.record('dxgxx', dxgxx, it, separate_iterations=False)
        self.record('dxgxy', dxgxy, it, separate_iterations=False)
        self.record('dxgxz', dxgxz, it, separate_iterations=False)
        self.record('dxgyy', dxgyy, it, separate_iterations=False)
        self.record('dxgyz', dxgyz, it, separate_iterations=False)
        self.record('dxgzz', dxgzz, it, separate_iterations=False)
        self.record('dxgdown4', np.array([[self.box0, self.box0, self.box0, self.box0],
                                          [self.box0, dxgxx, dxgxy, dxgxz],
                                          [self.box0, dxgxy, dxgyy, dxgyz],
                                          [self.box0, dxgxz, dxgyz, dxgzz]]), it)
    def calc_dygdown(self, it):
        gdown = self.read('gdown', it)
        dygxx = self.FD.D3y(gdown[0,0])
        dygxy = self.FD.D3y(gdown[0,1])
        dygxz = self.FD.D3y(gdown[0,2])
        dygyy = self.FD.D3y(gdown[1,1])
        dygyz = self.FD.D3y(gdown[1,2])
        dygzz = self.FD.D3y(gdown[2,2])
        self.record('dygxx', dygxx, it, separate_iterations=False)
        self.record('dygxy', dygxy, it, separate_iterations=False)
        self.record('dygxz', dygxz, it, separate_iterations=False)
        self.record('dygyy', dygyy, it, separate_iterations=False)
        self.record('dygyz', dygyz, it, separate_iterations=False)
        self.record('dygzz', dygzz, it, separate_iterations=False)
        self.record('dygdown4', np.array([[self.box0, self.box0, self.box0, self.box0],
                                          [self.box0, dygxx, dygxy, dygxz],
                                          [self.box0, dygxy, dygyy, dygyz],
                                          [self.box0, dygxz, dygyz, dygzz]]), it)
    def calc_dzgdown(self, it):
        gdown = self.read('gdown', it)
        dzgxx = self.FD.D3z(gdown[0,0])
        dzgxy = self.FD.D3z(gdown[0,1])
        dzgxz = self.FD.D3z(gdown[0,2])
        dzgyy = self.FD.D3z(gdown[1,1])
        dzgyz = self.FD.D3z(gdown[1,2])
        dzgzz = self.FD.D3z(gdown[2,2])
        self.record('dzgxx', dzgxx, it, separate_iterations=False)
        self.record('dzgxy', dzgxy, it, separate_iterations=False)
        self.record('dzgxz', dzgxz, it, separate_iterations=False)
        self.record('dzgyy', dzgyy, it, separate_iterations=False)
        self.record('dzgyz', dzgyz, it, separate_iterations=False)
        self.record('dzgzz', dzgzz, it, separate_iterations=False)
        self.record('dzgdown4', np.array([[self.box0, self.box0, self.box0, self.box0],
                                          [self.box0, dzgxx, dzgxy, dzgxz],
                                          [self.box0, dzgxy, dzgyy, dzgyz],
                                          [self.box0, dzgxz, dzgyz, dzgzz]]), it)
    def calc_dtdtgdown4(self):
        dtdtgxy = self.FD.Dt(self.read('dtgxy', 'all'), self.dt)
        dtdtgxz = self.FD.Dt(self.read('dtgxz', 'all'), self.dt)
        dtdtgyz = self.FD.Dt(self.read('dtgyz', 'all'), self.dt)
        self.record('dtdtgdown4', np.array([[self.Abox0, self.Abox0, self.Abox0, self.Abox0],
                               [self.Abox0, self.FD.Dt(self.read('dtgxx', 'all'), self.dt), dtdtgxy, dtdtgxz],
                               [self.Abox0, dtdtgxy, self.FD.Dt(self.read('dtgyy', 'all'), self.dt), dtdtgyz],
                               [self.Abox0, dtdtgxz, dtdtgyz, self.FD.Dt(self.read('dtgzz', 'all'), self.dt)]]), 'all')
    def calc_dtdxgdown4(self):
        dtdxgxy = self.FD.Dt(self.read('dxgxy', 'all'), self.dt)
        dtdxgxz = self.FD.Dt(self.read('dxgxz', 'all'), self.dt)
        dtdxgyz = self.FD.Dt(self.read('dxgyz', 'all'), self.dt)
        self.record('dtdxgdown4', np.array([[self.Abox0, self.Abox0, self.Abox0, self.Abox0],
                               [self.Abox0, self.FD.Dt(self.read('dxgxx', 'all'), self.dt), dtdxgxy, dtdxgxz],
                               [self.Abox0, dtdxgxy, self.FD.Dt(self.read('dxgyy', 'all'), self.dt), dtdxgyz],
                               [self.Abox0, dtdxgxz, dtdxgyz, self.FD.Dt(self.read('dxgzz', 'all'), self.dt)]]), 'all')
    def calc_dtdygdown4(self):
        dtdygxy = self.FD.Dt(self.read('dygxy', 'all'), self.dt)
        dtdygxz = self.FD.Dt(self.read('dygxz', 'all'), self.dt)
        dtdygyz = self.FD.Dt(self.read('dygyz', 'all'), self.dt)
        self.record('dtdygdown4', np.array([[self.Abox0, self.Abox0, self.Abox0, self.Abox0],
                               [self.Abox0, self.FD.Dt(self.read('dygxx', 'all'), self.dt), dtdygxy, dtdygxz],
                               [self.Abox0, dtdygxy, self.FD.Dt(self.read('dygyy', 'all'), self.dt), dtdygyz],
                               [self.Abox0, dtdygxz, dtdygyz, self.FD.Dt(self.read('dygzz', 'all'), self.dt)]]), 'all')
    def calc_dtdzgdown4(self):
        dtdzgxy = self.FD.Dt(self.read('dzgxy', 'all'), self.dt)
        dtdzgxz = self.FD.Dt(self.read('dzgxz', 'all'), self.dt)
        dtdzgyz = self.FD.Dt(self.read('dzgyz', 'all'), self.dt)
        self.record('dtdzgdown4', np.array([[self.Abox0, self.Abox0, self.Abox0, self.Abox0],
                               [self.Abox0, self.FD.Dt(self.read('dzgxx', 'all'), self.dt), dtdzgxy, dtdzgxz],
                               [self.Abox0, dtdzgxy, self.FD.Dt(self.read('dzgyy', 'all'), self.dt), dtdzgyz],
                               [self.Abox0, dtdzgxz, dtdzgyz, self.FD.Dt(self.read('dzgzz', 'all'), self.dt)]]), 'all')
        
    def calc_metric(self, it):
        gup, gdown = self.read('gup', it)
        gdown4 = np.array([[-self.box1, self.box0, self.box0, self.box0],
                           [self.box0, gdown[0,0], gdown[0,1], gdown[0,2]],
                           [self.box0, gdown[1,0], gdown[1,1], gdown[1,2]],
                           [self.box0, gdown[2,0], gdown[2,1], gdown[2,2]]])
        gup4 = np.array([[-self.box1, self.box0, self.box0, self.box0],
                         [self.box0, gup[0,0], gup[0,1], gup[0,2]],
                         [self.box0, gup[1,0], gup[1,1], gup[1,2]],
                         [self.box0, gup[2,0], gup[2,1], gup[2,2]]])
        return gdown4, gup4
        
    def calc_dmetric(self, it):
        Kdown = self.read('Kdown', it)
        gup, gdown = self.read('gup', it)
        Kup = np.einsum('ib...,ja...,ab... -> ij...', gup, gup, Kdown)
        del gdown, Kdown
        
        # -- 1st derivative, down
        dgdown4 = np.array([self.read('dtgdown4', it), self.read('dxgdown4', it), self.read('dygdown4', it), self.read('dzgdown4', it)])
        
        # -- 1st derivative, up
        dtgup4 = np.array([[self.box0, self.box0, self.box0, self.box0],
                           [self.box0, 2*Kup[0,0], 2*Kup[0,1], 2*Kup[0,2]],
                           [self.box0, 2*Kup[1,0], 2*Kup[1,1], 2*Kup[1,2]],
                           [self.box0, 2*Kup[2,0], 2*Kup[2,1], 2*Kup[2,2]]])
        dgxx = self.FD.D3_scalar(gup[0,0]) # [dxgxx, dygxx, dzgxx]
        dgxy = self.FD.D3_scalar(gup[0,1])
        dgxz = self.FD.D3_scalar(gup[0,2])
        dgyy = self.FD.D3_scalar(gup[1,1])
        dgyz = self.FD.D3_scalar(gup[1,2])
        dgzz = self.FD.D3_scalar(gup[2,2])
        dxgup4 = np.array([[self.box0, self.box0, self.box0, self.box0],
                             [self.box0, dgxx[0], dgxy[0], dgxz[0]],
                             [self.box0, dgxy[0], dgyy[0], dgyz[0]],
                             [self.box0, dgxz[0], dgyz[0], dgzz[0]]])
        dygup4 = np.array([[self.box0, self.box0, self.box0, self.box0],
                             [self.box0, dgxx[1], dgxy[1], dgxz[1]],
                             [self.box0, dgxy[1], dgyy[1], dgyz[1]],
                             [self.box0, dgxz[1], dgyz[1], dgzz[1]]])
        dzgup4 = np.array([[self.box0, self.box0, self.box0, self.box0],
                             [self.box0, dgxx[2], dgxy[2], dgxz[2]],
                             [self.box0, dgxy[2], dgyy[2], dgyz[2]],
                             [self.box0, dgxz[2], dgyz[2], dgzz[2]]])
        del dgxx, dgxy, dgxz, dgyy, dgyz, dgzz
        dgup4 = np.array([dtgup4, dxgup4, dygup4, dzgup4])
        del dtgup4, dxgup4, dygup4, dzgup4
        return dgdown4, dgup4
                                    
    def calc_ddmetric(self, dgdown4, it):        
        # dt
        dtdtgdown4 = self.read('dtdtgdown4', it)
        dtdxgdown4 = self.read('dtdxgdown4', it)
        dtdygdown4 = self.read('dtdygdown4', it)
        dtdzgdown4 = self.read('dtdzgdown4', it)
        # dx
        dxdxgxy = self.FD.D3x(dgdown4[1,1,2])
        dxdxgxz = self.FD.D3x(dgdown4[1,1,3])
        dxdxgyz = self.FD.D3x(dgdown4[1,2,3])
        dxdxgdown4 = np.array([[self.box0, self.box0, self.box0, self.box0],
                               [self.box0, self.FD.D3x(dgdown4[1,1,1]), dxdxgxy,              dxdxgxz],
                               [self.box0, dxdxgxy,              self.FD.D3x(dgdown4[1,2,2]), dxdxgyz],
                               [self.box0, dxdxgxz,              dxdxgyz,              self.FD.D3x(dgdown4[1,3,3])]])
        del dxdxgxy, dxdxgxz, dxdxgyz
        
        dxdygxy = self.FD.D3x(dgdown4[2,1,2])
        dxdygxz = self.FD.D3x(dgdown4[2,1,3])
        dxdygyz = self.FD.D3x(dgdown4[2,2,3])
        dxdygdown4 = np.array([[self.box0, self.box0, self.box0, self.box0],
                               [self.box0, self.FD.D3x(dgdown4[2,1,1]), dxdygxy,              dxdygxz],
                               [self.box0, dxdygxy,              self.FD.D3x(dgdown4[2,2,2]), dxdygyz],
                               [self.box0, dxdygxz,              dxdygyz,              self.FD.D3x(dgdown4[2,3,3])]])
        del dxdygxy, dxdygxz, dxdygyz
        
        dxdzgxy = self.FD.D3x(dgdown4[3,1,2])
        dxdzgxz = self.FD.D3x(dgdown4[3,1,3])
        dxdzgyz = self.FD.D3x(dgdown4[3,2,3])
        dxdzgdown4 = np.array([[self.box0, self.box0, self.box0, self.box0],
                               [self.box0, self.FD.D3x(dgdown4[3,1,1]), dxdzgxy,              dxdzgxz],
                               [self.box0, dxdzgxy,              self.FD.D3x(dgdown4[3,2,2]), dxdzgyz],
                               [self.box0, dxdzgxz,              dxdzgyz,              self.FD.D3x(dgdown4[3,3,3])]])
        del dxdzgxy, dxdzgxz, dxdzgyz
        
        # dy
        dydygxy = self.FD.D3y(dgdown4[2,1,2])
        dydygxz = self.FD.D3y(dgdown4[2,1,3])
        dydygyz = self.FD.D3y(dgdown4[2,2,3])
        dydygdown4 = np.array([[self.box0, self.box0, self.box0, self.box0],
                               [self.box0, self.FD.D3y(dgdown4[2,1,1]), dydygxy,              dydygxz],
                               [self.box0, dydygxy,              self.FD.D3y(dgdown4[2,2,2]), dydygyz],
                               [self.box0, dydygxz,              dydygyz,              self.FD.D3y(dgdown4[2,3,3])]])
        del dydygxy, dydygxz, dydygyz
        
        dydzgxy = self.FD.D3y(dgdown4[3,1,2])
        dydzgxz = self.FD.D3y(dgdown4[3,1,3])
        dydzgyz = self.FD.D3y(dgdown4[3,2,3])
        dydzgdown4 = np.array([[self.box0, self.box0, self.box0, self.box0],
                               [self.box0, self.FD.D3y(dgdown4[3,1,1]), dydzgxy,              dydzgxz],
                               [self.box0, dydzgxy,              self.FD.D3y(dgdown4[3,2,2]), dydzgyz],
                               [self.box0, dydzgxz,              dydzgyz,              self.FD.D3y(dgdown4[3,3,3])]])
        del dydzgxy, dydzgxz, dydzgyz
                                    
        # dz
        dzdzgxy = self.FD.D3z(dgdown4[3,1,2])
        dzdzgxz = self.FD.D3z(dgdown4[3,1,3])
        dzdzgyz = self.FD.D3z(dgdown4[3,2,3])
        dzdzgdown4 = np.array([[self.box0, self.box0, self.box0, self.box0],
                               [self.box0, self.FD.D3z(dgdown4[3,1,1]), dzdzgxy,              dzdzgxz],
                               [self.box0, dzdzgxy,              self.FD.D3z(dgdown4[3,2,2]), dzdzgyz],
                               [self.box0, dzdzgxz,              dzdzgyz,              self.FD.D3z(dgdown4[3,3,3])]])
        del dgdown4, dzdzgxy, dzdzgxz, dzdzgyz
                                    
        ddgdown4 = np.array([[dtdtgdown4, dtdxgdown4, dtdygdown4, dtdzgdown4],
                             [dtdxgdown4, dxdxgdown4, dxdygdown4, dxdzgdown4],
                             [dtdygdown4, dxdygdown4, dydygdown4, dydzgdown4],
                             [dtdzgdown4, dxdzgdown4, dydzgdown4, dzdzgdown4]])
        #hopefully symmetries allow me to do this
        return ddgdown4
        
    def calc_Weyl(self, it):
        # =========================================================================
        # Define metric and derivatives
        gdown4, gup4 = self.calc_metric(it)
        dgdown4, dgup4 = self.calc_dmetric(it)
        ddgdown4 = self.calc_ddmetric(dgdown4, it)
                                    
        # =========================================================================
        # Define Christoffel symbols and derivative
        Gudd4 =  np.einsum('ad...,bdc...->abc...', gup4, dgdown4)/2
        Gudd4 += np.einsum('ad...,cbd...->abc...', gup4, dgdown4)/2
        Gudd4 += -np.einsum('ad...,dbc...->abc...', gup4, dgdown4)/2
                  
        dGudd4 =  np.einsum('ead...,bdc...->eabc...', dgup4, dgdown4)/2 + np.einsum('ad...,ebdc...->eabc...', gup4, ddgdown4)/2
        dGudd4 += np.einsum('ead...,cbd...->eabc...', dgup4, dgdown4)/2 + np.einsum('ad...,ecbd...->eabc...', gup4, ddgdown4)/2
        dGudd4 += -np.einsum('ead...,dbc...->eabc...', dgup4, dgdown4)/2 -np.einsum('ad...,edbc...->eabc...', gup4, ddgdown4)/2
        del dgdown4, dgup4, ddgdown4
        
        # =========================================================================
        # Define Riemann and Ricci tensors
        RiemannTuddd4 =  np.einsum('cedb...->ebcd...', dGudd4)
        RiemannTuddd4 += -np.einsum('decb...->ebcd...', dGudd4)
        RiemannTuddd4 += np.einsum('ecf...,fdb...->ebcd...', Gudd4, Gudd4)
        RiemannTuddd4 += -np.einsum('edf...,fcb...->ebcd...', Gudd4, Gudd4)
        del dGudd4, Gudd4
                                    
        RiemannTdown4 = np.einsum('ae...,ebcd...->abcd...', gdown4, RiemannTuddd4)
        del RiemannTuddd4
        
        # Force symmetry
        if self.Rsym:
            for a in range(4):
                for b in range(4):
                    for c in range(4):
                        for d in range(4):
                            RiemannTdown4[b,a,c,d] = - RiemannTdown4[a,b,c,d]
                            RiemannTdown4[a,b,d,c] = - RiemannTdown4[a,b,c,d]
                            RiemannTdown4[c,d,a,b] =   RiemannTdown4[a,b,c,d]
        
        RicciTdown4 = np.einsum('ac...,abcd...->bd...', gup4, RiemannTdown4)
        RicciS4 = np.einsum('ab...,ab... -> ...', gup4, RicciTdown4)
        self.record('RicciS4', RicciS4, it)
        
        # =========================================================================
        # Define Weyl tensor
        Wterm1 = (-1/2)*(np.einsum('ac...,bd...->abcd...', gdown4, RicciTdown4)
                         -np.einsum('ad...,bc...->abcd...', gdown4, RicciTdown4)
                         -np.einsum('bc...,ad...->abcd...', gdown4, RicciTdown4)
                         +np.einsum('bd...,ac...->abcd...', gdown4, RicciTdown4))
        del RicciTdown4
        Wterm2 = RicciS4*(np.einsum('ac...,bd...->abcd...', gdown4, gdown4)-np.einsum('ad...,bc...->abcd...', gdown4, gdown4))/6
        del RicciS4
        Wdown = RiemannTdown4 + Wterm1 + Wterm2
        self.record('Bianchi_ID', RiemannTdown4[1,0,2,3] - RiemannTdown4[2,0,1,3] + RiemannTdown4[3,0,1,2], it)
        del RiemannTdown4, Wterm1, Wterm2
        
        
        # =========================================================================
        # Define electric and magnetic parts
        nup = np.array([self.box1,self.box0,self.box0,self.box0])
        Edown = np.einsum('abcd...,b...,d...->ac...', Wdown, nup, nup)
        Etrace = np.einsum('ij...,ij...->...', gup4, Edown)
        self.record('Etrace', Etrace, it)
        Etrace = np.nanmean(abs(Etrace))
        Eup    = np.einsum('ab...,cd...,ac... -> bd...', gup4, gup4, Edown)
        E2     = np.einsum('ab...,ab...->...', Eup, Edown)
        self.record('E2', E2, it)
        del Eup, E2, Edown
                         
        LCuudd = np.einsum('eg...,fh...,ghcd...->efcd...', gup4, gup4, self.LeviCivita(gdown4))
        del gdown4
        WSdown = np.einsum('abef...,efcd...->abcd...',Wdown, LCuudd)/2
        del Wdown, LCuudd
        Bdown  = np.einsum('abcd...,d...,b...->ac...', WSdown, nup, nup)
        del WSdown, nup
        Btrace = np.einsum('ac...,ac...->...', Bdown, gup4)
        self.record('Btrace', Btrace, it)
        Btrace = np.nanmean(abs(Btrace))
        Bup    = np.einsum('ab...,cd...,ac... -> bd...', gup4, gup4, Bdown)
        del gup4
        B2     = np.einsum('ab...,ab...->...', Bup, Bdown)
        self.record('B2', B2, it)
        del Bup, B2, Bdown
        
        global it_counter
        with it_counter.get_lock():
            it_counter.value += 1
        percent_done = it_counter.value*100/len(self.all_it[2:-2])
        print('progress: {:.2f}%, ETrace = {:.2e}, BTrace = {:.2e}'.format(percent_done, Etrace, Btrace))
            
    def LeviCivita(self, gdown4): # 4 indices down
        LC = np.zeros((4, 4, 4, 4, self.N, self.N, self.N))
        for i0 in [0,1,2,3]:
            for i1 in np.delete([0,1,2,3], i0):
                for i2 in np.delete([0,1,2,3], [i0, i1]):
                    for i3 in np.delete([0,1,2,3], [i0, i1, i2]):
                        LC[i0,i1,i2,i3,:,:,:] = float(((i1-i0)*(i2-i0)*(i3-i0)*(i2-i1)*(i3-i1)*(i3-i2))/(abs(i1-i0)*abs(i2-i0)*abs(i3-i0)*abs(i2-i1)*abs(i3-i1)*abs(i3-i2)))
        return LC*np.sqrt(abs(RRead.det4synch(gdown4)))
                         
if __name__ == "__main__":
    
    print("\n#################################")
    print("\n       Turn Around Radius        \n")
    print("#################################\n")
    
    
    OGpath = EBWeylparam.data_path
    all_hdf5it = list(np.arange(int(RRead.BASH('ls '+OGpath+'all_iterations/ | wc -l'))))
    funcs = MyClass(OGpath, all_hdf5it)
    it_counter = Value('i',0)
    pool = Pool(processes=1, initializer = init, initargs = (it_counter, ))
    
    RRead.MakeDir(OGpath+'EBWeyl')
    for it in all_hdf5it: funcs.calc_dtgdown(it)
    funcs.calc_dtdtgdown4()
    for var in ['dtgxx', 'dtgxy', 'dtgxz', 'dtgyy', 'dtgyz', 'dtgzz']:RRead.BASH('rm '+OGpath+var+'.hdf5')
    print('INFO: dtdtgdown4 calculated and recorded')
    
    for it in all_hdf5it: funcs.calc_dxgdown(it)
    funcs.calc_dtdxgdown4()
    for var in ['dxgxx', 'dxgxy', 'dxgxz', 'dxgyy', 'dxgyz', 'dxgzz']:RRead.BASH('rm '+OGpath+var+'.hdf5')
    print('INFO: dtdxgdown4 calculated and recorded')
    
    for it in all_hdf5it: funcs.calc_dygdown(it)
    funcs.calc_dtdygdown4()
    for var in ['dygxx', 'dygxy', 'dygxz', 'dygyy', 'dygyz', 'dygzz']:RRead.BASH('rm '+OGpath+var+'.hdf5')
    print('INFO: dtdygdown4 calculated and recorded')
    
    for it in all_hdf5it: funcs.calc_dzgdown(it)
    funcs.calc_dtdzgdown4()
    for var in ['dzgxx', 'dzgxy', 'dzgxz', 'dzgyy', 'dzgyz', 'dzgzz']:RRead.BASH('rm '+OGpath+var+'.hdf5')
    print('INFO: dtdzgdown4 calculated and recorded')
        
    for i in all_hdf5it[funcs.FD.mask*2:]: funcs.calc_Weyl(i)
    #pool.map(funcs.calc_Weyl,  all_hdf5it[4:])  #!!! dt FD here 2:-2 or 4:
    print('Done!') 
        
        
        
        
        
        
        
        
        
        
        
                         
                         
                         