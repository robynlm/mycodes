"""
Author: Robyn Munoz
Date: 9th June 2020

This Class finds the 3d index position of :
 - overdensity: findOD
 - underdensity: findUD
 - center: center
 of the simulation box.
 
 From the param.py file the following variables are needed:
 L, dx
"""

import numpy as np

class ODUDLocClass:
    def __init__(self, param, verbose=False):
        self.verbose = verbose
        self.param = param
        self.d1x = np.arange(param['xmin'], param['xmax'], param['dx'])
        self.d1y = np.arange(param['ymin'], param['ymax'], param['dy'])
        self.d1z = np.arange(param['zmin'], param['zmax'], param['dz'])
        self.d3x, self.d3y, self.d3z = np.meshgrid(self.d1x, self.d1y, 
                                                   self.d1z, indexing='ij')
    
    def findlocations(self):
        center = self.center()
        Amp_key = None
        for key in self.param.keys():
            if 'ICPertFLRW_Amp' in key:
                Amp_key = key
                break
            elif 'ICPertFLRW_GRH_Amp' in key:
                Amp_key = key
                break
        Rc_key = None
        for key in self.param.keys():
            if 'ICPertFLRW_Rcprofile' in key:
                Rc_key = key
                break
            elif 'ICPertFLRW_GRH_Rcprofile' in key:
                Rc_key = key
                break
        if Rc_key == None:
            Rc_key = 'Rcprofile'
            self.param[Rc_key] = 'sin'
                
        Vcenter = self.center()
        if Amp_key != None:
            if self.verbose: print(Amp_key, self.param[Amp_key], flush=True)
            if self.param[Amp_key] != 0:             
                if self.param[Rc_key]=='spin':
                    locs = [Vcenter]
                    loclab = ['_OD']
                else:
                    VOD, VUD = self.sinusoidal()
                    VmidOD = tuple([int(np.mean([VOD[i], Vcenter[i]])) 
                                    for i in range(3)])
                    VmidUD = tuple([int(np.mean([VUD[i], Vcenter[i]])) 
                                    for i in range(3)])
                    VF = (VOD[0], VOD[1], VUD[2])
                    locs = [VOD, VmidOD, Vcenter, VmidUD, VUD, VF]
                    loclab = ['_OD', '_midOD', '_cent', '_midUD', '_UD', '_F']
                return locs, loclab
            else:
                return [Vcenter], ['_cent']
        else:
            return [Vcenter], ['_cent']
        
    # find max and min of sinusoidal
    def sinusoidal(self):
        L_key = []
        for i in ['x', 'y', 'z']:
            for key in self.param.keys():
                if 'ICPertFLRW_lambda_'+i in key:
                    L_key += [key]
                    break
                elif 'ICPertFLRW_GRH_lambda_'+i in key:
                    L_key += [key]
                    break
            
        Rc = (np.sin(2 * np.pi * self.d3x / self.param[L_key[0]]) 
              + np.sin(2 * np.pi * self.d3y / self.param[L_key[1]]) 
              + np.sin(2 * np.pi * self.d3z / self.param[L_key[2]]))
        idxOD, idyOD, idzOD  = np.where(Rc == np.min(Rc))
        idxUD, idyUD, idzUD  = np.where(Rc == np.max(Rc))
        return (idxOD[0], idyOD[0], idzOD[0]), (idxUD[0], idyUD[0], idzUD[0])
    
    # Find the index clossest to the (0, 0, 0) position
    def center(self):
        return (np.argmin(abs(self.d1x)), 
                np.argmin(abs(self.d1y)), 
                np.argmin(abs(self.d1z)))
