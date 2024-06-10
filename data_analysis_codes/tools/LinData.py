import numpy as np
import sys
import pandas as pd
from . import ReadingTools as RRead
from . import LCDM
from . import EdS

class LinData_Class:
    def __init__(self, param):
        self.param = param
        self.t_initial = param['cctk_initial_time']
        self.delta_initial = {}
        try:
            path = (param['HorSpath'] + param['simname']
                    + '/output-0000/' + param['simname'] + '/')
            self.temporal_file = RRead.get_temporal_file(path)
            print(' read time file')
        except:
            self.temporal_file = pd.DataFrame({'A' : []})
            print(' could not read time file')
                
        if param['ICPertFLRW_Lambda']=='no':
            if 'ICPertFLRW_GRH_type_of_matter' in self.param.keys():
                self.evo = EdS.evo(matter = param['ICPertFLRW_GRH_type_of_matter'])
            elif 'ICPertFLRW_type_of_matter' in self.param.keys():
                self.evo = EdS.evo(matter = param['ICPertFLRW_type_of_matter'])
            else:
                self.evo = EdS.evo()
        else:
            self.evo = LCDM.evo()

        self.Amp_key = None
        for key in self.param.keys():
            if 'ICPertFLRW_Amp' in key:
                self.Amp_key = key
                break
            elif 'ICPertFLRW_GRH_Amp' in key:
                self.Amp_key = key
                break
        if self.Amp_key == None:
            self.Amp_key = 'ICPertFLRW_Amp'
            self.param['ICPertFLRW_Amp'] = 0.0
        
        speed_of_light = 299792458   # m.s^{-1}
        Grav_const = 6.67408e-11 # m^3.kg^{-1}.s^{-2}
        parsec = 3.0857e16   # m
        Megaparsecc = parsec*1e6  # m
        MassSun = 1.98847e30  # kg
        self.Massfac = ((Megaparsecc * speed_of_light**2) 
                        / (Grav_const * MassSun * self.evo.a_today**3))
        
        self.k_pert = 2 * np.pi / param['ICPertFLRW_lambda_pert1']
        self.d1x = np.arange(param['xmin'], param['xmax'], param['dx'])
        self.d1y = np.arange(param['ymin'], param['ymax'], param['dy'])
        self.d1z = np.arange(param['zmin'], param['zmax'], param['dz'])
        self.d3x, self.d3y, self.d3z = np.meshgrid(
            self.d1x, self.d1y, self.d1z, indexing='ij')
        
        # rescale dt
        if 'ICPertFLRW_time' in self.param.keys():
            if self.param['ICPertFLRW_time']=='proper':
                self.param['dt'] *= self.evo.a(self.param['cctk_initial_time'])
        elif 'ICPertFLRW_GRH_time' in self.param.keys():
            if self.param['ICPertFLRW_GRH_time']=='proper':
                self.param['dt'] *= self.evo.a(self.param['cctk_initial_time'])
        else:
            self.param['dt'] *= self.evo.a(self.param['cctk_initial_time'])
    
    def temp_from_temp(self, var_wanted_str, 
                       var_known_str, var_known_float):
        idx = np.argmin(abs(self.temporal_file[var_known_str] 
                            - var_known_float))
        return self.temporal_file[var_wanted_str][idx]            
        
    def F(self, torh5):
        return self.evo.fL(torh5) + (3/2)*self.evo.Omega_m(torh5)
        
    def xyz(self, loc):
        if type(loc)==str:
            loc_options = {'OD': - self.param['Lx'] / 4, 
                           'midOD': - self.param['Lx'] / 8, 
                           'cent':0.0, 
                           'midUD': self.param['Lx'] / 8, 
                           'UD': self.param['Lx'] / 4}
            return loc_options[loc]
        elif type(loc)==float:
            return loc

    ###############################
    # Background functions
    ###############################
    
    def gdet(self, torh5):
        return self.evo.a(torh5)**6   
    
    def an_initial(self, torh5):
        return self.evo.a(torh5)/self.evo.a(self.t_initial)
    
    def gdown(self, torh5):
        return np.identity(3) * self.evo.a(torh5)**2
    
    def gup(self, torh5):
        return np.identity(3) / self.evo.a(torh5)**2
        
    def Theta(self, torh5):
        return 3 * self.evo.Hprop(torh5)
    
    def rho(self, torh5):
        return self.evo.rho(torh5)
    def rho_u(self, torh5):
        return self.evo.rho(torh5)

    ###############################
    # 1st order functions
    ###############################
    
    def xyz(self, loc):
        if type(loc)==str:
            loc_options = {'OD' : - self.param['Lx'] / 4, 
                           'midOD' : - self.param['Lx'] / 8, 
                           'cent' : 0.0, 
                           'midUD' : self.param['Lx'] / 8, 
                           'UD' : self.param['Lx'] / 4}
            return loc_options[loc]
        elif type(loc)==float:
            return loc

    def Rc(self, loc):
        loc = loc.replace('_', '')
        UOD = {'OD': -1, 'midOD': -0.5, 'cent': 0, 'midUD': 0.5, 'UD': 1,
               'av': 0, 'L1': 0, 'var': 0.5, 'max': -1, 'min': 1}
        return 3 * UOD[loc] * self.param[self.Amp_key]
    
    def fac_delta(self, torh5):
        return -(self.k_pert**2)/((self.evo.a(torh5)**2) 
                                  * self.F(torh5) 
                                  *(self.evo.Hprop(torh5)**2))

    def delta_u(self, torh5, loc):
        try:
            return (self.delta_initial[loc] * self.evo.a(torh5) 
                    / self.evo.a(self.t_initial))
        except:
            print('WARNING: delta calculated at 1st order from Rc')
            print('You should define self.delta_initial[loc]')
            return self.fac_delta(torh5)*self.Rc(loc)
    
    
    def dTheta(self, torh5, loc):
        return - self.evo.fL(torh5) * self.delta_u(torh5, loc) / 3

    def drho(self, torh5, loc):
        return self.delta_u(torh5, loc)
    def drho_u(self, torh5, loc):
        return self.delta_u(torh5, loc)
    
    def dgdet(self, torh5, loc):
        return - 6 * ( (self.delta_u(torh5, loc) / 3) 
                      + self.Rc(loc))
    
    def dRicciS(self, torh5, loc):
        return (-4 * (self.k_pert**2) * self.Rc(loc) 
                / (self.evo.a(torh5)**2))
    
    def dRicciSconf(self, torh5, loc):
        return np.array([-4 * (self.k_pert**2) * self.Rc(loc)] 
                        * len(torh5))
