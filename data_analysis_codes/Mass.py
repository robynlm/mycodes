import numpy as np
import pandas as pd
from tools import ChooseParam
from tools import LinData
from tools import MassCalc
from tools import ODUDLoc

class MyClass:
    def __init__(self, p, Lin):
        self.p = p
        self.Lin = Lin
        locfinder = ODUDLoc.ODUDLocClass(p)
        self.defined_locs = locfinder.findlocations()
        self.MC = MassCalc.MassCalcClass(p, Lin, self.defined_locs[0])
        
    def getvals(self, Radius):
        gdet, rho = self.ini_data()
        print('gdet and rho calculated')
        NEWROW = []
        for i, r in enumerate(Radius):
            NEWROW += [self.MC.getvar(rho, gdet, r)]
            print('Radius = {}, PercentDone = {:.2f}%'.format(r, 100*(i+1.0)/len(Radius)), flush=True)
        return NEWROW
    
    def ini_data(self):
        ta2 = 2*(self.Lin.a_initial**2)
        
        T1 = (self.Lin.a_initial**2) * ( 1 - 2 * self.Lin.Rc() )
        gxx = T1 - 2 * self.Lin.dxdxRc() * self.Lin.iFH2
        gyy = T1 - 2 * self.Lin.dydyRc() * self.Lin.iFH2
        gzz = T1 - 2 * self.Lin.dzdzRc() * self.Lin.iFH2
        del T1
        
        dxgxx = -ta2 * self.Lin.dxRc() - 2 * self.Lin.dxdxdxRc() * self.Lin.iFH2
        dxgyy = -ta2 * self.Lin.dxRc()
        dxgzz = -ta2 * self.Lin.dxRc()
        dxgdet = dxgxx*gyy*gzz+gxx*dxgyy*gzz+gxx*gyy*dxgzz
        del dxgyy, dxgzz
        
        dygxx = -ta2 * self.Lin.dyRc()
        dygyy = -ta2 * self.Lin.dyRc() - 2 * self.Lin.dydydyRc() * self.Lin.iFH2
        dygzz = -ta2 * self.Lin.dyRc()
        dygdet = dygxx*gyy*gzz+gxx*dygyy*gzz+gxx*gyy*dygzz
        del dygxx, dygzz
        
        dzgxx = -ta2 * self.Lin.dzRc()
        dzgyy = -ta2 * self.Lin.dzRc()
        dzgzz = -ta2 * self.Lin.dzRc() - 2 * self.Lin.dzdzdzRc() * self.Lin.iFH2
        dzgdet = dzgxx*gyy*gzz+gxx*dzgyy*gzz+gxx*gyy*dzgzz
        del dzgxx, dzgyy
        
        gdet = gxx*gyy*gzz
        
        RicciS = ta2 * (( self.Lin.dydyRc() + self.Lin.dzdzRc() ) * gxx + ( self.Lin.dxdxRc() + self.Lin.dzdzRc() ) * gyy + ( self.Lin.dxdxRc() + self.Lin.dydyRc() ) * gzz) / gdet - 10 * (self.Lin.a_initial**4) * ( self.Lin.dxRc()**2 + self.Lin.dyRc()**2 + self.Lin.dzRc()**2 ) / gdet - (self.Lin.a_initial**2) * ( gxx**2 * ( gyy*(ta2*self.Lin.dyRc()**2 - dzgzz*self.Lin.dzRc()) + gzz*(ta2*self.Lin.dzRc()**2 - dygyy*self.Lin.dyRc())) + gyy**2 * ( gxx*(ta2*self.Lin.dxRc()**2 - dzgzz*self.Lin.dzRc()) + gzz*(ta2*self.Lin.dzRc()**2 - dxgxx*self.Lin.dxRc())) + gzz**2 * ( gxx*(ta2*self.Lin.dxRc()**2 - dygyy*self.Lin.dyRc()) + gyy*(ta2*self.Lin.dyRc()**2 - dxgxx*self.Lin.dxRc())) + 2 * ( dxgdet*self.Lin.dxRc()*(gyy+gzz) + dygdet*self.Lin.dyRc()*(gxx+gzz) + dzgdet*self.Lin.dzRc()*(gxx+gyy))) / gdet
        del ta2, dxgxx, dxgdet, dygyy, dygdet, dzgzz, dzgdet
        
        gxxu = 1 / gxx
        gyyu = 1 / gyy
        gzzu = 1 / gzz
        del gxx, gyy, gzz
        
        T2 = - (self.Lin.a_initial**2) * self.Lin.Hprop_initial * ( 1 - 2 * self.Lin.Rc() )
        kxx = T2 + ( 2 + self.Lin.fL_initial ) * self.Lin.dxdxRc() * self.Lin.iFH
        kyy = T2 + ( 2 + self.Lin.fL_initial ) * self.Lin.dydyRc() * self.Lin.iFH
        kzz = T2 + ( 2 + self.Lin.fL_initial ) * self.Lin.dzdzRc() * self.Lin.iFH
        del T2
        
        K = gxxu*kxx + gyyu*kyy + gzzu*kzz
        KijKji = gxxu*kxx*gxxu*kxx + gyyu*kyy*gyyu*kyy + gzzu*kzz*gzzu*kzz
        del kxx, kyy, kzz, gxxu, gyyu, gzzu
        
        rho = (RicciS + K**2 - KijKji - 2*self.Lin.Lambda) / (2 * self.Lin.kappa)
        del RicciS, K, KijKji
        return gdet, rho
        
if __name__ == "__main__":
    
    print("\n#################################")
    print("\n       Turn Around Radius        \n")
    print("#################################\n")
    
    sim, path, OGpath = ChooseParam.ChooseP('pflrw_d3e2_L1821_t1_N256_LCDM', 'home')
    Lin = LinData.LinData_Class(sim, path)
    
    funcs = MyClass(sim, Lin)
    Radius = np.array([332.90181249999995]) #np.linspace(np.sqrt(3)*sim.dx/2+0.01*sim.dx, sim.L/4, 20)
    data_values = funcs.getvals(Radius)
    data_header = ['Radius', 'Mass', 'Mass_MSun', 'delta_average']
    print(data_header)
    print(data_values)
    #pd.DataFrame(data_values).to_csv(path+'mass_in_sphere.csv', header=data_header, index=False)
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
            
