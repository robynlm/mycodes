import h5py
import sys
import numpy as np
import pandas as pd
from tools import LinData
from sphereint import sphereint
from tools import GetVars_Plot2d as GVar
from tools import ReadingTools as RRead

class MyClass:
    def __init__(self, param, Lin):
        self.param = param
        self.cell_vol = param['dx']**3
        self.Lin = Lin
        self.get_var = GVar.Get_var(param, self.Lin)
        self.it_file_name = param['h5datapath'] + param['simname']
        
        self.L = param['Lx']
        self.N = param['Nx']
        self.dx = param['dx']
        #self.tNd4 = int(3 * self.N / 4)
        centerloc = [int(self.N/4), int(self.N/4), int(self.N/4)]
        self.MC = sphereint.MassCalcClass(self.N, self.L, centerloc)
        
    #def shift(self, phi):
    #    phi = np.append(phi[self.tNd4:, :, :], phi[:self.tNd4, :, :], axis=0)
    #    phi = np.append(phi[:, self.tNd4:, :], phi[:, :self.tNd4, :], axis=1)
    #    phi = np.append(phi[:, :, self.tNd4:], phi[:, :, :self.tNd4], axis=2)
    #    return phi
    
    #def in_r(self, x, y, z, r):
    #    return np.sqrt(x**2 + y**2 + z**2) <= r
    
    #def get_weights(self, r):
    #    if r >= self.L/2:
    #        print('ERROR: that radius is too big for me sorry')
    #    else:
    #        weights = np.zeros((self.N, self.N, self.N))
    #        for ix in range(self.N):
    #            x = self.Lin.d3xyz[ix]
    #            for iy in range(self.N):
    #                y = self.Lin.d3xyz[iy]
    #                for iz in range(self.N):
    #                    z = self.Lin.d3xyz[iz]
    #                    in_rms = self.in_r(x, y, z, r - self.dx*np.sqrt(3))
    #                    #in_rps = self.in_r(x, y, z, r + self.dx*np.sqrt(3))
    #                    if in_rms: 
    #                        weights[ix,iy,iz] = 1.0
    #                    #elif in_rps:
    #                    #    dxmax = self.dx/2
    #                    #    p_in_r = [self.in_r(x+dxmax, y+dxmax, z+dxmax, r), 
    #                    #              self.in_r(x-dxmax, y+dxmax, z+dxmax, r), 
    #                    #              self.in_r(x+dxmax, y-dxmax, z+dxmax, r), 
    #                    #              self.in_r(x-dxmax, y-dxmax, z+dxmax, r), 
    #                    #              self.in_r(x+dxmax, y-dxmax, z-dxmax, r), 
    #                    #              self.in_r(x-dxmax, y-dxmax, z-dxmax, r), 
    #                    #              self.in_r(x+dxmax, y+dxmax, z-dxmax, r), 
    #                    #              self.in_r(x-dxmax, y+dxmax, z-dxmax, r)]
    #                    
    #                    #    if np.sum(p_in_r) == 8:
    #                    #        weights[ix,iy,iz] =  1.0
    #                    #    else:
    #                    #        weights[ix,iy,iz] =  0.0
    #                    else:
    #                        weights[ix,iy,iz] =  0.0
    #    return weights
        
    def getvals(self, Radius, all_hdf5it):
        data = []
        for r in Radius:
            weight = self.MC.get_box_weights(r)
            for i, it in enumerate(all_hdf5it):
                f = h5py.File('{}_it_{:06d}.hdf5'.format(self.it_file_name, 
                                                         int(it)), 'r')
                t = self.Lin.temp_from_temp('t', 'it', it)
            
                rho = self.get_var.get_the_rho(f, int(it))['rho']
                gdet = self.get_var.get_the_metric(f, int(it))['gdet']
        
                Mass_prop = np.sum(rho * weight * np.sqrt(gdet)) * self.cell_vol
                Mass = np.sum(rho * weight) * self.cell_vol
        
                if r == Radius[0]:
                    data += [[it, t, Mass_prop, Mass]]
                else:
                    data[i] += [Mass_prop, Mass]
            print('evo = {:.2f}'.format(100 
                                        * np.argmin(abs(np.array(Radius) - r)) 
                                        / len(Radius)), flush=True)
        return data
    
if __name__ == "__main__":
    
    print("\n#################################")
    print("\n         Mass in Sphere          \n")
    print("#################################\n")
    
    simname = sys.argv[1]
    param = RRead.read_parameters(simname)
    Lin = LinData.LinData_Class(param)
    funcs = MyClass(param, Lin)
    
    all_hdf5it = Lin.temporal_file['it'][0::param['IOHDF5::out_every']]
    nbrit = int(all_hdf5it[-1] / param['IOHDF5::out_every'])
    print(' -- number of iterations : ', nbrit, flush=True)
    
    radius_denominator = np.append(np.arange(3, 11, 1), np.arange(20, 128, 10))
    Radius = [param['Lx'] / i for i in radius_denominator]
    data_values = funcs.getvals(Radius, [0])
    
    data_header = ['it', 't'] + list(np.ravel([[var + str(i) 
                                                for var in ['Mass_prop_', 'Mass_']] 
                                               for i in radius_denominator]))
    pd.DataFrame(data_values).to_csv(param['datapath'] + 'new_mass_radius_evo.csv', 
                                     header = data_header, index = False)
    