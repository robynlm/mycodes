import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from data_analysis_codes.tools import LinData

def plotHamPert(p, it):
    Lin = LinData.LinData_Class(p, '/home/robynm/simulations/'+p.sim_name+'/')
    plt.figure()
    lw=3
    path = '/home/robynm/simulations/'+p.sim_name+'/DataSlice/Ham_contrib_it='
    f = pd.read_table(path+str(it)+'.csv', delimiter=',')
    x = np.arange(-p.L/2, p.L/2, p.dx)
    idx = np.argmin(abs(x))
    plt.semilogy(x, abs(f['HamPert']), linewidth=lw, color='C0', label=r'$|HamPert|$')
    #plt.semilogy(x, abs(f['Ham_ET']), linewidth=lw, color='C2', linestyle=':', label=r'$|Ham_{ET}|$')
    plt.semilogy(x, abs(f['RicciS']), linewidth=lw, color='C1', linestyle='--', label=r'$|{}^{(3)}R|$')
    #plt.semilogy(x, abs(f['Ricci_ET']), linewidth=lw, color='C1', linestyle=':', label=r'$|{}^{(3)}R_{ET}|$')
    plt.semilogy(x, abs(f['K2Pert']), linewidth=lw, color='C2', label=r'$|\frac{2}{3}(K^2-\bar{K}^2)|$')
    plt.semilogy(x, abs(f['A2']), linewidth=lw, label=r'$|-2A^2|$')
    plt.semilogy(x, abs(f['rhoPert']), linewidth=lw, color='C3', linestyle='--', label=r'$|-2\kappa(\rho-\bar{\rho})|$')

    plt.legend(bbox_to_anchor=(1,1))
    plt.grid()
    plt.xlabel('x [Mpc] (x=y=z)')
    plt.title('Hamiltonian Pert '+r'$a/a_{IN}=$'+' %.2f'%Lin.temp_from_temp('an', 'it', it))
    #plt.ylim(1e-15, 10000)
    
def plotHam(p, it):
    Lin = LinData.LinData_Class(p, '/home/robynm/simulations/'+p.sim_name+'/')
    plt.figure()
    lw=3
    path = '/home/robynm/simulations/'+p.sim_name+'/DataSlice/Ham_contrib_it='
    f = pd.read_table(path+str(it)+'.csv', delimiter=',')
    x = np.arange(-p.L/2, p.L/2, p.dx)
    plt.semilogy(x, abs(f['Ham']), linewidth=lw, label=r'$|Ham|$'+str(max(abs(f['Ham']))))
    plt.semilogy(x, abs(f['RicciS']), linewidth=lw, label=r'$|{}^{(3)}R|$'+str(max(abs(f['RicciS']))))
    plt.semilogy(x, abs(f['K2']), linewidth=lw, label=r'$|\frac{2}{3}K^2|$')
    plt.semilogy(x, abs(f['A2']), linewidth=lw, label=r'$|-2A^2|$')
    plt.semilogy(x, abs(f['rho']), linewidth=lw, linestyle='--', label=r'$|-2\kappa\rho|$')

    plt.legend(bbox_to_anchor=(1,1))
    plt.grid()
    plt.xlabel('x [Mpc] (x=y=z)')
    plt.title('Hamiltonian '+r'$a/a_{IN}=$'+' %.2f'%Lin.temp_from_temp('an', 'it', it))
    #plt.ylim(1e-15, 10000)
    
def plotRayPert(p, it):
    Lin = LinData.LinData_Class(p, '/home/robynm/simulations/'+p.sim_name+'/')
    plt.figure()
    lw=3
    path = '/home/robynm/simulations/'+p.sim_name+'/DataSlice/Ray_contrib_it='
    f = pd.read_table(path+str(it)+'.csv', delimiter=',')
    x = np.arange(-p.L/2, p.L/2, p.dx)
    plt.semilogy(x, abs(f['dKPert']), linewidth=lw, label=r'$|\dot{K}-\dot{\bar{K}}|$')
    plt.semilogy(x, abs(f['K2Pert']), linewidth=lw, label=r'$|\frac{1}{3}(K^2-\bar{K}^2)|$')
    plt.semilogy(x, abs(f['A2']), linewidth=lw, label=r'$|2A^2|$')
    plt.semilogy(x, abs(f['rhoPert']), linewidth=lw, label=r'$|\frac{\kappa}{2}(\rho-\bar{\rho})|$')

    plt.legend(bbox_to_anchor=(1,1))
    plt.grid()
    plt.xlabel('x [Mpc] (x=y=z)')
    plt.title('Raychaudhuri Pert '+r'$a/a_{IN}=$'+' %.2f'%Lin.temp_from_temp('an', 'it', it))
    plt.ylim(1e-15, 10000)
    
def plotRay(p, it):
    Lin = LinData.LinData_Class(p, '/home/robynm/simulations/'+p.sim_name)
    plt.figure()
    lw=3
    path = '/home/robynm/simulations/'+p.sim_name+'/DataSlice/Ray_contrib_it='
    f = pd.read_table(path+str(it)+'.csv', delimiter=',')
    x = np.arange(-p.L/2, p.L/2, p.dx)
    plt.semilogy(x, abs(f['dK']), linewidth=lw, label=r'$|\dot{K}|$')
    plt.semilogy(x, abs(f['K2']), linewidth=lw, label=r'$|\frac{1}{3}K^2|$')
    plt.semilogy(x, abs(f['A2']), linewidth=lw, label=r'$|2A^2|$')
    plt.semilogy(x, abs(f['rho']), linewidth=lw, label=r'$|\frac{\kappa}{2}\rho|$')

    plt.legend(bbox_to_anchor=(1,1))
    plt.grid()
    plt.xlabel('x [Mpc] (x=y=z)')
    plt.title('Raychaudhuri '+r'$a/a_{IN}=$'+' %.2f'%Lin.temp_from_temp('an', 'it', it))
    plt.ylim(1e-15, 10000)
    
def plotWeyl(p, it):
    Lin = LinData.LinData_Class(p, '/home/robynm/simulations/'+p.sim_name+'/')
    plt.figure()
    lw=3
    path = '/home/robynm/simulations/'+p.sim_name+'/DataSlice/WeylParts_it='
    f = pd.read_table(path+str(it)+'.csv', delimiter=',')
    x = np.arange(-p.L/2, p.L/2, p.dx)
    plt.semilogy(x, abs(f['E2']), linewidth=lw, label=r'$|E^2|$')
    plt.semilogy(x, abs(f['B2']), linewidth=lw, label=r'$|B^2|$')
    plt.semilogy(x, abs(f['L']), linewidth=lw, label=r'$|L|$')
    plt.semilogy(x, abs(f['M']), linewidth=lw, label=r'$|M|$')

    plt.legend(bbox_to_anchor=(1,1))
    plt.grid()
    plt.xlabel('x [Mpc] (x=y=z)')
    plt.title('Weyl Parts '+r'$a/a_{IN}=$'+' %.2f'%Lin.temp_from_temp('an', 'it', it))