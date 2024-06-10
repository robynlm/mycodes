"""
Author: Robyn Munoz
Date: 9th June 2020

Set of useful functions, see their respective descriptions.

"""

import numpy as np
import pandas as pd
import subprocess
import sys
import re
from . import LCDM
from . import EdS
import matplotlib.pyplot as plt
from . import Cstyle
from . import NumMethods
import time

" Submit bash command line "
def BASH(command):
    try:
        results = subprocess.check_output(command,shell=True)
        strresults = str(results, 'utf-8').strip()
    except:
        strresults = "ERROR: RRead.BASH('"+command+"') did not work"
    return strresults

def time_to_str(amount_of_time_sec):
    mins, sec = divmod(amount_of_time_sec, 60)
    hours, endmins = divmod(mins, 60)
    days, endhours = divmod(hours, 24)
    time_str = f"{int(days):02}d {int(endhours):02}h {int(endmins):02}m {sec:03.1f}s"
    return time_str

def progressbar(it):
    size = 30
    count = len(it)
    start = time.time() # time estimate start
    def show(j):
        x = int(size*j/count)
        # time estimate calculation and string
        remaining = ((time.time() - start) / j) * (count - j)
        percentage = j*100/count
        print(f"\r[{u'█'*x}{('.'*(size-x))}] {percentage:0.2f}% Est wait {time_to_str(remaining)}", 
              end="", file=sys.stdout, flush=True)
    show(0.1) # avoid div/0 
    for i, item in enumerate(it):
        yield item
        show(i+1)
    print("\n", flush=True, file=sys.stdout)
    print(f"total time used: {time_to_str(time.time() - start)}", flush=True, file=sys.stdout)
    print("\n", flush=True, file=sys.stdout)

def read_parameters(simname):
    # Find parameter file
    HorSloc = '/home/robynm/simulations/'
    path = HorSloc + simname + '/output-0000/' + simname + '.par'
    data_file = BASH('cat ' + path)
    if 'ERROR' in data_file:
        HorSloc = '/Users/robynmunoz/simulations/'
        path = HorSloc + simname + '/output-0000/' + simname + '.par'
        data_file = BASH('cat ' + path)
        if 'ERROR' in data_file:
            HorSloc = '/mnt/lustre2/ET_sims/'
            path = HorSloc + simname + '/output-0000/' + simname + '.par'
            data_file = BASH('cat ' + path)
            if 'ERROR' in data_file:
                print("Error in loading parameter file, check that it's in")
                print('/home/robynm/simulations/' + simname + '/output-0000/' 
                      + simname + '.par')
                print('Or')
                print('/Users/robynmunoz/simulations/' + simname + '/output-0000/'
                      + simname + '.par')
                print('Or')
                print('/mnt/lustre2/ET_sims/' + simname + '/output-0000/' 
                      + simname + '.par')
        
    # save paths of files
    parameters = {'simname':simname, 'HorSpath':HorSloc}
    parameters['datapath'] = (parameters['HorSpath'] + parameters['simname'] 
                                + '/output-0000/' + parameters['simname'] + '/')
    parameters['h5datapath'] = parameters['datapath'] + 'all_iterations/'
                
    # number of restarts
    files = BASH('ls ' + parameters['HorSpath'] 
                       + parameters['simname']).split('\n')
    parameters['nbr_restarts'] = len([fl for fl in files 
                                      if 'output' in fl and 'active' not in fl])
    
    # save all parameters
    for line in data_file.split('\n'):
        if '::' in line and '=' in line:
            thorn, variable, value = re.split('::|=', line)
            thorn = thorn.strip()
            variable = variable.strip()
            value = value.strip()
            
            # format value, float, int, string
            digitvalue = value.replace('.', '', 1).replace(
                '-', '', 1).replace('+', '', 1).replace('e', '', 1)
            if digitvalue.isdigit():
                if '.' in value or 'e' in value:
                    value = float(value)
                else:
                    value = int(value)
            else:
                if '"' in value:
                    value = value.split('"')[1]
                    
            if variable in ['out_every', 'one_file_per_group']:
                parameters[thorn + '::' + variable] = value
            else:
                parameters[variable] = value
                
    # calculate extra grid info
    for coord in ['x', 'y', 'z']:
        L = parameters[coord+'max']-parameters[coord+'min']
        N = int(L/parameters['d'+coord])
        parameters['L'+coord] = L
        parameters['N'+coord] = N
    if 'ICPertFLRW_lambda_pert1' not in parameters.keys():
        parameters['ICPertFLRW_lambda_pert1'] = parameters['Lx']
    
    if 'dtfac' not in parameters.keys():
        parameters['dtfac'] = 0.5
    parameters['dt'] = parameters['dtfac'] * parameters['dx']
        
    # is cosmological constant present
    if ('ICPertFLRW_Lambda' not in parameters.keys()):
        if ('ICPertFLRW_expansion' in parameters.keys()):
            if (parameters['ICPertFLRW_expansion'] == 'EdS'):
                parameters['ICPertFLRW_Lambda'] = 'no'
            else:
                parameters['ICPertFLRW_Lambda'] = 'yes'
        elif ('initial_hydro' in parameters.keys()):
            parameters['ICPertFLRW_Lambda'] = 'no'
        else:
            parameters['ICPertFLRW_Lambda'] = 'yes'
        
    # is it synchronous
    if ('harmonicF' in parameters.keys()):
        if parameters['harmonicF'] == 0.0:
            parameters['synchronous'] = True
        else:
            parameters['synchronous'] = False
    elif ('alphaF' in parameters.keys()):
        if parameters['alphaF'] == 0.0:
            parameters['synchronous'] = True
        else:
            parameters['synchronous'] = False
    elif ('alpha_func_of_rho' in parameters.keys()):
        if parameters['alpha_func_of_rho'] == 'no':
            parameters['synchronous'] = True
        else:
            parameters['synchronous'] = False
    else:
        print('I assume this is synchronous')
        parameters['synchronous'] = True
        
    return parameters

" Create a new directory "
def MakeDir(folder):
    try:
        BASH('mkdir '+folder)
        print("Created the folder: "+folder)
    except:
        print(folder+" already exists")

" This function extracts the components of h5py variable "
def read_xyz(f, keys, fixindexes=True):
    if fixindexes:
        xx = fixij(f[keys[0]])
        xy = fixij(f[keys[1]])
        xz = fixij(f[keys[2]])
        yy = fixij(f[keys[3]])
        yz = fixij(f[keys[4]])
        zz = fixij(f[keys[5]])
    else:
        xx = np.array(f[keys[0]])
        xy = np.array(f[keys[1]])
        xz = np.array(f[keys[2]])
        yy = np.array(f[keys[3]])
        yz = np.array(f[keys[4]])
        zz = np.array(f[keys[5]])
    phi = np.array([[xx, xy, xz],[xy, yy, yz],[xz, yz, zz]])
    return [xx, xy, xz, yy, yz, zz], phi

def fixij(f): # I do not need this for non simulations
    return np.transpose(np.array(f), (2, 1, 0))

" Set of functions to cut border around data array, eg: Ghost cells "
def cut0(data, i, N):
    if np.shape(data)[-1]==N:
        return data
    else:
        return data[i:-i, i:-i, i:-i]
def cut1(data, i, N):
    if np.shape(data)[-1]==N:
        return data
    else:
        return data[:, i:-i, i:-i, i:-i]
def cut2(data, i, N):
    if np.shape(data)[-1]==N:
        return data
    else:
        return data[:, :, i:-i, i:-i, i:-i]

def safe_division(a, b):
    if isinstance(a, float):
        c = 0.0 if b==0 else a/b
    else:
        c = np.divide(a, b, out=np.zeros_like(a), where=b!=0)
        c[np.where(np.logical_and(a==0, b==0))] = 1.0
    return c

# Simply calculate the norm of a given vector
def norm(f):
    x, y, z = f
    return np.sqrt(x**2+y**2+z**2)

def getcomponents3(f):
    if isinstance(f, list):
        xx, xy, xz, yy, yz, zz = f
    elif isinstance(f, np.ndarray):
        xx = f[0,0]
        xy = f[0,1]
        xz = f[0,2]
        yy = f[1,1]
        yz = f[1,2]
        zz = f[2,2]
    else:
        print("getcomponents3 doesn't understand f format")
    return xx, xy, xz, yy, yz, zz

def getcomponents4(f):
    if isinstance(f, list):
        tt, tx, ty, tz, xx, xy, xz, yy, yz, zz = f
    elif isinstance(f, np.ndarray):
        tt = f[0,0]
        tx = f[0,1]
        ty = f[0,2]
        tz = f[0,3]
        xx = f[1,1]
        xy = f[1,2]
        xz = f[1,3]
        yy = f[2,2]
        yz = f[2,3]
        zz = f[3,3]
    else:
        print("getcomponents4synch doesn't understand f format")
    return tt, tx, ty, tz, xx, xy, xz, yy, yz, zz
    
def det3(f):
    xx, xy, xz, yy, yz, zz = getcomponents3(f)
    return -xz*xz*yy + 2*xy*xz*yz - xx*yz*yz - xy*xy*zz + xx*yy*zz

def inv3(f):
    xx, xy, xz, yy, yz, zz = getcomponents3(f)
    fup = np.array([[   yy*zz-yz*yz, -(xy*zz-yz*xz),  xy*yz-yy*xz], 
                    [ -(xy*zz-xz*yz),  xx*zz-xz*xz, -(xx*yz-xy*xz)],
                    [   xy*yz-xz*yy, -(xx*yz-xz*xy),  xx*yy-xy*xy]])
    return fup / det3(f)

def det4(f):
    tt, tx, ty, tz, xx, xy, xz, yy, yz, zz = getcomponents4(f)
    det = (tz*tz*xy*xy - 2*ty*tz*xy*xz + ty*ty*xz*xz 
           - tz*tz*xx*yy + 2*tx*tz*xz*yy - tt*xz*xz*yy 
           + 2*ty*tz*xx*yz - 2*tx*tz*xy*yz - 2*tx*ty*xz*yz 
           + 2*tt*xy*xz*yz + tx*tx*yz*yz - tt*xx*yz*yz 
           - ty*ty*xx*zz + 2*tx*ty*xy*zz - tt*xy*xy*zz 
           - tx*tx*yy*zz + tt*xx*yy*zz)
    return det

def inv4(f):
    tt, tx, ty, tz, xx, xy, xz, yy, yz, zz = getcomponents4(f)
    fup = np.array([
        [-xz*xz*yy + 2*xy*xz*yz - xx*yz*yz - xy*xy*zz + xx*yy*zz, 
         tz*xz*yy - tz*xy*yz - ty*xz*yz + tx*yz*yz + ty*xy*zz - tx*yy*zz, 
         -tz*xy*xz + ty*xz*xz + tz*xx*yz - tx*xz*yz - ty*xx*zz + tx*xy*zz, 
         tz*xy*xy - ty*xy*xz - tz*xx*yy + tx*xz*yy + ty*xx*yz - tx*xy*yz], 
        [tz*xz*yy - tz*xy*yz - ty*xz*yz + tx*yz*yz + ty*xy*zz - tx*yy*zz, 
         -tz*tz*yy + 2*ty*tz*yz - tt*yz*yz - ty*ty*zz + tt*yy*zz, 
         tz*tz*xy - ty*tz*xz - tx*tz*yz + tt*xz*yz + tx*ty*zz - tt*xy*zz, 
         -ty*tz*xy + ty*ty*xz + tx*tz*yy - tt*xz*yy - tx*ty*yz + tt*xy*yz], 
        [-tz*xy*xz + ty*xz*xz + tz*xx*yz - tx*xz*yz - ty*xx*zz + tx*xy*zz, 
         tz*tz*xy - ty*tz*xz - tx*tz*yz + tt*xz*yz + tx*ty*zz - tt*xy*zz, 
         -tz*tz*xx + 2*tx*tz*xz - tt*xz*xz - tx*tx*zz + tt*xx*zz, 
         ty*tz*xx - tx*tz*xy - tx*ty*xz + tt*xy*xz + tx*tx*yz - tt*xx*yz], 
        [tz*xy*xy - ty*xy*xz - tz*xx*yy + tx*xz*yy + ty*xx*yz - tx*xy*yz, 
         -ty*tz*xy + ty*ty*xz + tx*tz*yy - tt*xz*yy - tx*ty*yz + tt*xy*yz, 
         ty*tz*xx - tx*tz*xy - tx*ty*xz + tt*xy*xz + tx*tx*yz - tt*xx*yz, 
         -ty*ty*xx + 2*tx*ty*xy - tt*xy*xy - tx*tx*yy + tt*xx*yy]])
    return fup / det4(f)

def normalise(f):
    return f/f[0]

def readasctable(name, colnames):
    try:
        data = pd.read_csv(name, header=6, names=colnames, delimiter=' ')
        return data
    except:
        try:
            data = pd.read_csv(name, header=9, names=colnames, delimiter=' ')
            return data
        except:
            print('Could not read table:', name)
            sys.exit()
        
def readasctable_scalar(name):
    colnames = ['iteration', 'time', 'data', 'nan']
    data = readasctable(name, colnames)
    it = np.array(data['iteration'])
    t = np.array(data['time'])
    S = np.array(data['data'])
    return it, t, S

def readasctable_vector(name):
    colnames = ['iteration', 'time', 'x', 'y', 'z', 'nan']
    data = readasctable(name, colnames)
    it = np.array(data['iteration'])
    t = np.array(data['time'])
    V = np.array([data['x'], data['y'], data['z']])
    return it, t, V

def readasctable_tensor(name):
    colnames = ['iteration', 'time', 'xx', 'xy', 'xz', 'yy', 'yz', 'zz', 'nan']
    data = readasctable(name, colnames)
    it = np.array(data['iteration'])
    t = np.array(data['time'])
    xx = np.array(data['xx'])
    xy = np.array(data['xy'])
    xz = np.array(data['xz'])
    yy = np.array(data['yy'])
    yz = np.array(data['yz'])
    zz = np.array(data['zz'])
    T = [xx, xy, xz, yy, yz, zz]
    Tdown = np.array([[xx, xy, xz],
                      [xy, yy, yz],
                      [xz, yz, zz]])
    return it, t, T, Tdown

def collect_iteration_and_time(param):
    files = re.split('\n', BASH('ls '+param['datapath']+'*.asc'))
    dic = {'admbase-curv':'T', 'admbase-metric':'T', 'admbase-shift':'V',
           'admbase-lapse':'S', 'ct_dust-ct_rho':'S', 'grid_coordinates':'No',
           'grid_structure':'No', 'ml_bssn-ml_ham':'S', 'ml_bssn-ml_mom':'V'}
            
    for file in files:
        try:
            var_type = dic[re.split(r"[/.]", file)[-3]]
        except:
            var_type = dic[re.split(r"[/.]", file)[-2]]

        if var_type=='T':
            it, t, dum, dummy = readasctable_tensor(file)
            del dum, dummy
            break
        elif var_type=='V':
            it, t, dum = readasctable_vector(file)
            del dum
            break
        elif var_type=='S':
            it, t, dum = readasctable_scalar(file)
            del dum
            break
        else: # var_type=='No':
            pass # This file did not record time
    try:
        return np.array(it), np.array(t)
    except:
        output = BASH('ls '+param['h5datapath'])
        last_it = int(output.split(param['simname']+'_it_')[-1].split('.hdf5')[0])
        all_it = np.arange(last_it+1)
        t = [param['cctk_initial_time']]
        for it in range(last_it):
            t += [t[-1] + param['dt']]
        all_t = np.array(t)
        return all_it, all_t

def collect_h5iteration(param):
    output = BASH('ls '+param['h5datapath'])
    output = output.split('\n')
    all_it = []
    for o in output:
        if param['simname']+'_it_' in o:
            itstr = o.replace(param['simname'] 
                              + '_it_', '').replace('.hdf5', '')
            all_it += [int(itstr)]
    return all_it

def get_temporal_file(path):
    filename = re.split('\n', BASH('ls '+path+'Time_*.csv'))[0]
    df = pd.read_csv(filename)
    return dict(zip(df.columns, df.values.T))


def plot_constraintsEScale(fc):
    plt.style.use(Cstyle.style1)
    plt.figure(figsize=(7, 10))
    labsize = 25
    alphafac = 0.1
    locs = ['maxabs', 'highQabs', 'medianabs', 'lowQabs', 'minabs']
    loclabs = ['max', 'first quartile', 'median', 'third quartile', 'min']
    cols = plt.cm.viridis(np.linspace(0,1,len(locs)))
    it = fc['it']

    ax = plt.subplot(211)
    for i, loc, in enumerate(locs):
        for n in range(3):
            varstr = 'Mom'+str(n+1)+'/MomEScale_'+loc
        var = np.array(fc[varstr])
        if n==0:
            ax.semilogy(it, var, color=cols[i], linestyle='-', label=loclabs[i])
        else:
            ax.semilogy(it, var, color=cols[i], linestyle='-')
    ax.grid()
    ax.legend(bbox_to_anchor=(1,1))
    ax.set_xticklabels([])
    ax.set_ylabel(r'$|\mathcal{M}_{i}/[\mathcal{M}]|$', fontsize=labsize)


    ax = plt.subplot(212)
    for i, loc, in enumerate(locs):
        varstr = 'Ham/HamEScale_'+loc
        var = np.array(fc[varstr])
        ax.semilogy(it, var, color=cols[i], linestyle='-')
    ax.grid()
    ax.set_xlabel(r'$it$', fontsize=labsize)
    ax.set_ylabel(r'$|\mathcal{H}/[\mathcal{H}]|$', fontsize=labsize)

    plt.subplots_adjust(hspace=0)
    
def plot_constraints_with_error(fc32, fc64, fc128):
    plt.style.use(Cstyle.style1)
    plt.figure(figsize=(7, 10))
    labsize = 25
    alphafac = 0.1
    locs = ['maxabs', 'highQabs', 'medianabs', 'lowQabs', 'minabs']
    loclabs = ['max', 'first quartile', 'median', 'third quartile', 'min']
    cols = plt.cm.viridis(np.linspace(0,1,len(locs)))
    it = fc128['it']

    ax = plt.subplot(211)
    for i, loc, in enumerate(locs):
        for n in range(3):
            varstr = 'Mom'+str(n+1)+'/MomEScale_'+loc
        var = np.array(fc128[varstr])
        err, c = NumMethods.get_error(np.array(fc32[varstr]), 
                                   np.array(fc64[varstr]), var)
        if n==0:
            ax.semilogy(it, var, color=cols[i], linestyle='-', label=loclabs[i])
        else:
            ax.semilogy(it, var, color=cols[i], linestyle='-')
        plt.fill_between(it, var-err, var+err, color=cols[i], alpha=alphafac)
    ax.grid()
    ax.legend(bbox_to_anchor=(1,1))
    ax.set_xticklabels([])
    ax.set_ylabel(r'$|\mathcal{M}_{i}/[\mathcal{M}]|$', fontsize=labsize)


    ax = plt.subplot(212)
    for i, loc, in enumerate(locs):
        varstr = 'Ham/HamEScale_'+loc
        var = np.array(fc128[varstr])
        err, c = NumMethods.get_error(np.array(fc32[varstr]), 
                                   np.array(fc64[varstr]), var)
        ax.semilogy(it, var, color=cols[i], linestyle='-')
        plt.fill_between(it, var-err, var+err, color=cols[i], alpha=alphafac)
    ax.grid()
    ax.set_xlabel(r'$it$', fontsize=labsize)
    ax.set_ylabel(r'$|\mathcal{H}/[\mathcal{H}]|$', fontsize=labsize)

    plt.subplots_adjust(hspace=0)
    
def plot_constraints(fc):
    plt.style.use(Cstyle.style1)
    plt.figure(figsize=(7, 10))
    labsize = 25
    alphafac = 0.1
    locs = ['maxabs', 'highQabs', 'medianabs', 'lowQabs', 'minabs']
    loclabs = ['max', 'first quartile', 'median', 'third quartile', 'min']
    cols = plt.cm.viridis(np.linspace(0,1,len(locs)))
    it = fc['it']

    ax = plt.subplot(211)
    for i, loc, in enumerate(locs):
        for n in range(3):
            varstr = 'Mom'+str(n+1)+'_'+loc
        var = np.array(fc[varstr])
        if n==0:
            ax.semilogy(it, var, color=cols[i], linestyle='-', label=loclabs[i])
        else:
            ax.semilogy(it, var, color=cols[i], linestyle='-')
    ax.grid()
    ax.legend(bbox_to_anchor=(1,1))
    ax.set_xticklabels([])
    ax.set_ylabel(r'$|\mathcal{M}_{i}|$', fontsize=labsize)


    ax = plt.subplot(212)
    for i, loc, in enumerate(locs):
        varstr = 'Ham_'+loc
        var = np.array(fc[varstr])
        ax.semilogy(it, var, color=cols[i], linestyle='-')
    ax.grid()
    ax.set_xlabel(r'$it$', fontsize=labsize)
    ax.set_ylabel(r'$|\mathcal{H}|$', fontsize=labsize)

    plt.subplots_adjust(hspace=0)
