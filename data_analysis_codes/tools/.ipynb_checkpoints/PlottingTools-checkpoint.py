import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from . import ReadingTools as RRead
import re

def read_data_file(p, file):
    path = '/home/robynm/simulations/'+p.sim_name
    if file=='h5':
        try:
            f = pd.read_table(path+'/h5_data.csv', delimiter=',')
        except:
            f = pd.read_table(path+'/asc_data.csv', delimiter=',')
    else:
        f = pd.read_table(path+'/asc_data.csv', delimiter=',')
    return f

def getdata(p, zevars=None):
    try:
        tfile = RRead.get_temporal_file('/home/robynm/simulations/'+p.sim_name+'/')
        for var in tfile.keys():
            tfile[var] = tfile[var][0::p.h5_every]
    except:
        tfile = {}
    
    dfile = read_data_file(p, 'h5')
    if zevars==None:
        zevars = dfile.keys()
    for var in zevars:
        tfile[var] = np.array(dfile[var])
    return tfile

def geterror_fromfile(p, filetype, zevar):
    path = '/home/robynm/simulations/'+p.sim_name
    try:
        file_err = pd.read_table(path+'/'+filetype+'_data_error.csv', delimiter=',')
        var_err = {}
        if 't' not in zevar:
            zevar+=['t']
        for var in zevar:
            var_err[var] = np.array(list(file_err[var+'_error']))
        return var_err
    except:
        return False
    
def ferr(x, var_error):
    return np.array(list(x)[0::2])[:len(var_error)]


def plot_locations(x, var, P, ap=0, locs='all', pt='plot', ax=None):
    locations_all = {'OD':['OD','C1'], 'midOD':['midOD','C9'], 'Center':['cent','C2'], 'midUD':['midUD','C8'], 'UD':['UD','C0'], 'Average':['av','C4']}
    locations = ['OD', 'midOD', 'Center', 'midUD', 'UD', 'Average'] if locs=='all' else locs
    for loc in locations:
        y = P['data'][var+'_'+locations_all[loc][0]]
        if ap!=0:
            L = 1821
            yap = y*(P['data']['a']**ap)*(L**2)
        else:
            yap = y
        if ax==None:
            if pt=='plot':
                plt.plot(x, yap, label=loc, color=locations_all[loc][1])
            elif pt=='semilogy':
                plt.semilogy(x, abs(yap), label=loc, color=locations_all[loc][1])
            elif pt=='loglog':
                plt.loglog(abs(x), abs(yap), label=loc, color=locations_all[loc][1])
        else:
            if pt=='plot':
                ax.plot(x, yap, label=loc, color=locations_all[loc][1])
            elif pt=='semilogy':
                ax.semilogy(x, abs(yap), label=loc, color=locations_all[loc][1])
            elif pt=='loglog':
                ax.loglog(abs(x), abs(yap), label=loc, color=locations_all[loc][1])
        
        if P['error']!=False:
            ey = P['error'][var+'_'+locations_all[loc][0]]
            yerr = ferr(y, ey)
            aerr = ferr(P['data']['a'], ey)
            L = 1821
            if pt in ['semilogy', 'loglog']:
                errortop = (abs(yerr)+abs(ey))*(aerr**ap)*(L**ap)
                errordow = (abs(yerr)-abs(ey))*(aerr**ap)*(L**ap)
            else:
                errortop = (yerr+ey)*(aerr**ap)*(L**ap)
                errordow = (yerr-ey)*(aerr**ap)*(L**ap)
            if ax==None:
                plt.fill_between(ferr(x, ey), errordow, errortop, facecolor=locations_all[loc][1], edgecolor=locations_all[loc][1], alpha=0.15)
            else:
                ax.fill_between(ferr(x, ey), errordow, errortop, facecolor=locations_all[loc][1], edgecolor=locations_all[loc][1], alpha=0.15)
    
                 
def plot_locations_th(x, t, f, ap=0, a=1, locs='all', pt='plot', ax=None):
    locations_all = {'OD':['OD','firebrick'], 'midOD':['midOD','C7'], 'Center':['cent','gold'], 'midUD':['midUD','C6'], 'UD':['UD','skyblue']}
    locations = ['OD', 'midOD', 'Center', 'midUD', 'UD'] if locs=='all' else locs
    for loc in locations:
        if ax==None:
            if pt=='plot':
                if ap==2:
                    L = 1821
                    data = (a**ap)*f(t, loc=locations_all[loc][0])*(L**2)
                else:
                    data = (a**ap)*f(t, loc=locations_all[loc][0])
                plt.plot(x, data, label=loc+r'$^{(1)}$', color=locations_all[loc][1], ls=':')
            elif pt=='semilogy':
                plt.semilogy(x, abs((a**ap)*f(t, loc=locations_all[loc][0])), label=loc+r'$^{(1)}$', color=locations_all[loc][1], ls=':')
            elif pt=='loglog':
                plt.loglog(abs(x), abs((a**ap)*f(t, loc=locations_all[loc][0])), label=loc+r'$^{(1)}$', color=locations_all[loc][1], ls=':')
        else:
            if pt=='plot':
                if ap==2:
                    L = 1821
                    data = (a**ap)*f(t, loc=locations_all[loc][0])*(L**2)
                else:
                    data = (a**ap)*f(t, loc=locations_all[loc][0])
                ax.plot(x, data, label=loc+r'$^{(1)}$', color=locations_all[loc][1], ls=':')
            elif pt=='semilogy':
                ax.semilogy(x, abs((a**ap)*f(t, loc=locations_all[loc][0])), label=loc+r'$^{(1)}$', color=locations_all[loc][1], ls=':')
            elif pt=='loglog':
                ax.loglog(abs(x), abs((a**ap)*f(t, loc=locations_all[loc][0])), label=loc+r'$^{(1)}$', color=locations_all[loc][1], ls=':')

def redzone():
    xmin, xmax = plt.gca().axes.get_xlim()
    x = np.arange(xmin-10, xmax+10, 10)
    plt.fill_between(x, [-10]*len(x), [-1]*len(x), facecolor='brown', edgecolor=None, alpha=0.2, label='Unphysical')
    plt.xlim(xmin, xmax)
    plt.text(1, -2.5, 'Unphysical \n region', c='brown')
    
def get_slice_iterations(p):
    path = p['HorSpath']+p['simname']+'/output-0000/'+p['simname']+'/DataSlice/Slice_'
    long = RRead.BASH('ls '+path+'*.csv').split('\n')
    strvar = re.split('Slice_|=|.csv', long[0])[1]
    var_all = np.sort([float(re.split('Slice_|=|.csv', l)[2]) for l in long])
    return strvar, var_all
    
def get_slice_data(p, str_var, var):
    path = p['HorSpath']+p['simname']+'/output-0000/'+p['simname']+'/DataSlice/Slice_'+str_var+'={:.2f}.csv'.format(var)
    f = pd.read_table(path, delimiter=',')
    return f

def geterror_fromdata(x128, data32, data64, data128, plot=False):
    d32 = np.ravel([[d]*2 for d in data32])
    d64 = np.array(data64)
    d128full = np.array(data128)
    d128 = np.array(data128)[0::2]
    maxlen = min([len(d32), len(d64), len(d128), len(x128)])
    d32 = d32[:maxlen]
    d64 = d64[:maxlen]
    d128 = d128[:maxlen]
    c = RRead.safe_division((d32-d64), (d64-d128))
    err = abs(RRead.safe_division((d64-d128),(c-1)))
    
    errextrap = [err[0]]
    for i in range(len(err)-1):
        errextrap +=[10**np.average([np.log10(err[i]), np.log10(err[i+1])]), err[i+1]]
    while len(x128)>len(errextrap):
        errextrap +=[err[-1]]
    errextrap = np.array(errextrap)
    if plot:
        plt.fill_between(x128, d128full+errextrap, d128full-errextrap, alpha=0.15)
    else:
        return x128, d128full+errextrap, d128full-errextrap