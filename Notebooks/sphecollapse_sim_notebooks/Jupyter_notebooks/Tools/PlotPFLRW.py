import re
import pandas as pd
import numpy as np
from data_analysis_codes.tools import LinData
import matplotlib.pyplot as plt
from data_analysis_codes.tools import ReadingTools as RRead

def geterror(p, filetype, var, a_lin):
    path = '/home/robynm/simulations/'+p.sim_name
    try:
        file_err = pd.read_table(path+'/'+filetype+'_data_error.csv', delimiter=',')
        var_err = np.array(list(file_err[var+'_error']))
        return var_err, True
    except:
        return np.zeros(len(a_lin)+2), False
    
def ferr(x, var_error):
    return np.array(list(x)[0::2])[:len(var_error)]

def plot(var, p, xval = 'a', file='h5', n='', nConf=0, plottype='', c='C0', ls='-', lw=3.5, lab='', err='no', ylab='', xlablog=True, ylablog=True, legin=False, leg=True, form=False):

    
    if err=='yes' and 'Th' not in n:
        var_error, have_error = geterror(p, file, var, a_lin)
        var_errortop = ferr(var_data, var_error)+var_error
        var_errordow = ferr(var_data, var_error)-var_error
    else:
        have_error=False
        
    if have_error:
        plt.fill_between(ferr(x, var_error), var_errordow, var_errortop, facecolor=c, edgecolor=c, alpha=0.2)