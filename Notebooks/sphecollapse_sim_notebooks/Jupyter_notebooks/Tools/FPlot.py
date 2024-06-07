from . import PlotPFLRW
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import gridspec
    
def plotHam(p):
    plt.figure(figsize=(20, 6))
    xlab = r'$a/a_{IN}$'
    
    plt.subplot(121)
    plot_4locations('Ham', p)
    plt.grid()
    plt.xlabel(xlab)
    plt.ylabel(r'$|Ham|$')
    
    plt.subplot(122)
    plot_4locations('Ham1', p)
    plt.grid()
    plt.xlabel(xlab)
    plt.ylabel(r'$|Ham^{(1)}|$')
    plt.legend(bbox_to_anchor=(1,1), fontsize=15)
    

    
def plotW(p):
    plt.figure(figsize=(20, 7.5))
    for nbr, var, title in zip([141, 142, 143, 144], 
                               ['E2', 'B2', 'L', 'M'], 
                               [r'$E^2$', r'$B^2$', r'$L$', r'$M$']):
        plt.subplot(nbr)
        plot_4locations(var, p)
        plt.grid()
        plt.title(title)
        plt.xlabel(r'$a/a_{IN}$')
    plt.subplots_adjust(wspace=0.25)
    plt.legend(bbox_to_anchor=(1,1), fontsize=15)