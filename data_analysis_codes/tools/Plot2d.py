import numpy as np
import matplotlib.pyplot as plt
from matplotlib import ticker
from matplotlib.colors import LogNorm
from . import Cstyle
from . import LinData as LinF

class MyPlot2dclass:
    def __init__(self, param, Lin, title='', 
                 filename = '', 
                 fixcolorbar = False, 
                 logcolor = False, 
                 CBarLimits = [0.0, 0.0]):
        self.param = param
        self.title = title
        self.filename = filename
        self.fixcolorbar = fixcolorbar
        self.logcolor = logcolor
        self.CBarLimits = CBarLimits

        self.options = ['xd', 'xnd', 
                        'yd', 'ynd', 
                        'zd', 'znd', 
                        'xy', 'xz', 'yz']
        
        self.Lin = Lin
        self.N = self.param['Nx']
        
        self.xmax = 0.5
        self.dmax = self.xmax*np.sqrt(2)
        self.coord = self.Lin.d1x*np.sqrt(2)/self.param['Lx']
        self.coorcoor = self.Lin.d1x/self.param['Lx']

        #Position at Over Density
        self.i = int(self.N/4)
        
    def formatplot(self, it, axes):
        plt.gca().set_aspect("equal")
        if self.logcolor:
            plt.colorbar(extend='both')
        else:
            plt.colorbar(format='%.1e', 
                         extend='both', 
                         ticks = ticker.MaxNLocator(nbins = 9))
        if self.fixcolorbar:
            plt.clim(self.CBarLimits[0], 
                     self.CBarLimits[1])

        plt.xlabel(axes[0]
                   + r'$\;\;[\lambda_{pert}]$')
        plt.ylabel(axes[1]
                   + r'$\;\;[\lambda_{pert}]$')
        try:
            an = self.Lin.temp_from_temp('an', 'it', it)
            if self.title!='':
                plt.title(self.title + '  ' + r'$a/a_{IN}$' + ' = %.0f'%an)
        except:
            t = self.Lin.temp_from_temp('t', 'it', it)
            if self.title!='':
                plt.title(self.title + '  ' + r'$t$' + ' = %.1f'%t)
        if self.filename!='':
            plt.savefig(self.filename + "_"+''.join(axes)+"_it=%06d.png" % it, bbox_inches='tight')
            plt.close()

    def iniplot(self):
        plt.style.use(Cstyle.style1)
        #plt.figure(figsize=(10, 10*np.sqrt(2)))

    def f_xd(self, data, it):
        self.iniplot()
        dataxd = np.array([data[:, j, j] 
                           for j in range(self.N)])
        self.plot(self.coorcoor, self.coord, dataxd)
        self.formatplot(it, ['x', 'd'])

    def f_xnd(self, data, it):
        self.iniplot()
        dataxd = np.array([data[:, j, self.N-1-j] 
                           for j in range(self.N)])
        self.plot(self.coorcoor, self.coord, dataxd)
        self.formatplot(it, ['x', 'nd'])

    def f_yd(self, data, it):
        self.iniplot()
        dataxd = np.array([data[j, :, j] 
                           for j in range(self.N)])
        self.plot(self.coorcoor, self.coord, dataxd)
        self.formatplot(it, ['y', 'd'])

    def f_ynd(self, data, it):
        self.iniplot()
        dataxd = np.array([data[j, :, self.N-1-j] 
                           for j in range(self.N)])
        self.plot(self.coorcoor, self.coord, dataxd)
        self.formatplot(it, ['y', 'nd'])

    def f_zd(self, data, it):
        self.iniplot()
        dataxd = np.array([data[j, j, :] 
                           for j in range(self.N)])
        self.plot(self.coorcoor, self.coord, dataxd)
        self.formatplot(it, ['z', 'd'])

    def f_znd(self, data, it):
        self.iniplot()
        dataxd = np.array([data[j, self.N-1-j,:] 
                           for j in range(self.N)])
        self.plot(self.coorcoor, self.coord, dataxd)
        self.formatplot(it, ['z', 'nd'])

    def f_xy(self, data, it, i=''):
        self.iniplot()
        if type(i)==str:
            i=self.i
        self.plot(self.coorcoor, self.coorcoor, data[:, :, i])
        self.formatplot(it, ['x', 'y'])

    def f_xz(self, data, it, i=''):
        self.iniplot()
        if type(i)==str:
            i=self.i
        self.plot(self.coorcoor, self.coorcoor, data[:, i, :])
        self.formatplot(it, ['x', 'z'])

    def f_yz(self, data, it, i=''):
        self.iniplot()
        if type(i)==str:
            i=self.i
        self.plot(self.coorcoor, self.coorcoor, data[i, :, :])
        self.formatplot(it, ['y', 'z'])
        
    def plot(self, x, y, z):
        colormap = plt.cm.get_cmap('inferno')
        if self.logcolor:
            plt.pcolormesh(x, y, abs(z), norm=LogNorm(vmin=abs(z).min(), vmax=abs(z).max()), cmap=colormap, shading='gouraud')
        else:
            plt.pcolormesh(x, y, z, cmap=colormap, shading='gouraud')
