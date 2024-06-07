import numpy as np

def dichotomy(y_wanted, function, lower_bound, upper_bound, tolerance):
    x_low = lower_bound
    x_upp = upper_bound
    x_mid = (x_low + x_upp) / 2
    y_low = function(x_low)
    y_upp = function(x_upp)
    y_mid = function(x_mid)
    while abs(y_wanted / y_mid - 1) > tolerance:
        if y_wanted > y_mid:
            y_low = y_mid
            x_low = x_mid
            x_mid = (x_low + x_upp) / 2
            y_mid = function(x_mid)
        else:
            y_upp = y_mid
            x_upp = x_mid
            x_mid = (x_low + x_upp) / 2
            y_mid = function(x_mid)
    return x_mid

def Lagrange_interp(x, y, xf, reduce_size=False):
    nbrpoints = 20
    yf = np.zeros(len(xf))
    for i in range(len(xf)):
        
        if reduce_size:
            ixmin = np.argmin(abs(x - xf[i]))
            if ixmin < int(nbrpoints/2):
                xconsidered = x[:ixmin + nbrpoints + 1]
                yconsidered = y[:ixmin + nbrpoints + 1]
            elif ixmin > len(x) - int(nbrpoints/2):
                xconsidered = x[ixmin - nbrpoints + 1:]
                yconsidered = y[ixmin - nbrpoints + 1:]
            else:
                xconsidered = x[ixmin - nbrpoints + 1:ixmin + nbrpoints + 1]
                yconsidered = y[ixmin - nbrpoints + 1:ixmin + nbrpoints + 1]
            k = len(xconsidered)
        else:
            xconsidered = x
            yconsidered = y
            k = len(x)
            
        for j in range(k):
            top = 1
            bot = 1
            for m in range(k):
                if j!=m:
                    top *= xf[i]-xconsidered[m]
                    bot *= xconsidered[j]-xconsidered[m]
            yf[i] += yconsidered[j]*top/bot
    return yf

def Midpoint_interp(fn, fnp1):
    return np.average([fn, fnp1])

def double_data(f):
    f_save = [f[0]]
    for i in range(len(f)-1):
        f_save += [Midpoint_interp(f[i], f[i+1]), f[i+1]]
    return np.array(f_save)

def extrapolate(f):
    return np.append(f, f[-1]+f[-1]-f[-2])

def get_error(f32, f64, f128):
    if len(f32)!=len(f128):
        f32 = double_data(double_data(f32))
        while len(f32)<len(f128):
            f32 = extrapolate(f32)
        while len(f32)>len(f128):
            f32 = f32[:-1]
    if len(f64)!=len(f128):
        f64 = double_data(f64)
        while len(f64)<len(f128):
            f64 = extrapolate(f64)
        while len(f64)>len(f128):
            f64 = f64[:-1]
    c = abs(f32-f64)/abs(f64-f128)
    err = abs((f64-f128)/(c-1))
    for ic in range(len(c)):
        if err[ic]>abs(f32[ic]-f64[ic]) or err[ic]>abs(f64[ic]-f128[ic]):
            err[ic] = np.max([abs(f32[ic]-f64[ic]), abs(f64[ic]-f128[ic])])
    return err, c