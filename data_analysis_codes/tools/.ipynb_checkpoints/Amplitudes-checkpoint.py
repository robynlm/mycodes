import numpy as np
from . import ReadingTools as RRead
from . import LCDM
from . import EdS

def f_Amp(ti, L, N, params, inhom='sin', expansion='LCDM'):
    output = {}
    
    # Type of expansion
    if expansion=='EdS':
        evo = EdS.evo()
    else:
        evo = LCDM.evo()
    a2 = evo.a(ti)**2
    ta2 = 2 * a2
    mta2 = - ta2
    F = evo.fL(ti) + (3/2) * evo.Omega_m(ti)
    iFH = 1 / ( F * evo.Hprop(ti) )
    iFH2 = 1 / ( F * evo.Hprop(ti)**2 )
    m2iFH2 = -2 * iFH2
    gdetflrw = (evo.a(ti))**6
    Kflrw = - 3 * evo.Hprop(ti)
    
    # Define the grid
    x, y, z = np.meshgrid(np.arange(-L[0]/2, L[0]/2, L[0]/N),
                          np.arange(-L[1]/2, L[1]/2, L[1]/N),
                          np.arange(-L[2]/2, L[2]/2, L[2]/N))
    cell_vol = (L[0]*L[1]*L[2])/(N**3)
    Box0 = np.zeros((N, N, N))
    
    
    #=========================================================================================
    # Define Rc
    #=========================================================================================
    if inhom=='sin':
        k_pertx = 2 * np.pi / params['lambda_pertx']
        k_perty = 2 * np.pi / params['lambda_perty']
        k_pertz = 2 * np.pi / params['lambda_pertz']
        
        Rc = (params['Amp_pertx'] * np.sin(k_pertx*x)
              + params['Amp_perty'] * np.sin(k_perty*y)
              + params['Amp_pertz'] * np.sin(k_pertz*z))
        dxRc = params['Amp_pertx'] * k_pertx * np.cos(k_pertx*x)
        dyRc = params['Amp_perty'] * k_perty * np.cos(k_perty*y)
        dzRc = params['Amp_pertz'] * k_pertz * np.cos(k_pertz*z)
        dxdxRc = - params['Amp_pertx'] * k_pertx**2 * np.sin(k_pertx*x)
        dxdyRc = Box0
        dxdzRc = Box0
        dydyRc = - params['Amp_perty'] * k_perty**2 * np.sin(k_perty*y)
        dydzRc = Box0
        dzdzRc = - params['Amp_pertz'] * k_pertz**2 * np.sin(k_pertz*z)
        dxdxdxRc = - params['Amp_pertx'] * k_pertx**3 * np.cos(k_pertx*x)
        dxdxdyRc = Box0
        dxdxdzRc = Box0
        dxdydyRc = Box0
        dxdydzRc = Box0
        dxdzdzRc = Box0
        dydydyRc = - params['Amp_perty'] * k_perty**3 * np.cos(k_perty*y)
        dydydzRc = Box0
        dydzdzRc = Box0
        dzdzdzRc = - params['Amp_pertz'] * k_pertz**3 * np.cos(k_pertz*z)
        dxdxdxdxRc = params['Amp_pertx'] * k_pertx**4 * np.sin(k_pertx*x)
        dxdxdxdyRc = Box0
        dxdxdxdzRc = Box0
        dxdxdydyRc = Box0
        dxdxdydzRc = Box0
        dxdxdzdzRc = Box0
        dxdydydyRc = Box0
        dxdydydzRc = Box0
        dxdydzdzRc = Box0
        dxdzdzdzRc = Box0
        dydydydyRc = params['Amp_perty'] * k_perty**4 * np.sin(k_perty*y)
        dydydydzRc = Box0
        dydydzdzRc = Box0
        dydzdzdzRc = Box0
        dzdzdzdzRc = params['Amp_pertz'] * k_pertz**4 * np.sin(k_pertz*z)
    elif inhom=='exp':
        print('exp')
        # Comoving curvature perturbation : Rc
        tsteepx = 2 * params['ICPertFLRW_steepness_x']
        tsteepy = 2 * params['ICPertFLRW_steepness_y']
        tsteepz = 2 * params['ICPertFLRW_steepness_z']
        expfacx = (- 1/2) * ((3**(1/4))*params['ICPertFLRW_variance_x']) ** (- tsteepx)
        expfacy = (- 1/2) * ((3**(1/4))*params['ICPertFLRW_variance_y']) ** (- tsteepy)
        expfacz = (- 1/2) * ((3**(1/4))*params['ICPertFLRW_variance_z']) ** (- tsteepz)
        xts = x ** tsteepx
        yts = y ** tsteepy
        zts = z ** tsteepz
        expx = np.exp( expfacx * xts)
        expy = np.exp( expfacy * yts)
        expz = np.exp( expfacz * zts)
        Ax = params['ICPertFLRW_exp_amplitude']
        Rc = Ax * expx * expy * expz

        # x^2s derivarives
        if (tsteepz >= 2):
            dxxts = tsteepx * (x ** (tsteepx - 1))
            dxdxxts = tsteepx * (tsteepx - 1) * (x ** (tsteepx - 2))
            if (tsteepz >= 4):
                dxdxdxxts = (tsteepx * (tsteepx - 1) * (tsteepx - 2)   
                               * (x ** (tsteepx - 3)))
                dxdxdxdxxts = (tsteepx * (tsteepx - 1) * (tsteepx - 2)  
                                 * (tsteepx - 3) * (x ** (tsteepx - 4)))
            elif (tsteepz < 4):
                dxdxdxxts = 0
                dxdxdxdxxts = 0
        elif (tsteepz < 2):
            dxxts = 0
            dxdxxts = 0
            dxdxdxxts = 0
            dxdxdxdxxts = 0

        # y^2s derivarives
        if (tsteepz >= 2):
            dyyts = tsteepy * (y ** (tsteepy - 1))
            dydyyts = tsteepy * (tsteepy - 1) * (y ** (tsteepy - 2))
            if (tsteepz >= 4):
                dydydyyts = (tsteepy * (tsteepy - 1) * (tsteepy - 2)    
                               * (y ** (tsteepy - 3)))
                dydydydyyts = (tsteepy * (tsteepy - 1) * (tsteepy - 2)  
                                 * (tsteepy - 3) * (y ** (tsteepy - 4)))
            elif (tsteepz < 4):
                dydydyyts = 0
                dydydydyyts = 0
        elif (tsteepz < 2):
            dyyts = 0
            dydyyts = 0
            dydydyyts = 0
            dydydydyyts = 0

        # z^2s derivarives
        if (tsteepz >= 2):
            dzzts = tsteepz * (z ** (tsteepz - 1))
            dzdzzts = tsteepz * (tsteepz - 1) * (z ** (tsteepz - 2))
            if (tsteepz >= 4):
                dzdzdzzts = (tsteepz * (tsteepz - 1) * (tsteepz - 2)    
                               * (z ** (tsteepz - 3)))
                dzdzdzdzzts = (tsteepz * (tsteepz - 1) * (tsteepz - 2)  
                                 * (tsteepz - 3) * (z ** (tsteepz - 4)))
            elif (tsteepz < 4):
                dzdzdzzts = 0
                dzdzdzdzzts = 0
        elif (tsteepz < 2):
            dzzts = 0
            dzdzzts = 0
            dzdzdzzts = 0
            dzdzdzdzzts = 0
          
        # Single derivatives
        dxexpx = expfacx * dxxts * expx
        dyexpy = expfacy * dyyts * expy
        dzexpz = expfacz * dzzts * expz
        dxRc = Ax * dxexpx * expy * expz
        dyRc = Ax * expx * dyexpy * expz
        dzRc = Ax * expx * expy * dzexpz

        # Double derivatives
        dxdxexpx = expfacx * (dxdxxts * expx + dxxts * dxexpx)
        dydyexpy = expfacy * (dydyyts * expy + dyyts * dyexpy)
        dzdzexpz = expfacz * (dzdzzts * expz + dzzts * dzexpz)
        dxdxRc = Ax * dxdxexpx * expy * expz
        dxdyRc = Ax * dxexpx * dyexpy * expz
        dxdzRc = Ax * dxexpx * expy * dzexpz
        dydyRc = Ax * expx * dydyexpy * expz
        dydzRc = Ax * expx * dyexpy * dzexpz
        dzdzRc = Ax * expx * expy * dzdzexpz

        # Triple derivatives
        dxdxdxexpx = expfacx * (dxdxdxxts * expx + 2 * dxdxxts * dxexpx + dxxts * dxdxexpx)
        dydydyexpy = expfacy * (dydydyyts * expy + 2 * dydyyts * dyexpy + dyyts * dydyexpy)
        dzdzdzexpz = expfacz * (dzdzdzzts * expz + 2 * dzdzzts * dzexpz + dzzts * dzdzexpz)
        dxdxdxRc = Ax * dxdxdxexpx * expy * expz
        dxdxdyRc = Ax * dxdxexpx * dyexpy * expz
        dxdxdzRc = Ax * dxdxexpx * expy * dzexpz
        dxdydyRc = Ax * dxexpx * dydyexpy * expz
        dxdydzRc = Ax * dxexpx * dyexpy * dzexpz
        dxdzdzRc = Ax * dxexpx * expy * dzdzexpz
        dydydyRc = Ax * expx * dydydyexpy * expz
        dydydzRc = Ax * expx * dydyexpy * dzexpz
        dydzdzRc = Ax * expx * dyexpy * dzdzexpz
        dzdzdzRc = Ax * expx * expy * dzdzdzexpz

        # Quadruple derivatives
        dxdxdxdxexpx = expfacx * (dxdxdxdxxts * expx + 3 * dxdxdxxts * dxexpx    
                                    + 3 * dxdxxts * dxdxexpx + dxxts * dxdxdxexpx)
        dydydydyexpy = expfacy * (dydydydyyts * expy + 3 * dydydyyts * dyexpy    
                                    + 3 * dydyyts * dydyexpy + dyyts * dydydyexpy)
        dzdzdzdzexpz = expfacz * (dzdzdzdzzts * expz + 3 * dzdzdzzts * dzexpz    
                                    + 3 * dzdzzts * dzdzexpz + dzzts * dzdzdzexpz)
        dxdxdxdxRc = Ax * dxdxdxdxexpx * expy * expz
        dxdxdxdyRc = Ax * dxdxdxexpx * dyexpy * expz
        dxdxdxdzRc = Ax * dxdxdxexpx * expy * dzexpz
        dxdxdydyRc = Ax * dxdxexpx * dydyexpy * expz
        dxdxdydzRc = Ax * dxdxexpx * dyexpy * dzexpz
        dxdxdzdzRc = Ax * dxdxexpx * expy * dzdzexpz
        dxdydydyRc = Ax * dxexpx * dydydyexpy * expz
        dxdydydzRc = Ax * dxexpx * dydyexpy * dzexpz
        dxdydzdzRc = Ax * dxexpx * dyexpy * dzdzexpz
        dxdzdzdzRc = Ax * dxexpx * expy * dzdzdzexpz
        dydydydyRc = Ax * expx * dydydydyexpy * expz
        dydydydzRc = Ax * expx * dydydyexpy * dzexpz
        dydydzdzRc = Ax * expx * dydyexpy * dzdzexpz
        dydzdzdzRc = Ax * expx * dyexpy * dzdzdzexpz
        dzdzdzdzRc = Ax * expx * expy * dzdzdzdzexpz
    ddRc  = (dxdxRc + dydyRc + dzdzRc) / a2
    
    #=========================================================================================
    # Metric
    #=========================================================================================
    gxx = a2 * ( 1 - 2 * Rc ) + m2iFH2 * dxdxRc
    gxy = m2iFH2 * dxdyRc
    gxz = m2iFH2 * dxdzRc
    gyy = a2 * ( 1 - 2 * Rc ) + m2iFH2 * dydyRc
    gyz = m2iFH2 * dydzRc
    gzz = a2 * ( 1 - 2 * Rc ) + m2iFH2 * dzdzRc
    gdown = np.array([[gxx, gxy, gxz], 
                      [gxy, gyy, gyz], 
                      [gxz, gyz, gzz]])
    
    # Derivatives
    dxgxx = mta2 * dxRc + m2iFH2 * dxdxdxRc
    dygxx = mta2 * dyRc + m2iFH2 * dxdxdyRc
    dzgxx = mta2 * dzRc + m2iFH2 * dxdxdzRc
    dxgxy = m2iFH2 * dxdxdyRc
    dygxy = m2iFH2 * dxdydyRc
    dzgxy = m2iFH2 * dxdydzRc
    dxgxz = m2iFH2 * dxdxdzRc
    dygxz = m2iFH2 * dxdydzRc
    dzgxz = m2iFH2 * dxdzdzRc
    dxgyy = mta2 * dxRc + m2iFH2 * dxdydyRc
    dygyy = mta2 * dyRc + m2iFH2 * dydydyRc
    dzgyy = mta2 * dzRc + m2iFH2 * dydydzRc
    dxgyz = m2iFH2 * dxdydzRc
    dygyz = m2iFH2 * dydydzRc
    dzgyz = m2iFH2 * dydzdzRc
    dxgzz = mta2 * dxRc + m2iFH2 * dxdzdzRc
    dygzz = mta2 * dyRc + m2iFH2 * dydzdzRc
    dzgzz = mta2 * dzRc + m2iFH2 * dzdzdzRc
    dg = np.array([[[dxgxx, dxgxy, dxgxz],
                    [dxgxy, dxgyy, dxgyz],
                    [dxgxz, dxgyz, dxgzz]],
                   [[dygxx, dygxy, dygxz],
                    [dygxy, dygyy, dygyz],
                    [dygxz, dygyz, dygzz]],
                   [[dzgxx, dzgxy, dzgxz],
                    [dzgxy, dzgyy, dzgyz],
                    [dzgxz, dzgyz, dzgzz]]])

    dxdxgxx = mta2 * dxdxRc + m2iFH2 * dxdxdxdxRc
    dxdygxx = mta2 * dxdyRc + m2iFH2 * dxdxdxdyRc
    dxdzgxx = mta2 * dxdzRc + m2iFH2 * dxdxdxdzRc
    dxdxgxy = m2iFH2 * dxdxdxdyRc
    dxdygxy = m2iFH2 * dxdxdydyRc
    dxdzgxy = m2iFH2 * dxdxdydzRc
    dxdxgxz = m2iFH2 * dxdxdxdzRc
    dxdygxz = m2iFH2 * dxdxdydzRc
    dxdzgxz = m2iFH2 * dxdxdzdzRc
    dxdxgyy = mta2 * dxdxRc + m2iFH2 * dxdxdydyRc
    dxdygyy = mta2 * dxdyRc + m2iFH2 * dxdydydyRc
    dxdzgyy = mta2 * dxdzRc + m2iFH2 * dxdydydzRc
    dxdxgyz = m2iFH2 * dxdxdydzRc
    dxdygyz = m2iFH2 * dxdydydzRc
    dxdzgyz = m2iFH2 * dxdydzdzRc
    dxdxgzz = mta2 * dxdxRc + m2iFH2 * dxdxdzdzRc
    dxdygzz = mta2 * dxdyRc + m2iFH2 * dxdydzdzRc
    dxdzgzz = mta2 * dxdzRc + m2iFH2 * dxdzdzdzRc

    dydygxx = mta2 * dydyRc + m2iFH2 * dxdxdydyRc
    dydzgxx = mta2 * dydzRc + m2iFH2 * dxdxdydzRc
    dydygxy = m2iFH2 * dxdydydyRc
    dydzgxy = m2iFH2 * dxdydydzRc
    dydygxz = m2iFH2 * dxdydydzRc
    dydzgxz = m2iFH2 * dxdydzdzRc
    dydygyy = mta2 * dydyRc + m2iFH2 * dydydydyRc
    dydzgyy = mta2 * dydzRc + m2iFH2 * dydydydzRc
    dydygyz = m2iFH2 * dydydydzRc
    dydzgyz = m2iFH2 * dydydzdzRc
    dydygzz = mta2 * dydyRc + m2iFH2 * dydydzdzRc
    dydzgzz = mta2 * dydzRc + m2iFH2 * dydzdzdzRc

    dzdzgxx = mta2 * dzdzRc + m2iFH2 * dxdxdzdzRc
    dzdzgxy = m2iFH2 * dxdydzdzRc
    dzdzgxz = m2iFH2 * dxdzdzdzRc
    dzdzgyy = mta2 * dzdzRc + m2iFH2 * dydydzdzRc
    dzdzgyz = m2iFH2 * dydzdzdzRc
    dzdzgzz = mta2 * dzdzRc + m2iFH2 * dzdzdzdzRc
    ddg = np.array([[[[dxdxgxx, dxdxgxy, dxdxgxz],
                      [dxdxgxy, dxdxgyy, dxdxgyz],
                      [dxdxgxz, dxdxgyz, dxdxgzz]],
                     [[dxdygxx, dxdygxy, dxdygxz],
                      [dxdygxy, dxdygyy, dxdygyz],
                      [dxdygxz, dxdygyz, dxdygzz]],
                     [[dxdzgxx, dxdzgxy, dxdzgxz],
                      [dxdzgxy, dxdzgyy, dxdzgyz],
                      [dxdzgxz, dxdzgyz, dxdzgzz]]],
                    [[[dxdygxx, dxdygxy, dxdygxz],
                      [dxdygxy, dxdygyy, dxdygyz],
                      [dxdygxz, dxdygyz, dxdygzz]],
                     [[dydygxx, dydygxy, dydygxz],
                      [dydygxy, dydygyy, dydygyz],
                      [dydygxz, dydygyz, dydygzz]],
                     [[dydzgxx, dydzgxy, dydzgxz],
                      [dydzgxy, dydzgyy, dydzgyz],
                      [dydzgxz, dydzgyz, dydzgzz]]],
                    [[[dxdzgxx, dxdzgxy, dxdzgxz],
                      [dxdzgxy, dxdzgyy, dxdzgyz],
                      [dxdzgxz, dxdzgyz, dxdzgzz]],
                     [[dydzgxx, dydzgxy, dydzgxz],
                      [dydzgxy, dydzgyy, dydzgyz],
                      [dydzgxz, dydzgyz, dydzgzz]],
                     [[dzdzgxx, dzdzgxy, dzdzgxz],
                      [dzdzgxy, dzdzgyy, dzdzgyz],
                      [dzdzgxz, dzdzgyz, dzdzgzz]]]])
    
    # Determinant
    gdet = (gxx * ( gyy*gzz - gyz*gyz )
            - gxy * ( gxy*gzz - gxz*gyz )
            + gxz * ( gxy*gyz - gxz*gyy ))

    dgdet = np.array([(dxgxx * ( gyy*gzz - gyz*gyz ) - dxgxy * ( gxy*gzz - gxz*gyz )
              + dxgxz * ( gxy*gyz - gxz*gyy ) + gxx * ( dxgyy*gzz - dxgyz*gyz )
              - gxy * ( dxgxy*gzz - dxgxz*gyz ) + gxz * ( dxgxy*gyz - dxgxz*gyy )
              + gxx * ( gyy*dxgzz - gyz*dxgyz )  - gxy * ( gxy*dxgzz - gxz*dxgyz )
              + gxz * ( gxy*dxgyz - gxz*dxgyy )),
              (dygxx * ( gyy*gzz - gyz*gyz ) - dygxy * ( gxy*gzz - gxz*gyz )
              + dygxz * ( gxy*gyz - gxz*gyy ) + gxx * ( dygyy*gzz - dygyz*gyz )
              - gxy * ( dygxy*gzz - dygxz*gyz ) + gxz * ( dygxy*gyz - dygxz*gyy )
              + gxx * ( gyy*dygzz - gyz*dygyz ) - gxy * ( gxy*dygzz - gxz*dygyz )
              + gxz * ( gxy*dygyz - gxz*dygyy )),
              (dzgxx * ( gyy*gzz - gyz*gyz ) - dzgxy * ( gxy*gzz - gxz*gyz )
              + dzgxz * ( gxy*gyz - gxz*gyy ) + gxx * ( dzgyy*gzz - dzgyz*gyz )
              - gxy * ( dzgxy*gzz - dzgxz*gyz ) + gxz * ( dzgxy*gyz - dzgxz*gyy )
              + gxx * ( gyy*dzgzz - gyz*dzgyz ) - gxy * ( gxy*dzgzz - gxz*dzgyz )
              + gxz * ( gxy*dzgyz - gxz*dzgyy ))])
    
    # indices up
    gxxu = ( gyy*gzz - gyz*gyz ) / gdet
    gxyu = ( gyz*gxz - gxy*gzz ) / gdet
    gxzu = ( gxy*gyz - gyy*gxz ) / gdet
    gyyu = ( gxx*gzz - gxz*gxz ) / gdet
    gyzu = ( gxy*gxz - gxx*gyz ) / gdet
    gzzu = ( gxx*gyy - gxy*gxy ) / gdet
    gup = np.array([[gxxu, gxyu, gxzu],
                    [gxyu, gyyu, gyzu],
                    [gxzu, gyzu, gzzu]])
                      
    dxgxxu = (( dxgyy*gzz - dxgyz*gyz + gyy*dxgzz - gyz*dxgyz )*gdet
              - ( gyy*gzz - gyz*gyz )*dgdet[0]) / (gdet*gdet)
    dxgxyu = (( dxgyz*gxz - dxgxy*gzz + gyz*dxgxz - gxy*dxgzz )*gdet
              - ( gyz*gxz - gxy*gzz )*dgdet[0]) / (gdet*gdet)
    dxgxzu = (( dxgxy*gyz - dxgyy*gxz + gxy*dxgyz - gyy*dxgxz )*gdet
              - ( gxy*gyz - gyy*gxz )*dgdet[0]) / (gdet*gdet)
    dxgyyu = (( dxgxx*gzz - dxgxz*gxz + gxx*dxgzz - gxz*dxgxz )*gdet
              - ( gxx*gzz - gxz*gxz )*dgdet[0]) / (gdet*gdet)
    dxgyzu = (( dxgxy*gxz - dxgxx*gyz + gxy*dxgxz - gxx*dxgyz )*gdet
              - ( gxy*gxz - gxx*gyz )*dgdet[0]) / (gdet*gdet)
    dxgzzu = (( dxgxx*gyy - dxgxy*gxy + gxx*dxgyy - gxy*dxgxy )*gdet
              - ( gxx*gyy - gxy*gxy )*dgdet[0]) / (gdet*gdet)

    dygxxu = (( dygyy*gzz - dygyz*gyz + gyy*dygzz - gyz*dygyz )*gdet
              - ( gyy*gzz - gyz*gyz )*dgdet[1]) / (gdet*gdet)
    dygxyu = (( dygyz*gxz - dygxy*gzz + gyz*dygxz - gxy*dygzz )*gdet
              - ( gyz*gxz - gxy*gzz )*dgdet[1]) / (gdet*gdet)
    dygxzu = (( dygxy*gyz - dygyy*gxz + gxy*dygyz - gyy*dygxz )*gdet
              - ( gxy*gyz - gyy*gxz )*dgdet[1]) / (gdet*gdet)
    dygyyu = (( dygxx*gzz - dygxz*gxz + gxx*dygzz - gxz*dygxz )*gdet
              - ( gxx*gzz - gxz*gxz )*dgdet[1]) / (gdet*gdet)
    dygyzu = (( dygxy*gxz - dygxx*gyz + gxy*dygxz - gxx*dygyz )*gdet
              - ( gxy*gxz - gxx*gyz )*dgdet[1]) / (gdet*gdet)
    dygzzu = (( dygxx*gyy - dygxy*gxy + gxx*dygyy - gxy*dygxy )*gdet
              - ( gxx*gyy - gxy*gxy )*dgdet[1]) / (gdet*gdet)

    dzgxxu = (( dzgyy*gzz - dzgyz*gyz + gyy*dzgzz - gyz*dzgyz )*gdet
              - ( gyy*gzz - gyz*gyz )*dgdet[2]) / (gdet*gdet)
    dzgxyu = (( dzgyz*gxz - dzgxy*gzz + gyz*dzgxz - gxy*dzgzz )*gdet
              - ( gyz*gxz - gxy*gzz )*dgdet[2]) / (gdet*gdet)
    dzgxzu = (( dzgxy*gyz - dzgyy*gxz + gxy*dzgyz - gyy*dzgxz )*gdet
              - ( gxy*gyz - gyy*gxz )*dgdet[2]) / (gdet*gdet)
    dzgyyu = (( dzgxx*gzz - dzgxz*gxz + gxx*dzgzz - gxz*dzgxz )*gdet
              - ( gxx*gzz - gxz*gxz )*dgdet[2]) / (gdet*gdet)
    dzgyzu = (( dzgxy*gxz - dzgxx*gyz + gxy*dzgxz - gxx*dzgyz )*gdet
              - ( gxy*gxz - gxx*gyz )*dgdet[2]) / (gdet*gdet)
    dzgzzu = (( dzgxx*gyy - dzgxy*gxy + gxx*dzgyy - gxy*dzgxy )*gdet
              - ( gxx*gyy - gxy*gxy )*dgdet[2]) / (gdet*gdet)
              
    dgu = np.array([[[dxgxxu, dxgxyu, dxgxzu],
                     [dxgxyu, dxgyyu, dxgyzu],
                     [dxgxzu, dxgyzu, dxgzzu]],
                    [[dygxxu, dygxyu, dygxzu],
                     [dygxyu, dygyyu, dygyzu],
                     [dygxzu, dygyzu, dygzzu]],
                    [[dzgxxu, dzgxyu, dzgxzu],
                     [dzgxyu, dzgyyu, dzgyzu],
                     [dzgxzu, dzgyzu, dzgzzu]]])
    
    #=========================================================================================
    # Christoffel symbols
    #=========================================================================================
    
    Gddd = np.zeros((3,3,3,N,N,N))
    for k in range(3):
        for i in range(3):
            for j in range(3):
                Gddd[k, i, j] = ( dg[i, k, j] + dg[j, i, k] - dg[k, i, j] ) / 2
    Gudd = np.einsum('ij..., jkl... -> ikl...', gup, Gddd)
        
    dGddd = np.zeros((3,3,3,3,N,N,N))
    for l in range(3):
        for k in range(3):
            for i in range(3):
                for j in range(3):
                    dGddd[l, k, i, j] = ( ddg[l, i, k, j] + ddg[l, j, i, k] - ddg[l, k, i, j] ) / 2
    
    #=========================================================================================
    # Ricci Scalar
    #=========================================================================================
    
    # R1
    dnGnlk = (np.einsum('iij..., jlk... -> lk...', dgu, Gddd)
              + np.einsum('ij..., ijlk... -> lk...', gup, dGddd))
    R1 = np.einsum('lk..., lk... -> ...', gup, dnGnlk)
    
    # R2
    dlGnkn = (np.einsum('lij..., ikj... -> lk...', dgu, Gddd)
              + np.einsum('ij..., likj... -> lk...', gup, dGddd))
    R2 = np.einsum('lk..., lk... -> ...', gup, dlGnkn)

    # R3
    R3 = np.einsum('ij..., aal..., lij... -> ...', gup, Gudd, Gudd)

    # R4
    R4 = np.einsum('ij..., alj..., lai... -> ...', gup, Gudd, Gudd)

    RicciScalar = R1 - R2 + R3 - R4
    RicciScalar_first = 4 * ddRc
    RicciS = (ta2 * (( dydyRc + dzdzRc ) * gxx + ( dxdxRc + dzdzRc ) * gyy + ( dxdxRc + dydyRc ) * gzz) / gdet 
              - 10 * a2 * a2 * ( dxRc**2 + dyRc**2 + dzRc**2) / gdet 
              - a2 * ( gxx**2 * ( gyy * ( ta2 * dyRc**2 - dzgzz * dzRc ) + gzz * ( ta2 * dzRc**2 - dygyy * dyRc )) 
                      + gyy**2 * ( gxx * ( ta2 * dxRc**2 - dzgzz * dzRc ) + gzz * ( ta2 * dzRc**2 - dxgxx * dxRc )) 
                      + gzz**2 * ( gxx * ( ta2 * dxRc**2 - dygyy * dyRc ) + gyy * ( ta2 * dyRc**2 - dxgxx * dxRc )) 
                      + 2 * ( dgdet[0] * dxRc * ( gyy + gzz ) 
                             + dgdet[1] * dyRc * ( gxx + gzz ) 
                             + dgdet[2] * dzRc * ( gxx + gyy ))) / gdet**2)
    
    #=========================================================================================
    # Extrinsic curvature
    #=========================================================================================

    kxx = - a2 * evo.Hprop(ti) * ( 1 - 2 * Rc ) + ( 2 + evo.fL(ti) ) * dxdxRc * iFH
    kxy = ( 2 + evo.fL(ti) ) * dxdyRc * iFH
    kxz = ( 2 + evo.fL(ti) ) * dxdzRc * iFH
    kyy = - a2 * evo.Hprop(ti) * ( 1 - 2 * Rc ) + ( 2 + evo.fL(ti) ) * dydyRc * iFH
    kyz = ( 2 + evo.fL(ti) ) * dydzRc * iFH
    kzz = - a2 * evo.Hprop(ti) * ( 1 - 2 * Rc ) + ( 2 + evo.fL(ti) ) * dzdzRc * iFH
    Kdown = np.array([[kxx, kxy, kxz],
                      [kxy, kyy, kyz],
                      [kxz, kyz, kzz]])
    """
    dxkxx = - a2 * Hprop * ( - 2 * dxRc ) + ( 2 + fL ) * dxdxdxRc * iFH
    dykxx = - a2 * Hprop * ( - 2 * dyRc )
    dzkxx = - a2 * Hprop * ( - 2 * dzRc )
    dxkyy = - a2 * Hprop * ( - 2 * dxRc )
    dykyy = - a2 * Hprop * ( - 2 * dyRc ) + ( 2 + fL ) * dydydyRc * iFH
    dzkyy = - a2 * Hprop * ( - 2 * dzRc )
    dxkzz = - a2 * Hprop * ( - 2 * dxRc )
    dykzz = - a2 * Hprop * ( - 2 * dyRc )
    dzkzz = - a2 * Hprop * ( - 2 * dzRc ) + ( 2 + fL ) * dzdzdzRc * iFH
    """
    K_L = np.einsum('lk..., lk... -> ...', gup, Kdown)
    Kup = np.einsum('ij..., kl..., jl... -> ik...', gup, gup, Kdown)
    KijKji = np.einsum('lk..., lk... -> ...', Kup, Kdown)

    #=========================================================================================
    # Energy density
    #=========================================================================================
    rho = (RicciScalar + K_L**2 - KijKji - 2*evo.Lambda) / (2 * evo.kappa)
    
    #=========================================================================================
    # Momentum constraint
    #=========================================================================================
    """
    DxK = (dxgxxu*kxx + dxgyyu*kyy + dxgzzu*kzz
           + gxxu*dxkxx + gyyu*dxkyy + gzzu*dxkzz)
    DyK = (dygxxu*kxx + dygyyu*kyy + dygzzu*kzz 
           + gxxu*dykxx + gyyu*dykyy + gzzu*dykzz)
    DzK = (dzgxxu*kxx + dzgyyu*kyy + dzgzzu*kzz 
           + gxxu*dzkxx + gyyu*dzkyy + gzzu*dzkzz)
    DdK = np.array([DxK, DyK, DzK])
    
    DiKix = (gxxu * (dxkxx - 2 * Gxuxx * kxx)
             + gyyu * ( - Gxuyy * kxx - Gyuxy * kyy)
             + gzzu * ( - Gxuzz * kxx - Gzuxz * kzz))
    DiKiy = (gxxu * ( - Gyuxx * kyy - Gxuxy * kxx)
             + gyyu * (dykyy - 2 * Gyuyy * kyy)
             + gzzu * ( - Gyuzz * kyy - Gzuyz * kzz))
    DiKiz = (gxxu * ( - Gzuxx * kzz - Gxuxz * kxx)
             + gyyu * ( - Gzuyy * kzz - Gyuyz * kyy)
             + gzzu * (dzkzz - 2 * Gzuzz * kzz))
    DdKud = np.array([DiKix, DiKiy, DiKiz])
    
    Mom = DdKud - DdK
    
    DdKud_DdKuu = np.einsum('a..., ad..., d... -> ...', DdKud, gup, DdKud)
    DdK_DuK = np.einsum('a..., ad..., d... -> ...', DdK, gup, DdK)
    Mom_energy_scale = np.sqrt(abs(DdKud_DdKuu + DdK_DuK))
    
    MidME = [abs(RRead.safe_division(Mom[0], Mom_energy_scale)), 
             abs(RRead.safe_division(Mom[1], Mom_energy_scale)), 
             abs(RRead.safe_division(Mom[2], Mom_energy_scale))]
    
    output['Mom_max'] = np.nanmax(MidME)
    output['Mom_Q3'] = np.nanpercentile(MidME, 75)
    output['Mom_med'] = np.nanpercentile(MidME, 50)
    output['Mom_Q1'] = np.nanpercentile(MidME, 25)
    output['Mom_min'] = np.nanmin(MidME)
    """
    
    # Inhomogeneity
    dgdet = gdet/gdetflrw - 1
    dK    = K_L/Kflrw - 1
    delta = rho/evo.rho(ti) - 1

    delta1 = ddRc * iFH2
    dgdet1 = -6*(delta1/3+Rc)
    dK1    = -evo.fL(ti)*delta1/3
    
    # Mass fac
    speed_of_light = 299792458   # m.s^{-1}
    Grav_const     = 6.67408e-11 # m^3.kg^{-1}.s^{-2}
    parsec         = 3.0857e16   # m
    Megaparsecc    = parsec*1e6  # m
    MassSun        = 1.98847e30  # kg
    Massfac = (Megaparsecc*speed_of_light**2)/(Grav_const*MassSun*evo.a_today**3)
    # Compute Mass
    V = np.sqrt(abs(gdet))*cell_vol
    Mass = rho * V * Massfac
    output['TotMass'] = np.sum(Mass)
    
    # data
    output['Rc'] = Rc
    output['gdet'] = gdet
    output['gxxu'] = gxxu
    output['gxyu'] = gxyu
    output['delta'] = delta
    output['delta1'] = delta1
    output['RicciScalar'] = RicciScalar
    output['RicciScalarDiag'] = RicciS
    output['R1'] = R1
    output['R2'] = R2
    output['R3'] = R3
    output['R4'] = R4
    output['K_L'] = K_L
    output['KijKji'] = KijKji
        
    # At peak of OD
    if inhom=='sin':
        idx = int(N/4)
    else:
        idx = int(N/2)
    output['deltaOD'] = delta[idx, idx, idx]
    output['delta1OD'] = delta1[idx, idx, idx]
    output['MassOD'] = Mass[idx, idx, idx]
    output['gxxOD'] = gxx[idx, idx, idx]
    output['gdetOD'] = gdet[idx, idx, idx]
    output['dgdetOD'] = dgdet[idx, idx, idx]
    output['dgdetOD1'] = dgdet1[idx, idx, idx]
    output['RicciSOD'] = RicciScalar[idx, idx, idx]
    output['RicciSDiagOD'] = RicciS[idx, idx, idx]
    output['RicciSOD1'] = RicciScalar_first[idx, idx, idx]
    output['KOD'] = K_L[idx, idx, idx]
    output['KxxOD'] = kxx[idx, idx, idx]
    output['dKOD'] = dK[idx, idx, idx]
    output['dKOD1'] = dK1[idx, idx, idx]
        
    return output
    
    
    
    
