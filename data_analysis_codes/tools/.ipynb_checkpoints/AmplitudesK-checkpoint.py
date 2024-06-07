import numpy as np

def f_Amp(ti, lambda_pert, Amp_pert, Amp_pertK):
    
    lambda_pertx, lambda_perty, lambda_pertz = lambda_pert
    Amp_pertx, Amp_perty, Amp_pertz = Amp_pert
    Amp_pertxK, Amp_pertyK, Amp_pertzK = Amp_pertK
    
    # Constants
    h = 0.6737
    expansion='LCDM'
    c = 1
    H0 = c*h/2997.9 # Units are Mpc
    G = 1
    kappa = 8*np.pi*G
    Omega_m0 = 0.3147
    Omega_l0 = 1 - Omega_m0
    zR = 0
    a0 = 1.0+zR
    t0_EdS = 2/(3*H0)
    
    if expansion=='EdS':
        aa = a0*(ti/t0_EdS)**(2/3)
        Hprop = 2 / ( 3 * ti )
        Omega_m = 1
        Lambda = 0
    else:
        aa = a0 * ( Omega_m0 / Omega_l0 )**(1/3) * np.sinh( np.sqrt(Omega_l0) * ti / t0_EdS )**(2/3)
        Hprop = H0 * np.sqrt( Omega_m0 * ( aa / a0 )**(-3) + Omega_l0 )
        Omega_m = Omega_m0 / ( Omega_m0 + Omega_l0 * ( aa / a0 )**3 )
        Omega_l = 1 - Omega_m
        Lambda = 3 * Omega_l0 * (H0**2) / (c**2)
    
    a2 = aa**2
    ta2 = 2 * a2
    
    fL   = Omega_m**(6/11)
    F    = fL + (3/2) * Omega_m
    iFH  = 1 / ( F * Hprop )
    iFH2  = 1 / ( F * Hprop**2 )

    k_pertx = 2 * np.pi / lambda_pertx
    k_perty = 2 * np.pi / lambda_perty
    k_pertz = 2 * np.pi / lambda_pertz

    gdetflrw = (aa)**6
    rhoflrw  = 3 * Hprop**2 * Omega_m / kappa
    Kflrw    = - 3 * Hprop

    N = 64
    x, y, z = np.meshgrid(np.arange(-lambda_pertx/2, 0, lambda_pertx/N),
                          np.arange(-lambda_perty/2, 0, lambda_perty/N), 
                          np.arange(-lambda_pertz/2, 0, lambda_pertz/N))
    cell_vol = (lambda_pertx*lambda_perty*lambda_perty)/N**3
    Rc = (Amp_pertx * np.sin(k_pertx*x)
          + Amp_perty * np.sin(k_perty*y)
          + Amp_pertz * np.sin(k_pertz*z))
    RcK = (Amp_pertxK * np.sin(k_pertx*x)
           + Amp_pertyK * np.sin(k_perty*y)
           + Amp_pertzK * np.sin(k_pertz*z))
    dxRc = Amp_pertx * k_pertx * np.cos(k_pertx*x)
    dyRc = Amp_perty * k_perty * np.cos(k_perty*y)
    dzRc = Amp_pertz * k_pertz * np.cos(k_pertz*z)
    dxdxRc = - Amp_pertx * k_pertx**2 * np.sin(k_pertx*x)
    dydyRc = - Amp_perty * k_perty**2 * np.sin(k_perty*y)
    dzdzRc = - Amp_pertz * k_pertz**2 * np.sin(k_pertz*z)
    dxdxRcK = - Amp_pertxK * k_pertx**2 * np.sin(k_pertx*x)
    dydyRcK = - Amp_pertyK * k_perty**2 * np.sin(k_perty*y)
    dzdzRcK = - Amp_pertzK * k_pertz**2 * np.sin(k_pertz*z)
    dxdxdxRc = - Amp_pertx * k_pertx**3 * np.cos(k_pertx*x)
    dydydyRc = - Amp_perty * k_perty**3 * np.cos(k_perty*y)
    dzdzdzRc = - Amp_pertz * k_pertz**3 * np.cos(k_pertz*z)
        
    ddRc     = (dxdxRc + dydyRc + dzdzRc) / a2
    gxx = a2 * ( 1 - 2 * Rc ) - 2 * dxdxRc * iFH2
    gyy = a2 * ( 1 - 2 * Rc ) - 2 * dydyRc * iFH2
    gzz = a2 * ( 1 - 2 * Rc ) - 2 * dzdzRc * iFH2
    gdet = gxx*gyy*gzz
    gxxu = ( gyy*gzz ) / gdet
    gyyu = ( gxx*gzz ) / gdet
    gzzu = ( gxx*gyy ) / gdet
    
    dxgxx = a2 * ( - 2 * dxRc ) - 2 * dxdxdxRc * iFH2
    dygxx = a2 * ( - 2 * dyRc )
    dzgxx = a2 * ( - 2 * dzRc )
    dxgyy = a2 * ( - 2 * dxRc )
    dygyy = a2 * ( - 2 * dyRc ) - 2 * dydydyRc * iFH2
    dzgyy = a2 * ( - 2 * dzRc )
    dxgzz = a2 * ( - 2 * dxRc )
    dygzz = a2 * ( - 2 * dyRc )
    dzgzz = a2 * ( - 2 * dzRc ) - 2 * dzdzdzRc * iFH2
    
    dxgdet = dxgxx*gyy*gzz + gxx*dxgyy*gzz + gxx*gyy*dxgzz
    dygdet = dygxx*gyy*gzz + gxx*dygyy*gzz + gxx*gyy*dygzz
    dzgdet = dzgxx*gyy*gzz + gxx*dzgyy*gzz + gxx*gyy*dzgzz

    RicciS = 4 * ddRc
    RicciSfull = ta2 * (( dydyRc + dzdzRc ) * gxx + ( dxdxRc + dzdzRc ) * gyy + ( dxdxRc + dydyRc ) * gzz) / gdet - 10 * a2 * a2 * ( dxRc**2 + dyRc**2 + dzRc**2) / gdet - a2 * ( gxx**2 * ( gyy * ( ta2 * dyRc**2 - dzgzz * dzRc ) + gzz * ( ta2 * dzRc**2 - dygyy * dyRc )) + gyy**2 * ( gxx * ( ta2 * dxRc**2 - dzgzz * dzRc ) + gzz * ( ta2 * dzRc**2 - dxgxx * dxRc )) + gzz**2 * ( gxx * ( ta2 * dxRc**2 - dygyy * dyRc ) + gyy * ( ta2 * dyRc**2 - dxgxx * dxRc )) + 2 * ( dxgdet * dxRc * ( gyy + gzz ) + dygdet * dyRc * ( gxx + gzz ) + dzgdet * dzRc * ( gxx + gyy ))) / gdet**2

    kxx = - a2 * Hprop * ( 1 - 2 * RcK ) + ( 2 + fL ) * dxdxRcK * iFH
    kyy = - a2 * Hprop * ( 1 - 2 * RcK ) + ( 2 + fL ) * dydyRcK * iFH
    kzz = - a2 * Hprop * ( 1 - 2 * RcK ) + ( 2 + fL ) * dzdzRcK * iFH
    K_L = gxxu*kxx + gyyu*kyy + gzzu*kzz
    KijKji = gxxu*kxx*gxxu*kxx + gyyu*kyy*gyyu*kyy + gzzu*kzz*gzzu*kzz

    rho = (RicciSfull + K_L**2 - KijKji - 2*Lambda) / (2 * kappa)

    dgdet = gdet/gdetflrw - 1
    dK    = K_L/Kflrw - 1
    delta = rho/rhoflrw - 1

    delta1 = ddRc * iFH2
    dgdet1 = -6*(delta1/3+Rc)
    dK1    = -fL*delta1/3
    
    # Mass fac
    speed_of_light = 299792458   # m.s^{-1}
    Grav_const     = 6.67408e-11 # m^3.kg^{-1}.s^{-2}
    parsec         = 3.0857e16   # m
    Megaparsecc    = parsec*1e6  # m
    MassSun        = 1.98847e30  # kg
    Massfac = (Megaparsecc*speed_of_light**2)/(Grav_const*MassSun*a0**3)
    # Compute Mass
    V = np.sqrt(abs(gdet))*cell_vol
    Mass = rho * V * Massfac
        
    idx = int(N/4)
    deltaOD = delta[idx, idx, idx]
    delta1OD = delta1[idx, idx, idx]
    MassOD = Mass[idx, idx, idx]
    gdetOD = gdet[idx, idx, idx]
    RicciSOD = RicciSfull[idx, idx, idx]
    KOD = K_L[idx, idx, idx]
        
    return deltaOD, delta1OD, np.sum(Mass), MassOD, gdetOD, RicciSOD, KOD
    
    
    
    