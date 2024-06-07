import numpy as np
from data_analysis_codes.ebweyl import ebweyl as ebw
from data_analysis_codes.tools import LCDM
from data_analysis_codes.tools import EdS

class ICPertFLRW():
    def __init__(self, Amp_pert, tIN, lambda_pert, N, expansion):
            
        if expansion=='EdS':
            evo = EdS.evo()
        else:
            evo = LCDM.evo()
    
        a2 = evo.a(tIN) **2
        ma2H = - a2 * evo.Hprop(tIN)
        F = evo.fL(tIN) + (3/2) * evo.Omega_m(tIN)
        m2iFH2 = - 2 / ( F * (evo.Hprop(tIN)**2) )
        tfiFH = ( 2 + evo.fL(tIN) ) / ( F * evo.Hprop(tIN) )

        k_pert = 2 * np.pi / lambda_pert

        dx = lambda_pert / N
        xyz = np.arange(-lambda_pert/2, lambda_pert/2, dx)
        x, y, z = np.meshgrid(xyz, xyz, xyz, indexing='ij')
        Box0 = np.zeros((N, N, N))
        Rc = Amp_pert * (np.sin(k_pert*x) + np.sin(k_pert*y) + np.sin(k_pert*z))
        dxdxRc = - Amp_pert * (k_pert**2) * np.sin(k_pert*x)
        dydyRc = - Amp_pert * (k_pert**2) * np.sin(k_pert*y)
        dzdzRc = - Amp_pert * (k_pert**2) * np.sin(k_pert*z)

        self.gdown4 = np.array([[Box0 - 1, Box0, Box0, Box0],
                                [Box0, a2*(1 - 2*Rc) + dxdxRc*m2iFH2, Box0, Box0],
                                [Box0, Box0, a2*(1 - 2*Rc) + dydyRc*m2iFH2, Box0],
                                [Box0, Box0, Box0, a2*(1 - 2*Rc) + dzdzRc*m2iFH2]])
        self.Kdown4 = np.array([[Box0, Box0, Box0, Box0],
                                [Box0, ma2H*(1 - 2*Rc) + dxdxRc*tfiFH, Box0, Box0],
                                [Box0, Box0, ma2H*(1 - 2*Rc) + dydyRc*tfiFH, Box0],
                                [Box0, Box0, Box0, ma2H*(1 - 2*Rc) + dzdzRc*tfiFH]])

        FD = ebw.FiniteDifference(dx, N)
        EBW = ebw.Weyl(FD, self.gdown4, self.Kdown4)
        #Gudd3 = EBW.christoffel_symbol_udd3()
        #RicciTdown3 = EBW.ricci_tensor_down3(Gudd3)
        #RicciS3 = EBW.trace_rank2tensor3(RicciTdown3)

        self.K = EBW.trace_rank2tensor3(self.Kdown4[1:,1:])
        #Kmixed3 = np.einsum('ij..., jk... -> ik...', EBW.gammaup3, EBW.Kdown3)
        #KijKji = np.einsum('ij..., ji... -> ...', Kmixed3, Kmixed3)

        #rho = (RicciS3 + self.K**2 - KijKji - 2*evo.Lambda) / (2 * evo.kappa)




