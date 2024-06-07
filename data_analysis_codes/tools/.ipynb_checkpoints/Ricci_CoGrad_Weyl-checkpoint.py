"""
Author: Robyn Munoz
Date: 4th Jan 2021

This class calculates:
 - The Christoffel symbols (Gudd)
 - The Ricci tensor with indices down (RicciTdown)
 - The Ricci scalar (RicciS)
 
Additional parameters can be calculated but the formulas used here are for dust in synchronous comoving gauge.
 - The comoving curvature, fractional density, and expansion gradients (Ca, Da, and Za) (M.Bruni et al. 1992)
 - The electric (E) and magnetic (B) part of the Weyl Tensor (Alcubierre's book)
 - Invariants from E and B (W.Bonnor 1995)

"""

import numpy as np
from . import ReadingTools as RRead

class Ricci_CoGrad_Weyl_Class():
    def __init__(self, FD):
        self.FD = FD
        
    def average(self, V, gdet, phi):
        try:
            return np.sum(phi*np.sqrt(gdet))*(self.FD.dx**3)/V
        except:
            return 0.0
        

    def Christoffel_symbol(self, gdown, gup):
        
        #--- Metric derivative
        dgxx = self.FD.D3_scalar(gdown[0,0]) # [dxgxx, dygxx, dzgxx]
        dgxy = self.FD.D3_scalar(gdown[0,1])
        dgxz = self.FD.D3_scalar(gdown[0,2])
        dgyy = self.FD.D3_scalar(gdown[1,1])
        dgyz = self.FD.D3_scalar(gdown[1,2])
        dgzz = self.FD.D3_scalar(gdown[2,2])
        
        #--- Christoffel Symbols down \Gamma_{jkl}
        Gxyz = dgxz[1]+dgxy[2]-dgyz[0]
        Gx = np.array([[dgxx[0], dgxx[1],           dgxx[2]],
                       [dgxx[1], 2*dgxy[1]-dgyy[0], Gxyz],
                       [dgxx[2], Gxyz,              2*dgxz[2]-dgzz[0]]])/2
        
        Gyxz = dgyz[0]+dgxy[2]-dgxz[1]
        Gy = np.array([[2*dgxy[0]-dgxx[1], dgyy[0], Gyxz],
                       [dgyy[0],           dgyy[1], dgyy[2]],
                       [Gyxz,              dgyy[2], 2*dgyz[2]-dgzz[1]]])/2
        
        Gzxy = dgyz[0]+dgxz[1]-dgxy[2]
        Gz = np.array([[2*dgxz[0]-dgxx[2], Gzxy,              dgzz[0]],
                       [Gzxy,              2*dgyz[1]-dgyy[2], dgzz[1]],
                       [dgzz[0],           dgzz[1],           dgzz[2]]])/2
        Gddd = np.array([Gx,Gy,Gz]) #\Gamma_{jkl}
        
        #--- Christoffel Symbols up \Gamma^{i}_{kl}
        Gudd = np.einsum('ij...,jkl... -> ikl...', gup, Gddd)  
        return Gudd, Gddd #\Gamma^{i}_{kl}
    

    def Christoffel_symbol4(self, alpha, dtalpha, Kdown, gup, Gudd3):
        alpha2 = alpha**2
        # not valid for beta =/= 0
        Gttt = -alpha*dtalpha
        Gtti = -alpha*self.FD.D3_scalar(alpha)
        Gtij =  alpha*Kdown
        
        N = len(alpha)
        Gutdtt = np.reshape(-Gttt/alpha2, (1, 1, 1, N, N, N))
        Gutdti = np.reshape(-Gtti/alpha2, (1, 1, 3, N, N, N))
        Gutdit = np.reshape(-Gtti/alpha2, (1, 3, 1, N, N, N))
        Gutdij = np.reshape(-Gtij/alpha2, (1, 3, 3, N, N, N))
        Gt = np.append(np.append(Gutdtt, Gutdti, axis = 2),
                       np.append(Gutdit, Gutdij, axis = 2), axis = 1)
        
        Guidtt = np.einsum('ij...,j...->i...', gup, -Gtti)
        Guidtj = np.einsum('ik...,kj...->ij...', gup, -Gtij)
        
        Guidtt = np.reshape(Guidtt, (3, 1, 1, N, N, N))
        Guidtj = np.reshape(Guidtj, (3, 1, 3, N, N, N))
        Guidjt = np.reshape(Guidtj, (3, 3, 1, N, N, N))
        Gi = np.append(np.append(Guidtt, Guidtj, axis = 2),
                       np.append(Guidjt, Gudd3, axis = 2), axis = 1)
        
        Gudd4 = np.append(Gt, Gi, axis = 0)
        return Gudd4 #\Gamma^{a}_{bc}
    
    def Christoffel_symbol4beta(self, alpha, dtalpha, betaup, dtbetaup, 
                                Kdown, gup, gdown, Gddd3):
        alpha2 = alpha**2
        digdown = self.FD.D3_tensor2(gdown)
        dibetaup = self.FD.D3_tensor1(betaup)
        dibetadown = (np.einsum('ijk...,k...->ij...', digdown, betaup)
                      + np.einsum('jk...,ik...->ij...', gdown, dibetaup))
        LieBgamma = (np.einsum('k...,kij...->ij...', betaup, digdown)
                     + np.einsum('ik...,kj...->ij...', dibetaup, gdown)
                     + np.einsum('jk...,ik...->ij...', dibetaup, gdown))
        dtgdown = LieBgamma - 2 * alpha * Kdown
        dtbetadown = (np.einsum('i...,ik...->k...', betaup, dtgdown) 
                      + np.einsum('i...,ik...->k...', dtbetaup, gdown))
        
        bsq = np.einsum('i...,j...,ij...->...', betaup, betaup, gdown)
        dtbsq = (np.einsum('i...,j...,ij...->...', dtbetaup, betaup, gdown) + 
                 np.einsum('i...,j...,ij...->...', betaup, dtbetaup, gdown) + 
                 np.einsum('i...,j...,ij...->...', betaup, betaup, dtgdown))
        dig00 = self.FD.D3_scalar( - alpha2 + bsq )
        dtg00 = - 2 * alpha * dtalpha + dtbsq
        
        g0i4up = betaup / alpha2
        gij4up = gup - (np.einsum('i...,j...->ij...', betaup, betaup)/alpha2)
        
        # Christoffel ddd
        Gttt = dtg00 / 2
        Gtti = dig00 / 2
        Gtij = self.symmetric_tensor(dibetadown) - (dtgdown / 2)
        Gktt = dtbetadown - Gtti
        Gkti = self.antisymmetric_tensor(dibetadown) + (dtgdown / 2)
        
        # Christoffel udd
        # G^{0}_{ab}
        N = len(alpha)
        Gutdtt = np.reshape((np.einsum('k...,k...->...', 
                                       betaup, Gktt) 
                             - Gttt)/alpha2, (1, 1, 1, N, N, N))
        Gutdti = np.reshape((np.einsum('k...,ki...->i...', 
                                       betaup, Gkti) 
                             - Gtti)/alpha2, (1, 1, 3, N, N, N))
        Gutdit = np.reshape(Gutdti, (1, 3, 1, N, N, N))
        Gutdij = np.reshape((np.einsum('k...,kij...->ij...', 
                                       betaup, Gddd3)
                             - Gtij)/alpha2, (1, 3, 3, N, N, N))
        Gt = np.append(np.append(Gutdtt, Gutdti, axis = 2),
                       np.append(Gutdit, Gutdij, axis = 2), 
                       axis = 1)
        
        # G^{i}_{ab}
        Gukdtt = np.reshape(g0i4up * Gttt 
                            + np.einsum('ij...,j...->i...', 
                                        gij4up, Gktt), 
                            (3, 1, 1, N, N, N))
        Gukdti = np.reshape(np.einsum('k...,i...->ki...', 
                                      g0i4up, Gtti)
                            + np.einsum('jk...,ji...->ki...', 
                                        gij4up, Gkti), 
                            (3, 1, 3, N, N, N))
        Gukdit = np.reshape(Gukdti, (3, 3, 1, N, N, N))
        Gukdij = np.reshape(np.einsum('k...,ij...->kij...', 
                                      g0i4up, Gtij)
                            + np.einsum('lk...,lij...->kij...', 
                                        gij4up, Gddd3), 
                            (3, 3, 3, N, N, N))
        Gi = np.append(np.append(Gukdtt, Gukdti, axis = 2),
                       np.append(Gukdit, Gukdij, axis = 2), 
                       axis = 1)
        
        Gudd4 = np.append(Gt, Gi, axis = 0)
        return Gudd4 #\Gamma^{a}_{bc}
    
    
    def Ricci_TandS(self, gup, Gudd):
        #--- Ricci Tensor
        Rterm0 = np.array([[self.FD.D3x(Gudd[0, j, k]) 
                            + self.FD.D3y(Gudd[1, j, k]) 
                            + self.FD.D3z(Gudd[2, j, k]) 
                            for k in range(3)] 
                           for j in range(3)])
        Gd     = np.einsum('iik... -> k...', Gudd)
        Rterm1 = np.array([self.FD.D3_scalar(Gd[j]) 
                           for j in range(3)])
        Rterm2 = np.einsum('iip...,pjk... -> jk...', Gudd, Gudd)
        Rterm3 = np.einsum('ijp...,pik... -> jk...', Gudd, Gudd)
        RicciTdown = Rterm0 - Rterm1 + Rterm2 - Rterm3   #R_{jk}
        
        # Force symmetry
        #for a in range(3):
        #    for b in range(3):
        #        RicciTdown[b,a] =   RicciTdown[a,b]
                
        #--- Ricci Scalar
        RicciS = np.einsum('jk...,jk... -> ...', gup, RicciTdown)
        
        return RicciTdown, RicciS
        
        
    # Covariant derivative of a scalar, a 2nd order tensor (indices down), 
    # and a 2nd order tensor (indices up), 
    # projected along the fluid flow (using the comoving gauge)
    def CovD3_scalar(self, f):
        return self.FD.D3_scalar(f)
    
    def CovD3_tensor1down(self, Gudd, fdown):
        df = self.FD.D3_tensor1(fdown)
        G1 = - np.einsum('abc..., c... -> ab...', Gudd, fdown)
        return df + G1
    
    def CovD4_tensor1down(self, Gudd, fdown, dtfdown):
        df = np.append(np.array([dtfdown]),
                       self.FD.D3_tensor1(fdown), axis = 0)
        G1 = - np.einsum('abc..., a... -> cb...', Gudd, fdown)
        return df + G1
    
    def CovD3_tensor1up(self, Gudd, fup):
        df = self.FD.D3_tensor1(fup)
        G1 = np.einsum('abc..., a... -> cb...', Gudd, fup)
        return df + G1

    def CovD3_tensor2down(self, Gudd, fdown):
        df = self.FD.D3_tensor2(fdown)
        G1 = - np.einsum('dca..., db... -> cab...', Gudd, fdown)
        G2 = - np.einsum('dcb..., ad... -> cab...', Gudd, fdown)
        return df + G1 + G2
    
    def CovD3_tensor2up(self, Gudd, fup):
        df = self.FD.D3_tensor2(fup)
        G1 = np.einsum('acd..., db... -> cab...', Gudd, fup)
        G2 = np.einsum('bcd..., ad... -> cab...', Gudd, fup)
        return df + G1 + G2
    
    def CovD3_tensor2mixed(self, Gudd, fmixed):
        df = self.FD.D3_tensor2(fmixed)
        G1 = np.einsum('acd..., db... -> cab...', Gudd, fmixed)
        G2 = - np.einsum('dcb..., ad... -> cab...', Gudd, fmixed)
        return df + G1 + G2
    
    # Comoving Gradients
    def CoGrad_Ca(self, gmixed, a, Riccis):
        return a**2 * self.PCovD_scalar(gmixed, Riccis)
    
    def CoGrad_Da(self, gmixed, a, rho):
        return self.PCovD_scalar(gmixed, rho) / rho
    
    def CoGrad_Za(self, gmixed, a, theta):
        return self.PCovD_scalar(gmixed, theta)
    
    
    # Electric and magnetic parts of Weyl tensor and invariants
    def Weyl_E(self, gdown, gup, LCuud3, Christoffeludd, RicciS, RicciTdown, Kdown):
        #, kappa, Tdown3):
        # In synchronous comoving gauge
        Kmixed = np.einsum('ij...,jk...->ik...', gup, Kdown)
        K = np.einsum('ab...,ab...->...', gup, Kdown)
        KKterm = np.einsum('im...,mj... -> ij...', Kdown, Kmixed)
        KKtermHam = np.einsum('ij...,ji... -> ...', Kmixed, Kmixed)
        
        #gmixed = np.einsum('ab...,bc...->ac...', gup, gdown)
        #Sdown  = np.einsum('ca...,db...,cd...->ab...', gmixed, gmixed, Tdown3)
        #Strace = np.einsum('ab...,ab...->...', Sdown, gup)
        
        #Edown  = RicciTdown + K*Kdown - KKterm 
        # - (1/3)*gdown*(kappa*rho+Lambda) - (kappa/2)*(Sdown - gdown*Strace/3)
        Edown  = (RicciTdown + K*Kdown - KKterm 
                  - (1/3)*gdown*(RicciS + K*K - KKtermHam))
        #- (kappa/2)*(Sdown - gdown*Strace/3)
        Emixed = np.einsum('ij...,jk...->ik...', gup, Edown)
        Eup    = np.einsum('ib...,ja...,ab... -> ij...', gup, gup, Edown)
        Etrace = np.einsum('ab...,ab...->...', gup, Edown)
        E2     = np.einsum('ab...,ab...->...', Eup, Edown)
        DdEdd  = self.CovD3_tensor2down(Christoffeludd, Edown)
        divE   = np.einsum('ab...,abc...->c...', gup, DdEdd)
        divE_norm = np.sqrt(np.einsum('ab...,a...,b...->...', gup, divE, divE))
        # TO DO: add lapse in curl calculation
        curlE  = self.symmetric_tensor(np.einsum('cda...,cbd...->ab...',LCuud3, DdEdd))
        curlE_norm = np.sqrt(np.einsum('ac...,bd...,ab...,cd...->...', gup, gup, curlE, curlE))
        
        return {'Edown':Edown, 'Emixed':Emixed, 'Eup':Eup, 
                'Etrace':Etrace, 'E2':E2, 
                'divE_norm':divE_norm, 'curlE_norm':curlE_norm}    
    
    def LeviCivita3symbol(self, Nx, Ny, Nz): # 3 indices
        LC = np.zeros((3, 3, 3, Nx, Ny, Nz))
        for i1 in [0,1,2]:
            for i2 in np.delete([0,1,2], i1):
                for i3 in np.delete([0,1,2], [i1, i2]):
                    num = (i2-i1) * (i3-i1) * (i3-i2)
                    den = (abs(i2-i1) * abs(i3-i1) * abs(i3-i2))
                    LC[i1,i2,i3,:,:,:] = float(num / den)
        return LC
    
    def LeviCivita4symbol(self, Nx, Ny, Nz): # 4 indices
        LC = np.zeros((4, 4, 4, 4, Nx, Ny, Nz))
        for i0 in [0,1,2,3]:
            for i1 in np.delete([0,1,2,3], i0):
                for i2 in np.delete([0,1,2,3], [i0, i1]):
                    for i3 in np.delete([0,1,2,3], [i0, i1, i2]):
                        num = (i1-i0) * (i2-i0) * (i3-i0) * (i2-i1) * (i3-i1) * (i3-i2)
                        den = (abs(i1-i0) * abs(i2-i0) * abs(i3-i0) 
                               * abs(i2-i1) * abs(i3-i1) * abs(i3-i2))
                        LC[i0,i1,i2,i3,:,:,:] = float(num / den)
        return LC
    
    def LeviCivita4(self, gdown4): # 4 indices down
        Nx, Ny, Nz = np.shape(gdown4[0,0])
        LC = self.LeviCivita4symbol(Nx, Ny, Nz)
        return LC*np.sqrt(abs(RRead.det4(gdown4)))
    
    def symmetric_tensor(self, Tdown):
        return (Tdown + np.einsum('ab...->ba...', Tdown))/2
    
    def antisymmetric_tensor(self, Tdown):
        return (Tdown - np.einsum('ab...->ba...', Tdown))/2
    
    def Weyl_B(self, gdown, gup, LCuud3, Christoffeludd, Kdown):        
        Ktrace  = np.einsum('ij...,ij...->...', gup, Kdown)
        Kmixed3 = np.einsum('ij...,jk...->ik...', gup, Kdown)
        Bterm2K = (self.CovD3_scalar(Ktrace) 
                   - np.einsum('ccb... -> b...', 
                               self.CovD3_tensor2mixed(Christoffeludd, Kmixed3)))
        Bterm2  = np.einsum('cdb...,ac...,d...->ab...', LCuud3, gdown, Bterm2K)/2
        
        DKdown = self.CovD3_tensor2down(Christoffeludd, Kdown)
        Bterm1 = np.einsum('cdb...,cda... -> ab...', LCuud3, DKdown)
        
        Bdown  = Bterm1 + Bterm2  
        Bmixed = np.einsum('ij...,jk...->ik...', gup, Bdown)
        Bup    = np.einsum('ib...,ja...,ab... -> ij...', gup, gup, Bdown)
        Btrace = np.einsum('ab...,ab...->...', gup, Bdown)
        B2     = np.einsum('ab...,ab...->...', Bup, Bdown)
        DdBdd  = self.CovD3_tensor2down(Christoffeludd, Bdown)
        divB   = np.einsum('ab...,abc...->c...', gup, DdBdd)
        divB_norm = np.sqrt(np.einsum('ab...,a...,b...->...', gup, divB, divB))
        # TO DO: add lapse in curl calculation
        curlB  = self.symmetric_tensor(np.einsum('cda...,cbd...->ab...',LCuud3, DdBdd))
        curlB_norm = np.sqrt(np.einsum('ac...,bd...,ab...,cd...->...', gup, gup, curlB, curlB))
        return {'Btrace':Btrace, 'B2':B2, 
                'Bdown':Bdown, 'Bmixed':Bmixed, 'Bup':Bup, 
                'divB_norm':divB_norm, 'curlB_norm':curlB_norm}
    
    def Weyl(self, gdown4, gup4, ndown, Edown4, Bdown4): # Alcubierre 2008 eq:8.3.13
        # not finished
        ldown = gdown4 + 2.0 * np.einsum('a...,b...->ab...', ndown, ndown)
        LCudd4 = np.einsum('bf...,ae...,e...,afcd...->bcd...', gup4, gup4, ndown, self.LeviCivita4(gdown4))
        
        Cdown = np.einsum('am...,nb...->abmn...', ldown, Edown4) - np.einsum('an...,mb...->abmn...', ldown, Edown4)
        Cdown -= np.einsum('bm...,na...->abmn...', ldown, Edown4) - np.einsum('bn...,ma...->abmn...', ldown, Edown4)
        Cdown -= np.einsum('mnl...,lab...->abmn...', (np.einsum('m...,nl...->mnl...',ndown,Bdown4) 
                                                      - np.einsum('n...,ml...->mnl...',ndown,Bdown4)),  LCudd4)
        Cdown -= np.einsum('abl...,lmn...->abmn...', (np.einsum('a...,bl...->abl...',ndown,Bdown4) 
                                                      - np.einsum('b...,al...->abl...',ndown,Bdown4)),  LCudd4)
        return Cdown
    
    def Weyl_LB(self, E2, B2): # Bonnor 1995
        return E2 - B2
    
    def Weyl_M(self, Eup, Bdown): # Bonnor 1995
        return np.einsum('ab...,ab...->...', Eup, Bdown)
    
    def Weyl_I(self, LB, M):  # McIntosh et al 1994
        return LB/2 + 1.0j*M
    
    def Weyl_J(self, Emixed, Bmixed):  # McIntosh et al 1994
        return ((1/6) * np.einsum('ab...,bc...,ca...->...', Emixed, Emixed, Emixed) 
                - (1/2) * np.einsum('ab...,bc...,ca...->...', Emixed, Bmixed, Bmixed) 
                - 1.0j * ((1/6) * np.einsum('ab...,bc...,ca...->...', Bmixed, Bmixed, Bmixed) 
                          - (1/2) * np.einsum('ab...,bc...,ca...->...', Emixed, Emixed, Bmixed)))
    
    def Weyl_S(self, I, J): # Baker, Campanelli 2000
        return 27.0*J*J/(I*I*I)
    
    def proj(self, gdown4, a, b):
        return (np.einsum('a...,b...,ab...->...', a, b, gdown4)
                * a / np.einsum('a...,b...,ab...->...', a, a, gdown4))
    
    def norm(self, gdown4, a):
        return np.sqrt(abs(np.einsum('a...,b...,ab...->...', a, a, gdown4)))
    
    def tetradbase(self, gdown4):
        Box_0 = np.zeros(np.shape(gdown4[0,0]))
        v0 = np.array([np.ones(np.shape(gdown4[0,0])), Box_0, Box_0, Box_0])
        v1 = np.array([Box_0, 1.0/np.sqrt(gdown4[1,1]), Box_0, Box_0])
        v2 = np.array([Box_0, Box_0, 1.0/np.sqrt(gdown4[2,2]), Box_0])
        v3 = np.array([Box_0, Box_0, Box_0, 1.0/np.sqrt(gdown4[3,3])])
        
        u0 = v0
        u1 = v1 - self.proj(gdown4, u0, v1)
        u2 = v2 - self.proj(gdown4, u0, v2) - self.proj(gdown4, u1, v2)
        u3 = v3 - self.proj(gdown4, u0, v3) - self.proj(gdown4, u1, v3) - self.proj(gdown4, u2, v3)
        
        e0 = u0 / self.norm(gdown4, u0)
        e1 = u1 / self.norm(gdown4, u1)
        e2 = u2 / self.norm(gdown4, u2)
        e3 = u3 / self.norm(gdown4, u3)
        
        return e0, e1, e2, e3
    
    def nullvectors(self, gdown4):
        e0, e1, e2, e3 = self.tetradbase(gdown4)
        kup  = ( e0 + e1 ) / np.sqrt(2)
        lup  = ( e0 - e1 ) / np.sqrt(2)
        mup  = ( e2 + 1j*e3) / np.sqrt(2)
        mbup = ( e2 - 1j*e3) / np.sqrt(2)
        return lup, kup, mup, mbup
    
    def Weyl_psi(self, gdown, Edict, Bdict):
        
        #with null vectors
        Box_0 = np.zeros(np.shape(Edict['E2']))
        Edown4 = np.array([[Box_0, Box_0, Box_0, Box_0],
                           [Box_0, Edict['Edown'][0,0], Edict['Edown'][0,1], Edict['Edown'][0,2]],
                           [Box_0, Edict['Edown'][1,0], Edict['Edown'][1,1], Edict['Edown'][1,2]],
                           [Box_0, Edict['Edown'][2,0], Edict['Edown'][2,1], Edict['Edown'][2,2]]])
        del Edict
        Bdown4 = np.array([[Box_0, Box_0, Box_0, Box_0],
                           [Box_0, Bdict['Bdown'][0,0], Bdict['Bdown'][0,1], Bdict['Bdown'][0,2]],
                           [Box_0, Bdict['Bdown'][1,0], Bdict['Bdown'][1,1], Bdict['Bdown'][1,2]],
                           [Box_0, Bdict['Bdown'][2,0], Bdict['Bdown'][2,1], Bdict['Bdown'][2,2]]])
        del Bdict
        Box_1 = np.ones(np.shape(Box_0))
        gdown4 = np.array([[-Box_1, Box_0, Box_0, Box_0],
                           [Box_0, gdown[0,0], gdown[0,1], gdown[0,2]],
                           [Box_0, gdown[1,0], gdown[1,1], gdown[1,2]],
                           [Box_0, gdown[2,0], gdown[2,1], gdown[2,2]]])
        del gdown
        gup4 = RRead.inv4(gdown4)
        lup, kup, mup, mbup = self.nullvectors(gdown4)
        ndown = np.array([-Box_1, Box_0, Box_0, Box_0])
        del Box_1, Box_0
        Cdown = self.Weyl(gdown4, gup4, ndown, Edown4, Bdown4)
        del gdown4, gup4, ndown, Edown4, Bdown4
        psi0 = np.einsum('abcd...,a...,b...,c...,d...->...', Cdown, kup, mup, kup, mup)
        psi1 = np.einsum('abcd...,a...,b...,c...,d...->...', Cdown, kup, lup, kup, mup)
        psi2 = np.einsum('abcd...,a...,b...,c...,d...->...', Cdown, kup, mup, mbup, lup)
        psi3 = np.einsum('abcd...,a...,b...,c...,d...->...', Cdown, kup, lup, mbup, lup)
        psi4 = np.einsum('abcd...,a...,b...,c...,d...->...', Cdown, mbup, lup, mbup, lup)
        
        mask = np.where(np.logical_and(abs(psi4)<1e-5, abs(psi0)>1e-5))
        psi0new = psi0
        psi0new[mask] = psi4[mask]
        psi4[mask] = psi0[mask]
        psi0 = psi0new
        psi1new = psi1
        psi1new[mask] = psi3[mask]
        psi3[mask] = psi1[mask]
        psi1 = psi1new
        return [psi0, psi1, psi2, psi3, psi4]
    
    def Weyl_LS(self, psi):
        return psi[2]*psi[4] - (psi[3]**2)
        
    def Weyl_K(self, psi):
        return psi[1]*(psi[4]**2) - 3*psi[4]*psi[3]*psi[2] + 2*(psi[3]**3)
     
    def Weyl_N(self, psi, LS, I):
        return 12*(LS**2) - (psi[4]**2)*I
   
    def Weyl_Invar(self, gdown, Edict, Bdict):
        LB = self.Weyl_LB(Edict['E2'], Bdict['B2'])
        M  = self.Weyl_M(Edict['Eup'], Bdict['Bdown'])
        J  = self.Weyl_J(Edict['Emixed'], Bdict['Bmixed']) #complex
        I  = self.Weyl_I(LB, M) #complex
        S  = self.Weyl_S(I, J) #complex
        psi = self.Weyl_psi(gdown, Edict, Bdict)
        LS = self.Weyl_LS(psi)
        K  = self.Weyl_K(psi)
        N  = self.Weyl_N(psi, LS, I)
        return {'LB':LB, 'M':M, 'J':J, 'I':I, 'S':S, 'LS':LS, 'K':K, 'N':N, 
                'E2':Edict['E2'], 'B2':Bdict['B2'], 
                'psi0':psi[0], 'psi1':psi[1], 'psi2':psi[2], 'psi3':psi[3], 'psi4':psi[4]}
    
    def Backreaction(self, gdet, K, A2):
        Volume = np.sum(np.sqrt(gdet)) * self.FD.dx**3
        return ((2/3) * ( self.average(Volume, gdet, K**2) 
                         - self.average(Volume, gdet, K)**2 ) 
                - 2*self.average(Volume, gdet, A2))
    
    
    
