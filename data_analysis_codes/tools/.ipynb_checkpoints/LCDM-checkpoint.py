import numpy as np

class evo:
    def __init__(self):
        self.h = 0.6737
        self.c = 1.0
        self.G = 1.0
        self.a_today = 1.0
        self.Omega_m_today = 0.3147

        self.kappa = 8 * np.pi * self.G
        self.Omega_l_today = 1 - self.Omega_m_today
        self.Hprop_today = (self.h * self.c) / 2997.9
        self.t_today_EdS = 2 / (3 * self.Hprop_today)
        self.Lambda = (self.Omega_l_today * 3 * (self.Hprop_today ** 2) 
                       / (self.c ** 2))

    def a(self, t):
        """Scale factor"""
        return (self.a_today 
                * (self.Omega_m_today / self.Omega_l_today) ** (1 / 3) 
                * np.sinh(np.sqrt(self.Omega_l_today) * t 
                          / self.t_today_EdS) ** (2 / 3))

    def Hprop(self, t):
        """Proper Hubble function"""
        return (self.Hprop_today 
                * np.sqrt( self.Omega_m_today / (self.an_today(t) ** 3) 
                          + self.Omega_l_today ))

    def Omega_m(self, t):
        """Omega_matter """
        return (self.Omega_m_today 
                / ( self.Omega_m_today 
                   + self.Omega_l_today * (self.an_today(t) ** 3) ))

    def an_today(self, t):
        """Scale factor normalised by a(z=0)"""
        return self.a(t) / self.a_today

    def z(self, t):
        """Redshift"""
        return -1 + self.a_today / self.a(t)

    def Hconf(self, t):
        """Conformal Hubble function"""
        return self.a(t) * self.Hprop(t)

    def fL(self, t):
        """Growth index"""
        return self.Omega_m(t) ** (6/11)

    def rho(self, t):
        """Energy density"""
        return (3 * self.Omega_m(t) * self.Hprop(t)**2) / self.kappa

