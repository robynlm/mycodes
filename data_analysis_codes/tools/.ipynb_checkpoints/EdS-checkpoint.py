import numpy as np

class evo:
    def __init__(self):

        self.h = 0.6737
        self.c = 1.0
        self.G = 1.0
        self.a_today = 1.0
        self.Omega_m_EdS = 1.0

        self.kappa = 8 * np.pi * self.G
        self.Omega_l_today = 0.0
        self.Hprop_today = (self.h * self.c) / 2997.9
        self.t_today = 2.0 / (3.0 * self.Hprop_today)
        self.Lambda = 0.0

    def a(self, t):
        """Scale factor"""
        return self.a_today * ((t / self.t_today) ** (2 / 3))
    
    def t_func_a(self, a):
        """Proper time from scale factor"""
        return self.t_today * ((a / self.a_today) ** (3 / 2))

    def Hprop(self, t):
        """Proper Hubble function"""
        return self.Hprop_today * self.t_today / t
    
    def t_func_Hprop(self, Hprop):
        """Proper time from Hubble scalar"""
        return self.Hprop_today * self.t_today / Hprop

    def Omega_m(self, t):
        """Omega_matter """
        return self.Omega_m_EdS

    def an_today(self, t):
        """Scale factor normalised by a(z=0)"""
        return self.a(t) / self.a_today

    def z(self, t):
        """Redshift"""
        return -1 + self.a_today / self.a(t)
    
    def a_func_z(self, z):
        """Scale factor from redshift"""
        return self.a_today / (1 + z)
    
    def t_func_z(self, z):
        """Proper time from redshift"""
        return self.t_func_a(self.a_func_z(z))

    def Hconf(self, t):
        """Conformal Hubble function"""
        return self.a(t) * self.Hprop(t)

    def fL(self, t):
        """Growth index"""
        return self.Omega_m(t) ** (6 / 11)

    def rho(self, t):
        """Energy density"""
        return (3 * self.Omega_m(t) * self.Hprop(t)**2) / self.kappa

