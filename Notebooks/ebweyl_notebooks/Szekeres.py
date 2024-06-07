import numpy as np
import scipy.special as sc
from data_analysis_codes.tools import ReadingTools as RRead

t = 1.5

N = 64
L = 20
dx = L/N
xyz = np.arange(-L/2, L/2, dx)
x, y, z = np.meshgrid(xyz, xyz, xyz, indexing='ij')
Box_zero = np.zeros([N, N, N])
Box_ones = np.ones([N, N, N])

h = 0.6737
c = 1
G = 1
kappa = 8*np.pi*G
a_today = 1
Hprop_today = (h*c)/2997.9
t_today_EdS = 2/(3*Hprop_today)
Omega_m_today = 0.3147
rho_flrw_today = (3*Omega_m_today*Hprop_today**2)/kappa
Omega_l_today = 1 - Omega_m_today
Lambda = Omega_l_today*3*(Hprop_today**2)/(c**2)
tauC = np.sqrt(3*Lambda/4)
tau = tauC*t

a = a_today*(Omega_m_today/Omega_l_today)**(1/3)*np.sinh(np.sqrt(Omega_l_today)*t/t_today_EdS)**(2/3)
Hprop =  Hprop_today * np.sqrt( Omega_m_today/((a/a_today)**3) + Omega_l_today )
dta = Hprop*a
dtHprop =  -3*Omega_m_today*dta*(a_today**3)*(a**(-4))*(Hprop_today**2)/(2*Hprop)
dtdta =  dta*Hprop+a*dtHprop#((Lambda/3)-(kappa*rho_flrw(t)/6))

atau = a_today*(Omega_m_today/Omega_l_today)**(1/3)*np.sinh(tau)**(2/3)
dtaua = (2/3)*atau*np.cosh(tau)/np.sinh(tau)

Omega_m = Omega_m_today / ( Omega_m_today + Omega_l_today*((a/a_today)**3) )
rho_flrw = rho_flrw_today/a**3

Amp = 1000
k = 2*np.pi/L
betaP = Amp*(1-np.sin(k*z))
dzbetaP = -k*Amp*np.cos(k*z)
betaM = 0
B = (3/4)*(Hprop_today**2)*(Omega_l_today*(Omega_m_today**2))**(1/3)
A = 1+B*betaP*(x**2+y**2)

fM = np.cosh(tau)/np.sinh(tau)
hyperthing = sc.hyp2f1(5/6, 3/2, 11/6, -np.sinh(tau)**2) 
integrated_part = (3/5)*np.sqrt(np.cosh(tau)**2)*hyperthing*(np.sinh(tau)**(5/3))/np.cosh(tau)
fP = fM*integrated_part
dtaufM = -1/np.sinh(tau)**2
dtfM = tauC*dtaufM
part_to_integrate = (np.sinh(tau)**(2/3))/np.cosh(tau)**2
dtaufP = dtaufM*integrated_part+fM*part_to_integrate
dtfP = tauC*dtaufP
dtaudtaufM = 2*np.cosh(tau)/np.sinh(tau)**3
dtdtfM = (tauC**2)*dtaudtaufM
dtau_part_to_integrate = (2/3)*(np.sinh(tau)**(-1/3))*(np.cosh(tau)**(-1))-2*(np.sinh(tau)**(5/3))*(np.cosh(tau)**(-3))
dtaudtaufP = dtaudtaufM*integrated_part+2*dtaufM*part_to_integrate+fM*dtau_part_to_integrate
dtdtfP = (tauC**2)*dtaudtaufP

F     = betaM*fM + betaP*fP
dzF   = dzbetaP*fP
Z     = F + A
dxZ   = 2*B*betaP*x
dxdxZ = 2*B*betaP
dyZ   = 2*B*betaP*y
dydyZ = 2*B*betaP
dzZ   = dzF + B*dzbetaP*(x**2+y**2)
dtZ   = betaM*dtfM + betaP*dtfP
dtdtZ = betaM*dtdtfM + betaP*dtdtfP

delta =  -F/Z
rho =  rho_flrw*(1+delta)

gdown = (a**2)*np.array([[Box_ones, Box_zero, Box_zero],
                         [Box_zero, Box_ones, Box_zero],
                         [Box_zero, Box_zero, Z**2]])
gup = RRead.inv3(gdown)
gdown4 = np.array([[-Box_ones, Box_zero, Box_zero, Box_zero],
                   [Box_zero, gdown[0,0], gdown[0,1], gdown[0,2]],
                   [Box_zero, gdown[0,1], gdown[1,1], gdown[1,2]],
                   [Box_zero, gdown[0,2], gdown[1,2], gdown[2,2]]])
gup4 = RRead.inv4(gdown4)

Kdown = -(a**2)*np.array([[Box_ones*Hprop, Box_zero, Box_zero],
                          [Box_zero, Box_ones*Hprop, Box_zero],
                          [Box_zero, Box_zero, (Z**2)*(Hprop+dtZ/Z)]])

udown = np.array([-Box_ones, Box_zero, Box_zero, Box_zero])
uup = np.einsum('ab...,b...->a...', gup4, udown)
Tdown4 = rho*np.einsum('a...,b...->ab...',udown,udown)
Tdown3 = Tdown4[1:,1:]