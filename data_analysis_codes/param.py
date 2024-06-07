import numpy as np
#self.dx_ref = p.dx/2**p.refinement_level

class BASE():
    def __init__(self):
        self.refinement_level = 0
        self.nbr_ghost        = 3
        self.h                = 0.6737
        self.c                = 1
        self.G                = 1
        self.zc               = 0
        self.Amp_pert         = 0

class BASELCDM():
    def __init__(self):
        BASE.__init__(self)
        self.expansion_ini    = 'LCDM'
        self.expansion_evo    = 'LCDM'
        self.Omega_m_today    = 0.3147

class BASEEdS():
    def __init__(self):
        BASE.__init__(self)
        self.expansion_ini    = 'EdS'
        self.expansion_evo    = 'EdS'
        self.Omega_m_initial  = 1

########################################################
########################################################
#
#                     PURE FLRW
#
########################################################
########################################################


class pflrw_A0_L1821_t1_N16_LCDM():
    def __init__(self):
        BASELCDM.__init__(self)
        self.sim_name      = "pflrw_A0_L1821_t1_N16_LCDM"
        self.parfile_name  = self.sim_name
        self.expansion_evo = 'LCDM'
        self.t_initial     = 1
        self.t_final       = 2600
        self.L             = 1821
        self.dx            = 113.8125
        self.dtfac         = 0.0006
        self.h5_every      = 100
        self.lambda_pert   = self.L
class pflrw_A0_L1821_t1_N32_LCDM():
    def __init__(self):
        pflrw_A0_L1821_t1_N16_LCDM.__init__(self)
        self.sim_name      = "pflrw_A0_L1821_t1_N32_LCDM"
        self.dx            = 56.90625
class pflrw_A0_L1821_t1_N64_LCDM():
    def __init__(self):
        pflrw_A0_L1821_t1_N16_LCDM.__init__(self)
        self.sim_name      = "pflrw_A0_L1821_t1_N64_LCDM"
        self.dx            = 28.453125
class pflrw_A0_L1821_t1_N128_LCDM():
    def __init__(self):
        pflrw_A0_L1821_t1_N16_LCDM.__init__(self)
        self.sim_name      = "pflrw_A0_L1821_t1_N128_LCDM"
        self.dx            = 14.2265625
class pflrw_A0_L1821_t1_N256_LCDM():
    def __init__(self):
        pflrw_A0_L1821_t1_N16_LCDM.__init__(self)
        self.sim_name      = "pflrw_A0_L1821_t1_N256_LCDM"
        self.dx            = 7.11328125
class pflrw_A0_L1821_t1_N512_LCDM():
    def __init__(self):
        pflrw_A0_L1821_t1_N16_LCDM.__init__(self)
        self.sim_name      = "pflrw_A0_L1821_t1_N512_LCDM"
        self.dx            = 3.556640625

        
class pflrw_A0_L1821_t1_N16_LCDM_CTDL():
    def __init__(self):
        BASELCDM.__init__(self)
        self.sim_name      = "pflrw_A0_L1821_t1_N16_LCDM_CTDL"
        self.parfile_name  = self.sim_name
        self.expansion_evo = 'EdS'
        self.t_initial     = 1
        self.t_final       = 2600
        self.L             = 1821
        self.dx            = 113.8125
        self.dtfac         = 0.0006
        self.h5_every      = 100
        self.lambda_pert   = self.L
class pflrw_A0_L1821_t1_N32_LCDM_CTDL():
    def __init__(self):
        pflrw_A0_L1821_t1_N16_LCDM_CTDL.__init__(self)
        self.sim_name      = "pflrw_A0_L1821_t1_N32_LCDM_CTDL"
        self.dx            = 56.90625
class pflrw_A0_L1821_t1_N64_LCDM_CTDL():
    def __init__(self):
        pflrw_A0_L1821_t1_N16_LCDM_CTDL.__init__(self)
        self.sim_name      = "pflrw_A0_L1821_t1_N64_LCDM_CTDL"
        self.dx            = 28.453125
class pflrw_A0_L1821_t1_N128_LCDM_CTDL():
    def __init__(self):
        pflrw_A0_L1821_t1_N16_LCDM_CTDL.__init__(self)
        self.sim_name      = "pflrw_A0_L1821_t1_N128_LCDM_CTDL"
        self.dx            = 14.2265625

        
class pflrw_A0_L1821_t1_N16_EdS():
    def __init__(self):
        BASEEdS.__init__(self)
        self.sim_name     = "pflrw_A0_L1821_t1_N16_EdS"
        self.parfile_name = self.sim_name
        self.t_initial    = 1
        self.t_final      = 2830
        self.L            = 1821
        self.dx           = 113.8125
        self.dtfac        = 0.0006
        self.h5_every     = 100
        self.lambda_pert  = self.L
class pflrw_A0_L1821_t1_N32_EdS():
    def __init__(self):
        pflrw_A0_L1821_t1_N16_EdS.__init__(self)
        self.sim_name     = "pflrw_A0_L1821_t1_N32_EdS"
        self.dx           = 56.90625
class pflrw_A0_L1821_t1_N64_EdS():
    def __init__(self):
        pflrw_A0_L1821_t1_N16_EdS.__init__(self)
        self.sim_name     = "pflrw_A0_L1821_t1_N64_EdS"
        self.dx           = 28.453125
class pflrw_A0_L1821_t1_N128_EdS():
    def __init__(self):
        pflrw_A0_L1821_t1_N16_EdS.__init__(self)
        self.sim_name     = "pflrw_A0_L1821_t1_N128_EdS"
        self.dx           = 14.2265625

########################################################
#
#            Small pert to check convergence
#
########################################################

class pflrw_d75e4_L1821_t1_N16_LCDM():
    def __init__(self):
        pflrw_A0_L1821_t1_N16_LCDM.__init__(self)
        self.sim_name     = "pflrw_d75e4_L1821_t1_N16_LCDM"
        self.Amp_pert     = 0.00260172735246902
class pflrw_d75e4_L1821_t1_N32_LCDM():
    def __init__(self):
        pflrw_A0_L1821_t1_N32_LCDM.__init__(self)
        self.sim_name     = "pflrw_d75e4_L1821_t1_N32_LCDM"
        self.Amp_pert     = 0.00260172735246902
class pflrw_d75e4_L1821_t1_N64_LCDM():
    def __init__(self):
        pflrw_A0_L1821_t1_N64_LCDM.__init__(self)
        self.sim_name     = "pflrw_d75e4_L1821_t1_N64_LCDM"
        self.Amp_pert     = 0.00260172735246902
class pflrw_d75e4_L1821_t1_N128_LCDM():
    def __init__(self):
        pflrw_A0_L1821_t1_N128_LCDM.__init__(self)
        self.sim_name     = "pflrw_d75e4_L1821_t1_N128_LCDM"
        self.Amp_pert     = 0.00260172735246902
        self.parfile_name = self.sim_name


class pflrw_A1e3_L1821_t1_N16_LCDM():
    def __init__(self):
        pflrw_A0_L1821_t1_N16_LCDM.__init__(self)
        self.sim_name     = "pflrw_A1e3_L1821_t1_N16_LCDM"
        self.Amp_pert     = 1e-3
class pflrw_A1e3_L1821_t1_N32_LCDM():
    def __init__(self):
        pflrw_A0_L1821_t1_N32_LCDM.__init__(self)
        self.sim_name     = "pflrw_A1e3_L1821_t1_N32_LCDM"
        self.Amp_pert     = 1e-3
class pflrw_A1e3_L1821_t1_N64_LCDM():
    def __init__(self):
        pflrw_A0_L1821_t1_N64_LCDM.__init__(self)
        self.sim_name     = "pflrw_A1e3_L1821_t1_N64_LCDM"
        self.Amp_pert     = 1e-3
class pflrw_A1e3_L1821_t1_N128_LCDM():
    def __init__(self):
        pflrw_A0_L1821_t1_N128_LCDM.__init__(self)
        self.sim_name     = "pflrw_A1e3_L1821_t1_N128_LCDM"
        self.Amp_pert     = 1e-3
        
        
class pflrw_A1e4_L1821_t1_N16_LCDM():
    def __init__(self):
        pflrw_A0_L1821_t1_N16_LCDM.__init__(self)
        self.sim_name     = "pflrw_A1e4_L1821_t1_N16_LCDM"
        self.Amp_pert     = 1e-4
class pflrw_A1e4_L1821_t1_N32_LCDM():
    def __init__(self):
        pflrw_A0_L1821_t1_N32_LCDM.__init__(self)
        self.sim_name     = "pflrw_A1e4_L1821_t1_N32_LCDM"
        self.Amp_pert     = 1e-4
class pflrw_A1e4_L1821_t1_N64_LCDM():
    def __init__(self):
        pflrw_A0_L1821_t1_N64_LCDM.__init__(self)
        self.sim_name     = "pflrw_A1e4_L1821_t1_N64_LCDM"
        self.Amp_pert     = 1e-4
class pflrw_A1e4_L1821_t1_N128_LCDM():
    def __init__(self):
        pflrw_A0_L1821_t1_N128_LCDM.__init__(self)
        self.sim_name     = "pflrw_A1e4_L1821_t1_N128_LCDM"
        self.Amp_pert     = 1e-4
        
        
class pflrw_A1e5_L1821_t1_N16_LCDM():
    def __init__(self):
        pflrw_A0_L1821_t1_N16_LCDM.__init__(self)
        self.sim_name      = "pflrw_A1e5_L1821_t1_N16_LCDM"
        self.Amp_pert     = 1e-5
class pflrw_A1e5_L1821_t1_N32_LCDM():
    def __init__(self):
        pflrw_A0_L1821_t1_N32_LCDM.__init__(self)
        self.sim_name     = "pflrw_A1e5_L1821_t1_N32_LCDM"
        self.Amp_pert     = 1e-5
class pflrw_A1e5_L1821_t1_N64_LCDM():
    def __init__(self):
        pflrw_A0_L1821_t1_N64_LCDM.__init__(self)
        self.sim_name     = "pflrw_A1e5_L1821_t1_N64_LCDM"
        self.Amp_pert     = 1e-5
class pflrw_A1e5_L1821_t1_N128_LCDM():
    def __init__(self):
        pflrw_A0_L1821_t1_N128_LCDM.__init__(self)
        self.sim_name     = "pflrw_A1e5_L1821_t1_N128_LCDM"
        self.Amp_pert     = 1e-5
        self.parfile_name = self.sim_name

        
class pflrw_A1e5_L1821_t1_N16_LCDM_rho1st():
    def __init__(self):
        pflrw_A1e5_L1821_t1_N16_LCDM.__init__(self)
        self.sim_name      = "pflrw_A1e5_L1821_t1_N16_LCDM_rho1st"
class pflrw_A1e5_L1821_t1_N32_LCDM_rho1st():
    def __init__(self):
        pflrw_A1e5_L1821_t1_N32_LCDM.__init__(self)
        self.sim_name     = "pflrw_A1e5_L1821_t1_N32_LCDM_rho1st"
class pflrw_A1e5_L1821_t1_N64_LCDM_rho1st():
    def __init__(self):
        pflrw_A1e5_L1821_t1_N64_LCDM.__init__(self)
        self.sim_name     = "pflrw_A1e5_L1821_t1_N64_LCDM_rho1st"
class pflrw_A1e5_L1821_t1_N128_LCDM_rho1st():
    def __init__(self):
        pflrw_A1e5_L1821_t1_N128_LCDM.__init__(self)
        self.sim_name     = "pflrw_A1e5_L1821_t1_N128_LCDM_rho1st"

        
class pflrw_A1e5_L1821_t1_N16_LCDM_rhoPHam():
    def __init__(self):
        pflrw_A1e5_L1821_t1_N16_LCDM.__init__(self)
        self.sim_name      = "pflrw_A1e5_L1821_t1_N16_LCDM_rhoPHam"
class pflrw_A1e5_L1821_t1_N32_LCDM_rhoPHam():
    def __init__(self):
        pflrw_A1e5_L1821_t1_N32_LCDM.__init__(self)
        self.sim_name     = "pflrw_A1e5_L1821_t1_N32_LCDM_rhoPHam"
class pflrw_A1e5_L1821_t1_N64_LCDM_rhoPHam():
    def __init__(self):
        pflrw_A1e5_L1821_t1_N64_LCDM.__init__(self)
        self.sim_name     = "pflrw_A1e5_L1821_t1_N64_LCDM_rhoPHam"
class pflrw_A1e5_L1821_t1_N128_LCDM_rhoPHam():
    def __init__(self):
        pflrw_A1e5_L1821_t1_N128_LCDM.__init__(self)
        self.sim_name     = "pflrw_A1e5_L1821_t1_N128_LCDM_rhoPHam"

########################################################
########################################################
#
#                     E.B. & M.B.(2016)
#
########################################################
########################################################

#Note, I forgot to remove CT_DUst::Lambda=...e-7 and I did for this set

class pflrw_d3e2_L1821_t1_N16_LCDM():
    def __init__(self):
        pflrw_A0_L1821_t1_N16_LCDM.__init__(self)
        self.sim_name     = "pflrw_d3e2_L1821_t1_N16_LCDM"
        self.Amp_pert     = 0.01138486133517004756
class pflrw_d3e2_L1821_t1_N32_LCDM():
    def __init__(self):
        pflrw_A0_L1821_t1_N32_LCDM.__init__(self)
        self.sim_name     = "pflrw_d3e2_L1821_t1_N32_LCDM"
        self.Amp_pert     = 0.01138486133517004756
class pflrw_d3e2_L1821_t1_N64_LCDM():
    def __init__(self):
        pflrw_A0_L1821_t1_N64_LCDM.__init__(self)
        self.sim_name     = "pflrw_d3e2_L1821_t1_N64_LCDM"
        self.Amp_pert     = 0.01138486133517004756
class pflrw_d3e2_L1821_t1_N128_LCDM():
    def __init__(self):
        pflrw_A0_L1821_t1_N128_LCDM.__init__(self)
        self.sim_name     = "pflrw_d3e2_L1821_t1_N128_LCDM"
        self.Amp_pert     = 0.01138486133517004756
class pflrw_d3e2_L1821_t1_N256_LCDM():
    def __init__(self):
        pflrw_A0_L1821_t1_N256_LCDM.__init__(self)
        self.sim_name     = "pflrw_d3e2_L1821_t1_N256_LCDM"
        self.Amp_pert     = 0.01138486133517004756
class pflrw_d3e2_L1821_t1_N512_LCDM():
    def __init__(self):
        pflrw_A0_L1821_t1_N512_LCDM.__init__(self)
        self.sim_name     = "pflrw_d3e2_L1821_t1_N512_LCDM"
        self.Amp_pert     = 0.01138486133517004756
        
class pflrw_d3e2_L1821_t1_N16_LCDM_TimeVar():
    def __init__(self):
        pflrw_d3e2_L1821_t1_N16_LCDM.__init__(self)
        self.sim_name     = "pflrw_d3e2_L1821_t1_N16_LCDM_TimeVar"
        self.parfile_name = self.sim_name
class pflrw_d3e2_L1821_t1_N32_LCDM_TimeVar():
    def __init__(self):
        pflrw_d3e2_L1821_t1_N32_LCDM.__init__(self)
        self.sim_name     = "pflrw_d3e2_L1821_t1_N32_LCDM_TimeVar"
        self.parfile_name = self.sim_name
class pflrw_d3e2_L1821_t1_N64_LCDM_TimeVar():
    def __init__(self):
        pflrw_d3e2_L1821_t1_N64_LCDM.__init__(self)
        self.sim_name     = "pflrw_d3e2_L1821_t1_N64_LCDM_TimeVar"
        self.parfile_name = self.sim_name
class pflrw_d3e2_L1821_t1_N128_LCDM_TimeVar():
    def __init__(self):
        pflrw_d3e2_L1821_t1_N128_LCDM.__init__(self)
        self.sim_name     = "pflrw_d3e2_L1821_t1_N128_LCDM_TimeVar"
        self.parfile_name = self.sim_name

        
class pflrw_d3e2_L1206_t1_N16_EdS():
    def __init__(self):
        pflrw_A0_L1821_t1_N16_EdS.__init__(self)
        self.sim_name     = "pflrw_d3e2_L1206_t1_N16_EdS"
        self.parfile_name = self.sim_name
        self.t_final      = 1000
        self.L            = 1206
        self.dx           = 75.375
        self.lambda_pert  = self.L
        self.Amp_pert     = 0.01068203445570103244
class pflrw_d3e2_L1206_t1_N32_EdS():
    def __init__(self):
        pflrw_d3e2_L1206_t1_N16_EdS.__init__(self)
        self.sim_name     = "pflrw_d3e2_L1206_t1_N32_EdS"
        self.dx           = 37.6875
class pflrw_d3e2_L1206_t1_N64_EdS():
    def __init__(self):
        pflrw_d3e2_L1206_t1_N16_EdS.__init__(self)
        self.sim_name     = "pflrw_d3e2_L1206_t1_N64_EdS"
        self.dx           = 18.84375
class pflrw_d3e2_L1206_t1_N128_EdS():
    def __init__(self):
        pflrw_d3e2_L1206_t1_N16_EdS.__init__(self)
        self.sim_name     = "pflrw_d3e2_L1206_t1_N128_EdS"
        self.dx           = 9.421875

########################################################
########################################################
#
#                  Linear collapse
#
########################################################
########################################################

class pflrw_d5e2_L50_z500_N16_LCDM():
    def __init__(self):
        BASELCDM.__init__(self)
        self.sim_name      = "pflrw_d5e2_L50_z50_N16_LCDM"
        self.parfile_name  = self.sim_name
        self.t_initial     = 0.4715789099
        self.t_final       = 400
        self.L             = 50
        self.dx            = 3.125
        self.dtfac         = 0.0004
        self.h5_every      = 1000
        self.lambda_pert   = self.L
        self.Amp_pert      = 1.097498765e-05
class pflrw_d5e2_L50_z500_N32_LCDM():
    def __init__(self):
        pflrw_d5e2_L50_z500_N16_LCDM.__init__(self)
        self.sim_name      = "pflrw_d5e2_L50_z500_N32_LCDM"
        self.dx            = 1.5625
class pflrw_d5e2_L50_z500_N64_LCDM():
    def __init__(self):
        pflrw_d5e2_L50_z500_N16_LCDM.__init__(self)
        self.sim_name      = "pflrw_d5e2_L50_z500_N64_LCDM"
        self.dx            = 0.78125
class pflrw_d5e2_L50_z500_N128_LCDM():
    def __init__(self):
        pflrw_d5e2_L50_z500_N16_LCDM.__init__(self)
        self.sim_name      = "pflrw_d5e2_L50_z500_N128_LCDM"
        self.dx            = 0.390625

class pflrw_d46e3_L50_z500_N16_LCDM():
    def __init__(self):
        BASELCDM.__init__(self)
        self.sim_name      = "pflrw_d46e3_L50_z50_N16_LCDM"
        self.parfile_name  = self.sim_name
        self.t_initial     = 0.4715789099
        self.t_final       = 400
        self.L             = 50
        self.dx            = 3.125
        self.dtfac         = 0.0003
        self.h5_every      = 1000
        self.lambda_pert   = self.L
        self.Amp_pert      = 1.788649529e-05
class pflrw_d46e3_L50_z500_N32_LCDM():
    def __init__(self):
        pflrw_d46e3_L50_z500_N16_LCDM.__init__(self)
        self.sim_name      = "pflrw_d46e3_L50_z500_N32_LCDM"
        self.dx            = 1.5625
class pflrw_d46e3_L50_z500_N64_LCDM():
    def __init__(self):
        pflrw_d46e3_L50_z500_N16_LCDM.__init__(self)
        self.sim_name      = "pflrw_d46e3_L50_z500_N64_LCDM"
        self.dx            = 0.78125
class pflrw_d46e3_L50_z500_N128_LCDM():
    def __init__(self):
        pflrw_d46e3_L50_z500_N16_LCDM.__init__(self)
        self.sim_name      = "pflrw_d46e3_L50_z500_N128_LCDM"
        self.dx            = 0.390625











