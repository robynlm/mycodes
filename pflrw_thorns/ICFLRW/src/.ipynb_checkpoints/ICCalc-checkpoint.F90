#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine ICFLRW_ICCalc (CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS
  integer :: i,j,k
  integer, parameter :: dp = 8
  real(dp), parameter :: pi = 4._dp*atan(1._dp)
  CCTK_REAL :: H0, a0, t0_EdS, Omega_lambda0
  CCTK_REAL :: t, aa, a2, kappa, Hprop
  CCTK_REAL :: gdet
  CCTK_REAL :: RicciScalar
  CCTK_REAL :: K_Loc, KijKji

  logical :: want_flrw, want_EdS

  want_flrw   = CCTK_EQUALS (my_initial_data, "ICFLRW")
  if (want_flrw) then
    want_EdS     = CCTK_EQUALS (ICFLRW_expansion, "EdS")

    H0 = ICFLRW_h * ICFLRW_c / 2997.9_dp ! Units are Mpc
    a0 = 1._dp + ICFLRW_z_comoving_ref
    t0_EdS = 2._dp / ( 3._dp * H0 )
    Omega_lambda0 = 1._dp - ICFLRW_Omega_matter0

    t = cctk_time
    if (want_EdS) then
      aa = a0 * ( t / t0_EdS )**(2._dp/3._dp)
      Hprop = 2._dp / ( 3._dp * t )
      Lambda = 0._dp
    else
      aa = a0 * ( ICFLRW_Omega_matter0 / Omega_lambda0 )**(1._dp/3._dp) &
           * sinh( sqrt(Omega_lambda0) * t / t0_EdS )**(2._dp/3._dp)
      Hprop = H0 * sqrt( ICFLRW_Omega_matter0 * ( aa / a0 )**(-3._dp) + Omega_lambda0 )
      Lambda = 3._dp * Omega_lambda0 * H0**2._dp / ICFLRW_c**2._dp
    endif

    a2 = aa**2._dp
    kappa = 8._dp * pi * ICFLRW_G

    do k = 1, cctk_lsh(3)
      do j = 1, cctk_lsh(2)
        do i = 1, cctk_lsh(1)

          !============================================================================================
          ! Metric
          !============================================================================================

          gxx(i,j,k) = a2
          gxy(i,j,k) = 0._dp
          gxz(i,j,k) = 0._dp
          gyy(i,j,k) = a2
          gyz(i,j,k) = 0._dp
          gzz(i,j,k) = a2

          ! Determinant
          gdet = a2 * a2 * a2

          !============================================================================================
          ! Ricci Scalar
          !============================================================================================

          RicciScalar = 0._dp

          !============================================================================================
          ! Extrinsic Curvature
          !============================================================================================

          kxx(i,j,k) = - a2 * Hprop
          kxy(i,j,k) = 0._dp
          kxz(i,j,k) = 0._dp
          kyy(i,j,k) = - a2 * Hprop
          kyz(i,j,k) = 0._dp
          kzz(i,j,k) = - a2 * Hprop

          K_Loc = - 3._dp * Hprop

          KijKji = 3._dp * Hprop * Hprop

          !============================================================================================
          ! Density
          !============================================================================================

          rho(i,j,k) = (RicciScalar + K_Loc**2._dp - KijKji - 2._dp*Lambda) / (2._dp * kappa)

          rhodp(i,j,k) = rho(i,j,k) * gdet ** (1._dp/2._dp)
          prs(i,j,k) = 0._dp
          eps(i,j,k) = 0._dp
          u1(i,j,k) = 0._dp
          u2(i,j,k) = 0._dp
          u3(i,j,k) = 0._dp
           
        enddo  !i
      enddo  !j
    enddo  !k
  endif  !if you want flrw
end subroutine ICFLRW_ICCalc




