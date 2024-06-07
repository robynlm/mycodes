#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine ICCalcPertGRHydro (CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS
  integer :: i,j,k
  integer, parameter :: dp = 8
  real(dp), parameter :: pi = 4._dp*atan(1._dp)
  CCTK_REAL :: t, aa, a2, Hprop, kflrw, kappa, rhoflrw
  logical :: want_data, want_hydro

  want_data  = CCTK_EQUALS (initial_data,  "ICPertFLRW_GRHydro")
  want_hydro = CCTK_EQUALS (initial_hydro, "ICPertFLRW_GRHydro")

  t = cctk_time  !proper time
  aa = ai * ( t / ti )**(2._dp/3._dp)
  a2 = aa**2._dp
  kappa = 8._dp * pi * G

  Hprop = 2._dp / ( 3._dp * ti )
  kflrw = -a2 * Hprop
  rhoflrw = 3._dp * Hprop**2._dp * Omega_mi / kappa

  do k = 1, cctk_lsh(3)
     do j = 1, cctk_lsh(2)
        do i = 1, cctk_lsh(1)

           if (want_data) then
              gxx(i,j,k) = a2
              gxy(i,j,k) = 0._dp
              gxz(i,j,k) = 0._dp
              gyy(i,j,k) = a2
              gyz(i,j,k) = 0._dp
              gzz(i,j,k) = a2

              kxx(i,j,k) = kflrw
              kxy(i,j,k) = 0._dp
              kxz(i,j,k) = 0._dp
              kyy(i,j,k) = kflrw
              kyz(i,j,k) = 0._dp
              kzz(i,j,k) = kflrw
           endif

           if (want_hydro) then
              rho(i,j,k) = rhoflrw
              press(i,j,k) = 0._dp
              eps(i,j,k) = 0._dp
              vel(i,j,k,:) = 0._dp
           endif

        enddo
     enddo
  enddo
end subroutine ICCalcPertGRHydro
