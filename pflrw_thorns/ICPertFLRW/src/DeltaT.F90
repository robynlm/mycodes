#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

subroutine ICPertFLRW_DeltaT (CCTK_ARGUMENTS)
  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_FUNCTIONS
  DECLARE_CCTK_PARAMETERS
  integer :: i,j,k
  integer, parameter :: dp = 8
  real(dp), parameter :: pi = 4._dp*atan(1._dp)
  CCTK_REAL :: H0, a0, tau0_EdS, Omega_lambda0
  CCTK_REAL :: tau, aa, aa0, redshift
  logical :: want_pflrw, want_Lambda, want_radiation, want_proper_time

  character :: msg*1000

  want_pflrw = CCTK_EQUALS (my_initial_data, "ICPertFLRW")
  want_Lambda = CCTK_EQUALS (ICPertFLRW_Lambda, "yes")
  want_radiation = CCTK_EQUALS (ICPertFLRW_type_of_matter, "radiation")
  want_proper_time = CCTK_EQUALS (ICPertFLRW_time, "proper")

  if (ICPertFLRW_need_to_fix_dtfac) then
    ! This condition is so this calculation is
    ! not repeated when using mesh refinement
    if (want_proper_time) then
      if (want_pflrw) then

        ! Define the pressure to energy density ratio
        if (want_radiation) then
          eosw = 1._dp / 3._dp
        else
          eosw = 0._dp
        endif

        ! Cosmological measurements today
        H0 = ICPertFLRW_h * ICPertFLRW_c / 2997.9_dp ! Units are Mpc
        tau0_EdS = 2._dp / ( 3._dp * H0 * (1._dp + eosw) )
        Omega_lambda0 = 1._dp - ICPertFLRW_Omega_matter0

        ! Define scale factor in terms of the reference redshift
        if (ICPertFLRW_z_comoving_ref .LT. 0._dp) then
          aa = 1._dp
        else
          a0 = 1._dp + ICPertFLRW_z_comoving_ref
        endif

        tau = cctk_time
        if (want_Lambda) then ! Lambda CDM model
          aa0 = (( ICPertFLRW_Omega_matter0 / Omega_lambda0 )**(1._dp/3._dp) &
                 * sinh( sqrt(Omega_lambda0) * tau / tau0_EdS )**(2._dp/3._dp))
          redshift = (1._dp / aa0) - 1._dp
        else ! Einstein - de Sitter model
          aa0 = ( tau / tau0_EdS )**(2._dp/(3._dp * (1._dp + eosw)))
          redshift = (1._dp / aa0) - 1._dp
        endif


        if (ICPertFLRW_z_comoving_ref .LT. 0._dp) then
          a0 = aa * (1._dp + redshift)
        else
          aa = a0 / (1._dp + redshift)
        endif

        ! This overwrites the Time::dtfac parameter
        call CCTK_INFO("Define time spacing: dtfac = dtfac * scale factor")
        dtfac = dtfac * aa
        write (msg, '("The initial scale factor is ",g16.6)') aa
        call CCTK_INFO (msg)
        write (msg, '("dtfac is ",g16.6)') dtfac
        call CCTK_INFO (msg)

        ICPertFLRW_need_to_fix_dtfac = .false.        

      endif  !if you want pflrw
    endif  !if you want proper time
  endif !if ICPertFLRW_need_to_fix_dtfac
end subroutine ICPertFLRW_DeltaT




