#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

module DefineRcGRH

  implicit none
  contains

    subroutine ICPertFLRW_GRH_DefineRc ( x, y, z, &
                                     Rc, dxRc, dyRc, dzRc, &
                                     dxdxRc, dxdyRc, dxdzRc, dydyRc, dydzRc, dzdzRc, &
                                     dxdxdxRc, dxdxdyRc, dxdxdzRc, dxdydyRc, dxdydzRc, dxdzdzRc, &
                                     dydydyRc, dydydzRc, dydzdzRc, dzdzdzRc, &
                                     dxdxdxdxRc, dxdxdxdyRc, dxdxdxdzRc, dxdxdydyRc, dxdxdydzRc, &
                                     dxdxdzdzRc, dxdydydyRc, dxdydydzRc, dxdydzdzRc, dxdzdzdzRc, &
                                     dydydydyRc, dydydydzRc, dydydzdzRc, dydzdzdzRc, dzdzdzdzRc )
      DECLARE_CCTK_PARAMETERS
      DECLARE_CCTK_FUNCTIONS
      ! Input
      CCTK_REAL :: x, y, z

      ! Output
      CCTK_REAL :: Rc, dxRc, dyRc, dzRc
      CCTK_REAL :: dxdxRc, dxdyRc, dxdzRc, dydyRc, dydzRc, dzdzRc
      CCTK_REAL :: dxdxdxRc, dxdxdyRc, dxdxdzRc, dxdydyRc, dxdydzRc, dxdzdzRc
      CCTK_REAL :: dydydyRc, dydydzRc, dydzdzRc, dzdzdzRc
      CCTK_REAL :: dxdxdxdxRc, dxdxdxdyRc, dxdxdxdzRc, dxdxdydyRc, dxdxdydzRc
      CCTK_REAL :: dxdxdzdzRc, dxdydydyRc, dxdydydzRc, dxdydzdzRc, dxdzdzdzRc
      CCTK_REAL :: dydydydyRc, dydydydzRc, dydydzdzRc, dydzdzdzRc, dzdzdzdzRc

      ! Local variables
      integer, parameter :: dp = 8
      integer :: m
      real(dp), parameter :: pi = 4._dp*atan(1._dp)
      logical :: SinRc, SpinRc
      CCTK_REAL :: twopi, kx, ky, kz
      CCTK_REAL :: sinx, siny, sinz
      CCTK_REAL :: cosx, cosy, cosz
      twopi = 2._dp * pi

      SinRc = CCTK_EQUALS (ICPertFLRW_GRH_Rcprofile, "sin")
      SpinRc = CCTK_EQUALS (ICPertFLRW_GRH_Rcprofile, "spin")

      ! Initialise them to zero
      Rc = 0._dp
      ! Single derivatives
      dxRc = 0._dp
      dyRc = 0._dp
      dzRc = 0._dp
      ! Double derivatives
      dxdxRc = 0._dp
      dxdyRc = 0._dp
      dxdzRc = 0._dp
      dydyRc = 0._dp
      dydzRc = 0._dp
      dzdzRc = 0._dp
      ! Triple derivatives
      dxdxdxRc = 0._dp
      dxdxdyRc = 0._dp
      dxdxdzRc = 0._dp
      dxdydyRc = 0._dp
      dxdydzRc = 0._dp
      dxdzdzRc = 0._dp
      dydydyRc = 0._dp
      dydydzRc = 0._dp
      dydzdzRc = 0._dp
      dzdzdzRc = 0._dp
      ! Quadruple derivatives
      dxdxdxdxRc = 0._dp
      dxdxdxdyRc = 0._dp
      dxdxdxdzRc = 0._dp
      dxdxdydyRc = 0._dp
      dxdxdydzRc = 0._dp
      dxdxdzdzRc = 0._dp
      dxdydydyRc = 0._dp
      dxdydydzRc = 0._dp
      dxdydzdzRc = 0._dp
      dxdzdzdzRc = 0._dp
      dydydydyRc = 0._dp
      dydydydzRc = 0._dp
      dydydzdzRc = 0._dp
      dydzdzdzRc = 0._dp
      dzdzdzdzRc = 0._dp

      if (SinRc) then
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !           Sinusoidal
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do m = 1, 20        
          kx = twopi / ICPertFLRW_GRH_lambda_x(m)
          ky = twopi / ICPertFLRW_GRH_lambda_y(m)
          kz = twopi / ICPertFLRW_GRH_lambda_z(m)

          sinx = sin(x * kx + ICPertFLRW_GRH_phi_x(m))
          siny = sin(y * ky + ICPertFLRW_GRH_phi_y(m))
          sinz = sin(z * kz + ICPertFLRW_GRH_phi_z(m))
          cosx = cos(x * kx + ICPertFLRW_GRH_phi_x(m))
          cosy = cos(y * ky + ICPertFLRW_GRH_phi_y(m))
          cosz = cos(z * kz + ICPertFLRW_GRH_phi_z(m))
  
          Rc = Rc + ICPertFLRW_GRH_Amp_x(m) * sinx &
                  + ICPertFLRW_GRH_Amp_y(m) * siny &
                  + ICPertFLRW_GRH_Amp_z(m) * sinz

          ! Single derivatives
          dxRc = dxRc + ICPertFLRW_GRH_Amp_x(m) * kx * cosx
          dyRc = dyRc + ICPertFLRW_GRH_Amp_y(m) * ky * cosy
          dzRc = dzRc + ICPertFLRW_GRH_Amp_z(m) * kz * cosz

          ! Double derivatives
          dxdxRc = dxdxRc - ICPertFLRW_GRH_Amp_x(m) * (kx**2._dp) * sinx
          dydyRc = dydyRc - ICPertFLRW_GRH_Amp_y(m) * (ky**2._dp) * siny
          dzdzRc = dzdzRc - ICPertFLRW_GRH_Amp_z(m) * (kz**2._dp) * sinz

          ! Triple derivatives
          dxdxdxRc = dxdxdxRc - ICPertFLRW_GRH_Amp_x(m) * (kx**3._dp) * cosx
          dydydyRc = dydydyRc - ICPertFLRW_GRH_Amp_y(m) * (ky**3._dp) * cosy
          dzdzdzRc = dzdzdzRc - ICPertFLRW_GRH_Amp_z(m) * (kz**3._dp) * cosz
 
          ! Quadruple derivatives
          dxdxdxdxRc = dxdxdxdxRc + ICPertFLRW_GRH_Amp_x(m) * (kx**4._dp) * sinx
          dydydydyRc = dydydydyRc + ICPertFLRW_GRH_Amp_y(m) * (ky**4._dp) * siny
          dzdzdzdzRc = dzdzdzdzRc + ICPertFLRW_GRH_Amp_z(m) * (kz**4._dp) * sinz
        enddo
      else if (SpinRc) then
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !           Cosine Spin
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        kx = twopi / ICPertFLRW_GRH_lambda_x(1)
        ky = twopi / ICPertFLRW_GRH_lambda_y(1)
        kz = twopi / ICPertFLRW_GRH_lambda_z(1)

        cosx = cos(x * kx)
        cosy = cos((x - y) * ky)
        cosz = cos(z * kz)
        sinx = sin(x * kx)
        siny = sin((x - y) * ky)
        sinz = sin(z * kz)

        Rc = ICPertFLRW_GRH_Amp_x(1) * cosx &
             + ICPertFLRW_GRH_Amp_y(1) * cosy &
             + ICPertFLRW_GRH_Amp_z(1) * cosz

        ! Single derivatives
        dxRc = - ICPertFLRW_GRH_Amp_x(1) * kx * sinx &
               - ICPertFLRW_GRH_Amp_y(1) * ky * siny
        dyRc = ICPertFLRW_GRH_Amp_y(1) * ky * siny
        dzRc = - ICPertFLRW_GRH_Amp_z(1) * kz * sinz

        ! Double derivatives
        dxdxRc = - ICPertFLRW_GRH_Amp_x(1) * (kx**2._dp) * cosx &
                 - ICPertFLRW_GRH_Amp_y(1) * (ky**2._dp) * cosy
        dxdyRc = ICPertFLRW_GRH_Amp_y(1) * (ky**2._dp) * cosy
        dydyRc = - ICPertFLRW_GRH_Amp_y(1) * (ky**2._dp) * cosy
        dzdzRc = - ICPertFLRW_GRH_Amp_z(1) * (kz**2._dp) * cosz

        ! Triple derivatives
        dxdxdxRc = ICPertFLRW_GRH_Amp_x(1) * (kx**3._dp) * sinx &
                   + ICPertFLRW_GRH_Amp_y(1) * (ky**3._dp) * siny
        dxdxdyRc = - ICPertFLRW_GRH_Amp_y(1) * (ky**3._dp) * siny
        dxdydyRc = ICPertFLRW_GRH_Amp_y(1) * (ky**3._dp) * siny
        dydydyRc = - ICPertFLRW_GRH_Amp_y(1) * (ky**3._dp) * siny
        dzdzdzRc = ICPertFLRW_GRH_Amp_z(1) * (kz**3._dp) * sinz

        ! Quadruple derivatives
        dxdxdxdxRc = ICPertFLRW_GRH_Amp_x(1) * (kx**4._dp) * cosx &
                     + ICPertFLRW_GRH_Amp_y(1) * (ky**4._dp) * cosy
        dxdxdxdyRc = - ICPertFLRW_GRH_Amp_y(1) * (ky**4._dp) * cosy
        dxdxdydyRc = ICPertFLRW_GRH_Amp_y(1) * (ky**4._dp) * cosy
        dxdydydyRc = - ICPertFLRW_GRH_Amp_y(1) * (ky**4._dp) * cosy
        dydydydyRc = ICPertFLRW_GRH_Amp_y(1) * (ky**4._dp) * cosy
        dzdzdzdzRc = ICPertFLRW_GRH_Amp_z(1) * (kz**4._dp) * cosz
      endif
  end subroutine ICPertFLRW_GRH_DefineRc
end module DefineRcGRH
