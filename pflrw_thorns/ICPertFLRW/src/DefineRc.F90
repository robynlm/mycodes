#include "cctk.h"
!#include "cctk_Arguments.h"
#include "cctk_Functions.h"
#include "cctk_Parameters.h"

module DefineRc

  implicit none
  contains

    subroutine ICPertFLRW_DefineRc ( x, y, z, &
                                     Rc, dxRc, dyRc, dzRc, &
                                     dxdxRc, dxdyRc, dxdzRc, dydyRc, dydzRc, dzdzRc, &
                                     dxdxdxRc, dxdxdyRc, dxdxdzRc, dxdydyRc, dxdydzRc, dxdzdzRc, &
                                     dydydyRc, dydydzRc, dydzdzRc, dzdzdzRc, &
                                     dxdxdxdxRc, dxdxdxdyRc, dxdxdxdzRc, dxdxdydyRc, dxdxdydzRc, &
                                     dxdxdzdzRc, dxdydydyRc, dxdydydzRc, dxdydzdzRc, dxdzdzdzRc, &
                                     dydydydyRc, dydydydzRc, dydydzdzRc, dydzdzdzRc, dzdzdzdzRc )
      !DECLARE_CCTK_ARGUMENTS
      DECLARE_CCTK_FUNCTIONS
      DECLARE_CCTK_PARAMETERS
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
      logical :: SinRc, SpinRc, ExpRc
      ! sin
      integer :: m
      real(dp), parameter :: pi = 4._dp*atan(1._dp)
      CCTK_REAL :: twopi = 2._dp * pi
      CCTK_REAL :: kx, ky, kz
      CCTK_REAL :: sinx, siny, sinz
      CCTK_REAL :: cosx, cosy, cosz
      ! exp
      CCTK_REAL :: tsteepx, tsteepy, tsteepz
      CCTK_REAL :: expfacx, expfacy, expfacz
      CCTK_REAL :: xts, yts, zts
      CCTK_REAL :: dxxts, dyyts, dzzts
      CCTK_REAL :: dxdxxts, dydyyts, dzdzzts
      CCTK_REAL :: dxdxdxxts, dydydyyts, dzdzdzzts
      CCTK_REAL :: dxdxdxdxxts, dydydydyyts, dzdzdzdzzts
      CCTK_REAL :: expx, expy, expz
      CCTK_REAL :: dxexpx, dyexpy, dzexpz
      CCTK_REAL :: dxdxexpx, dydyexpy, dzdzexpz
      CCTK_REAL :: dxdxdxexpx, dydydyexpy, dzdzdzexpz
      CCTK_REAL :: dxdxdxdxexpx, dydydydyexpy, dzdzdzdzexpz

      SinRc = CCTK_EQUALS (ICPertFLRW_Rcprofile, "sin")
      SpinRc = CCTK_EQUALS (ICPertFLRW_Rcprofile, "spin")
      ExpRc = CCTK_EQUALS (ICPertFLRW_Rcprofile, "exp")

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
          kx = twopi / ICPertFLRW_lambda_x(m)
          ky = twopi / ICPertFLRW_lambda_y(m)
          kz = twopi / ICPertFLRW_lambda_z(m)

          sinx = sin(x * kx + ICPertFLRW_phi_x(m))
          siny = sin(y * ky + ICPertFLRW_phi_y(m))
          sinz = sin(z * kz + ICPertFLRW_phi_z(m))
          cosx = cos(x * kx + ICPertFLRW_phi_x(m))
          cosy = cos(y * ky + ICPertFLRW_phi_y(m))
          cosz = cos(z * kz + ICPertFLRW_phi_z(m))

          Rc = Rc + ICPertFLRW_Amp_x(m) * sinx &
                  + ICPertFLRW_Amp_y(m) * siny &
                  + ICPertFLRW_Amp_z(m) * sinz

          ! Single derivatives
          dxRc = dxRc + ICPertFLRW_Amp_x(m) * kx * cosx
          dyRc = dyRc + ICPertFLRW_Amp_y(m) * ky * cosy
          dzRc = dzRc + ICPertFLRW_Amp_z(m) * kz * cosz

          ! Double derivatives
          dxdxRc = dxdxRc - ICPertFLRW_Amp_x(m) * (kx**2._dp) * sinx
          dydyRc = dydyRc - ICPertFLRW_Amp_y(m) * (ky**2._dp) * siny
          dzdzRc = dzdzRc - ICPertFLRW_Amp_z(m) * (kz**2._dp) * sinz

          ! Triple derivatives
          dxdxdxRc = dxdxdxRc - ICPertFLRW_Amp_x(m) * (kx**3._dp) * cosx
          dydydyRc = dydydyRc - ICPertFLRW_Amp_y(m) * (ky**3._dp) * cosy
          dzdzdzRc = dzdzdzRc - ICPertFLRW_Amp_z(m) * (kz**3._dp) * cosz

          ! Quadruple derivatives
          dxdxdxdxRc = dxdxdxdxRc + ICPertFLRW_Amp_x(m) * (kx**4._dp) * sinx
          dydydydyRc = dydydydyRc + ICPertFLRW_Amp_y(m) * (ky**4._dp) * siny
          dzdzdzdzRc = dzdzdzdzRc + ICPertFLRW_Amp_z(m) * (kz**4._dp) * sinz
        enddo
      else if (SpinRc) then
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !           Cosine Spin
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          kx = twopi / ICPertFLRW_lambda_x(1)
          ky = twopi / ICPertFLRW_lambda_y(1)
          kz = twopi / ICPertFLRW_lambda_z(1)

          cosx = cos(x * kx)
          cosy = cos((x - y) * ky)
          cosz = cos(z * kz)
          sinx = sin(x * kx)
          siny = sin((x - y) * ky)
          sinz = sin(z * kz)

          Rc = ICPertFLRW_Amp_x(1) * cosx &
               + ICPertFLRW_Amp_y(1) * cosy &
               + ICPertFLRW_Amp_z(1) * cosz

          ! Single derivatives
          dxRc = - ICPertFLRW_Amp_x(1) * kx * sinx &
                 - ICPertFLRW_Amp_y(1) * ky * siny
          dyRc = ICPertFLRW_Amp_y(1) * ky * siny
          dzRc = - ICPertFLRW_Amp_z(1) * kz * sinz

          ! Double derivatives
          dxdxRc = - ICPertFLRW_Amp_x(1) * (kx**2._dp) * cosx &
                   - ICPertFLRW_Amp_y(1) * (ky**2._dp) * cosy
          dxdyRc = ICPertFLRW_Amp_y(1) * (ky**2._dp) * cosy
          dydyRc = - ICPertFLRW_Amp_y(1) * (ky**2._dp) * cosy
          dzdzRc = - ICPertFLRW_Amp_z(1) * (kz**2._dp) * cosz

          ! Triple derivatives
          dxdxdxRc = ICPertFLRW_Amp_x(1) * (kx**3._dp) * sinx &
                     + ICPertFLRW_Amp_y(1) * (ky**3._dp) * siny
          dxdxdyRc = - ICPertFLRW_Amp_y(1) * (ky**3._dp) * siny
          dxdydyRc = ICPertFLRW_Amp_y(1) * (ky**3._dp) * siny
          dydydyRc = - ICPertFLRW_Amp_y(1) * (ky**3._dp) * siny
          dzdzdzRc = ICPertFLRW_Amp_z(1) * (kz**3._dp) * sinz

          ! Quadruple derivatives
          dxdxdxdxRc = ICPertFLRW_Amp_x(1) * (kx**4._dp) * cosx &
                       + ICPertFLRW_Amp_y(1) * (ky**4._dp) * cosy
          dxdxdxdyRc = - ICPertFLRW_Amp_y(1) * (ky**4._dp) * cosy
          dxdxdydyRc = ICPertFLRW_Amp_y(1) * (ky**4._dp) * cosy
          dxdydydyRc = - ICPertFLRW_Amp_y(1) * (ky**4._dp) * cosy
          dydydydyRc = ICPertFLRW_Amp_y(1) * (ky**4._dp) * cosy
          dzdzdzdzRc = ICPertFLRW_Amp_z(1) * (kz**4._dp) * cosz

      else if (ExpRc) then
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          !           Exponential
          !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          ! Comoving curvature perturbation : Rc
          tsteepx = 2._dp * ICPertFLRW_steepness_x
          tsteepy = 2._dp * ICPertFLRW_steepness_y
          tsteepz = 2._dp * ICPertFLRW_steepness_z
          expfacx = (- 1._dp / 2._dp) * (ICPertFLRW_variance_x) ** (- tsteepx)
          expfacy = (- 1._dp / 2._dp) * (ICPertFLRW_variance_y) ** (- tsteepy)
          expfacz = (- 1._dp / 2._dp) * (ICPertFLRW_variance_z) ** (- tsteepz)
          xts = x ** tsteepx
          yts = y ** tsteepy
          zts = z ** tsteepz
          expx = exp( expfacx * xts)
          expy = exp( expfacy * yts)
          expz = exp( expfacz * zts)
          Rc = ICPertFLRW_exp_amplitude * expx * expy * expz

          ! x^2s derivarives
          if (tsteepz >= 2._dp) then
              dxxts = tsteepx * (x ** (tsteepx - 1._dp))
              dxdxxts = tsteepx * (tsteepx - 1._dp) * (x ** (tsteepx - 2._dp))
              if (tsteepz >= 4._dp) then
                  dxdxdxxts = (tsteepx * (tsteepx - 1._dp) * (tsteepx - 2._dp)    &
                               * (x ** (tsteepx - 3._dp)))
                  dxdxdxdxxts = (tsteepx * (tsteepx - 1._dp) * (tsteepx - 2._dp)  &
                                 * (tsteepx - 3._dp) * (x ** (tsteepx - 4._dp)))
              else if (tsteepz < 4._dp) then
                  dxdxdxxts = 0._dp
                  dxdxdxdxxts = 0._dp
              endif
          else if (tsteepz < 2._dp) then
              dxxts = 0._dp
              dxdxxts = 0._dp
              dxdxdxxts = 0._dp
              dxdxdxdxxts = 0._dp
          endif

          ! y^2s derivarives
          if (tsteepz >= 2._dp) then
              dyyts = tsteepy * (y ** (tsteepy - 1._dp))
              dydyyts = tsteepy * (tsteepy - 1._dp) * (y ** (tsteepy - 2._dp))
              if (tsteepz >= 4._dp) then
                  dydydyyts = (tsteepy * (tsteepy - 1._dp) * (tsteepy - 2._dp)    &
                               * (y ** (tsteepy - 3._dp)))
                  dydydydyyts = (tsteepy * (tsteepy - 1._dp) * (tsteepy - 2._dp)  &
                                 * (tsteepy - 3._dp) * (y ** (tsteepy - 4._dp)))
              else if (tsteepz < 4._dp) then
                  dydydyyts = 0._dp
                  dydydydyyts = 0._dp
              endif
          else if (tsteepz < 2._dp) then
              dyyts = 0._dp
              dydyyts = 0._dp
              dydydyyts = 0._dp
              dydydydyyts = 0._dp
          endif

          ! z^2s derivarives
          if (tsteepz >= 2._dp) then
              dzzts = tsteepz * (z ** (tsteepz - 1._dp))
              dzdzzts = tsteepz * (tsteepz - 1._dp) * (z ** (tsteepz - 2._dp))
              if (tsteepz >= 4._dp) then
                  dzdzdzzts = (tsteepz * (tsteepz - 1._dp) * (tsteepz - 2._dp)    &
                               * (z ** (tsteepz - 3._dp)))
                  dzdzdzdzzts = (tsteepz * (tsteepz - 1._dp) * (tsteepz - 2._dp)  &
                                 * (tsteepz - 3._dp) * (z ** (tsteepz - 4._dp)))
              else if (tsteepz < 4._dp) then
                  dzdzdzzts = 0._dp
                  dzdzdzdzzts = 0._dp
              endif
          else if (tsteepz < 2._dp) then
              dzzts = 0._dp
              dzdzzts = 0._dp
              dzdzdzzts = 0._dp
              dzdzdzdzzts = 0._dp
          endif
          
          ! Single derivatives
          dxexpx = expfacx * dxxts * expx
          dyexpy = expfacy * dyyts * expy
          dzexpz = expfacz * dzzts * expz
          dxRc = ICPertFLRW_exp_amplitude * dxexpx * expy * expz
          dyRc = ICPertFLRW_exp_amplitude * expx * dyexpy * expz
          dzRc = ICPertFLRW_exp_amplitude * expx * expy * dzexpz

          ! Double derivatives
          dxdxexpx = expfacx * (dxdxxts * expx + dxxts * dxexpx)
          dydyexpy = expfacy * (dydyyts * expy + dyyts * dyexpy)
          dzdzexpz = expfacz * (dzdzzts * expz + dzzts * dzexpz)
          dxdxRc = ICPertFLRW_exp_amplitude * dxdxexpx * expy * expz
          dxdyRc = ICPertFLRW_exp_amplitude * dxexpx * dyexpy * expz
          dxdzRc = ICPertFLRW_exp_amplitude * dxexpx * expy * dzexpz
          dydyRc = ICPertFLRW_exp_amplitude * expx * dydyexpy * expz
          dydzRc = ICPertFLRW_exp_amplitude * expx * dyexpy * dzexpz
          dzdzRc = ICPertFLRW_exp_amplitude * expx * expy * dzdzexpz

          ! Triple derivatives
          dxdxdxexpx = expfacx * (dxdxdxxts * expx + 2._dp * dxdxxts * dxexpx + dxxts * dxdxexpx)
          dydydyexpy = expfacy * (dydydyyts * expy + 2._dp * dydyyts * dyexpy + dyyts * dydyexpy)
          dzdzdzexpz = expfacz * (dzdzdzzts * expz + 2._dp * dzdzzts * dzexpz + dzzts * dzdzexpz)
          dxdxdxRc = ICPertFLRW_exp_amplitude * dxdxdxexpx * expy * expz
          dxdxdyRc = ICPertFLRW_exp_amplitude * dxdxexpx * dyexpy * expz
          dxdxdzRc = ICPertFLRW_exp_amplitude * dxdxexpx * expy * dzexpz
          dxdydyRc = ICPertFLRW_exp_amplitude * dxexpx * dydyexpy * expz
          dxdydzRc = ICPertFLRW_exp_amplitude * dxexpx * dyexpy * dzexpz
          dxdzdzRc = ICPertFLRW_exp_amplitude * dxexpx * expy * dzdzexpz
          dydydyRc = ICPertFLRW_exp_amplitude * expx * dydydyexpy * expz
          dydydzRc = ICPertFLRW_exp_amplitude * expx * dydyexpy * dzexpz
          dydzdzRc = ICPertFLRW_exp_amplitude * expx * dyexpy * dzdzexpz
          dzdzdzRc = ICPertFLRW_exp_amplitude * expx * expy * dzdzdzexpz

          ! Quadruple derivatives
          dxdxdxdxexpx = expfacx * (dxdxdxdxxts * expx + 3._dp * dxdxdxxts * dxexpx    &
                                    + 3._dp * dxdxxts * dxdxexpx + dxxts * dxdxdxexpx)
          dydydydyexpy = expfacy * (dydydydyyts * expy + 3._dp * dydydyyts * dyexpy    &
                                    + 3._dp * dydyyts * dydyexpy + dyyts * dydydyexpy)
          dzdzdzdzexpz = expfacz * (dzdzdzdzzts * expz + 3._dp * dzdzdzzts * dzexpz    &
                                    + 3._dp * dzdzzts * dzdzexpz + dzzts * dzdzdzexpz)
          dxdxdxdxRc = ICPertFLRW_exp_amplitude * dxdxdxdxexpx * expy * expz
          dxdxdxdyRc = ICPertFLRW_exp_amplitude * dxdxdxexpx * dyexpy * expz
          dxdxdxdzRc = ICPertFLRW_exp_amplitude * dxdxdxexpx * expy * dzexpz
          dxdxdydyRc = ICPertFLRW_exp_amplitude * dxdxexpx * dydyexpy * expz
          dxdxdydzRc = ICPertFLRW_exp_amplitude * dxdxexpx * dyexpy * dzexpz
          dxdxdzdzRc = ICPertFLRW_exp_amplitude * dxdxexpx * expy * dzdzexpz
          dxdydydyRc = ICPertFLRW_exp_amplitude * dxexpx * dydydyexpy * expz
          dxdydydzRc = ICPertFLRW_exp_amplitude * dxexpx * dydyexpy * dzexpz
          dxdydzdzRc = ICPertFLRW_exp_amplitude * dxexpx * dyexpy * dzdzexpz
          dxdzdzdzRc = ICPertFLRW_exp_amplitude * dxexpx * expy * dzdzdzexpz
          dydydydyRc = ICPertFLRW_exp_amplitude * expx * dydydydyexpy * expz
          dydydydzRc = ICPertFLRW_exp_amplitude * expx * dydydyexpy * dzexpz
          dydydzdzRc = ICPertFLRW_exp_amplitude * expx * dydyexpy * dzdzexpz
          dydzdzdzRc = ICPertFLRW_exp_amplitude * expx * dyexpy * dzdzdzexpz
          dzdzdzdzRc = ICPertFLRW_exp_amplitude * expx * expy * dzdzdzdzexpz

      endif
    end subroutine ICPertFLRW_DefineRc
end module DefineRc
