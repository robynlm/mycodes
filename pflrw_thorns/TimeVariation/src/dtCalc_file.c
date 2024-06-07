#include <math.h>
#include <stdlib.h>
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

void dtCalc(CCTK_ARGUMENTS) {
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS;
  
  CCTK_REAL H0 = ICPertFLRW_h * ICPertFLRW_c / 2997.9; // Units are Mpc
  CCTK_REAL a0 = 1.0;
  CCTK_REAL t0EdS = 2.0 / ( 3.0 * H0 );
    
  if (ICPertFLRW_Lambda == "yes") {
    CCTK_REAL a_varying = a0 * pow(ICPertFLRW_Omega_matter0/(1.0 - ICPertFLRW_Omega_matter0), 1.0/3.0) * pow(sinh( sqrt(1.0 - ICPertFLRW_Omega_matter0) * cctkGH->cctk_time/t0EdS), 2.0/3.0);
    CCTK_REAL a_constant = a0 * pow(ICPertFLRW_Omega_matter0/(1.0 - ICPertFLRW_Omega_matter0), 1.0/3.0) * pow(sinh( sqrt(1.0 - ICPertFLRW_Omega_matter0) * dtfac/t0EdS), 2.0/3.0);
  }  else  {
    CCTK_REAL a_varying = a0 * pow(cctkGH->cctk_time / t0EdS, 2.0/3.0);
    CCTK_REAL a_constant = a0 * pow(dtfac / t0EdS, 2.0/3.0);
  }
      
  CCTK_REAL courant_min_time_varying   = 3.0 * a_varying * cctk_delta_space[0];
  CCTK_REAL dt_varying  = courant_fac * courant_min_time_varying / sqrt(3.0);

  CCTK_REAL courant_min_time_constant  = 3.0 * a_constant * cctk_delta_space[0];
  CCTK_REAL dt_constant = courant_fac * courant_min_time_constant / sqrt(3.0);
  
  if (verbose) {
    CCTK_VInfo(CCTK_THORNSTRING, "dt_varying  %g", dt_varying);
    CCTK_VInfo(CCTK_THORNSTRING, "dt_constant %g", dt_constant);
  }
  

  if (dt_varying < dt_constant) {
    *courant_min_time = courant_min_time_varying;
  }  else  {
    *courant_min_time = courant_min_time_constant;
  }
}
