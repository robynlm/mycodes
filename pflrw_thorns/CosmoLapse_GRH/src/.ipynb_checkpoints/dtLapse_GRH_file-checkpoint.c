#include <math.h>
#include <stdlib.h>
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

void dtLapse_GRH(CCTK_ARGUMENTS) {
  DECLARE_CCTK_PARAMETERS;
  DECLARE_CCTK_ARGUMENTS;
 
  CCTK_REAL Kterm, dtKterm; 

    if (CCTK_EQUALS(CosmoLapse_GRH_Kexpression, "initial_background_K")) {
        Kterm = - 2.0 * pow(cctk_initial_time, -1);
        dtKterm = 0.0;
    } else {
        Kterm = 0.0;
        dtKterm = 0.0;
    }
    
    for (int i = 0; i < cctk_ash[0] * cctk_ash[1] * cctk_ash[2]; i++) {
        if (CCTK_EQUALS(CosmoLapse_GRH_Kexpression, "background_K")) {
            Kterm = - 2.0 * pow(tau[i], -1);
            dtKterm = 2.0 * pow(tau[i], -2);
        }

        dtalp[i] = - harmonicF * pow(alpha[i], harmonicN) * (trK[i] - Kterm);
        A[i] = dtalp[i];
        alpharhs[i] = dtalp[i];
        Arhs[i] = - harmonicF * ( harmonicN * dtalp[i] * pow(alpha[i], harmonicN - 1.0) * (trK[i] - Kterm)
                                  + pow(alpha[i], harmonicN) * (trKrhs[i] - dtKterm) );
    }
}
