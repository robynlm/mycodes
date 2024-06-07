#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

void CosmoLapse_dtLapse(CCTK_ARGUMENTS) 
  {
    DECLARE_CCTK_PARAMETERS;
    DECLARE_CCTK_ARGUMENTS;
    
    //CCTK_INFO("Activate CosmoLapse");
    
    for (int i = 0; i < cctk_ash[0] * cctk_ash[1] * cctk_ash[2]; i++)
      {
        dtalp[i] = - CosmoLapse_f * pow(alpha[i], CosmoLapse_n) * trK[i];
        A[i] = dtalp[i];
        alpharhs[i] = dtalp[i];
        Arhs[i] = - CosmoLapse_f * ( CosmoLapse_n * dtalp[i] * pow(alpha[i], CosmoLapse_n - 1.0) * trK[i] 
                                     + pow(alpha[i], CosmoLapse_n) * trKrhs[i] );
      }
  }
