/*  File produced by Kranc */

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

extern "C" void CT_Dust_RegisterVars(CCTK_ARGUMENTS)
{
  #ifdef DECLARE_CCTK_ARGUMENTS_CT_Dust_RegisterVars
  DECLARE_CCTK_ARGUMENTS_CHECKED(CT_Dust_RegisterVars);
  #else
  DECLARE_CCTK_ARGUMENTS;
  #endif
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT ierr CCTK_ATTRIBUTE_UNUSED = 0;
  /* Register all the evolved grid functions with MoL */
  ierr += MoLRegisterEvolved(CCTK_VarIndex("CT_Dust::DD"),  CCTK_VarIndex("CT_Dust::DDrhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("CT_Dust::EE"),  CCTK_VarIndex("CT_Dust::EErhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("CT_Dust::SS1"),  CCTK_VarIndex("CT_Dust::SS1rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("CT_Dust::SS2"),  CCTK_VarIndex("CT_Dust::SS2rhs"));
  ierr += MoLRegisterEvolved(CCTK_VarIndex("CT_Dust::SS3"),  CCTK_VarIndex("CT_Dust::SS3rhs"));
  /* Register all the evolved Array functions with MoL */
  return;
}
