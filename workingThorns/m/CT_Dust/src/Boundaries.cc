/*  File produced by Kranc */

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Faces.h"
#include "util_Table.h"
#include "Symmetry.h"
#include "Kranc.hh"


/* the boundary treatment is split into 3 steps:    */
/* 1. excision                                      */
/* 2. symmetries                                    */
/* 3. "other" boundary conditions, e.g. radiative */

/* to simplify scheduling and testing, the 3 steps  */
/* are currently applied in separate functions      */


extern "C" void CT_Dust_CheckBoundaries(CCTK_ARGUMENTS)
{
#ifdef DECLARE_CCTK_ARGUMENTS_CT_Dust_CheckBoundaries
  DECLARE_CCTK_ARGUMENTS_CHECKED(CT_Dust_CheckBoundaries);
#else
  DECLARE_CCTK_ARGUMENTS;
#endif
  DECLARE_CCTK_PARAMETERS;
  
  return;
}

extern "C" void CT_Dust_SelectBoundConds(CCTK_ARGUMENTS)
{
#ifdef DECLARE_CCTK_ARGUMENTS_CT_Dust_SelectBoundConds
  DECLARE_CCTK_ARGUMENTS_CHECKED(CT_Dust_SelectBoundConds);
#else
  DECLARE_CCTK_ARGUMENTS;
#endif
  DECLARE_CCTK_PARAMETERS;
  
  CCTK_INT ierr CCTK_ATTRIBUTE_UNUSED = 0;
  
  if (CCTK_EQUALS(CT_D_bound, "none"  ) ||
      CCTK_EQUALS(CT_D_bound, "static") ||
      CCTK_EQUALS(CT_D_bound, "flat"  ) ||
      CCTK_EQUALS(CT_D_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "CT_Dust::CT_D", CT_D_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register CT_D_bound BC for CT_Dust::CT_D!");
  }
  
  if (CCTK_EQUALS(CT_E_bound, "none"  ) ||
      CCTK_EQUALS(CT_E_bound, "static") ||
      CCTK_EQUALS(CT_E_bound, "flat"  ) ||
      CCTK_EQUALS(CT_E_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "CT_Dust::CT_E", CT_E_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register CT_E_bound BC for CT_Dust::CT_E!");
  }
  
  if (CCTK_EQUALS(CT_S_bound, "none"  ) ||
      CCTK_EQUALS(CT_S_bound, "static") ||
      CCTK_EQUALS(CT_S_bound, "flat"  ) ||
      CCTK_EQUALS(CT_S_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "CT_Dust::CT_S", CT_S_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register CT_S_bound BC for CT_Dust::CT_S!");
  }
  
  if (CCTK_EQUALS(DD_bound, "none"  ) ||
      CCTK_EQUALS(DD_bound, "static") ||
      CCTK_EQUALS(DD_bound, "flat"  ) ||
      CCTK_EQUALS(DD_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "CT_Dust::DD", DD_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register DD_bound BC for CT_Dust::DD!");
  }
  
  if (CCTK_EQUALS(EE_bound, "none"  ) ||
      CCTK_EQUALS(EE_bound, "static") ||
      CCTK_EQUALS(EE_bound, "flat"  ) ||
      CCTK_EQUALS(EE_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "CT_Dust::EE", EE_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register EE_bound BC for CT_Dust::EE!");
  }
  
  if (CCTK_EQUALS(SS1_bound, "none"  ) ||
      CCTK_EQUALS(SS1_bound, "static") ||
      CCTK_EQUALS(SS1_bound, "flat"  ) ||
      CCTK_EQUALS(SS1_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "CT_Dust::SS1", SS1_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register SS1_bound BC for CT_Dust::SS1!");
  }
  
  if (CCTK_EQUALS(SS2_bound, "none"  ) ||
      CCTK_EQUALS(SS2_bound, "static") ||
      CCTK_EQUALS(SS2_bound, "flat"  ) ||
      CCTK_EQUALS(SS2_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "CT_Dust::SS2", SS2_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register SS2_bound BC for CT_Dust::SS2!");
  }
  
  if (CCTK_EQUALS(SS3_bound, "none"  ) ||
      CCTK_EQUALS(SS3_bound, "static") ||
      CCTK_EQUALS(SS3_bound, "flat"  ) ||
      CCTK_EQUALS(SS3_bound, "zero"  ) )
  {
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, -1,
                      "CT_Dust::SS3", SS3_bound);
    if (ierr < 0)
       CCTK_ERROR("Failed to register SS3_bound BC for CT_Dust::SS3!");
  }
  
  if (CCTK_EQUALS(CT_D_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_CT_D_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_CT_D_bound < 0) handle_CT_D_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_CT_D_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_CT_D_bound , CT_D_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_CT_D_bound ,CT_D_bound_speed, "SPEED") < 0)
       CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_CT_D_bound, 
                      "CT_Dust::CT_D", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for CT_Dust::CT_D!");
  
  }
  
  if (CCTK_EQUALS(CT_E_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_CT_E_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_CT_E_bound < 0) handle_CT_E_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_CT_E_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_CT_E_bound , CT_E_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_CT_E_bound ,CT_E_bound_speed, "SPEED") < 0)
       CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_CT_E_bound, 
                      "CT_Dust::CT_E", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for CT_Dust::CT_E!");
  
  }
  
  if (CCTK_EQUALS(CT_S_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_CT_S_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_CT_S_bound < 0) handle_CT_S_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_CT_S_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_CT_S_bound , CT_S_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_CT_S_bound ,CT_S_bound_speed, "SPEED") < 0)
       CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_CT_S_bound, 
                      "CT_Dust::CT_S", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for CT_Dust::CT_S!");
  
  }
  
  if (CCTK_EQUALS(DD_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_DD_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_DD_bound < 0) handle_DD_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_DD_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_DD_bound , DD_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_DD_bound ,DD_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_DD_bound, 
                      "CT_Dust::DD", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for CT_Dust::DD!");
  
  }
  
  if (CCTK_EQUALS(EE_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_EE_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_EE_bound < 0) handle_EE_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_EE_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_EE_bound , EE_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_EE_bound ,EE_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_EE_bound, 
                      "CT_Dust::EE", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for CT_Dust::EE!");
  
  }
  
  if (CCTK_EQUALS(SS1_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_SS1_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_SS1_bound < 0) handle_SS1_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_SS1_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_SS1_bound , SS1_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_SS1_bound ,SS1_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_SS1_bound, 
                      "CT_Dust::SS1", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for CT_Dust::SS1!");
  
  }
  
  if (CCTK_EQUALS(SS2_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_SS2_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_SS2_bound < 0) handle_SS2_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_SS2_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_SS2_bound , SS2_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_SS2_bound ,SS2_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_SS2_bound, 
                      "CT_Dust::SS2", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for CT_Dust::SS2!");
  
  }
  
  if (CCTK_EQUALS(SS3_bound, "radiative"))
  {
   /* select radiation boundary condition */
    static CCTK_INT handle_SS3_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_SS3_bound < 0) handle_SS3_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_SS3_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_SS3_bound , SS3_bound_limit, "LIMIT") < 0)
       CCTK_ERROR("could not set LIMIT value in table!");
    if (Util_TableSetReal(handle_SS3_bound ,SS3_bound_speed, "SPEED") < 0)
        CCTK_ERROR("could not set SPEED value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_SS3_bound, 
                      "CT_Dust::SS3", "Radiation");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Radiation BC for CT_Dust::SS3!");
  
  }
  
  if (CCTK_EQUALS(CT_D_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_CT_D_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_CT_D_bound < 0) handle_CT_D_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_CT_D_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_CT_D_bound ,CT_D_bound_scalar, "SCALAR") < 0)
        CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_CT_D_bound, 
                      "CT_Dust::CT_D", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Scalar BC for CT_Dust::CT_D!");
  
  }
  
  if (CCTK_EQUALS(CT_E_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_CT_E_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_CT_E_bound < 0) handle_CT_E_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_CT_E_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_CT_E_bound ,CT_E_bound_scalar, "SCALAR") < 0)
        CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_CT_E_bound, 
                      "CT_Dust::CT_E", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Scalar BC for CT_Dust::CT_E!");
  
  }
  
  if (CCTK_EQUALS(CT_S_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_CT_S_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_CT_S_bound < 0) handle_CT_S_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_CT_S_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_CT_S_bound ,CT_S_bound_scalar, "SCALAR") < 0)
        CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, 1, handle_CT_S_bound, 
                      "CT_Dust::CT_S", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Failed to register Scalar BC for CT_Dust::CT_S!");
  
  }
  
  if (CCTK_EQUALS(DD_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_DD_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_DD_bound < 0) handle_DD_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_DD_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_DD_bound ,DD_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_DD_bound, 
                      "CT_Dust::DD", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for CT_Dust::DD!");
  
  }
  
  if (CCTK_EQUALS(EE_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_EE_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_EE_bound < 0) handle_EE_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_EE_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_EE_bound ,EE_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_EE_bound, 
                      "CT_Dust::EE", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for CT_Dust::EE!");
  
  }
  
  if (CCTK_EQUALS(SS1_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_SS1_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_SS1_bound < 0) handle_SS1_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_SS1_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_SS1_bound ,SS1_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_SS1_bound, 
                      "CT_Dust::SS1", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for CT_Dust::SS1!");
  
  }
  
  if (CCTK_EQUALS(SS2_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_SS2_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_SS2_bound < 0) handle_SS2_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_SS2_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_SS2_bound ,SS2_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_SS2_bound, 
                      "CT_Dust::SS2", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for CT_Dust::SS2!");
  
  }
  
  if (CCTK_EQUALS(SS3_bound, "scalar"))
  {
   /* select scalar boundary condition */
    static CCTK_INT handle_SS3_bound CCTK_ATTRIBUTE_UNUSED = -1;
    if (handle_SS3_bound < 0) handle_SS3_bound = Util_TableCreate(UTIL_TABLE_FLAGS_CASE_INSENSITIVE);
    if (handle_SS3_bound < 0) CCTK_ERROR("could not create table!");
    if (Util_TableSetReal(handle_SS3_bound ,SS3_bound_scalar, "SCALAR") < 0)
      CCTK_ERROR("could not set SCALAR value in table!");
  
    ierr = KrancBdy_SelectVarForBC(cctkGH, CCTK_ALL_FACES, 1, handle_SS3_bound, 
                      "CT_Dust::SS3", "scalar");
  
    if (ierr < 0)
       CCTK_ERROR("Error in registering Scalar BC for CT_Dust::SS3!");
  
  }
  return;
}



/* template for entries in parameter file:
#$bound$#CT_Dust::CT_D_bound       = "skip"
#$bound$#CT_Dust::CT_D_bound_speed = 1.0
#$bound$#CT_Dust::CT_D_bound_limit = 0.0
#$bound$#CT_Dust::CT_D_bound_scalar = 0.0

#$bound$#CT_Dust::CT_E_bound       = "skip"
#$bound$#CT_Dust::CT_E_bound_speed = 1.0
#$bound$#CT_Dust::CT_E_bound_limit = 0.0
#$bound$#CT_Dust::CT_E_bound_scalar = 0.0

#$bound$#CT_Dust::CT_S_bound       = "skip"
#$bound$#CT_Dust::CT_S_bound_speed = 1.0
#$bound$#CT_Dust::CT_S_bound_limit = 0.0
#$bound$#CT_Dust::CT_S_bound_scalar = 0.0

#$bound$#CT_Dust::DD_bound       = "skip"
#$bound$#CT_Dust::DD_bound_speed = 1.0
#$bound$#CT_Dust::DD_bound_limit = 0.0
#$bound$#CT_Dust::DD_bound_scalar = 0.0

#$bound$#CT_Dust::EE_bound       = "skip"
#$bound$#CT_Dust::EE_bound_speed = 1.0
#$bound$#CT_Dust::EE_bound_limit = 0.0
#$bound$#CT_Dust::EE_bound_scalar = 0.0

#$bound$#CT_Dust::SS1_bound       = "skip"
#$bound$#CT_Dust::SS1_bound_speed = 1.0
#$bound$#CT_Dust::SS1_bound_limit = 0.0
#$bound$#CT_Dust::SS1_bound_scalar = 0.0

#$bound$#CT_Dust::SS2_bound       = "skip"
#$bound$#CT_Dust::SS2_bound_speed = 1.0
#$bound$#CT_Dust::SS2_bound_limit = 0.0
#$bound$#CT_Dust::SS2_bound_scalar = 0.0

#$bound$#CT_Dust::SS3_bound       = "skip"
#$bound$#CT_Dust::SS3_bound_speed = 1.0
#$bound$#CT_Dust::SS3_bound_limit = 0.0
#$bound$#CT_Dust::SS3_bound_scalar = 0.0

*/

