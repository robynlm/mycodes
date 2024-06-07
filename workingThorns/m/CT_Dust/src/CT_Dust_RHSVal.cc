/*  File produced by Kranc */

#define KRANC_C

#include <algorithm>
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "Kranc.hh"
#include "Differencing.h"
#include "loopcontrol.h"

namespace CT_Dust {

extern "C" void CT_Dust_RHSVal_SelectBCs(CCTK_ARGUMENTS)
{
  #ifdef DECLARE_CCTK_ARGUMENTS_CT_Dust_RHSVal_SelectBCs
  DECLARE_CCTK_ARGUMENTS_CHECKED(CT_Dust_RHSVal_SelectBCs);
  #else
  DECLARE_CCTK_ARGUMENTS;
  #endif
  DECLARE_CCTK_PARAMETERS;
  
  if (cctk_iteration % CT_Dust_RHSVal_calc_every != CT_Dust_RHSVal_calc_offset)
    return;
  CCTK_INT ierr CCTK_ATTRIBUTE_UNUSED = 0;
  ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "CT_Dust::CT_Drhs","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for CT_Dust::CT_Drhs.");
  ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "CT_Dust::CT_Erhs","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for CT_Dust::CT_Erhs.");
  ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "CT_Dust::CT_Srhs","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for CT_Dust::CT_Srhs.");
  return;
}

static void CT_Dust_RHSVal_Body(const cGH* restrict const cctkGH, const int dir, const int face, const CCTK_REAL normal[3], const CCTK_REAL tangentA[3], const CCTK_REAL tangentB[3], const int imin[3], const int imax[3], const int n_subblock_gfs, CCTK_REAL* restrict const subblock_gfs[])
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  /* Include user-supplied include files */
  /* Initialise finite differencing variables */
  const ptrdiff_t di CCTK_ATTRIBUTE_UNUSED = 1;
  const ptrdiff_t dj CCTK_ATTRIBUTE_UNUSED = 
    CCTK_GFINDEX3D(cctkGH,0,1,0) - CCTK_GFINDEX3D(cctkGH,0,0,0);
  const ptrdiff_t dk CCTK_ATTRIBUTE_UNUSED = 
    CCTK_GFINDEX3D(cctkGH,0,0,1) - CCTK_GFINDEX3D(cctkGH,0,0,0);
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL) * di;
  const ptrdiff_t cdj CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL) * dj;
  const ptrdiff_t cdk CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL) * dk;
  const ptrdiff_t cctkLbnd1 CCTK_ATTRIBUTE_UNUSED = cctk_lbnd[0];
  const ptrdiff_t cctkLbnd2 CCTK_ATTRIBUTE_UNUSED = cctk_lbnd[1];
  const ptrdiff_t cctkLbnd3 CCTK_ATTRIBUTE_UNUSED = cctk_lbnd[2];
  const CCTK_REAL t CCTK_ATTRIBUTE_UNUSED = cctk_time;
  const CCTK_REAL cctkOriginSpace1 CCTK_ATTRIBUTE_UNUSED = 
    CCTK_ORIGIN_SPACE(0);
  const CCTK_REAL cctkOriginSpace2 CCTK_ATTRIBUTE_UNUSED = 
    CCTK_ORIGIN_SPACE(1);
  const CCTK_REAL cctkOriginSpace3 CCTK_ATTRIBUTE_UNUSED = 
    CCTK_ORIGIN_SPACE(2);
  const CCTK_REAL dt CCTK_ATTRIBUTE_UNUSED = CCTK_DELTA_TIME;
  const CCTK_REAL dx CCTK_ATTRIBUTE_UNUSED = CCTK_DELTA_SPACE(0);
  const CCTK_REAL dy CCTK_ATTRIBUTE_UNUSED = CCTK_DELTA_SPACE(1);
  const CCTK_REAL dz CCTK_ATTRIBUTE_UNUSED = CCTK_DELTA_SPACE(2);
  const CCTK_REAL dxi CCTK_ATTRIBUTE_UNUSED = pow(dx,-1);
  const CCTK_REAL dyi CCTK_ATTRIBUTE_UNUSED = pow(dy,-1);
  const CCTK_REAL dzi CCTK_ATTRIBUTE_UNUSED = pow(dz,-1);
  const CCTK_REAL khalf CCTK_ATTRIBUTE_UNUSED = 0.5;
  const CCTK_REAL kthird CCTK_ATTRIBUTE_UNUSED = 
    0.333333333333333333333333333333;
  const CCTK_REAL ktwothird CCTK_ATTRIBUTE_UNUSED = 
    0.666666666666666666666666666667;
  const CCTK_REAL kfourthird CCTK_ATTRIBUTE_UNUSED = 
    1.33333333333333333333333333333;
  const CCTK_REAL hdxi CCTK_ATTRIBUTE_UNUSED = 0.5*dxi;
  const CCTK_REAL hdyi CCTK_ATTRIBUTE_UNUSED = 0.5*dyi;
  const CCTK_REAL hdzi CCTK_ATTRIBUTE_UNUSED = 0.5*dzi;
  /* Initialize predefined quantities */
  const CCTK_REAL p1o1024dx CCTK_ATTRIBUTE_UNUSED = 0.0009765625*pow(dx,-1);
  const CCTK_REAL p1o1024dy CCTK_ATTRIBUTE_UNUSED = 0.0009765625*pow(dy,-1);
  const CCTK_REAL p1o1024dz CCTK_ATTRIBUTE_UNUSED = 0.0009765625*pow(dz,-1);
  const CCTK_REAL p1o120dx CCTK_ATTRIBUTE_UNUSED = 0.00833333333333333333333333333333*pow(dx,-1);
  const CCTK_REAL p1o120dy CCTK_ATTRIBUTE_UNUSED = 0.00833333333333333333333333333333*pow(dy,-1);
  const CCTK_REAL p1o120dz CCTK_ATTRIBUTE_UNUSED = 0.00833333333333333333333333333333*pow(dz,-1);
  const CCTK_REAL p1o12dx CCTK_ATTRIBUTE_UNUSED = 0.0833333333333333333333333333333*pow(dx,-1);
  const CCTK_REAL p1o12dy CCTK_ATTRIBUTE_UNUSED = 0.0833333333333333333333333333333*pow(dy,-1);
  const CCTK_REAL p1o12dz CCTK_ATTRIBUTE_UNUSED = 0.0833333333333333333333333333333*pow(dz,-1);
  const CCTK_REAL p1o144dxdy CCTK_ATTRIBUTE_UNUSED = 0.00694444444444444444444444444444*pow(dx,-1)*pow(dy,-1);
  const CCTK_REAL p1o144dxdz CCTK_ATTRIBUTE_UNUSED = 0.00694444444444444444444444444444*pow(dx,-1)*pow(dz,-1);
  const CCTK_REAL p1o144dydz CCTK_ATTRIBUTE_UNUSED = 0.00694444444444444444444444444444*pow(dy,-1)*pow(dz,-1);
  const CCTK_REAL p1o1680dx CCTK_ATTRIBUTE_UNUSED = 0.000595238095238095238095238095238*pow(dx,-1);
  const CCTK_REAL p1o1680dy CCTK_ATTRIBUTE_UNUSED = 0.000595238095238095238095238095238*pow(dy,-1);
  const CCTK_REAL p1o1680dz CCTK_ATTRIBUTE_UNUSED = 0.000595238095238095238095238095238*pow(dz,-1);
  const CCTK_REAL p1o16dx CCTK_ATTRIBUTE_UNUSED = 0.0625*pow(dx,-1);
  const CCTK_REAL p1o16dy CCTK_ATTRIBUTE_UNUSED = 0.0625*pow(dy,-1);
  const CCTK_REAL p1o16dz CCTK_ATTRIBUTE_UNUSED = 0.0625*pow(dz,-1);
  const CCTK_REAL p1o180dx2 CCTK_ATTRIBUTE_UNUSED = 0.00555555555555555555555555555556*pow(dx,-2);
  const CCTK_REAL p1o180dy2 CCTK_ATTRIBUTE_UNUSED = 0.00555555555555555555555555555556*pow(dy,-2);
  const CCTK_REAL p1o180dz2 CCTK_ATTRIBUTE_UNUSED = 0.00555555555555555555555555555556*pow(dz,-2);
  const CCTK_REAL p1o24dx CCTK_ATTRIBUTE_UNUSED = 0.0416666666666666666666666666667*pow(dx,-1);
  const CCTK_REAL p1o24dy CCTK_ATTRIBUTE_UNUSED = 0.0416666666666666666666666666667*pow(dy,-1);
  const CCTK_REAL p1o24dz CCTK_ATTRIBUTE_UNUSED = 0.0416666666666666666666666666667*pow(dz,-1);
  const CCTK_REAL p1o256dx CCTK_ATTRIBUTE_UNUSED = 0.00390625*pow(dx,-1);
  const CCTK_REAL p1o256dy CCTK_ATTRIBUTE_UNUSED = 0.00390625*pow(dy,-1);
  const CCTK_REAL p1o256dz CCTK_ATTRIBUTE_UNUSED = 0.00390625*pow(dz,-1);
  const CCTK_REAL p1o2dx CCTK_ATTRIBUTE_UNUSED = 0.5*pow(dx,-1);
  const CCTK_REAL p1o2dy CCTK_ATTRIBUTE_UNUSED = 0.5*pow(dy,-1);
  const CCTK_REAL p1o2dz CCTK_ATTRIBUTE_UNUSED = 0.5*pow(dz,-1);
  const CCTK_REAL p1o3600dxdy CCTK_ATTRIBUTE_UNUSED = 0.000277777777777777777777777777778*pow(dx,-1)*pow(dy,-1);
  const CCTK_REAL p1o3600dxdz CCTK_ATTRIBUTE_UNUSED = 0.000277777777777777777777777777778*pow(dx,-1)*pow(dz,-1);
  const CCTK_REAL p1o3600dydz CCTK_ATTRIBUTE_UNUSED = 0.000277777777777777777777777777778*pow(dy,-1)*pow(dz,-1);
  const CCTK_REAL p1o4dx CCTK_ATTRIBUTE_UNUSED = 0.25*pow(dx,-1);
  const CCTK_REAL p1o4dxdy CCTK_ATTRIBUTE_UNUSED = 0.25*pow(dx,-1)*pow(dy,-1);
  const CCTK_REAL p1o4dxdz CCTK_ATTRIBUTE_UNUSED = 0.25*pow(dx,-1)*pow(dz,-1);
  const CCTK_REAL p1o4dy CCTK_ATTRIBUTE_UNUSED = 0.25*pow(dy,-1);
  const CCTK_REAL p1o4dydz CCTK_ATTRIBUTE_UNUSED = 0.25*pow(dy,-1)*pow(dz,-1);
  const CCTK_REAL p1o4dz CCTK_ATTRIBUTE_UNUSED = 0.25*pow(dz,-1);
  const CCTK_REAL p1o5040dx2 CCTK_ATTRIBUTE_UNUSED = 0.000198412698412698412698412698413*pow(dx,-2);
  const CCTK_REAL p1o5040dy2 CCTK_ATTRIBUTE_UNUSED = 0.000198412698412698412698412698413*pow(dy,-2);
  const CCTK_REAL p1o5040dz2 CCTK_ATTRIBUTE_UNUSED = 0.000198412698412698412698412698413*pow(dz,-2);
  const CCTK_REAL p1o560dx CCTK_ATTRIBUTE_UNUSED = 0.00178571428571428571428571428571*pow(dx,-1);
  const CCTK_REAL p1o560dy CCTK_ATTRIBUTE_UNUSED = 0.00178571428571428571428571428571*pow(dy,-1);
  const CCTK_REAL p1o560dz CCTK_ATTRIBUTE_UNUSED = 0.00178571428571428571428571428571*pow(dz,-1);
  const CCTK_REAL p1o60dx CCTK_ATTRIBUTE_UNUSED = 0.0166666666666666666666666666667*pow(dx,-1);
  const CCTK_REAL p1o60dy CCTK_ATTRIBUTE_UNUSED = 0.0166666666666666666666666666667*pow(dy,-1);
  const CCTK_REAL p1o60dz CCTK_ATTRIBUTE_UNUSED = 0.0166666666666666666666666666667*pow(dz,-1);
  const CCTK_REAL p1o64dx CCTK_ATTRIBUTE_UNUSED = 0.015625*pow(dx,-1);
  const CCTK_REAL p1o64dy CCTK_ATTRIBUTE_UNUSED = 0.015625*pow(dy,-1);
  const CCTK_REAL p1o64dz CCTK_ATTRIBUTE_UNUSED = 0.015625*pow(dz,-1);
  const CCTK_REAL p1o705600dxdy CCTK_ATTRIBUTE_UNUSED = 1.41723356009070294784580498866e-6*pow(dx,-1)*pow(dy,-1);
  const CCTK_REAL p1o705600dxdz CCTK_ATTRIBUTE_UNUSED = 1.41723356009070294784580498866e-6*pow(dx,-1)*pow(dz,-1);
  const CCTK_REAL p1o705600dydz CCTK_ATTRIBUTE_UNUSED = 1.41723356009070294784580498866e-6*pow(dy,-1)*pow(dz,-1);
  const CCTK_REAL p1o840dx CCTK_ATTRIBUTE_UNUSED = 0.00119047619047619047619047619048*pow(dx,-1);
  const CCTK_REAL p1o840dy CCTK_ATTRIBUTE_UNUSED = 0.00119047619047619047619047619048*pow(dy,-1);
  const CCTK_REAL p1o840dz CCTK_ATTRIBUTE_UNUSED = 0.00119047619047619047619047619048*pow(dz,-1);
  const CCTK_REAL p1odx CCTK_ATTRIBUTE_UNUSED = pow(dx,-1);
  const CCTK_REAL p1odx2 CCTK_ATTRIBUTE_UNUSED = pow(dx,-2);
  const CCTK_REAL p1ody CCTK_ATTRIBUTE_UNUSED = pow(dy,-1);
  const CCTK_REAL p1ody2 CCTK_ATTRIBUTE_UNUSED = pow(dy,-2);
  const CCTK_REAL p1odz CCTK_ATTRIBUTE_UNUSED = pow(dz,-1);
  const CCTK_REAL p1odz2 CCTK_ATTRIBUTE_UNUSED = pow(dz,-2);
  const CCTK_REAL pm1o120dx CCTK_ATTRIBUTE_UNUSED = -0.00833333333333333333333333333333*pow(dx,-1);
  const CCTK_REAL pm1o120dy CCTK_ATTRIBUTE_UNUSED = -0.00833333333333333333333333333333*pow(dy,-1);
  const CCTK_REAL pm1o120dz CCTK_ATTRIBUTE_UNUSED = -0.00833333333333333333333333333333*pow(dz,-1);
  const CCTK_REAL pm1o12dx2 CCTK_ATTRIBUTE_UNUSED = -0.0833333333333333333333333333333*pow(dx,-2);
  const CCTK_REAL pm1o12dy2 CCTK_ATTRIBUTE_UNUSED = -0.0833333333333333333333333333333*pow(dy,-2);
  const CCTK_REAL pm1o12dz2 CCTK_ATTRIBUTE_UNUSED = -0.0833333333333333333333333333333*pow(dz,-2);
  const CCTK_REAL pm1o2dx CCTK_ATTRIBUTE_UNUSED = -0.5*pow(dx,-1);
  const CCTK_REAL pm1o2dy CCTK_ATTRIBUTE_UNUSED = -0.5*pow(dy,-1);
  const CCTK_REAL pm1o2dz CCTK_ATTRIBUTE_UNUSED = -0.5*pow(dz,-1);
  const CCTK_REAL pm1o4dx CCTK_ATTRIBUTE_UNUSED = -0.25*pow(dx,-1);
  const CCTK_REAL pm1o4dy CCTK_ATTRIBUTE_UNUSED = -0.25*pow(dy,-1);
  const CCTK_REAL pm1o4dz CCTK_ATTRIBUTE_UNUSED = -0.25*pow(dz,-1);
  const CCTK_REAL pm1o60dx CCTK_ATTRIBUTE_UNUSED = -0.0166666666666666666666666666667*pow(dx,-1);
  const CCTK_REAL pm1o60dy CCTK_ATTRIBUTE_UNUSED = -0.0166666666666666666666666666667*pow(dy,-1);
  const CCTK_REAL pm1o60dz CCTK_ATTRIBUTE_UNUSED = -0.0166666666666666666666666666667*pow(dz,-1);
  const CCTK_REAL pm1o840dx CCTK_ATTRIBUTE_UNUSED = -0.00119047619047619047619047619048*pow(dx,-1);
  const CCTK_REAL pm1o840dy CCTK_ATTRIBUTE_UNUSED = -0.00119047619047619047619047619048*pow(dy,-1);
  const CCTK_REAL pm1o840dz CCTK_ATTRIBUTE_UNUSED = -0.00119047619047619047619047619048*pow(dz,-1);
  /* Jacobian variable pointers */
  const bool usejacobian1 = (!CCTK_IsFunctionAliased("MultiPatch_GetMap") || MultiPatch_GetMap(cctkGH) != jacobian_identity_map)
                        && strlen(jacobian_group) > 0;
  const bool usejacobian = assume_use_jacobian>=0 ? assume_use_jacobian : usejacobian1;
  if (usejacobian && (strlen(jacobian_derivative_group) == 0))
  {
    CCTK_WARN(CCTK_WARN_ALERT, "GenericFD::jacobian_group and GenericFD::jacobian_derivative_group must both be set to valid group names");
  }
  
  const CCTK_REAL* restrict jacobian_ptrs[9];
  if (usejacobian) GroupDataPointers(cctkGH, jacobian_group,
                                                9, jacobian_ptrs);
  
  const CCTK_REAL* restrict const J11 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_ptrs[0] : 0;
  const CCTK_REAL* restrict const J12 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_ptrs[1] : 0;
  const CCTK_REAL* restrict const J13 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_ptrs[2] : 0;
  const CCTK_REAL* restrict const J21 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_ptrs[3] : 0;
  const CCTK_REAL* restrict const J22 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_ptrs[4] : 0;
  const CCTK_REAL* restrict const J23 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_ptrs[5] : 0;
  const CCTK_REAL* restrict const J31 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_ptrs[6] : 0;
  const CCTK_REAL* restrict const J32 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_ptrs[7] : 0;
  const CCTK_REAL* restrict const J33 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_ptrs[8] : 0;
  
  const CCTK_REAL* restrict jacobian_determinant_ptrs[1] CCTK_ATTRIBUTE_UNUSED;
  if (usejacobian && strlen(jacobian_determinant_group) > 0) GroupDataPointers(cctkGH, jacobian_determinant_group,
                                                1, jacobian_determinant_ptrs);
  
  const CCTK_REAL* restrict const detJ CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_determinant_ptrs[0] : 0;
  
  const CCTK_REAL* restrict jacobian_inverse_ptrs[9] CCTK_ATTRIBUTE_UNUSED;
  if (usejacobian && strlen(jacobian_inverse_group) > 0) GroupDataPointers(cctkGH, jacobian_inverse_group,
                                                9, jacobian_inverse_ptrs);
  
  const CCTK_REAL* restrict const iJ11 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_inverse_ptrs[0] : 0;
  const CCTK_REAL* restrict const iJ12 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_inverse_ptrs[1] : 0;
  const CCTK_REAL* restrict const iJ13 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_inverse_ptrs[2] : 0;
  const CCTK_REAL* restrict const iJ21 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_inverse_ptrs[3] : 0;
  const CCTK_REAL* restrict const iJ22 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_inverse_ptrs[4] : 0;
  const CCTK_REAL* restrict const iJ23 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_inverse_ptrs[5] : 0;
  const CCTK_REAL* restrict const iJ31 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_inverse_ptrs[6] : 0;
  const CCTK_REAL* restrict const iJ32 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_inverse_ptrs[7] : 0;
  const CCTK_REAL* restrict const iJ33 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_inverse_ptrs[8] : 0;
  
  const CCTK_REAL* restrict jacobian_derivative_ptrs[18] CCTK_ATTRIBUTE_UNUSED;
  if (usejacobian) GroupDataPointers(cctkGH, jacobian_derivative_group,
                                      18, jacobian_derivative_ptrs);
  
  const CCTK_REAL* restrict const dJ111 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[0] : 0;
  const CCTK_REAL* restrict const dJ112 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[1] : 0;
  const CCTK_REAL* restrict const dJ113 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[2] : 0;
  const CCTK_REAL* restrict const dJ122 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[3] : 0;
  const CCTK_REAL* restrict const dJ123 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[4] : 0;
  const CCTK_REAL* restrict const dJ133 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[5] : 0;
  const CCTK_REAL* restrict const dJ211 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[6] : 0;
  const CCTK_REAL* restrict const dJ212 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[7] : 0;
  const CCTK_REAL* restrict const dJ213 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[8] : 0;
  const CCTK_REAL* restrict const dJ222 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[9] : 0;
  const CCTK_REAL* restrict const dJ223 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[10] : 0;
  const CCTK_REAL* restrict const dJ233 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[11] : 0;
  const CCTK_REAL* restrict const dJ311 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[12] : 0;
  const CCTK_REAL* restrict const dJ312 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[13] : 0;
  const CCTK_REAL* restrict const dJ313 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[14] : 0;
  const CCTK_REAL* restrict const dJ322 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[15] : 0;
  const CCTK_REAL* restrict const dJ323 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[16] : 0;
  const CCTK_REAL* restrict const dJ333 CCTK_ATTRIBUTE_UNUSED = usejacobian ? jacobian_derivative_ptrs[17] : 0;
  /* Assign local copies of arrays functions */
  
  
  /* Calculate temporaries and arrays functions */
  /* Copy local copies back to grid functions */
  /* Loop over the grid points */
  const int imin0=imin[0];
  const int imin1=imin[1];
  const int imin2=imin[2];
  const int imax0=imax[0];
  const int imax1=imax[1];
  const int imax2=imax[2];
  const int N = cctk_ash[0]*cctk_ash[1]*cctk_ash[2];
  #pragma omp parallel
  CCTK_LOOP3(CT_Dust_RHSVal,
    i,j,k, imin0,imin1,imin2, imax0,imax1,imax2,
    cctk_ash[0],cctk_ash[1],cctk_ash[2])
  {
    const ptrdiff_t index CCTK_ATTRIBUTE_UNUSED = di*i + dj*j + dk*k;
    /* Assign local copies of grid functions */
    
    CCTK_REAL alpL CCTK_ATTRIBUTE_UNUSED = alp[index];
    CCTK_REAL betaxL CCTK_ATTRIBUTE_UNUSED = betax[index];
    CCTK_REAL betayL CCTK_ATTRIBUTE_UNUSED = betay[index];
    CCTK_REAL betazL CCTK_ATTRIBUTE_UNUSED = betaz[index];
    CCTK_REAL EEL CCTK_ATTRIBUTE_UNUSED = EE[index];
    CCTK_REAL epsL CCTK_ATTRIBUTE_UNUSED = eps[index];
    CCTK_REAL gxxL CCTK_ATTRIBUTE_UNUSED = gxx[index];
    CCTK_REAL gxyL CCTK_ATTRIBUTE_UNUSED = gxy[index];
    CCTK_REAL gxzL CCTK_ATTRIBUTE_UNUSED = gxz[index];
    CCTK_REAL gyyL CCTK_ATTRIBUTE_UNUSED = gyy[index];
    CCTK_REAL gyzL CCTK_ATTRIBUTE_UNUSED = gyz[index];
    CCTK_REAL gzzL CCTK_ATTRIBUTE_UNUSED = gzz[index];
    CCTK_REAL kxxL CCTK_ATTRIBUTE_UNUSED = kxx[index];
    CCTK_REAL kxyL CCTK_ATTRIBUTE_UNUSED = kxy[index];
    CCTK_REAL kxzL CCTK_ATTRIBUTE_UNUSED = kxz[index];
    CCTK_REAL kyyL CCTK_ATTRIBUTE_UNUSED = kyy[index];
    CCTK_REAL kyzL CCTK_ATTRIBUTE_UNUSED = kyz[index];
    CCTK_REAL kzzL CCTK_ATTRIBUTE_UNUSED = kzz[index];
    CCTK_REAL pressL CCTK_ATTRIBUTE_UNUSED = press[index];
    CCTK_REAL rhoL CCTK_ATTRIBUTE_UNUSED = rho[index];
    CCTK_REAL SS1L CCTK_ATTRIBUTE_UNUSED = SS1[index];
    CCTK_REAL SS2L CCTK_ATTRIBUTE_UNUSED = SS2[index];
    CCTK_REAL SS3L CCTK_ATTRIBUTE_UNUSED = SS3[index];
    CCTK_REAL velxL CCTK_ATTRIBUTE_UNUSED = vel[index];
    CCTK_REAL velyL CCTK_ATTRIBUTE_UNUSED = vel[index+N];
    CCTK_REAL velzL CCTK_ATTRIBUTE_UNUSED = vel[index+2*N];
    CCTK_REAL w_lorentzL CCTK_ATTRIBUTE_UNUSED = w_lorentz[index];
    
    
    CCTK_REAL J11L, J12L, J13L, J21L, J22L, J23L, J31L, J32L, J33L CCTK_ATTRIBUTE_UNUSED ;
    
    if (usejacobian)
    {
      J11L = J11[index];
      J12L = J12[index];
      J13L = J13[index];
      J21L = J21[index];
      J22L = J22[index];
      J23L = J23[index];
      J31L = J31[index];
      J32L = J32[index];
      J33L = J33[index];
    }
    /* Include user supplied include files */
    /* Precompute derivatives */
    CCTK_REAL PDstandardNth1alp CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2alp CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3alp CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth1betax CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2betax CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3betax CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth1betay CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2betay CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3betay CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth1betaz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2betaz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3betaz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth1EE CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2EE CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3EE CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth1eps CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2eps CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3eps CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth1gxx CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2gxx CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3gxx CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth1gxy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2gxy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3gxy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth1gxz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2gxz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3gxz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth1gyy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2gyy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3gyy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth1gyz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2gyz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3gyz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth1gzz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2gzz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3gzz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth1press CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2press CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3press CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth1rho CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2rho CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3rho CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth1SS1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2SS1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3SS1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth1SS2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2SS2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3SS2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth1SS3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2SS3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3SS3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth1velx CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2velx CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3velx CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth1vely CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2vely CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3vely CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth1velz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2velz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3velz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth1w_lorentz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2w_lorentz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3w_lorentz CCTK_ATTRIBUTE_UNUSED;
    
    switch (fdOrder)
    {
      case 2:
      {
        PDstandardNth1alp = PDstandardNthfdOrder21(&alp[index]);
        PDstandardNth2alp = PDstandardNthfdOrder22(&alp[index]);
        PDstandardNth3alp = PDstandardNthfdOrder23(&alp[index]);
        PDstandardNth1betax = PDstandardNthfdOrder21(&betax[index]);
        PDstandardNth2betax = PDstandardNthfdOrder22(&betax[index]);
        PDstandardNth3betax = PDstandardNthfdOrder23(&betax[index]);
        PDstandardNth1betay = PDstandardNthfdOrder21(&betay[index]);
        PDstandardNth2betay = PDstandardNthfdOrder22(&betay[index]);
        PDstandardNth3betay = PDstandardNthfdOrder23(&betay[index]);
        PDstandardNth1betaz = PDstandardNthfdOrder21(&betaz[index]);
        PDstandardNth2betaz = PDstandardNthfdOrder22(&betaz[index]);
        PDstandardNth3betaz = PDstandardNthfdOrder23(&betaz[index]);
        PDstandardNth1EE = PDstandardNthfdOrder21(&EE[index]);
        PDstandardNth2EE = PDstandardNthfdOrder22(&EE[index]);
        PDstandardNth3EE = PDstandardNthfdOrder23(&EE[index]);
        PDstandardNth1eps = PDstandardNthfdOrder21(&eps[index]);
        PDstandardNth2eps = PDstandardNthfdOrder22(&eps[index]);
        PDstandardNth3eps = PDstandardNthfdOrder23(&eps[index]);
        PDstandardNth1gxx = PDstandardNthfdOrder21(&gxx[index]);
        PDstandardNth2gxx = PDstandardNthfdOrder22(&gxx[index]);
        PDstandardNth3gxx = PDstandardNthfdOrder23(&gxx[index]);
        PDstandardNth1gxy = PDstandardNthfdOrder21(&gxy[index]);
        PDstandardNth2gxy = PDstandardNthfdOrder22(&gxy[index]);
        PDstandardNth3gxy = PDstandardNthfdOrder23(&gxy[index]);
        PDstandardNth1gxz = PDstandardNthfdOrder21(&gxz[index]);
        PDstandardNth2gxz = PDstandardNthfdOrder22(&gxz[index]);
        PDstandardNth3gxz = PDstandardNthfdOrder23(&gxz[index]);
        PDstandardNth1gyy = PDstandardNthfdOrder21(&gyy[index]);
        PDstandardNth2gyy = PDstandardNthfdOrder22(&gyy[index]);
        PDstandardNth3gyy = PDstandardNthfdOrder23(&gyy[index]);
        PDstandardNth1gyz = PDstandardNthfdOrder21(&gyz[index]);
        PDstandardNth2gyz = PDstandardNthfdOrder22(&gyz[index]);
        PDstandardNth3gyz = PDstandardNthfdOrder23(&gyz[index]);
        PDstandardNth1gzz = PDstandardNthfdOrder21(&gzz[index]);
        PDstandardNth2gzz = PDstandardNthfdOrder22(&gzz[index]);
        PDstandardNth3gzz = PDstandardNthfdOrder23(&gzz[index]);
        PDstandardNth1press = PDstandardNthfdOrder21(&press[index]);
        PDstandardNth2press = PDstandardNthfdOrder22(&press[index]);
        PDstandardNth3press = PDstandardNthfdOrder23(&press[index]);
        PDstandardNth1rho = PDstandardNthfdOrder21(&rho[index]);
        PDstandardNth2rho = PDstandardNthfdOrder22(&rho[index]);
        PDstandardNth3rho = PDstandardNthfdOrder23(&rho[index]);
        PDstandardNth1SS1 = PDstandardNthfdOrder21(&SS1[index]);
        PDstandardNth2SS1 = PDstandardNthfdOrder22(&SS1[index]);
        PDstandardNth3SS1 = PDstandardNthfdOrder23(&SS1[index]);
        PDstandardNth1SS2 = PDstandardNthfdOrder21(&SS2[index]);
        PDstandardNth2SS2 = PDstandardNthfdOrder22(&SS2[index]);
        PDstandardNth3SS2 = PDstandardNthfdOrder23(&SS2[index]);
        PDstandardNth1SS3 = PDstandardNthfdOrder21(&SS3[index]);
        PDstandardNth2SS3 = PDstandardNthfdOrder22(&SS3[index]);
        PDstandardNth3SS3 = PDstandardNthfdOrder23(&SS3[index]);
        PDstandardNth1velx = PDstandardNthfdOrder21(&vel[index]);
        PDstandardNth2velx = PDstandardNthfdOrder22(&vel[index]);
        PDstandardNth3velx = PDstandardNthfdOrder23(&vel[index]);
        PDstandardNth1vely = PDstandardNthfdOrder21(&vel[index+N]);
        PDstandardNth2vely = PDstandardNthfdOrder22(&vel[index+N]);
        PDstandardNth3vely = PDstandardNthfdOrder23(&vel[index+N]);
        PDstandardNth1velz = PDstandardNthfdOrder21(&vel[index+2*N]);
        PDstandardNth2velz = PDstandardNthfdOrder22(&vel[index+2*N]);
        PDstandardNth3velz = PDstandardNthfdOrder23(&vel[index+2*N]);
        PDstandardNth1w_lorentz = PDstandardNthfdOrder21(&w_lorentz[index]);
        PDstandardNth2w_lorentz = PDstandardNthfdOrder22(&w_lorentz[index]);
        PDstandardNth3w_lorentz = PDstandardNthfdOrder23(&w_lorentz[index]);
        break;
      }
      
      case 4:
      {
        PDstandardNth1alp = PDstandardNthfdOrder41(&alp[index]);
        PDstandardNth2alp = PDstandardNthfdOrder42(&alp[index]);
        PDstandardNth3alp = PDstandardNthfdOrder43(&alp[index]);
        PDstandardNth1betax = PDstandardNthfdOrder41(&betax[index]);
        PDstandardNth2betax = PDstandardNthfdOrder42(&betax[index]);
        PDstandardNth3betax = PDstandardNthfdOrder43(&betax[index]);
        PDstandardNth1betay = PDstandardNthfdOrder41(&betay[index]);
        PDstandardNth2betay = PDstandardNthfdOrder42(&betay[index]);
        PDstandardNth3betay = PDstandardNthfdOrder43(&betay[index]);
        PDstandardNth1betaz = PDstandardNthfdOrder41(&betaz[index]);
        PDstandardNth2betaz = PDstandardNthfdOrder42(&betaz[index]);
        PDstandardNth3betaz = PDstandardNthfdOrder43(&betaz[index]);
        PDstandardNth1EE = PDstandardNthfdOrder41(&EE[index]);
        PDstandardNth2EE = PDstandardNthfdOrder42(&EE[index]);
        PDstandardNth3EE = PDstandardNthfdOrder43(&EE[index]);
        PDstandardNth1eps = PDstandardNthfdOrder41(&eps[index]);
        PDstandardNth2eps = PDstandardNthfdOrder42(&eps[index]);
        PDstandardNth3eps = PDstandardNthfdOrder43(&eps[index]);
        PDstandardNth1gxx = PDstandardNthfdOrder41(&gxx[index]);
        PDstandardNth2gxx = PDstandardNthfdOrder42(&gxx[index]);
        PDstandardNth3gxx = PDstandardNthfdOrder43(&gxx[index]);
        PDstandardNth1gxy = PDstandardNthfdOrder41(&gxy[index]);
        PDstandardNth2gxy = PDstandardNthfdOrder42(&gxy[index]);
        PDstandardNth3gxy = PDstandardNthfdOrder43(&gxy[index]);
        PDstandardNth1gxz = PDstandardNthfdOrder41(&gxz[index]);
        PDstandardNth2gxz = PDstandardNthfdOrder42(&gxz[index]);
        PDstandardNth3gxz = PDstandardNthfdOrder43(&gxz[index]);
        PDstandardNth1gyy = PDstandardNthfdOrder41(&gyy[index]);
        PDstandardNth2gyy = PDstandardNthfdOrder42(&gyy[index]);
        PDstandardNth3gyy = PDstandardNthfdOrder43(&gyy[index]);
        PDstandardNth1gyz = PDstandardNthfdOrder41(&gyz[index]);
        PDstandardNth2gyz = PDstandardNthfdOrder42(&gyz[index]);
        PDstandardNth3gyz = PDstandardNthfdOrder43(&gyz[index]);
        PDstandardNth1gzz = PDstandardNthfdOrder41(&gzz[index]);
        PDstandardNth2gzz = PDstandardNthfdOrder42(&gzz[index]);
        PDstandardNth3gzz = PDstandardNthfdOrder43(&gzz[index]);
        PDstandardNth1press = PDstandardNthfdOrder41(&press[index]);
        PDstandardNth2press = PDstandardNthfdOrder42(&press[index]);
        PDstandardNth3press = PDstandardNthfdOrder43(&press[index]);
        PDstandardNth1rho = PDstandardNthfdOrder41(&rho[index]);
        PDstandardNth2rho = PDstandardNthfdOrder42(&rho[index]);
        PDstandardNth3rho = PDstandardNthfdOrder43(&rho[index]);
        PDstandardNth1SS1 = PDstandardNthfdOrder41(&SS1[index]);
        PDstandardNth2SS1 = PDstandardNthfdOrder42(&SS1[index]);
        PDstandardNth3SS1 = PDstandardNthfdOrder43(&SS1[index]);
        PDstandardNth1SS2 = PDstandardNthfdOrder41(&SS2[index]);
        PDstandardNth2SS2 = PDstandardNthfdOrder42(&SS2[index]);
        PDstandardNth3SS2 = PDstandardNthfdOrder43(&SS2[index]);
        PDstandardNth1SS3 = PDstandardNthfdOrder41(&SS3[index]);
        PDstandardNth2SS3 = PDstandardNthfdOrder42(&SS3[index]);
        PDstandardNth3SS3 = PDstandardNthfdOrder43(&SS3[index]);
        PDstandardNth1velx = PDstandardNthfdOrder41(&vel[index]);
        PDstandardNth2velx = PDstandardNthfdOrder42(&vel[index]);
        PDstandardNth3velx = PDstandardNthfdOrder43(&vel[index]);
        PDstandardNth1vely = PDstandardNthfdOrder41(&vel[index+N]);
        PDstandardNth2vely = PDstandardNthfdOrder42(&vel[index+N]);
        PDstandardNth3vely = PDstandardNthfdOrder43(&vel[index+N]);
        PDstandardNth1velz = PDstandardNthfdOrder41(&vel[index+2*N]);
        PDstandardNth2velz = PDstandardNthfdOrder42(&vel[index+2*N]);
        PDstandardNth3velz = PDstandardNthfdOrder43(&vel[index+2*N]);
        PDstandardNth1w_lorentz = PDstandardNthfdOrder41(&w_lorentz[index]);
        PDstandardNth2w_lorentz = PDstandardNthfdOrder42(&w_lorentz[index]);
        PDstandardNth3w_lorentz = PDstandardNthfdOrder43(&w_lorentz[index]);
        break;
      }
      
      case 6:
      {
        PDstandardNth1alp = PDstandardNthfdOrder61(&alp[index]);
        PDstandardNth2alp = PDstandardNthfdOrder62(&alp[index]);
        PDstandardNth3alp = PDstandardNthfdOrder63(&alp[index]);
        PDstandardNth1betax = PDstandardNthfdOrder61(&betax[index]);
        PDstandardNth2betax = PDstandardNthfdOrder62(&betax[index]);
        PDstandardNth3betax = PDstandardNthfdOrder63(&betax[index]);
        PDstandardNth1betay = PDstandardNthfdOrder61(&betay[index]);
        PDstandardNth2betay = PDstandardNthfdOrder62(&betay[index]);
        PDstandardNth3betay = PDstandardNthfdOrder63(&betay[index]);
        PDstandardNth1betaz = PDstandardNthfdOrder61(&betaz[index]);
        PDstandardNth2betaz = PDstandardNthfdOrder62(&betaz[index]);
        PDstandardNth3betaz = PDstandardNthfdOrder63(&betaz[index]);
        PDstandardNth1EE = PDstandardNthfdOrder61(&EE[index]);
        PDstandardNth2EE = PDstandardNthfdOrder62(&EE[index]);
        PDstandardNth3EE = PDstandardNthfdOrder63(&EE[index]);
        PDstandardNth1eps = PDstandardNthfdOrder61(&eps[index]);
        PDstandardNth2eps = PDstandardNthfdOrder62(&eps[index]);
        PDstandardNth3eps = PDstandardNthfdOrder63(&eps[index]);
        PDstandardNth1gxx = PDstandardNthfdOrder61(&gxx[index]);
        PDstandardNth2gxx = PDstandardNthfdOrder62(&gxx[index]);
        PDstandardNth3gxx = PDstandardNthfdOrder63(&gxx[index]);
        PDstandardNth1gxy = PDstandardNthfdOrder61(&gxy[index]);
        PDstandardNth2gxy = PDstandardNthfdOrder62(&gxy[index]);
        PDstandardNth3gxy = PDstandardNthfdOrder63(&gxy[index]);
        PDstandardNth1gxz = PDstandardNthfdOrder61(&gxz[index]);
        PDstandardNth2gxz = PDstandardNthfdOrder62(&gxz[index]);
        PDstandardNth3gxz = PDstandardNthfdOrder63(&gxz[index]);
        PDstandardNth1gyy = PDstandardNthfdOrder61(&gyy[index]);
        PDstandardNth2gyy = PDstandardNthfdOrder62(&gyy[index]);
        PDstandardNth3gyy = PDstandardNthfdOrder63(&gyy[index]);
        PDstandardNth1gyz = PDstandardNthfdOrder61(&gyz[index]);
        PDstandardNth2gyz = PDstandardNthfdOrder62(&gyz[index]);
        PDstandardNth3gyz = PDstandardNthfdOrder63(&gyz[index]);
        PDstandardNth1gzz = PDstandardNthfdOrder61(&gzz[index]);
        PDstandardNth2gzz = PDstandardNthfdOrder62(&gzz[index]);
        PDstandardNth3gzz = PDstandardNthfdOrder63(&gzz[index]);
        PDstandardNth1press = PDstandardNthfdOrder61(&press[index]);
        PDstandardNth2press = PDstandardNthfdOrder62(&press[index]);
        PDstandardNth3press = PDstandardNthfdOrder63(&press[index]);
        PDstandardNth1rho = PDstandardNthfdOrder61(&rho[index]);
        PDstandardNth2rho = PDstandardNthfdOrder62(&rho[index]);
        PDstandardNth3rho = PDstandardNthfdOrder63(&rho[index]);
        PDstandardNth1SS1 = PDstandardNthfdOrder61(&SS1[index]);
        PDstandardNth2SS1 = PDstandardNthfdOrder62(&SS1[index]);
        PDstandardNth3SS1 = PDstandardNthfdOrder63(&SS1[index]);
        PDstandardNth1SS2 = PDstandardNthfdOrder61(&SS2[index]);
        PDstandardNth2SS2 = PDstandardNthfdOrder62(&SS2[index]);
        PDstandardNth3SS2 = PDstandardNthfdOrder63(&SS2[index]);
        PDstandardNth1SS3 = PDstandardNthfdOrder61(&SS3[index]);
        PDstandardNth2SS3 = PDstandardNthfdOrder62(&SS3[index]);
        PDstandardNth3SS3 = PDstandardNthfdOrder63(&SS3[index]);
        PDstandardNth1velx = PDstandardNthfdOrder61(&vel[index]);
        PDstandardNth2velx = PDstandardNthfdOrder62(&vel[index]);
        PDstandardNth3velx = PDstandardNthfdOrder63(&vel[index]);
        PDstandardNth1vely = PDstandardNthfdOrder61(&vel[index+N]);
        PDstandardNth2vely = PDstandardNthfdOrder62(&vel[index+N]);
        PDstandardNth3vely = PDstandardNthfdOrder63(&vel[index+N]);
        PDstandardNth1velz = PDstandardNthfdOrder61(&vel[index+2*N]);
        PDstandardNth2velz = PDstandardNthfdOrder62(&vel[index+2*N]);
        PDstandardNth3velz = PDstandardNthfdOrder63(&vel[index+2*N]);
        PDstandardNth1w_lorentz = PDstandardNthfdOrder61(&w_lorentz[index]);
        PDstandardNth2w_lorentz = PDstandardNthfdOrder62(&w_lorentz[index]);
        PDstandardNth3w_lorentz = PDstandardNthfdOrder63(&w_lorentz[index]);
        break;
      }
      
      case 8:
      {
        PDstandardNth1alp = PDstandardNthfdOrder81(&alp[index]);
        PDstandardNth2alp = PDstandardNthfdOrder82(&alp[index]);
        PDstandardNth3alp = PDstandardNthfdOrder83(&alp[index]);
        PDstandardNth1betax = PDstandardNthfdOrder81(&betax[index]);
        PDstandardNth2betax = PDstandardNthfdOrder82(&betax[index]);
        PDstandardNth3betax = PDstandardNthfdOrder83(&betax[index]);
        PDstandardNth1betay = PDstandardNthfdOrder81(&betay[index]);
        PDstandardNth2betay = PDstandardNthfdOrder82(&betay[index]);
        PDstandardNth3betay = PDstandardNthfdOrder83(&betay[index]);
        PDstandardNth1betaz = PDstandardNthfdOrder81(&betaz[index]);
        PDstandardNth2betaz = PDstandardNthfdOrder82(&betaz[index]);
        PDstandardNth3betaz = PDstandardNthfdOrder83(&betaz[index]);
        PDstandardNth1EE = PDstandardNthfdOrder81(&EE[index]);
        PDstandardNth2EE = PDstandardNthfdOrder82(&EE[index]);
        PDstandardNth3EE = PDstandardNthfdOrder83(&EE[index]);
        PDstandardNth1eps = PDstandardNthfdOrder81(&eps[index]);
        PDstandardNth2eps = PDstandardNthfdOrder82(&eps[index]);
        PDstandardNth3eps = PDstandardNthfdOrder83(&eps[index]);
        PDstandardNth1gxx = PDstandardNthfdOrder81(&gxx[index]);
        PDstandardNth2gxx = PDstandardNthfdOrder82(&gxx[index]);
        PDstandardNth3gxx = PDstandardNthfdOrder83(&gxx[index]);
        PDstandardNth1gxy = PDstandardNthfdOrder81(&gxy[index]);
        PDstandardNth2gxy = PDstandardNthfdOrder82(&gxy[index]);
        PDstandardNth3gxy = PDstandardNthfdOrder83(&gxy[index]);
        PDstandardNth1gxz = PDstandardNthfdOrder81(&gxz[index]);
        PDstandardNth2gxz = PDstandardNthfdOrder82(&gxz[index]);
        PDstandardNth3gxz = PDstandardNthfdOrder83(&gxz[index]);
        PDstandardNth1gyy = PDstandardNthfdOrder81(&gyy[index]);
        PDstandardNth2gyy = PDstandardNthfdOrder82(&gyy[index]);
        PDstandardNth3gyy = PDstandardNthfdOrder83(&gyy[index]);
        PDstandardNth1gyz = PDstandardNthfdOrder81(&gyz[index]);
        PDstandardNth2gyz = PDstandardNthfdOrder82(&gyz[index]);
        PDstandardNth3gyz = PDstandardNthfdOrder83(&gyz[index]);
        PDstandardNth1gzz = PDstandardNthfdOrder81(&gzz[index]);
        PDstandardNth2gzz = PDstandardNthfdOrder82(&gzz[index]);
        PDstandardNth3gzz = PDstandardNthfdOrder83(&gzz[index]);
        PDstandardNth1press = PDstandardNthfdOrder81(&press[index]);
        PDstandardNth2press = PDstandardNthfdOrder82(&press[index]);
        PDstandardNth3press = PDstandardNthfdOrder83(&press[index]);
        PDstandardNth1rho = PDstandardNthfdOrder81(&rho[index]);
        PDstandardNth2rho = PDstandardNthfdOrder82(&rho[index]);
        PDstandardNth3rho = PDstandardNthfdOrder83(&rho[index]);
        PDstandardNth1SS1 = PDstandardNthfdOrder81(&SS1[index]);
        PDstandardNth2SS1 = PDstandardNthfdOrder82(&SS1[index]);
        PDstandardNth3SS1 = PDstandardNthfdOrder83(&SS1[index]);
        PDstandardNth1SS2 = PDstandardNthfdOrder81(&SS2[index]);
        PDstandardNth2SS2 = PDstandardNthfdOrder82(&SS2[index]);
        PDstandardNth3SS2 = PDstandardNthfdOrder83(&SS2[index]);
        PDstandardNth1SS3 = PDstandardNthfdOrder81(&SS3[index]);
        PDstandardNth2SS3 = PDstandardNthfdOrder82(&SS3[index]);
        PDstandardNth3SS3 = PDstandardNthfdOrder83(&SS3[index]);
        PDstandardNth1velx = PDstandardNthfdOrder81(&vel[index]);
        PDstandardNth2velx = PDstandardNthfdOrder82(&vel[index]);
        PDstandardNth3velx = PDstandardNthfdOrder83(&vel[index]);
        PDstandardNth1vely = PDstandardNthfdOrder81(&vel[index+N]);
        PDstandardNth2vely = PDstandardNthfdOrder82(&vel[index+N]);
        PDstandardNth3vely = PDstandardNthfdOrder83(&vel[index+N]);
        PDstandardNth1velz = PDstandardNthfdOrder81(&vel[index+2*N]);
        PDstandardNth2velz = PDstandardNthfdOrder82(&vel[index+2*N]);
        PDstandardNth3velz = PDstandardNthfdOrder83(&vel[index+2*N]);
        PDstandardNth1w_lorentz = PDstandardNthfdOrder81(&w_lorentz[index]);
        PDstandardNth2w_lorentz = PDstandardNthfdOrder82(&w_lorentz[index]);
        PDstandardNth3w_lorentz = PDstandardNthfdOrder83(&w_lorentz[index]);
        break;
      }
      default:
        CCTK_BUILTIN_UNREACHABLE();
    }
    /* Calculate temporaries and grid functions */
    ptrdiff_t dir1 CCTK_ATTRIBUTE_UNUSED = isgn(betaxL);
    
    ptrdiff_t dir2 CCTK_ATTRIBUTE_UNUSED = isgn(betayL);
    
    ptrdiff_t dir3 CCTK_ATTRIBUTE_UNUSED = isgn(betazL);
    
    CCTK_REAL JacPDstandardNth1alp CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1betax CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1betay CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1betaz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1EE CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1eps CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1gxx CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1gxy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1gxz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1gyy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1gyz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1gzz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1press CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1rho CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1SS1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1SS2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1SS3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1velx CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1vely CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1velz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1w_lorentz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2alp CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2betax CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2betay CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2betaz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2EE CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2eps CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2gxx CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2gxy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2gxz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2gyy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2gyz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2gzz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2press CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2rho CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2SS1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2SS2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2SS3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2velx CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2vely CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2velz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2w_lorentz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3alp CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3betax CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3betay CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3betaz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3EE CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3eps CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3gxx CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3gxy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3gxz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3gyy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3gyz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3gzz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3press CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3rho CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3SS1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3SS2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3SS3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3velx CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3vely CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3velz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3w_lorentz CCTK_ATTRIBUTE_UNUSED;
    
    if (usejacobian)
    {
      JacPDstandardNth1alp = J11L*PDstandardNth1alp + J21L*PDstandardNth2alp 
        + J31L*PDstandardNth3alp;
      
      JacPDstandardNth1betax = J11L*PDstandardNth1betax + 
        J21L*PDstandardNth2betax + J31L*PDstandardNth3betax;
      
      JacPDstandardNth1betay = J11L*PDstandardNth1betay + 
        J21L*PDstandardNth2betay + J31L*PDstandardNth3betay;
      
      JacPDstandardNth1betaz = J11L*PDstandardNth1betaz + 
        J21L*PDstandardNth2betaz + J31L*PDstandardNth3betaz;
      
      JacPDstandardNth1EE = J11L*PDstandardNth1EE + J21L*PDstandardNth2EE + 
        J31L*PDstandardNth3EE;
      
      JacPDstandardNth1eps = J11L*PDstandardNth1eps + J21L*PDstandardNth2eps 
        + J31L*PDstandardNth3eps;
      
      JacPDstandardNth1gxx = J11L*PDstandardNth1gxx + J21L*PDstandardNth2gxx 
        + J31L*PDstandardNth3gxx;
      
      JacPDstandardNth1gxy = J11L*PDstandardNth1gxy + J21L*PDstandardNth2gxy 
        + J31L*PDstandardNth3gxy;
      
      JacPDstandardNth1gxz = J11L*PDstandardNth1gxz + J21L*PDstandardNth2gxz 
        + J31L*PDstandardNth3gxz;
      
      JacPDstandardNth1gyy = J11L*PDstandardNth1gyy + J21L*PDstandardNth2gyy 
        + J31L*PDstandardNth3gyy;
      
      JacPDstandardNth1gyz = J11L*PDstandardNth1gyz + J21L*PDstandardNth2gyz 
        + J31L*PDstandardNth3gyz;
      
      JacPDstandardNth1gzz = J11L*PDstandardNth1gzz + J21L*PDstandardNth2gzz 
        + J31L*PDstandardNth3gzz;
      
      JacPDstandardNth1press = J11L*PDstandardNth1press + 
        J21L*PDstandardNth2press + J31L*PDstandardNth3press;
      
      JacPDstandardNth1rho = J11L*PDstandardNth1rho + J21L*PDstandardNth2rho 
        + J31L*PDstandardNth3rho;
      
      JacPDstandardNth1SS1 = J11L*PDstandardNth1SS1 + J21L*PDstandardNth2SS1 
        + J31L*PDstandardNth3SS1;
      
      JacPDstandardNth1SS2 = J11L*PDstandardNth1SS2 + J21L*PDstandardNth2SS2 
        + J31L*PDstandardNth3SS2;
      
      JacPDstandardNth1SS3 = J11L*PDstandardNth1SS3 + J21L*PDstandardNth2SS3 
        + J31L*PDstandardNth3SS3;
      
      JacPDstandardNth1velx = J11L*PDstandardNth1velx + 
        J21L*PDstandardNth2velx + J31L*PDstandardNth3velx;
      
      JacPDstandardNth1vely = J11L*PDstandardNth1vely + 
        J21L*PDstandardNth2vely + J31L*PDstandardNth3vely;
      
      JacPDstandardNth1velz = J11L*PDstandardNth1velz + 
        J21L*PDstandardNth2velz + J31L*PDstandardNth3velz;
      
      JacPDstandardNth1w_lorentz = J11L*PDstandardNth1w_lorentz + 
        J21L*PDstandardNth2w_lorentz + J31L*PDstandardNth3w_lorentz;
      
      JacPDstandardNth2alp = J12L*PDstandardNth1alp + J22L*PDstandardNth2alp 
        + J32L*PDstandardNth3alp;
      
      JacPDstandardNth2betax = J12L*PDstandardNth1betax + 
        J22L*PDstandardNth2betax + J32L*PDstandardNth3betax;
      
      JacPDstandardNth2betay = J12L*PDstandardNth1betay + 
        J22L*PDstandardNth2betay + J32L*PDstandardNth3betay;
      
      JacPDstandardNth2betaz = J12L*PDstandardNth1betaz + 
        J22L*PDstandardNth2betaz + J32L*PDstandardNth3betaz;
      
      JacPDstandardNth2EE = J12L*PDstandardNth1EE + J22L*PDstandardNth2EE + 
        J32L*PDstandardNth3EE;
      
      JacPDstandardNth2eps = J12L*PDstandardNth1eps + J22L*PDstandardNth2eps 
        + J32L*PDstandardNth3eps;
      
      JacPDstandardNth2gxx = J12L*PDstandardNth1gxx + J22L*PDstandardNth2gxx 
        + J32L*PDstandardNth3gxx;
      
      JacPDstandardNth2gxy = J12L*PDstandardNth1gxy + J22L*PDstandardNth2gxy 
        + J32L*PDstandardNth3gxy;
      
      JacPDstandardNth2gxz = J12L*PDstandardNth1gxz + J22L*PDstandardNth2gxz 
        + J32L*PDstandardNth3gxz;
      
      JacPDstandardNth2gyy = J12L*PDstandardNth1gyy + J22L*PDstandardNth2gyy 
        + J32L*PDstandardNth3gyy;
      
      JacPDstandardNth2gyz = J12L*PDstandardNth1gyz + J22L*PDstandardNth2gyz 
        + J32L*PDstandardNth3gyz;
      
      JacPDstandardNth2gzz = J12L*PDstandardNth1gzz + J22L*PDstandardNth2gzz 
        + J32L*PDstandardNth3gzz;
      
      JacPDstandardNth2press = J12L*PDstandardNth1press + 
        J22L*PDstandardNth2press + J32L*PDstandardNth3press;
      
      JacPDstandardNth2rho = J12L*PDstandardNth1rho + J22L*PDstandardNth2rho 
        + J32L*PDstandardNth3rho;
      
      JacPDstandardNth2SS1 = J12L*PDstandardNth1SS1 + J22L*PDstandardNth2SS1 
        + J32L*PDstandardNth3SS1;
      
      JacPDstandardNth2SS2 = J12L*PDstandardNth1SS2 + J22L*PDstandardNth2SS2 
        + J32L*PDstandardNth3SS2;
      
      JacPDstandardNth2SS3 = J12L*PDstandardNth1SS3 + J22L*PDstandardNth2SS3 
        + J32L*PDstandardNth3SS3;
      
      JacPDstandardNth2velx = J12L*PDstandardNth1velx + 
        J22L*PDstandardNth2velx + J32L*PDstandardNth3velx;
      
      JacPDstandardNth2vely = J12L*PDstandardNth1vely + 
        J22L*PDstandardNth2vely + J32L*PDstandardNth3vely;
      
      JacPDstandardNth2velz = J12L*PDstandardNth1velz + 
        J22L*PDstandardNth2velz + J32L*PDstandardNth3velz;
      
      JacPDstandardNth2w_lorentz = J12L*PDstandardNth1w_lorentz + 
        J22L*PDstandardNth2w_lorentz + J32L*PDstandardNth3w_lorentz;
      
      JacPDstandardNth3alp = J13L*PDstandardNth1alp + J23L*PDstandardNth2alp 
        + J33L*PDstandardNth3alp;
      
      JacPDstandardNth3betax = J13L*PDstandardNth1betax + 
        J23L*PDstandardNth2betax + J33L*PDstandardNth3betax;
      
      JacPDstandardNth3betay = J13L*PDstandardNth1betay + 
        J23L*PDstandardNth2betay + J33L*PDstandardNth3betay;
      
      JacPDstandardNth3betaz = J13L*PDstandardNth1betaz + 
        J23L*PDstandardNth2betaz + J33L*PDstandardNth3betaz;
      
      JacPDstandardNth3EE = J13L*PDstandardNth1EE + J23L*PDstandardNth2EE + 
        J33L*PDstandardNth3EE;
      
      JacPDstandardNth3eps = J13L*PDstandardNth1eps + J23L*PDstandardNth2eps 
        + J33L*PDstandardNth3eps;
      
      JacPDstandardNth3gxx = J13L*PDstandardNth1gxx + J23L*PDstandardNth2gxx 
        + J33L*PDstandardNth3gxx;
      
      JacPDstandardNth3gxy = J13L*PDstandardNth1gxy + J23L*PDstandardNth2gxy 
        + J33L*PDstandardNth3gxy;
      
      JacPDstandardNth3gxz = J13L*PDstandardNth1gxz + J23L*PDstandardNth2gxz 
        + J33L*PDstandardNth3gxz;
      
      JacPDstandardNth3gyy = J13L*PDstandardNth1gyy + J23L*PDstandardNth2gyy 
        + J33L*PDstandardNth3gyy;
      
      JacPDstandardNth3gyz = J13L*PDstandardNth1gyz + J23L*PDstandardNth2gyz 
        + J33L*PDstandardNth3gyz;
      
      JacPDstandardNth3gzz = J13L*PDstandardNth1gzz + J23L*PDstandardNth2gzz 
        + J33L*PDstandardNth3gzz;
      
      JacPDstandardNth3press = J13L*PDstandardNth1press + 
        J23L*PDstandardNth2press + J33L*PDstandardNth3press;
      
      JacPDstandardNth3rho = J13L*PDstandardNth1rho + J23L*PDstandardNth2rho 
        + J33L*PDstandardNth3rho;
      
      JacPDstandardNth3SS1 = J13L*PDstandardNth1SS1 + J23L*PDstandardNth2SS1 
        + J33L*PDstandardNth3SS1;
      
      JacPDstandardNth3SS2 = J13L*PDstandardNth1SS2 + J23L*PDstandardNth2SS2 
        + J33L*PDstandardNth3SS2;
      
      JacPDstandardNth3SS3 = J13L*PDstandardNth1SS3 + J23L*PDstandardNth2SS3 
        + J33L*PDstandardNth3SS3;
      
      JacPDstandardNth3velx = J13L*PDstandardNth1velx + 
        J23L*PDstandardNth2velx + J33L*PDstandardNth3velx;
      
      JacPDstandardNth3vely = J13L*PDstandardNth1vely + 
        J23L*PDstandardNth2vely + J33L*PDstandardNth3vely;
      
      JacPDstandardNth3velz = J13L*PDstandardNth1velz + 
        J23L*PDstandardNth2velz + J33L*PDstandardNth3velz;
      
      JacPDstandardNth3w_lorentz = J13L*PDstandardNth1w_lorentz + 
        J23L*PDstandardNth2w_lorentz + J33L*PDstandardNth3w_lorentz;
    }
    else
    {
      JacPDstandardNth1alp = PDstandardNth1alp;
      
      JacPDstandardNth1betax = PDstandardNth1betax;
      
      JacPDstandardNth1betay = PDstandardNth1betay;
      
      JacPDstandardNth1betaz = PDstandardNth1betaz;
      
      JacPDstandardNth1EE = PDstandardNth1EE;
      
      JacPDstandardNth1eps = PDstandardNth1eps;
      
      JacPDstandardNth1gxx = PDstandardNth1gxx;
      
      JacPDstandardNth1gxy = PDstandardNth1gxy;
      
      JacPDstandardNth1gxz = PDstandardNth1gxz;
      
      JacPDstandardNth1gyy = PDstandardNth1gyy;
      
      JacPDstandardNth1gyz = PDstandardNth1gyz;
      
      JacPDstandardNth1gzz = PDstandardNth1gzz;
      
      JacPDstandardNth1press = PDstandardNth1press;
      
      JacPDstandardNth1rho = PDstandardNth1rho;
      
      JacPDstandardNth1SS1 = PDstandardNth1SS1;
      
      JacPDstandardNth1SS2 = PDstandardNth1SS2;
      
      JacPDstandardNth1SS3 = PDstandardNth1SS3;
      
      JacPDstandardNth1velx = PDstandardNth1velx;
      
      JacPDstandardNth1vely = PDstandardNth1vely;
      
      JacPDstandardNth1velz = PDstandardNth1velz;
      
      JacPDstandardNth1w_lorentz = PDstandardNth1w_lorentz;
      
      JacPDstandardNth2alp = PDstandardNth2alp;
      
      JacPDstandardNth2betax = PDstandardNth2betax;
      
      JacPDstandardNth2betay = PDstandardNth2betay;
      
      JacPDstandardNth2betaz = PDstandardNth2betaz;
      
      JacPDstandardNth2EE = PDstandardNth2EE;
      
      JacPDstandardNth2eps = PDstandardNth2eps;
      
      JacPDstandardNth2gxx = PDstandardNth2gxx;
      
      JacPDstandardNth2gxy = PDstandardNth2gxy;
      
      JacPDstandardNth2gxz = PDstandardNth2gxz;
      
      JacPDstandardNth2gyy = PDstandardNth2gyy;
      
      JacPDstandardNth2gyz = PDstandardNth2gyz;
      
      JacPDstandardNth2gzz = PDstandardNth2gzz;
      
      JacPDstandardNth2press = PDstandardNth2press;
      
      JacPDstandardNth2rho = PDstandardNth2rho;
      
      JacPDstandardNth2SS1 = PDstandardNth2SS1;
      
      JacPDstandardNth2SS2 = PDstandardNth2SS2;
      
      JacPDstandardNth2SS3 = PDstandardNth2SS3;
      
      JacPDstandardNth2velx = PDstandardNth2velx;
      
      JacPDstandardNth2vely = PDstandardNth2vely;
      
      JacPDstandardNth2velz = PDstandardNth2velz;
      
      JacPDstandardNth2w_lorentz = PDstandardNth2w_lorentz;
      
      JacPDstandardNth3alp = PDstandardNth3alp;
      
      JacPDstandardNth3betax = PDstandardNth3betax;
      
      JacPDstandardNth3betay = PDstandardNth3betay;
      
      JacPDstandardNth3betaz = PDstandardNth3betaz;
      
      JacPDstandardNth3EE = PDstandardNth3EE;
      
      JacPDstandardNth3eps = PDstandardNth3eps;
      
      JacPDstandardNth3gxx = PDstandardNth3gxx;
      
      JacPDstandardNth3gxy = PDstandardNth3gxy;
      
      JacPDstandardNth3gxz = PDstandardNth3gxz;
      
      JacPDstandardNth3gyy = PDstandardNth3gyy;
      
      JacPDstandardNth3gyz = PDstandardNth3gyz;
      
      JacPDstandardNth3gzz = PDstandardNth3gzz;
      
      JacPDstandardNth3press = PDstandardNth3press;
      
      JacPDstandardNth3rho = PDstandardNth3rho;
      
      JacPDstandardNth3SS1 = PDstandardNth3SS1;
      
      JacPDstandardNth3SS2 = PDstandardNth3SS2;
      
      JacPDstandardNth3SS3 = PDstandardNth3SS3;
      
      JacPDstandardNth3velx = PDstandardNth3velx;
      
      JacPDstandardNth3vely = PDstandardNth3vely;
      
      JacPDstandardNth3velz = PDstandardNth3velz;
      
      JacPDstandardNth3w_lorentz = PDstandardNth3w_lorentz;
    }
    
    CCTK_REAL detg CCTK_ATTRIBUTE_UNUSED = 2*gxyL*gxzL*gyzL - 
      gzzL*pow(gxyL,2) + gyyL*(gxxL*gzzL - pow(gxzL,2)) - gxxL*pow(gyzL,2);
    
    CCTK_REAL gu11 CCTK_ATTRIBUTE_UNUSED = (gyyL*gzzL - 
      pow(gyzL,2))*pow(detg,-1);
    
    CCTK_REAL gu12 CCTK_ATTRIBUTE_UNUSED = (gxzL*gyzL - 
      gxyL*gzzL)*pow(detg,-1);
    
    CCTK_REAL gu13 CCTK_ATTRIBUTE_UNUSED = (-(gxzL*gyyL) + 
      gxyL*gyzL)*pow(detg,-1);
    
    CCTK_REAL gu22 CCTK_ATTRIBUTE_UNUSED = (gxxL*gzzL - 
      pow(gxzL,2))*pow(detg,-1);
    
    CCTK_REAL gu23 CCTK_ATTRIBUTE_UNUSED = (gxyL*gxzL - 
      gxxL*gyzL)*pow(detg,-1);
    
    CCTK_REAL gu33 CCTK_ATTRIBUTE_UNUSED = (gxxL*gyyL - 
      pow(gxyL,2))*pow(detg,-1);
    
    CCTK_REAL mu CCTK_ATTRIBUTE_UNUSED = (1 + epsL)*rhoL;
    
    CCTK_REAL Sstress11 CCTK_ATTRIBUTE_UNUSED = pressL*gu11 + (pressL + 
      mu)*pow(velxL,2)*pow(w_lorentzL,2);
    
    CCTK_REAL Sstress12 CCTK_ATTRIBUTE_UNUSED = pressL*gu12 + 
      velxL*velyL*(pressL + mu)*pow(w_lorentzL,2);
    
    CCTK_REAL Sstress13 CCTK_ATTRIBUTE_UNUSED = pressL*gu13 + 
      velxL*velzL*(pressL + mu)*pow(w_lorentzL,2);
    
    CCTK_REAL Sstress21 CCTK_ATTRIBUTE_UNUSED = pressL*gu12 + 
      velxL*velyL*(pressL + mu)*pow(w_lorentzL,2);
    
    CCTK_REAL Sstress22 CCTK_ATTRIBUTE_UNUSED = pressL*gu22 + (pressL + 
      mu)*pow(velyL,2)*pow(w_lorentzL,2);
    
    CCTK_REAL Sstress23 CCTK_ATTRIBUTE_UNUSED = pressL*gu23 + 
      velyL*velzL*(pressL + mu)*pow(w_lorentzL,2);
    
    CCTK_REAL Sstress31 CCTK_ATTRIBUTE_UNUSED = pressL*gu13 + 
      velxL*velzL*(pressL + mu)*pow(w_lorentzL,2);
    
    CCTK_REAL Sstress32 CCTK_ATTRIBUTE_UNUSED = pressL*gu23 + 
      velyL*velzL*(pressL + mu)*pow(w_lorentzL,2);
    
    CCTK_REAL Sstress33 CCTK_ATTRIBUTE_UNUSED = pressL*gu33 + (pressL + 
      mu)*pow(velzL,2)*pow(w_lorentzL,2);
    
    CCTK_REAL DDrhsL CCTK_ATTRIBUTE_UNUSED = 0.;
    
    CCTK_REAL SS1rhsL CCTK_ATTRIBUTE_UNUSED = -(EEL*JacPDstandardNth1alp) 
      + SS1L*JacPDstandardNth1betax + SS2L*JacPDstandardNth1betay + 
      SS3L*JacPDstandardNth1betaz + betaxL*JacPDstandardNth1SS1 + 
      betayL*JacPDstandardNth2SS1 + SS1L*(JacPDstandardNth1betax + 
      JacPDstandardNth2betay + JacPDstandardNth3betaz) + 
      betazL*JacPDstandardNth3SS1 + (-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + pow(gyzL,2)))*((1 + 
      epsL)*rhoL*velxL*(gxxL*velxL + gxyL*velyL + 
      gxzL*velzL)*JacPDstandardNth1alp*pow(w_lorentzL,2) + (pressL + rhoL + 
      epsL*rhoL)*(gxxL*velxL + gxyL*velyL + 
      gxzL*velzL)*(velyL*JacPDstandardNth2alp + 
      velzL*JacPDstandardNth3alp)*pow(w_lorentzL,2) + 
      pressL*JacPDstandardNth1alp*(1 + gxyL*velxL*velyL*pow(w_lorentzL,2) + 
      gxzL*velxL*velzL*pow(w_lorentzL,2) + 
      gxxL*pow(velxL,2)*pow(w_lorentzL,2))) + alpL*((pressL + rhoL + 
      epsL*rhoL)*velyL*(gxxL*velxL + gxyL*velyL + 
      gxzL*velzL)*(gzzL*(2*gxyL*JacPDstandardNth2gxy - 
      gxxL*JacPDstandardNth2gyy) - 2*gxyL*gxzL*JacPDstandardNth2gyz - 
      2*gyzL*(gxzL*JacPDstandardNth2gxy + gxyL*JacPDstandardNth2gxz - 
      gxxL*JacPDstandardNth2gyz) - gyyL*(gzzL*JacPDstandardNth2gxx - 
      2*gxzL*JacPDstandardNth2gxz + gxxL*JacPDstandardNth2gzz) + 
      JacPDstandardNth2gzz*pow(gxyL,2) + JacPDstandardNth2gyy*pow(gxzL,2) + 
      JacPDstandardNth2gxx*pow(gyzL,2))*pow(w_lorentzL,2) + (pressL + rhoL + 
      epsL*rhoL)*velzL*(gxxL*velxL + gxyL*velyL + 
      gxzL*velzL)*(gzzL*(2*gxyL*JacPDstandardNth3gxy - 
      gxxL*JacPDstandardNth3gyy) - 2*gxyL*gxzL*JacPDstandardNth3gyz - 
      2*gyzL*(gxzL*JacPDstandardNth3gxy + gxyL*JacPDstandardNth3gxz - 
      gxxL*JacPDstandardNth3gyz) - gyyL*(gzzL*JacPDstandardNth3gxx - 
      2*gxzL*JacPDstandardNth3gxz + gxxL*JacPDstandardNth3gzz) + 
      JacPDstandardNth3gzz*pow(gxyL,2) + JacPDstandardNth3gyy*pow(gxzL,2) + 
      JacPDstandardNth3gxx*pow(gyzL,2))*pow(w_lorentzL,2) + 
      (-JacPDstandardNth1gxy - JacPDstandardNth2gxx)*((gxzL*gyzL - 
      gxyL*gzzL)*pressL - (pressL + rhoL + 
      epsL*rhoL)*velxL*velyL*(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + pow(gyzL,2)))*pow(w_lorentzL,2)) 
      + JacPDstandardNth1gxz*((gxzL*gyyL - gxyL*gyzL)*pressL + (pressL + rhoL 
      + epsL*rhoL)*velxL*velzL*(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + pow(gyzL,2)))*pow(w_lorentzL,2)) 
      + JacPDstandardNth3gxx*((gxzL*gyyL - gxyL*gyzL)*pressL + (pressL + rhoL 
      + epsL*rhoL)*velxL*velzL*(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + pow(gyzL,2)))*pow(w_lorentzL,2)) 
      + (-JacPDstandardNth2gxz - JacPDstandardNth3gxy)*((gxyL*gxzL - 
      gxxL*gyzL)*pressL - (pressL + rhoL + 
      epsL*rhoL)*velyL*velzL*(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + pow(gyzL,2)))*pow(w_lorentzL,2)) 
      + JacPDstandardNth1gxx*(pressL*(-(gyyL*gzzL) + pow(gyzL,2)) + (pressL + 
      rhoL + epsL*rhoL)*(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + 
      pow(gyzL,2)))*pow(velxL,2)*pow(w_lorentzL,2)) + 
      JacPDstandardNth2gxy*(pressL*(-(gxxL*gzzL) + pow(gxzL,2)) + (pressL + 
      rhoL + epsL*rhoL)*(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + 
      pow(gyzL,2)))*pow(velyL,2)*pow(w_lorentzL,2)) + 
      JacPDstandardNth3gxz*(pressL*(-(gxxL*gyyL) + pow(gxyL,2)) + (pressL + 
      rhoL + epsL*rhoL)*(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + 
      pow(gyzL,2)))*pow(velzL,2)*pow(w_lorentzL,2)) + 
      (gzzL*(2*gxyL*JacPDstandardNth1gxy - gxxL*JacPDstandardNth1gyy) - 
      2*gxyL*gxzL*JacPDstandardNth1gyz - 2*gyzL*(gxzL*JacPDstandardNth1gxy + 
      gxyL*JacPDstandardNth1gxz - gxxL*JacPDstandardNth1gyz) - 
      gyyL*(gzzL*JacPDstandardNth1gxx - 2*gxzL*JacPDstandardNth1gxz + 
      gxxL*JacPDstandardNth1gzz) + JacPDstandardNth1gzz*pow(gxyL,2) + 
      JacPDstandardNth1gyy*pow(gxzL,2) + 
      JacPDstandardNth1gxx*pow(gyzL,2))*((1 + epsL)*rhoL*velxL*(gxxL*velxL + 
      gxyL*velyL + gxzL*velzL)*pow(w_lorentzL,2) + pressL*(1 + 
      gxyL*velxL*velyL*pow(w_lorentzL,2) + gxzL*velxL*velzL*pow(w_lorentzL,2) + 
      gxxL*pow(velxL,2)*pow(w_lorentzL,2))) + 
      0.5*(JacPDstandardNth1gxx*Sstress11 + JacPDstandardNth1gxz*Sstress13 + 
      JacPDstandardNth1gxy*(Sstress12 + Sstress21) + 
      JacPDstandardNth1gyy*Sstress22 + JacPDstandardNth1gyz*Sstress23 + 
      JacPDstandardNth1gxz*Sstress31 + JacPDstandardNth1gyz*Sstress32 + 
      JacPDstandardNth1gzz*Sstress33)*pow(detg,0.5) - gxzL*((-(gxzL*gyyL) + 
      gxyL*gyzL)*JacPDstandardNth1press + (gxyL*gxzL - 
      gxxL*gyzL)*JacPDstandardNth2press + JacPDstandardNth3press*(gxxL*gyyL - 
      pow(gxyL,2)) + (-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + gyyL*pow(gxzL,2) 
      + gxxL*(-(gyyL*gzzL) + pow(gyzL,2)))*(-((pressL + rhoL + 
      epsL*rhoL)*w_lorentzL*(w_lorentzL*(velxL*JacPDstandardNth1velz + 
      velyL*JacPDstandardNth2velz + velzL*(JacPDstandardNth1velx + 
      JacPDstandardNth2vely + 2*JacPDstandardNth3velz)) + 
      2*velzL*(velxL*JacPDstandardNth1w_lorentz + 
      velyL*JacPDstandardNth2w_lorentz + velzL*JacPDstandardNth3w_lorentz))) - 
      velzL*(velxL*(rhoL*JacPDstandardNth1eps + JacPDstandardNth1press + 
      JacPDstandardNth1rho + epsL*JacPDstandardNth1rho) + 
      velyL*(rhoL*JacPDstandardNth2eps + JacPDstandardNth2press + 
      JacPDstandardNth2rho + epsL*JacPDstandardNth2rho) + 
      velzL*(rhoL*JacPDstandardNth3eps + JacPDstandardNth3press + 
      JacPDstandardNth3rho + epsL*JacPDstandardNth3rho))*pow(w_lorentzL,2)) + 
      pressL*(gxyL*JacPDstandardNth1gyz + gyzL*(JacPDstandardNth1gxy - 
      JacPDstandardNth2gxx) + gxzL*(-JacPDstandardNth1gyy + 
      JacPDstandardNth2gxy) + gxyL*JacPDstandardNth2gxz + 
      gyyL*(-JacPDstandardNth1gxz + JacPDstandardNth3gxx) - 
      2*gxyL*JacPDstandardNth3gxy + gxxL*(-JacPDstandardNth2gyz + 
      JacPDstandardNth3gyy) + (-(gxxL*gyyL) + 
      pow(gxyL,2))*(2*gxyL*gzzL*JacPDstandardNth3gxy + 
      2*gxzL*gyyL*JacPDstandardNth3gxz - 2*(gxzL*gyzL*JacPDstandardNth3gxy + 
      gxyL*gyzL*JacPDstandardNth3gxz + gxyL*gxzL*JacPDstandardNth3gyz) + 
      gxxL*(-(gzzL*JacPDstandardNth3gyy) + 2*gyzL*JacPDstandardNth3gyz - 
      gyyL*JacPDstandardNth3gzz) + JacPDstandardNth3gzz*pow(gxyL,2) + 
      JacPDstandardNth3gyy*pow(gxzL,2) + JacPDstandardNth3gxx*(-(gyyL*gzzL) + 
      pow(gyzL,2)))*pow(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + pow(gyzL,2)),-1) + ((gxzL*gyyL 
      - gxyL*gyzL)*(2*gxyL*gzzL*JacPDstandardNth1gxy + 
      2*gxzL*gyyL*JacPDstandardNth1gxz - 2*(gxzL*gyzL*JacPDstandardNth1gxy + 
      gxyL*gyzL*JacPDstandardNth1gxz + gxyL*gxzL*JacPDstandardNth1gyz) + 
      gxxL*(-(gzzL*JacPDstandardNth1gyy) + 2*gyzL*JacPDstandardNth1gyz - 
      gyyL*JacPDstandardNth1gzz) + JacPDstandardNth1gzz*pow(gxyL,2) + 
      JacPDstandardNth1gyy*pow(gxzL,2) + JacPDstandardNth1gxx*(-(gyyL*gzzL) + 
      pow(gyzL,2))) - (gxyL*gxzL - 
      gxxL*gyzL)*(2*gxyL*gzzL*JacPDstandardNth2gxy + 
      2*gxzL*gyyL*JacPDstandardNth2gxz - 2*(gxzL*gyzL*JacPDstandardNth2gxy + 
      gxyL*gyzL*JacPDstandardNth2gxz + gxyL*gxzL*JacPDstandardNth2gyz) + 
      gxxL*(-(gzzL*JacPDstandardNth2gyy) + 2*gyzL*JacPDstandardNth2gyz - 
      gyyL*JacPDstandardNth2gzz) + JacPDstandardNth2gzz*pow(gxyL,2) + 
      JacPDstandardNth2gyy*pow(gxzL,2) + JacPDstandardNth2gxx*(-(gyyL*gzzL) + 
      pow(gyzL,2))))*pow(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + pow(gyzL,2)),-1))) - 
      gxxL*((gxzL*gyzL - gxyL*gzzL)*JacPDstandardNth2press + (-(gxzL*gyyL) + 
      gxyL*gyzL)*JacPDstandardNth3press + JacPDstandardNth1press*(gyyL*gzzL - 
      pow(gyzL,2)) + (-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + gyyL*pow(gxzL,2) 
      + gxxL*(-(gyyL*gzzL) + pow(gyzL,2)))*(-((pressL + rhoL + 
      epsL*rhoL)*w_lorentzL*(w_lorentzL*(velxL*JacPDstandardNth1velx + 
      velyL*JacPDstandardNth2velx + velzL*JacPDstandardNth3velx + 
      velxL*(JacPDstandardNth1velx + JacPDstandardNth2vely + 
      JacPDstandardNth3velz)) + 2*velxL*(velxL*JacPDstandardNth1w_lorentz + 
      velyL*JacPDstandardNth2w_lorentz + velzL*JacPDstandardNth3w_lorentz))) - 
      velxL*(velxL*(rhoL*JacPDstandardNth1eps + JacPDstandardNth1press + 
      JacPDstandardNth1rho + epsL*JacPDstandardNth1rho) + 
      velyL*(rhoL*JacPDstandardNth2eps + JacPDstandardNth2press + 
      JacPDstandardNth2rho + epsL*JacPDstandardNth2rho) + 
      velzL*(rhoL*JacPDstandardNth3eps + JacPDstandardNth3press + 
      JacPDstandardNth3rho + epsL*JacPDstandardNth3rho))*pow(w_lorentzL,2)) + 
      pressL*(-2*gyzL*JacPDstandardNth1gyz + gzzL*(JacPDstandardNth1gyy - 
      JacPDstandardNth2gxy) + gyzL*JacPDstandardNth2gxz + 
      gyzL*JacPDstandardNth3gxy + gyyL*(JacPDstandardNth1gzz - 
      JacPDstandardNth3gxz) + gxzL*(JacPDstandardNth2gyz - 
      JacPDstandardNth3gyy) + gxyL*(-JacPDstandardNth2gzz + 
      JacPDstandardNth3gyz) + (gxzL*gyyL - 
      gxyL*gyzL)*(2*gxyL*gzzL*JacPDstandardNth3gxy + 
      2*gxzL*gyyL*JacPDstandardNth3gxz - 2*(gxzL*gyzL*JacPDstandardNth3gxy + 
      gxyL*gyzL*JacPDstandardNth3gxz + gxyL*gxzL*JacPDstandardNth3gyz) + 
      gxxL*(-(gzzL*JacPDstandardNth3gyy) + 2*gyzL*JacPDstandardNth3gyz - 
      gyyL*JacPDstandardNth3gzz) + JacPDstandardNth3gzz*pow(gxyL,2) + 
      JacPDstandardNth3gyy*pow(gxzL,2) + JacPDstandardNth3gxx*(-(gyyL*gzzL) + 
      pow(gyzL,2)))*pow(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + pow(gyzL,2)),-1) + 
      ((-(gyyL*gzzL) + pow(gyzL,2))*(2*gxyL*gzzL*JacPDstandardNth1gxy + 
      2*gxzL*gyyL*JacPDstandardNth1gxz - 2*(gxzL*gyzL*JacPDstandardNth1gxy + 
      gxyL*gyzL*JacPDstandardNth1gxz + gxyL*gxzL*JacPDstandardNth1gyz) + 
      gxxL*(-(gzzL*JacPDstandardNth1gyy) + 2*gyzL*JacPDstandardNth1gyz - 
      gyyL*JacPDstandardNth1gzz) + JacPDstandardNth1gzz*pow(gxyL,2) + 
      JacPDstandardNth1gyy*pow(gxzL,2) + JacPDstandardNth1gxx*(-(gyyL*gzzL) + 
      pow(gyzL,2))) - (gxzL*gyzL - 
      gxyL*gzzL)*(2*gxyL*gzzL*JacPDstandardNth2gxy + 
      2*gxzL*gyyL*JacPDstandardNth2gxz - 2*(gxzL*gyzL*JacPDstandardNth2gxy + 
      gxyL*gyzL*JacPDstandardNth2gxz + gxyL*gxzL*JacPDstandardNth2gyz) + 
      gxxL*(-(gzzL*JacPDstandardNth2gyy) + 2*gyzL*JacPDstandardNth2gyz - 
      gyyL*JacPDstandardNth2gzz) + JacPDstandardNth2gzz*pow(gxyL,2) + 
      JacPDstandardNth2gyy*pow(gxzL,2) + JacPDstandardNth2gxx*(-(gyyL*gzzL) + 
      pow(gyzL,2))))*pow(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + pow(gyzL,2)),-1))) - 
      gxyL*((gxzL*gyzL - gxyL*gzzL)*JacPDstandardNth1press + (gxyL*gxzL - 
      gxxL*gyzL)*JacPDstandardNth3press + JacPDstandardNth2press*(gxxL*gzzL - 
      pow(gxzL,2)) + (-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + gyyL*pow(gxzL,2) 
      + gxxL*(-(gyyL*gzzL) + pow(gyzL,2)))*(-((pressL + rhoL + 
      epsL*rhoL)*w_lorentzL*(w_lorentzL*(velxL*JacPDstandardNth1vely + 
      velyL*JacPDstandardNth2vely + velzL*JacPDstandardNth3vely + 
      velyL*(JacPDstandardNth1velx + JacPDstandardNth2vely + 
      JacPDstandardNth3velz)) + 2*velyL*(velxL*JacPDstandardNth1w_lorentz + 
      velyL*JacPDstandardNth2w_lorentz + velzL*JacPDstandardNth3w_lorentz))) - 
      velyL*(velxL*(rhoL*JacPDstandardNth1eps + JacPDstandardNth1press + 
      JacPDstandardNth1rho + epsL*JacPDstandardNth1rho) + 
      velyL*(rhoL*JacPDstandardNth2eps + JacPDstandardNth2press + 
      JacPDstandardNth2rho + epsL*JacPDstandardNth2rho) + 
      velzL*(rhoL*JacPDstandardNth3eps + JacPDstandardNth3press + 
      JacPDstandardNth3rho + epsL*JacPDstandardNth3rho))*pow(w_lorentzL,2)) + 
      pressL*(gxzL*JacPDstandardNth1gyz + gzzL*(-JacPDstandardNth1gxy + 
      JacPDstandardNth2gxx) - 2*gxzL*JacPDstandardNth2gxz + 
      gyzL*(JacPDstandardNth1gxz - JacPDstandardNth3gxx) + 
      gxzL*JacPDstandardNth3gxy + gxyL*(-JacPDstandardNth1gzz + 
      JacPDstandardNth3gxz) + gxxL*(JacPDstandardNth2gzz - 
      JacPDstandardNth3gyz) + (-((gxzL*gyzL - 
      gxyL*gzzL)*(2*gxyL*gzzL*JacPDstandardNth1gxy + 
      2*gxzL*gyyL*JacPDstandardNth1gxz - 2*(gxzL*gyzL*JacPDstandardNth1gxy + 
      gxyL*gyzL*JacPDstandardNth1gxz + gxyL*gxzL*JacPDstandardNth1gyz) + 
      gxxL*(-(gzzL*JacPDstandardNth1gyy) + 2*gyzL*JacPDstandardNth1gyz - 
      gyyL*JacPDstandardNth1gzz) + JacPDstandardNth1gzz*pow(gxyL,2) + 
      JacPDstandardNth1gyy*pow(gxzL,2) + JacPDstandardNth1gxx*(-(gyyL*gzzL) + 
      pow(gyzL,2)))) + (-(gxxL*gzzL) + 
      pow(gxzL,2))*(2*gxyL*gzzL*JacPDstandardNth2gxy + 
      2*gxzL*gyyL*JacPDstandardNth2gxz - 2*(gxzL*gyzL*JacPDstandardNth2gxy + 
      gxyL*gyzL*JacPDstandardNth2gxz + gxyL*gxzL*JacPDstandardNth2gyz) + 
      gxxL*(-(gzzL*JacPDstandardNth2gyy) + 2*gyzL*JacPDstandardNth2gyz - 
      gyyL*JacPDstandardNth2gzz) + JacPDstandardNth2gzz*pow(gxyL,2) + 
      JacPDstandardNth2gyy*pow(gxzL,2) + JacPDstandardNth2gxx*(-(gyyL*gzzL) + 
      pow(gyzL,2))) - (gxyL*gxzL - 
      gxxL*gyzL)*(2*gxyL*gzzL*JacPDstandardNth3gxy + 
      2*gxzL*gyyL*JacPDstandardNth3gxz - 2*(gxzL*gyzL*JacPDstandardNth3gxy + 
      gxyL*gyzL*JacPDstandardNth3gxz + gxyL*gxzL*JacPDstandardNth3gyz) + 
      gxxL*(-(gzzL*JacPDstandardNth3gyy) + 2*gyzL*JacPDstandardNth3gyz - 
      gyyL*JacPDstandardNth3gzz) + JacPDstandardNth3gzz*pow(gxyL,2) + 
      JacPDstandardNth3gyy*pow(gxzL,2) + JacPDstandardNth3gxx*(-(gyyL*gzzL) + 
      pow(gyzL,2))))*pow(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + pow(gyzL,2)),-1))));
    
    CCTK_REAL SS2rhsL CCTK_ATTRIBUTE_UNUSED = betaxL*JacPDstandardNth1SS2 
      - EEL*JacPDstandardNth2alp + SS1L*JacPDstandardNth2betax + 
      SS2L*JacPDstandardNth2betay + SS3L*JacPDstandardNth2betaz + 
      betayL*JacPDstandardNth2SS2 + SS2L*(JacPDstandardNth1betax + 
      JacPDstandardNth2betay + JacPDstandardNth3betaz) + 
      betazL*JacPDstandardNth3SS2 + (-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + pow(gyzL,2)))*((1 + 
      epsL)*rhoL*velyL*(gxyL*velxL + gyyL*velyL + 
      gyzL*velzL)*JacPDstandardNth2alp*pow(w_lorentzL,2) + (pressL + rhoL + 
      epsL*rhoL)*(gxyL*velxL + gyyL*velyL + 
      gyzL*velzL)*(velxL*JacPDstandardNth1alp + 
      velzL*JacPDstandardNth3alp)*pow(w_lorentzL,2) + 
      pressL*JacPDstandardNth2alp*(1 + gxyL*velxL*velyL*pow(w_lorentzL,2) + 
      gyzL*velyL*velzL*pow(w_lorentzL,2) + 
      gyyL*pow(velyL,2)*pow(w_lorentzL,2))) + alpL*((pressL + rhoL + 
      epsL*rhoL)*velxL*(gxyL*velxL + gyyL*velyL + 
      gyzL*velzL)*(gzzL*(2*gxyL*JacPDstandardNth1gxy - 
      gxxL*JacPDstandardNth1gyy) - 2*gxyL*gxzL*JacPDstandardNth1gyz - 
      2*gyzL*(gxzL*JacPDstandardNth1gxy + gxyL*JacPDstandardNth1gxz - 
      gxxL*JacPDstandardNth1gyz) - gyyL*(gzzL*JacPDstandardNth1gxx - 
      2*gxzL*JacPDstandardNth1gxz + gxxL*JacPDstandardNth1gzz) + 
      JacPDstandardNth1gzz*pow(gxyL,2) + JacPDstandardNth1gyy*pow(gxzL,2) + 
      JacPDstandardNth1gxx*pow(gyzL,2))*pow(w_lorentzL,2) + (pressL + rhoL + 
      epsL*rhoL)*velzL*(gxyL*velxL + gyyL*velyL + 
      gyzL*velzL)*(gzzL*(2*gxyL*JacPDstandardNth3gxy - 
      gxxL*JacPDstandardNth3gyy) - 2*gxyL*gxzL*JacPDstandardNth3gyz - 
      2*gyzL*(gxzL*JacPDstandardNth3gxy + gxyL*JacPDstandardNth3gxz - 
      gxxL*JacPDstandardNth3gyz) - gyyL*(gzzL*JacPDstandardNth3gxx - 
      2*gxzL*JacPDstandardNth3gxz + gxxL*JacPDstandardNth3gzz) + 
      JacPDstandardNth3gzz*pow(gxyL,2) + JacPDstandardNth3gyy*pow(gxzL,2) + 
      JacPDstandardNth3gxx*pow(gyzL,2))*pow(w_lorentzL,2) + 
      (-JacPDstandardNth1gyy - JacPDstandardNth2gxy)*((gxzL*gyzL - 
      gxyL*gzzL)*pressL - (pressL + rhoL + 
      epsL*rhoL)*velxL*velyL*(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + pow(gyzL,2)))*pow(w_lorentzL,2)) 
      + JacPDstandardNth1gyz*((gxzL*gyyL - gxyL*gyzL)*pressL + (pressL + rhoL 
      + epsL*rhoL)*velxL*velzL*(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + pow(gyzL,2)))*pow(w_lorentzL,2)) 
      + JacPDstandardNth3gxy*((gxzL*gyyL - gxyL*gyzL)*pressL + (pressL + rhoL 
      + epsL*rhoL)*velxL*velzL*(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + pow(gyzL,2)))*pow(w_lorentzL,2)) 
      + (-JacPDstandardNth2gyz - JacPDstandardNth3gyy)*((gxyL*gxzL - 
      gxxL*gyzL)*pressL - (pressL + rhoL + 
      epsL*rhoL)*velyL*velzL*(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + pow(gyzL,2)))*pow(w_lorentzL,2)) 
      + JacPDstandardNth1gxy*(pressL*(-(gyyL*gzzL) + pow(gyzL,2)) + (pressL + 
      rhoL + epsL*rhoL)*(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + 
      pow(gyzL,2)))*pow(velxL,2)*pow(w_lorentzL,2)) + 
      JacPDstandardNth2gyy*(pressL*(-(gxxL*gzzL) + pow(gxzL,2)) + (pressL + 
      rhoL + epsL*rhoL)*(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + 
      pow(gyzL,2)))*pow(velyL,2)*pow(w_lorentzL,2)) + 
      JacPDstandardNth3gyz*(pressL*(-(gxxL*gyyL) + pow(gxyL,2)) + (pressL + 
      rhoL + epsL*rhoL)*(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + 
      pow(gyzL,2)))*pow(velzL,2)*pow(w_lorentzL,2)) + 
      (gzzL*(2*gxyL*JacPDstandardNth2gxy - gxxL*JacPDstandardNth2gyy) - 
      2*gxyL*gxzL*JacPDstandardNth2gyz - 2*gyzL*(gxzL*JacPDstandardNth2gxy + 
      gxyL*JacPDstandardNth2gxz - gxxL*JacPDstandardNth2gyz) - 
      gyyL*(gzzL*JacPDstandardNth2gxx - 2*gxzL*JacPDstandardNth2gxz + 
      gxxL*JacPDstandardNth2gzz) + JacPDstandardNth2gzz*pow(gxyL,2) + 
      JacPDstandardNth2gyy*pow(gxzL,2) + 
      JacPDstandardNth2gxx*pow(gyzL,2))*((1 + epsL)*rhoL*velyL*(gxyL*velxL + 
      gyyL*velyL + gyzL*velzL)*pow(w_lorentzL,2) + pressL*(1 + 
      gxyL*velxL*velyL*pow(w_lorentzL,2) + gyzL*velyL*velzL*pow(w_lorentzL,2) + 
      gyyL*pow(velyL,2)*pow(w_lorentzL,2))) + 
      0.5*(JacPDstandardNth2gxx*Sstress11 + JacPDstandardNth2gxz*Sstress13 + 
      JacPDstandardNth2gxy*(Sstress12 + Sstress21) + 
      JacPDstandardNth2gyy*Sstress22 + JacPDstandardNth2gyz*Sstress23 + 
      JacPDstandardNth2gxz*Sstress31 + JacPDstandardNth2gyz*Sstress32 + 
      JacPDstandardNth2gzz*Sstress33)*pow(detg,0.5) - gyzL*((-(gxzL*gyyL) + 
      gxyL*gyzL)*JacPDstandardNth1press + (gxyL*gxzL - 
      gxxL*gyzL)*JacPDstandardNth2press + JacPDstandardNth3press*(gxxL*gyyL - 
      pow(gxyL,2)) + (-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + gyyL*pow(gxzL,2) 
      + gxxL*(-(gyyL*gzzL) + pow(gyzL,2)))*(-((pressL + rhoL + 
      epsL*rhoL)*w_lorentzL*(w_lorentzL*(velxL*JacPDstandardNth1velz + 
      velyL*JacPDstandardNth2velz + velzL*(JacPDstandardNth1velx + 
      JacPDstandardNth2vely + 2*JacPDstandardNth3velz)) + 
      2*velzL*(velxL*JacPDstandardNth1w_lorentz + 
      velyL*JacPDstandardNth2w_lorentz + velzL*JacPDstandardNth3w_lorentz))) - 
      velzL*(velxL*(rhoL*JacPDstandardNth1eps + JacPDstandardNth1press + 
      JacPDstandardNth1rho + epsL*JacPDstandardNth1rho) + 
      velyL*(rhoL*JacPDstandardNth2eps + JacPDstandardNth2press + 
      JacPDstandardNth2rho + epsL*JacPDstandardNth2rho) + 
      velzL*(rhoL*JacPDstandardNth3eps + JacPDstandardNth3press + 
      JacPDstandardNth3rho + epsL*JacPDstandardNth3rho))*pow(w_lorentzL,2)) + 
      pressL*(gxyL*JacPDstandardNth1gyz + gyzL*(JacPDstandardNth1gxy - 
      JacPDstandardNth2gxx) + gxzL*(-JacPDstandardNth1gyy + 
      JacPDstandardNth2gxy) + gxyL*JacPDstandardNth2gxz + 
      gyyL*(-JacPDstandardNth1gxz + JacPDstandardNth3gxx) - 
      2*gxyL*JacPDstandardNth3gxy + gxxL*(-JacPDstandardNth2gyz + 
      JacPDstandardNth3gyy) + (-(gxxL*gyyL) + 
      pow(gxyL,2))*(2*gxyL*gzzL*JacPDstandardNth3gxy + 
      2*gxzL*gyyL*JacPDstandardNth3gxz - 2*(gxzL*gyzL*JacPDstandardNth3gxy + 
      gxyL*gyzL*JacPDstandardNth3gxz + gxyL*gxzL*JacPDstandardNth3gyz) + 
      gxxL*(-(gzzL*JacPDstandardNth3gyy) + 2*gyzL*JacPDstandardNth3gyz - 
      gyyL*JacPDstandardNth3gzz) + JacPDstandardNth3gzz*pow(gxyL,2) + 
      JacPDstandardNth3gyy*pow(gxzL,2) + JacPDstandardNth3gxx*(-(gyyL*gzzL) + 
      pow(gyzL,2)))*pow(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + pow(gyzL,2)),-1) + ((gxzL*gyyL 
      - gxyL*gyzL)*(2*gxyL*gzzL*JacPDstandardNth1gxy + 
      2*gxzL*gyyL*JacPDstandardNth1gxz - 2*(gxzL*gyzL*JacPDstandardNth1gxy + 
      gxyL*gyzL*JacPDstandardNth1gxz + gxyL*gxzL*JacPDstandardNth1gyz) + 
      gxxL*(-(gzzL*JacPDstandardNth1gyy) + 2*gyzL*JacPDstandardNth1gyz - 
      gyyL*JacPDstandardNth1gzz) + JacPDstandardNth1gzz*pow(gxyL,2) + 
      JacPDstandardNth1gyy*pow(gxzL,2) + JacPDstandardNth1gxx*(-(gyyL*gzzL) + 
      pow(gyzL,2))) - (gxyL*gxzL - 
      gxxL*gyzL)*(2*gxyL*gzzL*JacPDstandardNth2gxy + 
      2*gxzL*gyyL*JacPDstandardNth2gxz - 2*(gxzL*gyzL*JacPDstandardNth2gxy + 
      gxyL*gyzL*JacPDstandardNth2gxz + gxyL*gxzL*JacPDstandardNth2gyz) + 
      gxxL*(-(gzzL*JacPDstandardNth2gyy) + 2*gyzL*JacPDstandardNth2gyz - 
      gyyL*JacPDstandardNth2gzz) + JacPDstandardNth2gzz*pow(gxyL,2) + 
      JacPDstandardNth2gyy*pow(gxzL,2) + JacPDstandardNth2gxx*(-(gyyL*gzzL) + 
      pow(gyzL,2))))*pow(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + pow(gyzL,2)),-1))) - 
      gxyL*((gxzL*gyzL - gxyL*gzzL)*JacPDstandardNth2press + (-(gxzL*gyyL) + 
      gxyL*gyzL)*JacPDstandardNth3press + JacPDstandardNth1press*(gyyL*gzzL - 
      pow(gyzL,2)) + (-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + gyyL*pow(gxzL,2) 
      + gxxL*(-(gyyL*gzzL) + pow(gyzL,2)))*(-((pressL + rhoL + 
      epsL*rhoL)*w_lorentzL*(w_lorentzL*(velxL*JacPDstandardNth1velx + 
      velyL*JacPDstandardNth2velx + velzL*JacPDstandardNth3velx + 
      velxL*(JacPDstandardNth1velx + JacPDstandardNth2vely + 
      JacPDstandardNth3velz)) + 2*velxL*(velxL*JacPDstandardNth1w_lorentz + 
      velyL*JacPDstandardNth2w_lorentz + velzL*JacPDstandardNth3w_lorentz))) - 
      velxL*(velxL*(rhoL*JacPDstandardNth1eps + JacPDstandardNth1press + 
      JacPDstandardNth1rho + epsL*JacPDstandardNth1rho) + 
      velyL*(rhoL*JacPDstandardNth2eps + JacPDstandardNth2press + 
      JacPDstandardNth2rho + epsL*JacPDstandardNth2rho) + 
      velzL*(rhoL*JacPDstandardNth3eps + JacPDstandardNth3press + 
      JacPDstandardNth3rho + epsL*JacPDstandardNth3rho))*pow(w_lorentzL,2)) + 
      pressL*(-2*gyzL*JacPDstandardNth1gyz + gzzL*(JacPDstandardNth1gyy - 
      JacPDstandardNth2gxy) + gyzL*JacPDstandardNth2gxz + 
      gyzL*JacPDstandardNth3gxy + gyyL*(JacPDstandardNth1gzz - 
      JacPDstandardNth3gxz) + gxzL*(JacPDstandardNth2gyz - 
      JacPDstandardNth3gyy) + gxyL*(-JacPDstandardNth2gzz + 
      JacPDstandardNth3gyz) + (gxzL*gyyL - 
      gxyL*gyzL)*(2*gxyL*gzzL*JacPDstandardNth3gxy + 
      2*gxzL*gyyL*JacPDstandardNth3gxz - 2*(gxzL*gyzL*JacPDstandardNth3gxy + 
      gxyL*gyzL*JacPDstandardNth3gxz + gxyL*gxzL*JacPDstandardNth3gyz) + 
      gxxL*(-(gzzL*JacPDstandardNth3gyy) + 2*gyzL*JacPDstandardNth3gyz - 
      gyyL*JacPDstandardNth3gzz) + JacPDstandardNth3gzz*pow(gxyL,2) + 
      JacPDstandardNth3gyy*pow(gxzL,2) + JacPDstandardNth3gxx*(-(gyyL*gzzL) + 
      pow(gyzL,2)))*pow(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + pow(gyzL,2)),-1) + 
      ((-(gyyL*gzzL) + pow(gyzL,2))*(2*gxyL*gzzL*JacPDstandardNth1gxy + 
      2*gxzL*gyyL*JacPDstandardNth1gxz - 2*(gxzL*gyzL*JacPDstandardNth1gxy + 
      gxyL*gyzL*JacPDstandardNth1gxz + gxyL*gxzL*JacPDstandardNth1gyz) + 
      gxxL*(-(gzzL*JacPDstandardNth1gyy) + 2*gyzL*JacPDstandardNth1gyz - 
      gyyL*JacPDstandardNth1gzz) + JacPDstandardNth1gzz*pow(gxyL,2) + 
      JacPDstandardNth1gyy*pow(gxzL,2) + JacPDstandardNth1gxx*(-(gyyL*gzzL) + 
      pow(gyzL,2))) - (gxzL*gyzL - 
      gxyL*gzzL)*(2*gxyL*gzzL*JacPDstandardNth2gxy + 
      2*gxzL*gyyL*JacPDstandardNth2gxz - 2*(gxzL*gyzL*JacPDstandardNth2gxy + 
      gxyL*gyzL*JacPDstandardNth2gxz + gxyL*gxzL*JacPDstandardNth2gyz) + 
      gxxL*(-(gzzL*JacPDstandardNth2gyy) + 2*gyzL*JacPDstandardNth2gyz - 
      gyyL*JacPDstandardNth2gzz) + JacPDstandardNth2gzz*pow(gxyL,2) + 
      JacPDstandardNth2gyy*pow(gxzL,2) + JacPDstandardNth2gxx*(-(gyyL*gzzL) + 
      pow(gyzL,2))))*pow(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + pow(gyzL,2)),-1))) - 
      gyyL*((gxzL*gyzL - gxyL*gzzL)*JacPDstandardNth1press + (gxyL*gxzL - 
      gxxL*gyzL)*JacPDstandardNth3press + JacPDstandardNth2press*(gxxL*gzzL - 
      pow(gxzL,2)) + (-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + gyyL*pow(gxzL,2) 
      + gxxL*(-(gyyL*gzzL) + pow(gyzL,2)))*(-((pressL + rhoL + 
      epsL*rhoL)*w_lorentzL*(w_lorentzL*(velxL*JacPDstandardNth1vely + 
      velyL*JacPDstandardNth2vely + velzL*JacPDstandardNth3vely + 
      velyL*(JacPDstandardNth1velx + JacPDstandardNth2vely + 
      JacPDstandardNth3velz)) + 2*velyL*(velxL*JacPDstandardNth1w_lorentz + 
      velyL*JacPDstandardNth2w_lorentz + velzL*JacPDstandardNth3w_lorentz))) - 
      velyL*(velxL*(rhoL*JacPDstandardNth1eps + JacPDstandardNth1press + 
      JacPDstandardNth1rho + epsL*JacPDstandardNth1rho) + 
      velyL*(rhoL*JacPDstandardNth2eps + JacPDstandardNth2press + 
      JacPDstandardNth2rho + epsL*JacPDstandardNth2rho) + 
      velzL*(rhoL*JacPDstandardNth3eps + JacPDstandardNth3press + 
      JacPDstandardNth3rho + epsL*JacPDstandardNth3rho))*pow(w_lorentzL,2)) + 
      pressL*(gxzL*JacPDstandardNth1gyz + gzzL*(-JacPDstandardNth1gxy + 
      JacPDstandardNth2gxx) - 2*gxzL*JacPDstandardNth2gxz + 
      gyzL*(JacPDstandardNth1gxz - JacPDstandardNth3gxx) + 
      gxzL*JacPDstandardNth3gxy + gxyL*(-JacPDstandardNth1gzz + 
      JacPDstandardNth3gxz) + gxxL*(JacPDstandardNth2gzz - 
      JacPDstandardNth3gyz) + (-((gxzL*gyzL - 
      gxyL*gzzL)*(2*gxyL*gzzL*JacPDstandardNth1gxy + 
      2*gxzL*gyyL*JacPDstandardNth1gxz - 2*(gxzL*gyzL*JacPDstandardNth1gxy + 
      gxyL*gyzL*JacPDstandardNth1gxz + gxyL*gxzL*JacPDstandardNth1gyz) + 
      gxxL*(-(gzzL*JacPDstandardNth1gyy) + 2*gyzL*JacPDstandardNth1gyz - 
      gyyL*JacPDstandardNth1gzz) + JacPDstandardNth1gzz*pow(gxyL,2) + 
      JacPDstandardNth1gyy*pow(gxzL,2) + JacPDstandardNth1gxx*(-(gyyL*gzzL) + 
      pow(gyzL,2)))) + (-(gxxL*gzzL) + 
      pow(gxzL,2))*(2*gxyL*gzzL*JacPDstandardNth2gxy + 
      2*gxzL*gyyL*JacPDstandardNth2gxz - 2*(gxzL*gyzL*JacPDstandardNth2gxy + 
      gxyL*gyzL*JacPDstandardNth2gxz + gxyL*gxzL*JacPDstandardNth2gyz) + 
      gxxL*(-(gzzL*JacPDstandardNth2gyy) + 2*gyzL*JacPDstandardNth2gyz - 
      gyyL*JacPDstandardNth2gzz) + JacPDstandardNth2gzz*pow(gxyL,2) + 
      JacPDstandardNth2gyy*pow(gxzL,2) + JacPDstandardNth2gxx*(-(gyyL*gzzL) + 
      pow(gyzL,2))) - (gxyL*gxzL - 
      gxxL*gyzL)*(2*gxyL*gzzL*JacPDstandardNth3gxy + 
      2*gxzL*gyyL*JacPDstandardNth3gxz - 2*(gxzL*gyzL*JacPDstandardNth3gxy + 
      gxyL*gyzL*JacPDstandardNth3gxz + gxyL*gxzL*JacPDstandardNth3gyz) + 
      gxxL*(-(gzzL*JacPDstandardNth3gyy) + 2*gyzL*JacPDstandardNth3gyz - 
      gyyL*JacPDstandardNth3gzz) + JacPDstandardNth3gzz*pow(gxyL,2) + 
      JacPDstandardNth3gyy*pow(gxzL,2) + JacPDstandardNth3gxx*(-(gyyL*gzzL) + 
      pow(gyzL,2))))*pow(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + pow(gyzL,2)),-1))));
    
    CCTK_REAL SS3rhsL CCTK_ATTRIBUTE_UNUSED = betaxL*JacPDstandardNth1SS3 
      + betayL*JacPDstandardNth2SS3 - EEL*JacPDstandardNth3alp + 
      SS1L*JacPDstandardNth3betax + SS2L*JacPDstandardNth3betay + 
      SS3L*JacPDstandardNth3betaz + SS3L*(JacPDstandardNth1betax + 
      JacPDstandardNth2betay + JacPDstandardNth3betaz) + 
      betazL*JacPDstandardNth3SS3 + (-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + pow(gyzL,2)))*((pressL + rhoL + 
      epsL*rhoL)*(gxzL*velxL + gyzL*velyL + 
      gzzL*velzL)*(velxL*JacPDstandardNth1alp + 
      velyL*JacPDstandardNth2alp)*pow(w_lorentzL,2) + (1 + 
      epsL)*rhoL*velzL*(gxzL*velxL + gyzL*velyL + 
      gzzL*velzL)*JacPDstandardNth3alp*pow(w_lorentzL,2) + 
      pressL*JacPDstandardNth3alp*(1 + gxzL*velxL*velzL*pow(w_lorentzL,2) + 
      gyzL*velyL*velzL*pow(w_lorentzL,2) + 
      gzzL*pow(velzL,2)*pow(w_lorentzL,2))) + alpL*((pressL + rhoL + 
      epsL*rhoL)*velxL*(gxzL*velxL + gyzL*velyL + 
      gzzL*velzL)*(gzzL*(2*gxyL*JacPDstandardNth1gxy - 
      gxxL*JacPDstandardNth1gyy) - 2*gxyL*gxzL*JacPDstandardNth1gyz - 
      2*gyzL*(gxzL*JacPDstandardNth1gxy + gxyL*JacPDstandardNth1gxz - 
      gxxL*JacPDstandardNth1gyz) - gyyL*(gzzL*JacPDstandardNth1gxx - 
      2*gxzL*JacPDstandardNth1gxz + gxxL*JacPDstandardNth1gzz) + 
      JacPDstandardNth1gzz*pow(gxyL,2) + JacPDstandardNth1gyy*pow(gxzL,2) + 
      JacPDstandardNth1gxx*pow(gyzL,2))*pow(w_lorentzL,2) + (pressL + rhoL + 
      epsL*rhoL)*velyL*(gxzL*velxL + gyzL*velyL + 
      gzzL*velzL)*(gzzL*(2*gxyL*JacPDstandardNth2gxy - 
      gxxL*JacPDstandardNth2gyy) - 2*gxyL*gxzL*JacPDstandardNth2gyz - 
      2*gyzL*(gxzL*JacPDstandardNth2gxy + gxyL*JacPDstandardNth2gxz - 
      gxxL*JacPDstandardNth2gyz) - gyyL*(gzzL*JacPDstandardNth2gxx - 
      2*gxzL*JacPDstandardNth2gxz + gxxL*JacPDstandardNth2gzz) + 
      JacPDstandardNth2gzz*pow(gxyL,2) + JacPDstandardNth2gyy*pow(gxzL,2) + 
      JacPDstandardNth2gxx*pow(gyzL,2))*pow(w_lorentzL,2) + 
      (-JacPDstandardNth1gyz - JacPDstandardNth2gxz)*((gxzL*gyzL - 
      gxyL*gzzL)*pressL - (pressL + rhoL + 
      epsL*rhoL)*velxL*velyL*(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + pow(gyzL,2)))*pow(w_lorentzL,2)) 
      + JacPDstandardNth1gzz*((gxzL*gyyL - gxyL*gyzL)*pressL + (pressL + rhoL 
      + epsL*rhoL)*velxL*velzL*(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + pow(gyzL,2)))*pow(w_lorentzL,2)) 
      + JacPDstandardNth3gxz*((gxzL*gyyL - gxyL*gyzL)*pressL + (pressL + rhoL 
      + epsL*rhoL)*velxL*velzL*(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + pow(gyzL,2)))*pow(w_lorentzL,2)) 
      + (-JacPDstandardNth2gzz - JacPDstandardNth3gyz)*((gxyL*gxzL - 
      gxxL*gyzL)*pressL - (pressL + rhoL + 
      epsL*rhoL)*velyL*velzL*(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + pow(gyzL,2)))*pow(w_lorentzL,2)) 
      + JacPDstandardNth1gxz*(pressL*(-(gyyL*gzzL) + pow(gyzL,2)) + (pressL + 
      rhoL + epsL*rhoL)*(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + 
      pow(gyzL,2)))*pow(velxL,2)*pow(w_lorentzL,2)) + 
      JacPDstandardNth2gyz*(pressL*(-(gxxL*gzzL) + pow(gxzL,2)) + (pressL + 
      rhoL + epsL*rhoL)*(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + 
      pow(gyzL,2)))*pow(velyL,2)*pow(w_lorentzL,2)) + 
      JacPDstandardNth3gzz*(pressL*(-(gxxL*gyyL) + pow(gxyL,2)) + (pressL + 
      rhoL + epsL*rhoL)*(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + 
      pow(gyzL,2)))*pow(velzL,2)*pow(w_lorentzL,2)) + 
      (gzzL*(2*gxyL*JacPDstandardNth3gxy - gxxL*JacPDstandardNth3gyy) - 
      2*gxyL*gxzL*JacPDstandardNth3gyz - 2*gyzL*(gxzL*JacPDstandardNth3gxy + 
      gxyL*JacPDstandardNth3gxz - gxxL*JacPDstandardNth3gyz) - 
      gyyL*(gzzL*JacPDstandardNth3gxx - 2*gxzL*JacPDstandardNth3gxz + 
      gxxL*JacPDstandardNth3gzz) + JacPDstandardNth3gzz*pow(gxyL,2) + 
      JacPDstandardNth3gyy*pow(gxzL,2) + 
      JacPDstandardNth3gxx*pow(gyzL,2))*((1 + epsL)*rhoL*velzL*(gxzL*velxL + 
      gyzL*velyL + gzzL*velzL)*pow(w_lorentzL,2) + pressL*(1 + 
      gxzL*velxL*velzL*pow(w_lorentzL,2) + gyzL*velyL*velzL*pow(w_lorentzL,2) + 
      gzzL*pow(velzL,2)*pow(w_lorentzL,2))) + 
      0.5*(JacPDstandardNth3gxx*Sstress11 + JacPDstandardNth3gxz*Sstress13 + 
      JacPDstandardNth3gxy*(Sstress12 + Sstress21) + 
      JacPDstandardNth3gyy*Sstress22 + JacPDstandardNth3gyz*Sstress23 + 
      JacPDstandardNth3gxz*Sstress31 + JacPDstandardNth3gyz*Sstress32 + 
      JacPDstandardNth3gzz*Sstress33)*pow(detg,0.5) - gzzL*((-(gxzL*gyyL) + 
      gxyL*gyzL)*JacPDstandardNth1press + (gxyL*gxzL - 
      gxxL*gyzL)*JacPDstandardNth2press + JacPDstandardNth3press*(gxxL*gyyL - 
      pow(gxyL,2)) + (-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + gyyL*pow(gxzL,2) 
      + gxxL*(-(gyyL*gzzL) + pow(gyzL,2)))*(-((pressL + rhoL + 
      epsL*rhoL)*w_lorentzL*(w_lorentzL*(velxL*JacPDstandardNth1velz + 
      velyL*JacPDstandardNth2velz + velzL*(JacPDstandardNth1velx + 
      JacPDstandardNth2vely + 2*JacPDstandardNth3velz)) + 
      2*velzL*(velxL*JacPDstandardNth1w_lorentz + 
      velyL*JacPDstandardNth2w_lorentz + velzL*JacPDstandardNth3w_lorentz))) - 
      velzL*(velxL*(rhoL*JacPDstandardNth1eps + JacPDstandardNth1press + 
      JacPDstandardNth1rho + epsL*JacPDstandardNth1rho) + 
      velyL*(rhoL*JacPDstandardNth2eps + JacPDstandardNth2press + 
      JacPDstandardNth2rho + epsL*JacPDstandardNth2rho) + 
      velzL*(rhoL*JacPDstandardNth3eps + JacPDstandardNth3press + 
      JacPDstandardNth3rho + epsL*JacPDstandardNth3rho))*pow(w_lorentzL,2)) + 
      pressL*(gxyL*JacPDstandardNth1gyz + gyzL*(JacPDstandardNth1gxy - 
      JacPDstandardNth2gxx) + gxzL*(-JacPDstandardNth1gyy + 
      JacPDstandardNth2gxy) + gxyL*JacPDstandardNth2gxz + 
      gyyL*(-JacPDstandardNth1gxz + JacPDstandardNth3gxx) - 
      2*gxyL*JacPDstandardNth3gxy + gxxL*(-JacPDstandardNth2gyz + 
      JacPDstandardNth3gyy) + (-(gxxL*gyyL) + 
      pow(gxyL,2))*(2*gxyL*gzzL*JacPDstandardNth3gxy + 
      2*gxzL*gyyL*JacPDstandardNth3gxz - 2*(gxzL*gyzL*JacPDstandardNth3gxy + 
      gxyL*gyzL*JacPDstandardNth3gxz + gxyL*gxzL*JacPDstandardNth3gyz) + 
      gxxL*(-(gzzL*JacPDstandardNth3gyy) + 2*gyzL*JacPDstandardNth3gyz - 
      gyyL*JacPDstandardNth3gzz) + JacPDstandardNth3gzz*pow(gxyL,2) + 
      JacPDstandardNth3gyy*pow(gxzL,2) + JacPDstandardNth3gxx*(-(gyyL*gzzL) + 
      pow(gyzL,2)))*pow(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + pow(gyzL,2)),-1) + ((gxzL*gyyL 
      - gxyL*gyzL)*(2*gxyL*gzzL*JacPDstandardNth1gxy + 
      2*gxzL*gyyL*JacPDstandardNth1gxz - 2*(gxzL*gyzL*JacPDstandardNth1gxy + 
      gxyL*gyzL*JacPDstandardNth1gxz + gxyL*gxzL*JacPDstandardNth1gyz) + 
      gxxL*(-(gzzL*JacPDstandardNth1gyy) + 2*gyzL*JacPDstandardNth1gyz - 
      gyyL*JacPDstandardNth1gzz) + JacPDstandardNth1gzz*pow(gxyL,2) + 
      JacPDstandardNth1gyy*pow(gxzL,2) + JacPDstandardNth1gxx*(-(gyyL*gzzL) + 
      pow(gyzL,2))) - (gxyL*gxzL - 
      gxxL*gyzL)*(2*gxyL*gzzL*JacPDstandardNth2gxy + 
      2*gxzL*gyyL*JacPDstandardNth2gxz - 2*(gxzL*gyzL*JacPDstandardNth2gxy + 
      gxyL*gyzL*JacPDstandardNth2gxz + gxyL*gxzL*JacPDstandardNth2gyz) + 
      gxxL*(-(gzzL*JacPDstandardNth2gyy) + 2*gyzL*JacPDstandardNth2gyz - 
      gyyL*JacPDstandardNth2gzz) + JacPDstandardNth2gzz*pow(gxyL,2) + 
      JacPDstandardNth2gyy*pow(gxzL,2) + JacPDstandardNth2gxx*(-(gyyL*gzzL) + 
      pow(gyzL,2))))*pow(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + pow(gyzL,2)),-1))) - 
      gxzL*((gxzL*gyzL - gxyL*gzzL)*JacPDstandardNth2press + (-(gxzL*gyyL) + 
      gxyL*gyzL)*JacPDstandardNth3press + JacPDstandardNth1press*(gyyL*gzzL - 
      pow(gyzL,2)) + (-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + gyyL*pow(gxzL,2) 
      + gxxL*(-(gyyL*gzzL) + pow(gyzL,2)))*(-((pressL + rhoL + 
      epsL*rhoL)*w_lorentzL*(w_lorentzL*(velxL*JacPDstandardNth1velx + 
      velyL*JacPDstandardNth2velx + velzL*JacPDstandardNth3velx + 
      velxL*(JacPDstandardNth1velx + JacPDstandardNth2vely + 
      JacPDstandardNth3velz)) + 2*velxL*(velxL*JacPDstandardNth1w_lorentz + 
      velyL*JacPDstandardNth2w_lorentz + velzL*JacPDstandardNth3w_lorentz))) - 
      velxL*(velxL*(rhoL*JacPDstandardNth1eps + JacPDstandardNth1press + 
      JacPDstandardNth1rho + epsL*JacPDstandardNth1rho) + 
      velyL*(rhoL*JacPDstandardNth2eps + JacPDstandardNth2press + 
      JacPDstandardNth2rho + epsL*JacPDstandardNth2rho) + 
      velzL*(rhoL*JacPDstandardNth3eps + JacPDstandardNth3press + 
      JacPDstandardNth3rho + epsL*JacPDstandardNth3rho))*pow(w_lorentzL,2)) + 
      pressL*(-2*gyzL*JacPDstandardNth1gyz + gzzL*(JacPDstandardNth1gyy - 
      JacPDstandardNth2gxy) + gyzL*JacPDstandardNth2gxz + 
      gyzL*JacPDstandardNth3gxy + gyyL*(JacPDstandardNth1gzz - 
      JacPDstandardNth3gxz) + gxzL*(JacPDstandardNth2gyz - 
      JacPDstandardNth3gyy) + gxyL*(-JacPDstandardNth2gzz + 
      JacPDstandardNth3gyz) + (gxzL*gyyL - 
      gxyL*gyzL)*(2*gxyL*gzzL*JacPDstandardNth3gxy + 
      2*gxzL*gyyL*JacPDstandardNth3gxz - 2*(gxzL*gyzL*JacPDstandardNth3gxy + 
      gxyL*gyzL*JacPDstandardNth3gxz + gxyL*gxzL*JacPDstandardNth3gyz) + 
      gxxL*(-(gzzL*JacPDstandardNth3gyy) + 2*gyzL*JacPDstandardNth3gyz - 
      gyyL*JacPDstandardNth3gzz) + JacPDstandardNth3gzz*pow(gxyL,2) + 
      JacPDstandardNth3gyy*pow(gxzL,2) + JacPDstandardNth3gxx*(-(gyyL*gzzL) + 
      pow(gyzL,2)))*pow(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + pow(gyzL,2)),-1) + 
      ((-(gyyL*gzzL) + pow(gyzL,2))*(2*gxyL*gzzL*JacPDstandardNth1gxy + 
      2*gxzL*gyyL*JacPDstandardNth1gxz - 2*(gxzL*gyzL*JacPDstandardNth1gxy + 
      gxyL*gyzL*JacPDstandardNth1gxz + gxyL*gxzL*JacPDstandardNth1gyz) + 
      gxxL*(-(gzzL*JacPDstandardNth1gyy) + 2*gyzL*JacPDstandardNth1gyz - 
      gyyL*JacPDstandardNth1gzz) + JacPDstandardNth1gzz*pow(gxyL,2) + 
      JacPDstandardNth1gyy*pow(gxzL,2) + JacPDstandardNth1gxx*(-(gyyL*gzzL) + 
      pow(gyzL,2))) - (gxzL*gyzL - 
      gxyL*gzzL)*(2*gxyL*gzzL*JacPDstandardNth2gxy + 
      2*gxzL*gyyL*JacPDstandardNth2gxz - 2*(gxzL*gyzL*JacPDstandardNth2gxy + 
      gxyL*gyzL*JacPDstandardNth2gxz + gxyL*gxzL*JacPDstandardNth2gyz) + 
      gxxL*(-(gzzL*JacPDstandardNth2gyy) + 2*gyzL*JacPDstandardNth2gyz - 
      gyyL*JacPDstandardNth2gzz) + JacPDstandardNth2gzz*pow(gxyL,2) + 
      JacPDstandardNth2gyy*pow(gxzL,2) + JacPDstandardNth2gxx*(-(gyyL*gzzL) + 
      pow(gyzL,2))))*pow(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + pow(gyzL,2)),-1))) - 
      gyzL*((gxzL*gyzL - gxyL*gzzL)*JacPDstandardNth1press + (gxyL*gxzL - 
      gxxL*gyzL)*JacPDstandardNth3press + JacPDstandardNth2press*(gxxL*gzzL - 
      pow(gxzL,2)) + (-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + gyyL*pow(gxzL,2) 
      + gxxL*(-(gyyL*gzzL) + pow(gyzL,2)))*(-((pressL + rhoL + 
      epsL*rhoL)*w_lorentzL*(w_lorentzL*(velxL*JacPDstandardNth1vely + 
      velyL*JacPDstandardNth2vely + velzL*JacPDstandardNth3vely + 
      velyL*(JacPDstandardNth1velx + JacPDstandardNth2vely + 
      JacPDstandardNth3velz)) + 2*velyL*(velxL*JacPDstandardNth1w_lorentz + 
      velyL*JacPDstandardNth2w_lorentz + velzL*JacPDstandardNth3w_lorentz))) - 
      velyL*(velxL*(rhoL*JacPDstandardNth1eps + JacPDstandardNth1press + 
      JacPDstandardNth1rho + epsL*JacPDstandardNth1rho) + 
      velyL*(rhoL*JacPDstandardNth2eps + JacPDstandardNth2press + 
      JacPDstandardNth2rho + epsL*JacPDstandardNth2rho) + 
      velzL*(rhoL*JacPDstandardNth3eps + JacPDstandardNth3press + 
      JacPDstandardNth3rho + epsL*JacPDstandardNth3rho))*pow(w_lorentzL,2)) + 
      pressL*(gxzL*JacPDstandardNth1gyz + gzzL*(-JacPDstandardNth1gxy + 
      JacPDstandardNth2gxx) - 2*gxzL*JacPDstandardNth2gxz + 
      gyzL*(JacPDstandardNth1gxz - JacPDstandardNth3gxx) + 
      gxzL*JacPDstandardNth3gxy + gxyL*(-JacPDstandardNth1gzz + 
      JacPDstandardNth3gxz) + gxxL*(JacPDstandardNth2gzz - 
      JacPDstandardNth3gyz) + (-((gxzL*gyzL - 
      gxyL*gzzL)*(2*gxyL*gzzL*JacPDstandardNth1gxy + 
      2*gxzL*gyyL*JacPDstandardNth1gxz - 2*(gxzL*gyzL*JacPDstandardNth1gxy + 
      gxyL*gyzL*JacPDstandardNth1gxz + gxyL*gxzL*JacPDstandardNth1gyz) + 
      gxxL*(-(gzzL*JacPDstandardNth1gyy) + 2*gyzL*JacPDstandardNth1gyz - 
      gyyL*JacPDstandardNth1gzz) + JacPDstandardNth1gzz*pow(gxyL,2) + 
      JacPDstandardNth1gyy*pow(gxzL,2) + JacPDstandardNth1gxx*(-(gyyL*gzzL) + 
      pow(gyzL,2)))) + (-(gxxL*gzzL) + 
      pow(gxzL,2))*(2*gxyL*gzzL*JacPDstandardNth2gxy + 
      2*gxzL*gyyL*JacPDstandardNth2gxz - 2*(gxzL*gyzL*JacPDstandardNth2gxy + 
      gxyL*gyzL*JacPDstandardNth2gxz + gxyL*gxzL*JacPDstandardNth2gyz) + 
      gxxL*(-(gzzL*JacPDstandardNth2gyy) + 2*gyzL*JacPDstandardNth2gyz - 
      gyyL*JacPDstandardNth2gzz) + JacPDstandardNth2gzz*pow(gxyL,2) + 
      JacPDstandardNth2gyy*pow(gxzL,2) + JacPDstandardNth2gxx*(-(gyyL*gzzL) + 
      pow(gyzL,2))) - (gxyL*gxzL - 
      gxxL*gyzL)*(2*gxyL*gzzL*JacPDstandardNth3gxy + 
      2*gxzL*gyyL*JacPDstandardNth3gxz - 2*(gxzL*gyzL*JacPDstandardNth3gxy + 
      gxyL*gyzL*JacPDstandardNth3gxz + gxyL*gxzL*JacPDstandardNth3gyz) + 
      gxxL*(-(gzzL*JacPDstandardNth3gyy) + 2*gyzL*JacPDstandardNth3gyz - 
      gyyL*JacPDstandardNth3gzz) + JacPDstandardNth3gzz*pow(gxyL,2) + 
      JacPDstandardNth3gyy*pow(gxzL,2) + JacPDstandardNth3gxx*(-(gyyL*gzzL) + 
      pow(gyzL,2))))*pow(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + pow(gyzL,2)),-1))));
    
    CCTK_REAL EErhsL CCTK_ATTRIBUTE_UNUSED = betaxL*JacPDstandardNth1EE + 
      betayL*JacPDstandardNth2EE + EEL*(JacPDstandardNth1betax + 
      JacPDstandardNth2betay + JacPDstandardNth3betaz) + 
      betazL*JacPDstandardNth3EE + alpL*((kxxL*Sstress11 + kxzL*Sstress13 + 
      kxyL*(Sstress12 + Sstress21) + kyyL*Sstress22 + kyzL*Sstress23 + 
      kxzL*Sstress31 + kyzL*Sstress32 + kzzL*Sstress33)*pow(detg,0.5) - 
      (((gxzL*gyyL - gxyL*gyzL)*JacPDstandardNth1SS3 - (gxyL*gxzL - 
      gxxL*gyzL)*JacPDstandardNth2SS3)*(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) 
      + gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + pow(gyzL,2))) + (gxzL*gyyL - 
      gxyL*gyzL)*JacPDstandardNth3SS1*(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + pow(gyzL,2))) + 
      JacPDstandardNth3SS3*(-(gxxL*gyyL) + pow(gxyL,2))*(-2*gxyL*gxzL*gyzL + 
      gzzL*pow(gxyL,2) + gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + 
      pow(gyzL,2))) + JacPDstandardNth2SS2*(-(gxxL*gzzL) + 
      pow(gxzL,2))*(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + gyyL*pow(gxzL,2) + 
      gxxL*(-(gyyL*gzzL) + pow(gyzL,2))) + (-2*gxyL*gxzL*gyzL + 
      gzzL*pow(gxyL,2) + gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + 
      pow(gyzL,2)))*((gxzL*gyzL - gxyL*gzzL)*(-JacPDstandardNth1SS2 - 
      JacPDstandardNth2SS1) - (gxyL*gxzL - gxxL*gyzL)*JacPDstandardNth3SS2 + 
      JacPDstandardNth1SS1*(-(gyyL*gzzL) + pow(gyzL,2))) + 
      SS2L*((-(gyzL*JacPDstandardNth1gxz) + gxyL*JacPDstandardNth1gzz + 
      gzzL*(JacPDstandardNth1gxy - JacPDstandardNth2gxx) + 
      gxzL*(-JacPDstandardNth1gyz + 2*JacPDstandardNth2gxz) - 
      gxxL*JacPDstandardNth2gzz)*(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + pow(gyzL,2))) + 
      (gyzL*JacPDstandardNth3gxx - gxzL*JacPDstandardNth3gxy - 
      gxyL*JacPDstandardNth3gxz + 
      gxxL*JacPDstandardNth3gyz)*(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + pow(gyzL,2))) + (gxzL*gyzL - 
      gxyL*gzzL)*(2*gxyL*gzzL*JacPDstandardNth1gxy + 
      2*gxzL*gyyL*JacPDstandardNth1gxz - 2*(gxzL*gyzL*JacPDstandardNth1gxy + 
      gxyL*gyzL*JacPDstandardNth1gxz + gxyL*gxzL*JacPDstandardNth1gyz) + 
      gxxL*(-(gzzL*JacPDstandardNth1gyy) + 2*gyzL*JacPDstandardNth1gyz - 
      gyyL*JacPDstandardNth1gzz) + JacPDstandardNth1gzz*pow(gxyL,2) + 
      JacPDstandardNth1gyy*pow(gxzL,2) + JacPDstandardNth1gxx*(-(gyyL*gzzL) + 
      pow(gyzL,2))) - (-(gxxL*gzzL) + 
      pow(gxzL,2))*(2*gxyL*gzzL*JacPDstandardNth2gxy + 
      2*gxzL*gyyL*JacPDstandardNth2gxz - 2*(gxzL*gyzL*JacPDstandardNth2gxy + 
      gxyL*gyzL*JacPDstandardNth2gxz + gxyL*gxzL*JacPDstandardNth2gyz) + 
      gxxL*(-(gzzL*JacPDstandardNth2gyy) + 2*gyzL*JacPDstandardNth2gyz - 
      gyyL*JacPDstandardNth2gzz) + JacPDstandardNth2gzz*pow(gxyL,2) + 
      JacPDstandardNth2gyy*pow(gxzL,2) + JacPDstandardNth2gxx*(-(gyyL*gzzL) + 
      pow(gyzL,2))) + (gxyL*gxzL - 
      gxxL*gyzL)*(2*gxyL*gzzL*JacPDstandardNth3gxy + 
      2*gxzL*gyyL*JacPDstandardNth3gxz - 2*(gxzL*gyzL*JacPDstandardNth3gxy + 
      gxyL*gyzL*JacPDstandardNth3gxz + gxyL*gxzL*JacPDstandardNth3gyz) + 
      gxxL*(-(gzzL*JacPDstandardNth3gyy) + 2*gyzL*JacPDstandardNth3gyz - 
      gyyL*JacPDstandardNth3gzz) + JacPDstandardNth3gzz*pow(gxyL,2) + 
      JacPDstandardNth3gyy*pow(gxzL,2) + JacPDstandardNth3gxx*(-(gyyL*gzzL) + 
      pow(gyzL,2)))) - SS1L*((gzzL*JacPDstandardNth1gyy - 
      2*gyzL*JacPDstandardNth1gyz + 
      gyyL*JacPDstandardNth1gzz)*(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + pow(gyzL,2))) + 
      (-(gzzL*JacPDstandardNth2gxy) + gyzL*JacPDstandardNth2gxz + 
      gxzL*JacPDstandardNth2gyz - 
      gxyL*JacPDstandardNth2gzz)*(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + pow(gyzL,2))) + 
      (gyzL*JacPDstandardNth3gxy - gyyL*JacPDstandardNth3gxz - 
      gxzL*JacPDstandardNth3gyy + 
      gxyL*JacPDstandardNth3gyz)*(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + pow(gyzL,2))) + (-(gyyL*gzzL) + 
      pow(gyzL,2))*(2*gxyL*gzzL*JacPDstandardNth1gxy + 
      2*gxzL*gyyL*JacPDstandardNth1gxz - 2*(gxzL*gyzL*JacPDstandardNth1gxy + 
      gxyL*gyzL*JacPDstandardNth1gxz + gxyL*gxzL*JacPDstandardNth1gyz) + 
      gxxL*(-(gzzL*JacPDstandardNth1gyy) + 2*gyzL*JacPDstandardNth1gyz - 
      gyyL*JacPDstandardNth1gzz) + JacPDstandardNth1gzz*pow(gxyL,2) + 
      JacPDstandardNth1gyy*pow(gxzL,2) + JacPDstandardNth1gxx*(-(gyyL*gzzL) + 
      pow(gyzL,2))) - (gxzL*gyzL - 
      gxyL*gzzL)*(2*gxyL*gzzL*JacPDstandardNth2gxy + 
      2*gxzL*gyyL*JacPDstandardNth2gxz - 2*(gxzL*gyzL*JacPDstandardNth2gxy + 
      gxyL*gyzL*JacPDstandardNth2gxz + gxyL*gxzL*JacPDstandardNth2gyz) + 
      gxxL*(-(gzzL*JacPDstandardNth2gyy) + 2*gyzL*JacPDstandardNth2gyz - 
      gyyL*JacPDstandardNth2gzz) + JacPDstandardNth2gzz*pow(gxyL,2) + 
      JacPDstandardNth2gyy*pow(gxzL,2) + JacPDstandardNth2gxx*(-(gyyL*gzzL) + 
      pow(gyzL,2))) + (gxzL*gyyL - 
      gxyL*gyzL)*(2*gxyL*gzzL*JacPDstandardNth3gxy + 
      2*gxzL*gyyL*JacPDstandardNth3gxz - 2*(gxzL*gyzL*JacPDstandardNth3gxy + 
      gxyL*gyzL*JacPDstandardNth3gxz + gxyL*gxzL*JacPDstandardNth3gyz) + 
      gxxL*(-(gzzL*JacPDstandardNth3gyy) + 2*gyzL*JacPDstandardNth3gyz - 
      gyyL*JacPDstandardNth3gzz) + JacPDstandardNth3gzz*pow(gxyL,2) + 
      JacPDstandardNth3gyy*pow(gxzL,2) + JacPDstandardNth3gxx*(-(gyyL*gzzL) + 
      pow(gyzL,2)))) - SS3L*((gyzL*JacPDstandardNth1gxy - 
      gyyL*JacPDstandardNth1gxz - gxzL*JacPDstandardNth1gyy + 
      gxyL*JacPDstandardNth1gyz)*(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + pow(gyzL,2))) + 
      (-(gyzL*JacPDstandardNth2gxx) + gxzL*JacPDstandardNth2gxy + 
      gxyL*JacPDstandardNth2gxz - 
      gxxL*JacPDstandardNth2gyz)*(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + pow(gyzL,2))) + 
      (gyyL*JacPDstandardNth3gxx - 2*gxyL*JacPDstandardNth3gxy + 
      gxxL*JacPDstandardNth3gyy)*(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + pow(gyzL,2))) + (gxzL*gyyL - 
      gxyL*gyzL)*(2*gxyL*gzzL*JacPDstandardNth1gxy + 
      2*gxzL*gyyL*JacPDstandardNth1gxz - 2*(gxzL*gyzL*JacPDstandardNth1gxy + 
      gxyL*gyzL*JacPDstandardNth1gxz + gxyL*gxzL*JacPDstandardNth1gyz) + 
      gxxL*(-(gzzL*JacPDstandardNth1gyy) + 2*gyzL*JacPDstandardNth1gyz - 
      gyyL*JacPDstandardNth1gzz) + JacPDstandardNth1gzz*pow(gxyL,2) + 
      JacPDstandardNth1gyy*pow(gxzL,2) + JacPDstandardNth1gxx*(-(gyyL*gzzL) + 
      pow(gyzL,2))) - (gxyL*gxzL - 
      gxxL*gyzL)*(2*gxyL*gzzL*JacPDstandardNth2gxy + 
      2*gxzL*gyyL*JacPDstandardNth2gxz - 2*(gxzL*gyzL*JacPDstandardNth2gxy + 
      gxyL*gyzL*JacPDstandardNth2gxz + gxyL*gxzL*JacPDstandardNth2gyz) + 
      gxxL*(-(gzzL*JacPDstandardNth2gyy) + 2*gyzL*JacPDstandardNth2gyz - 
      gyyL*JacPDstandardNth2gzz) + JacPDstandardNth2gzz*pow(gxyL,2) + 
      JacPDstandardNth2gyy*pow(gxzL,2) + JacPDstandardNth2gxx*(-(gyyL*gzzL) + 
      pow(gyzL,2))) + (-(gxxL*gyyL) + 
      pow(gxyL,2))*(2*gxyL*gzzL*JacPDstandardNth3gxy + 
      2*gxzL*gyyL*JacPDstandardNth3gxz - 2*(gxzL*gyzL*JacPDstandardNth3gxy + 
      gxyL*gyzL*JacPDstandardNth3gxz + gxyL*gxzL*JacPDstandardNth3gyz) + 
      gxxL*(-(gzzL*JacPDstandardNth3gyy) + 2*gyzL*JacPDstandardNth3gyz - 
      gyyL*JacPDstandardNth3gzz) + JacPDstandardNth3gzz*pow(gxyL,2) + 
      JacPDstandardNth3gyy*pow(gxzL,2) + JacPDstandardNth3gxx*(-(gyyL*gzzL) + 
      pow(gyzL,2)))))*pow(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + pow(gyzL,2)),-2)) + 
      SS3L*(-(gu13*JacPDstandardNth1alp) - gu23*JacPDstandardNth2alp - 
      gu33*JacPDstandardNth3alp - ((gxzL*gyyL - 
      gxyL*gyzL)*JacPDstandardNth1alp + (-(gxyL*gxzL) + 
      gxxL*gyzL)*JacPDstandardNth2alp + JacPDstandardNth3alp*(-(gxxL*gyyL) + 
      pow(gxyL,2)))*pow(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + pow(gyzL,2)),-1)) + 
      SS2L*(-(gu12*JacPDstandardNth1alp) - gu22*JacPDstandardNth2alp - 
      gu23*JacPDstandardNth3alp - (gzzL*(gxyL*JacPDstandardNth1alp - 
      gxxL*JacPDstandardNth2alp) + gxxL*gyzL*JacPDstandardNth3alp - 
      gxzL*(gyzL*JacPDstandardNth1alp + gxyL*JacPDstandardNth3alp) + 
      JacPDstandardNth2alp*pow(gxzL,2))*pow(-2*gxyL*gxzL*gyzL + 
      gzzL*pow(gxyL,2) + gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + 
      pow(gyzL,2)),-1)) + SS1L*(-(gu11*JacPDstandardNth1alp) - 
      gu12*JacPDstandardNth2alp - gu13*JacPDstandardNth3alp - 
      (gxyL*gzzL*JacPDstandardNth2alp + gxzL*gyyL*JacPDstandardNth3alp - 
      gyzL*(gxzL*JacPDstandardNth2alp + gxyL*JacPDstandardNth3alp) + 
      JacPDstandardNth1alp*(-(gyyL*gzzL) + 
      pow(gyzL,2)))*pow(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + pow(gyzL,2)),-1));
    /* Copy local copies back to grid functions */
    DDrhs[index] = DDrhsL;
    EErhs[index] = EErhsL;
    SS1rhs[index] = SS1rhsL;
    SS2rhs[index] = SS2rhsL;
    SS3rhs[index] = SS3rhsL;
  }
  CCTK_ENDLOOP3(CT_Dust_RHSVal);
}
extern "C" void CT_Dust_RHSVal(CCTK_ARGUMENTS)
{
  #ifdef DECLARE_CCTK_ARGUMENTS_CT_Dust_RHSVal
  DECLARE_CCTK_ARGUMENTS_CHECKED(CT_Dust_RHSVal);
  #else
  DECLARE_CCTK_ARGUMENTS;
  #endif
  DECLARE_CCTK_PARAMETERS;
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering CT_Dust_RHSVal_Body");
  }
  if (cctk_iteration % CT_Dust_RHSVal_calc_every != CT_Dust_RHSVal_calc_offset)
  {
    return;
  }
  
  const char* const groups[] = {
    "ADMBase::curv",
    "ADMBase::lapse",
    "ADMBase::metric",
    "ADMBase::shift",
    "CT_Dust::CT_Drhs",
    "CT_Dust::CT_E",
    "CT_Dust::CT_Erhs",
    "CT_Dust::CT_S",
    "CT_Dust::CT_Srhs",
    "HydroBase::eps",
    "HydroBase::press",
    "HydroBase::rho",
    "HydroBase::vel",
    "HydroBase::w_lorentz"};
  AssertGroupStorage(cctkGH, "CT_Dust_RHSVal", 14, groups);
  
  switch (fdOrder)
  {
    case 2:
    {
      EnsureStencilFits(cctkGH, "CT_Dust_RHSVal", 1, 1, 1);
      break;
    }
    
    case 4:
    {
      EnsureStencilFits(cctkGH, "CT_Dust_RHSVal", 2, 2, 2);
      break;
    }
    
    case 6:
    {
      EnsureStencilFits(cctkGH, "CT_Dust_RHSVal", 3, 3, 3);
      break;
    }
    
    case 8:
    {
      EnsureStencilFits(cctkGH, "CT_Dust_RHSVal", 4, 4, 4);
      break;
    }
    default:
      CCTK_BUILTIN_UNREACHABLE();
  }
  
  LoopOverInterior(cctkGH, CT_Dust_RHSVal_Body);
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Leaving CT_Dust_RHSVal_Body");
  }
}

} // namespace CT_Dust
