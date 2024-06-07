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

extern "C" void CT_Dust_Expansion_SelectBCs(CCTK_ARGUMENTS)
{
  #ifdef DECLARE_CCTK_ARGUMENTS_CT_Dust_Expansion_SelectBCs
  DECLARE_CCTK_ARGUMENTS_CHECKED(CT_Dust_Expansion_SelectBCs);
  #else
  DECLARE_CCTK_ARGUMENTS;
  #endif
  DECLARE_CCTK_PARAMETERS;
  
  if (cctk_iteration % CT_Dust_Expansion_calc_every != CT_Dust_Expansion_calc_offset)
    return;
  CCTK_INT ierr CCTK_ATTRIBUTE_UNUSED = 0;
  ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "CT_Dust::CT_sigmass","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for CT_Dust::CT_sigmass.");
  ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "CT_Dust::CT_sigmats","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for CT_Dust::CT_sigmats.");
  ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "CT_Dust::CT_sigmatt","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for CT_Dust::CT_sigmatt.");
  ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "CT_Dust::CT_Theta","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for CT_Dust::CT_Theta.");
  return;
}

static void CT_Dust_Expansion_Body(const cGH* restrict const cctkGH, const int dir, const int face, const CCTK_REAL normal[3], const CCTK_REAL tangentA[3], const CCTK_REAL tangentB[3], const int imin[3], const int imax[3], const int n_subblock_gfs, CCTK_REAL* restrict const subblock_gfs[])
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
  CCTK_LOOP3(CT_Dust_Expansion,
    i,j,k, imin0,imin1,imin2, imax0,imax1,imax2,
    cctk_ash[0],cctk_ash[1],cctk_ash[2])
  {
    const ptrdiff_t index CCTK_ATTRIBUTE_UNUSED = di*i + dj*j + dk*k;
    /* Assign local copies of grid functions */
    
    CCTK_REAL alpL CCTK_ATTRIBUTE_UNUSED = alp[index];
    CCTK_REAL betaxL CCTK_ATTRIBUTE_UNUSED = betax[index];
    CCTK_REAL betayL CCTK_ATTRIBUTE_UNUSED = betay[index];
    CCTK_REAL betazL CCTK_ATTRIBUTE_UNUSED = betaz[index];
    CCTK_REAL DDL CCTK_ATTRIBUTE_UNUSED = DD[index];
    CCTK_REAL DDrhsL CCTK_ATTRIBUTE_UNUSED = DDrhs[index];
    CCTK_REAL dtalpL CCTK_ATTRIBUTE_UNUSED = dtalp[index];
    CCTK_REAL dtbetaxL CCTK_ATTRIBUTE_UNUSED = dtbetax[index];
    CCTK_REAL dtbetayL CCTK_ATTRIBUTE_UNUSED = dtbetay[index];
    CCTK_REAL dtbetazL CCTK_ATTRIBUTE_UNUSED = dtbetaz[index];
    CCTK_REAL EEL CCTK_ATTRIBUTE_UNUSED = EE[index];
    CCTK_REAL EErhsL CCTK_ATTRIBUTE_UNUSED = EErhs[index];
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
    CCTK_REAL SS1rhsL CCTK_ATTRIBUTE_UNUSED = SS1rhs[index];
    CCTK_REAL SS2rhsL CCTK_ATTRIBUTE_UNUSED = SS2rhs[index];
    CCTK_REAL SS3rhsL CCTK_ATTRIBUTE_UNUSED = SS3rhs[index];
    CCTK_REAL ThetaL CCTK_ATTRIBUTE_UNUSED = Theta[index];
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
    CCTK_REAL JacPDstandardNth1gxx CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1gxy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1gxz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1gyy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1gyz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1gzz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1velx CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1vely CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1velz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1w_lorentz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2alp CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2betax CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2betay CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2betaz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2gxx CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2gxy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2gxz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2gyy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2gyz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2gzz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2velx CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2vely CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2velz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2w_lorentz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3alp CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3betax CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3betay CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3betaz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3gxx CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3gxy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3gxz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3gyy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3gyz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3gzz CCTK_ATTRIBUTE_UNUSED;
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
      
      JacPDstandardNth1gxx = PDstandardNth1gxx;
      
      JacPDstandardNth1gxy = PDstandardNth1gxy;
      
      JacPDstandardNth1gxz = PDstandardNth1gxz;
      
      JacPDstandardNth1gyy = PDstandardNth1gyy;
      
      JacPDstandardNth1gyz = PDstandardNth1gyz;
      
      JacPDstandardNth1gzz = PDstandardNth1gzz;
      
      JacPDstandardNth1velx = PDstandardNth1velx;
      
      JacPDstandardNth1vely = PDstandardNth1vely;
      
      JacPDstandardNth1velz = PDstandardNth1velz;
      
      JacPDstandardNth1w_lorentz = PDstandardNth1w_lorentz;
      
      JacPDstandardNth2alp = PDstandardNth2alp;
      
      JacPDstandardNth2betax = PDstandardNth2betax;
      
      JacPDstandardNth2betay = PDstandardNth2betay;
      
      JacPDstandardNth2betaz = PDstandardNth2betaz;
      
      JacPDstandardNth2gxx = PDstandardNth2gxx;
      
      JacPDstandardNth2gxy = PDstandardNth2gxy;
      
      JacPDstandardNth2gxz = PDstandardNth2gxz;
      
      JacPDstandardNth2gyy = PDstandardNth2gyy;
      
      JacPDstandardNth2gyz = PDstandardNth2gyz;
      
      JacPDstandardNth2gzz = PDstandardNth2gzz;
      
      JacPDstandardNth2velx = PDstandardNth2velx;
      
      JacPDstandardNth2vely = PDstandardNth2vely;
      
      JacPDstandardNth2velz = PDstandardNth2velz;
      
      JacPDstandardNth2w_lorentz = PDstandardNth2w_lorentz;
      
      JacPDstandardNth3alp = PDstandardNth3alp;
      
      JacPDstandardNth3betax = PDstandardNth3betax;
      
      JacPDstandardNth3betay = PDstandardNth3betay;
      
      JacPDstandardNth3betaz = PDstandardNth3betaz;
      
      JacPDstandardNth3gxx = PDstandardNth3gxx;
      
      JacPDstandardNth3gxy = PDstandardNth3gxy;
      
      JacPDstandardNth3gxz = PDstandardNth3gxz;
      
      JacPDstandardNth3gyy = PDstandardNth3gyy;
      
      JacPDstandardNth3gyz = PDstandardNth3gyz;
      
      JacPDstandardNth3gzz = PDstandardNth3gzz;
      
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
    
    CCTK_REAL a2 CCTK_ATTRIBUTE_UNUSED = pow(alpL,2);
    
    CCTK_REAL bl1 CCTK_ATTRIBUTE_UNUSED = betaxL*gxxL + betayL*gxyL + 
      betazL*gxzL;
    
    CCTK_REAL bl2 CCTK_ATTRIBUTE_UNUSED = betaxL*gxyL + betayL*gyyL + 
      betazL*gyzL;
    
    CCTK_REAL bl3 CCTK_ATTRIBUTE_UNUSED = betaxL*gxzL + betayL*gyzL + 
      betazL*gzzL;
    
    CCTK_REAL g4ss11 CCTK_ATTRIBUTE_UNUSED = a2*gu11 - pow(betaxL,2);
    
    CCTK_REAL g4ss12 CCTK_ATTRIBUTE_UNUSED = -(betaxL*betayL) + a2*gu12;
    
    CCTK_REAL g4ss13 CCTK_ATTRIBUTE_UNUSED = -(betaxL*betazL) + a2*gu13;
    
    CCTK_REAL g4ss21 CCTK_ATTRIBUTE_UNUSED = -(betaxL*betayL) + a2*gu12;
    
    CCTK_REAL g4ss22 CCTK_ATTRIBUTE_UNUSED = a2*gu22 - pow(betayL,2);
    
    CCTK_REAL g4ss23 CCTK_ATTRIBUTE_UNUSED = -(betayL*betazL) + a2*gu23;
    
    CCTK_REAL g4ss31 CCTK_ATTRIBUTE_UNUSED = -(betaxL*betazL) + a2*gu13;
    
    CCTK_REAL g4ss32 CCTK_ATTRIBUTE_UNUSED = -(betayL*betazL) + a2*gu23;
    
    CCTK_REAL g4ss33 CCTK_ATTRIBUTE_UNUSED = a2*gu33 - pow(betazL,2);
    
    CCTK_REAL ku11 CCTK_ATTRIBUTE_UNUSED = 2*(kyzL*gu12*gu13 + 
      gu11*(kxyL*gu12 + kxzL*gu13)) + kxxL*pow(gu11,2) + kyyL*pow(gu12,2) + 
      kzzL*pow(gu13,2);
    
    CCTK_REAL ku12 CCTK_ATTRIBUTE_UNUSED = gu12*(kxxL*gu11 + kxyL*gu12 + 
      kxzL*gu13) + (kxyL*gu11 + kyyL*gu12 + kyzL*gu13)*gu22 + (kxzL*gu11 + 
      kyzL*gu12 + kzzL*gu13)*gu23;
    
    CCTK_REAL ku13 CCTK_ATTRIBUTE_UNUSED = gu11*(kxxL*gu13 + kxyL*gu23 + 
      kxzL*gu33) + gu12*(kxyL*gu13 + kyyL*gu23 + kyzL*gu33) + gu13*(kxzL*gu13 
      + kyzL*gu23 + kzzL*gu33);
    
    CCTK_REAL ku21 CCTK_ATTRIBUTE_UNUSED = gu11*(kxxL*gu12 + kxyL*gu22 + 
      kxzL*gu23) + gu12*(kxyL*gu12 + kyyL*gu22 + kyzL*gu23) + gu13*(kxzL*gu12 
      + kyzL*gu22 + kzzL*gu23);
    
    CCTK_REAL ku22 CCTK_ATTRIBUTE_UNUSED = 2*(kyzL*gu22*gu23 + 
      gu12*(kxyL*gu22 + kxzL*gu23)) + kxxL*pow(gu12,2) + kyyL*pow(gu22,2) + 
      kzzL*pow(gu23,2);
    
    CCTK_REAL ku23 CCTK_ATTRIBUTE_UNUSED = gu13*(kxxL*gu12 + kxyL*gu22 + 
      kxzL*gu23) + gu23*(kxyL*gu12 + kyyL*gu22 + kyzL*gu23) + (kxzL*gu12 + 
      kyzL*gu22 + kzzL*gu23)*gu33;
    
    CCTK_REAL ku31 CCTK_ATTRIBUTE_UNUSED = gu11*(kxxL*gu13 + kxyL*gu23 + 
      kxzL*gu33) + gu12*(kxyL*gu13 + kyyL*gu23 + kyzL*gu33) + gu13*(kxzL*gu13 
      + kyzL*gu23 + kzzL*gu33);
    
    CCTK_REAL ku32 CCTK_ATTRIBUTE_UNUSED = gu12*(kxxL*gu13 + kxyL*gu23 + 
      kxzL*gu33) + gu22*(kxyL*gu13 + kyyL*gu23 + kyzL*gu33) + gu23*(kxzL*gu13 
      + kyzL*gu23 + kzzL*gu33);
    
    CCTK_REAL ku33 CCTK_ATTRIBUTE_UNUSED = 2*(kyzL*gu23*gu33 + 
      gu13*(kxyL*gu23 + kxzL*gu33)) + kxxL*pow(gu13,2) + kyyL*pow(gu23,2) + 
      kzzL*pow(gu33,2);
    
    CCTK_REAL dsa1 CCTK_ATTRIBUTE_UNUSED = JacPDstandardNth1alp;
    
    CCTK_REAL dsa2 CCTK_ATTRIBUTE_UNUSED = JacPDstandardNth2alp;
    
    CCTK_REAL dsa3 CCTK_ATTRIBUTE_UNUSED = JacPDstandardNth3alp;
    
    CCTK_REAL dsb11 CCTK_ATTRIBUTE_UNUSED = JacPDstandardNth1betax;
    
    CCTK_REAL dsb12 CCTK_ATTRIBUTE_UNUSED = JacPDstandardNth1betay;
    
    CCTK_REAL dsb13 CCTK_ATTRIBUTE_UNUSED = JacPDstandardNth1betaz;
    
    CCTK_REAL dsb21 CCTK_ATTRIBUTE_UNUSED = JacPDstandardNth2betax;
    
    CCTK_REAL dsb22 CCTK_ATTRIBUTE_UNUSED = JacPDstandardNth2betay;
    
    CCTK_REAL dsb23 CCTK_ATTRIBUTE_UNUSED = JacPDstandardNth2betaz;
    
    CCTK_REAL dsb31 CCTK_ATTRIBUTE_UNUSED = JacPDstandardNth3betax;
    
    CCTK_REAL dsb32 CCTK_ATTRIBUTE_UNUSED = JacPDstandardNth3betay;
    
    CCTK_REAL dsb33 CCTK_ATTRIBUTE_UNUSED = JacPDstandardNth3betaz;
    
    CCTK_REAL dsg111 CCTK_ATTRIBUTE_UNUSED = JacPDstandardNth1gxx;
    
    CCTK_REAL dsg211 CCTK_ATTRIBUTE_UNUSED = JacPDstandardNth2gxx;
    
    CCTK_REAL dsg311 CCTK_ATTRIBUTE_UNUSED = JacPDstandardNth3gxx;
    
    CCTK_REAL dsg112 CCTK_ATTRIBUTE_UNUSED = JacPDstandardNth1gxy;
    
    CCTK_REAL dsg212 CCTK_ATTRIBUTE_UNUSED = JacPDstandardNth2gxy;
    
    CCTK_REAL dsg312 CCTK_ATTRIBUTE_UNUSED = JacPDstandardNth3gxy;
    
    CCTK_REAL dsg113 CCTK_ATTRIBUTE_UNUSED = JacPDstandardNth1gxz;
    
    CCTK_REAL dsg213 CCTK_ATTRIBUTE_UNUSED = JacPDstandardNth2gxz;
    
    CCTK_REAL dsg313 CCTK_ATTRIBUTE_UNUSED = JacPDstandardNth3gxz;
    
    CCTK_REAL dsg121 CCTK_ATTRIBUTE_UNUSED = JacPDstandardNth1gxy;
    
    CCTK_REAL dsg221 CCTK_ATTRIBUTE_UNUSED = JacPDstandardNth2gxy;
    
    CCTK_REAL dsg321 CCTK_ATTRIBUTE_UNUSED = JacPDstandardNth3gxy;
    
    CCTK_REAL dsg122 CCTK_ATTRIBUTE_UNUSED = JacPDstandardNth1gyy;
    
    CCTK_REAL dsg222 CCTK_ATTRIBUTE_UNUSED = JacPDstandardNth2gyy;
    
    CCTK_REAL dsg322 CCTK_ATTRIBUTE_UNUSED = JacPDstandardNth3gyy;
    
    CCTK_REAL dsg123 CCTK_ATTRIBUTE_UNUSED = JacPDstandardNth1gyz;
    
    CCTK_REAL dsg223 CCTK_ATTRIBUTE_UNUSED = JacPDstandardNth2gyz;
    
    CCTK_REAL dsg323 CCTK_ATTRIBUTE_UNUSED = JacPDstandardNth3gyz;
    
    CCTK_REAL dsg131 CCTK_ATTRIBUTE_UNUSED = JacPDstandardNth1gxz;
    
    CCTK_REAL dsg231 CCTK_ATTRIBUTE_UNUSED = JacPDstandardNth2gxz;
    
    CCTK_REAL dsg331 CCTK_ATTRIBUTE_UNUSED = JacPDstandardNth3gxz;
    
    CCTK_REAL dsg132 CCTK_ATTRIBUTE_UNUSED = JacPDstandardNth1gyz;
    
    CCTK_REAL dsg232 CCTK_ATTRIBUTE_UNUSED = JacPDstandardNth2gyz;
    
    CCTK_REAL dsg332 CCTK_ATTRIBUTE_UNUSED = JacPDstandardNth3gyz;
    
    CCTK_REAL dsg133 CCTK_ATTRIBUTE_UNUSED = JacPDstandardNth1gzz;
    
    CCTK_REAL dsg233 CCTK_ATTRIBUTE_UNUSED = JacPDstandardNth2gzz;
    
    CCTK_REAL dsg333 CCTK_ATTRIBUTE_UNUSED = JacPDstandardNth3gzz;
    
    CCTK_REAL dsgu111 CCTK_ATTRIBUTE_UNUSED = -(((-(gyyL*gzzL) + 
      pow(gyzL,2))*(gzzL*(2*gxyL*JacPDstandardNth1gxy - 
      gxxL*JacPDstandardNth1gyy) - 2*(gxyL*gxzL*JacPDstandardNth1gyz + 
      gyzL*(gxzL*JacPDstandardNth1gxy + gxyL*JacPDstandardNth1gxz - 
      gxxL*JacPDstandardNth1gyz)) - gyyL*(gzzL*JacPDstandardNth1gxx - 
      2*gxzL*JacPDstandardNth1gxz + gxxL*JacPDstandardNth1gzz) + 
      JacPDstandardNth1gzz*pow(gxyL,2) + JacPDstandardNth1gyy*pow(gxzL,2) + 
      JacPDstandardNth1gxx*pow(gyzL,2)) + (gzzL*JacPDstandardNth1gyy - 
      2*gyzL*JacPDstandardNth1gyz + 
      gyyL*JacPDstandardNth1gzz)*(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + 
      pow(gyzL,2))))*pow(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + pow(gyzL,2)),-2));
    
    CCTK_REAL dsgu112 CCTK_ATTRIBUTE_UNUSED = ((gxzL*gyzL - 
      gxyL*gzzL)*(gzzL*(2*gxyL*JacPDstandardNth1gxy - 
      gxxL*JacPDstandardNth1gyy) - 2*(gxyL*gxzL*JacPDstandardNth1gyz + 
      gyzL*(gxzL*JacPDstandardNth1gxy + gxyL*JacPDstandardNth1gxz - 
      gxxL*JacPDstandardNth1gyz)) - gyyL*(gzzL*JacPDstandardNth1gxx - 
      2*gxzL*JacPDstandardNth1gxz + gxxL*JacPDstandardNth1gzz) + 
      JacPDstandardNth1gzz*pow(gxyL,2) + JacPDstandardNth1gyy*pow(gxzL,2) + 
      JacPDstandardNth1gxx*pow(gyzL,2)) - (-(gzzL*JacPDstandardNth1gxy) + 
      gyzL*JacPDstandardNth1gxz + gxzL*JacPDstandardNth1gyz - 
      gxyL*JacPDstandardNth1gzz)*(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + 
      pow(gyzL,2))))*pow(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + pow(gyzL,2)),-2);
    
    CCTK_REAL dsgu113 CCTK_ATTRIBUTE_UNUSED = (-((gxzL*gyyL - 
      gxyL*gyzL)*(gzzL*(2*gxyL*JacPDstandardNth1gxy - 
      gxxL*JacPDstandardNth1gyy) - 2*(gxyL*gxzL*JacPDstandardNth1gyz + 
      gyzL*(gxzL*JacPDstandardNth1gxy + gxyL*JacPDstandardNth1gxz - 
      gxxL*JacPDstandardNth1gyz)) - gyyL*(gzzL*JacPDstandardNth1gxx - 
      2*gxzL*JacPDstandardNth1gxz + gxxL*JacPDstandardNth1gzz) + 
      JacPDstandardNth1gzz*pow(gxyL,2) + JacPDstandardNth1gyy*pow(gxzL,2) + 
      JacPDstandardNth1gxx*pow(gyzL,2))) + (-(gyzL*JacPDstandardNth1gxy) + 
      gyyL*JacPDstandardNth1gxz + gxzL*JacPDstandardNth1gyy - 
      gxyL*JacPDstandardNth1gyz)*(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + 
      pow(gyzL,2))))*pow(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + pow(gyzL,2)),-2);
    
    CCTK_REAL dsgu121 CCTK_ATTRIBUTE_UNUSED = ((gxzL*gyzL - 
      gxyL*gzzL)*(gzzL*(2*gxyL*JacPDstandardNth1gxy - 
      gxxL*JacPDstandardNth1gyy) - 2*(gxyL*gxzL*JacPDstandardNth1gyz + 
      gyzL*(gxzL*JacPDstandardNth1gxy + gxyL*JacPDstandardNth1gxz - 
      gxxL*JacPDstandardNth1gyz)) - gyyL*(gzzL*JacPDstandardNth1gxx - 
      2*gxzL*JacPDstandardNth1gxz + gxxL*JacPDstandardNth1gzz) + 
      JacPDstandardNth1gzz*pow(gxyL,2) + JacPDstandardNth1gyy*pow(gxzL,2) + 
      JacPDstandardNth1gxx*pow(gyzL,2)) - (-(gzzL*JacPDstandardNth1gxy) + 
      gyzL*JacPDstandardNth1gxz + gxzL*JacPDstandardNth1gyz - 
      gxyL*JacPDstandardNth1gzz)*(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + 
      pow(gyzL,2))))*pow(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + pow(gyzL,2)),-2);
    
    CCTK_REAL dsgu122 CCTK_ATTRIBUTE_UNUSED = (-((-(gxxL*gzzL) + 
      pow(gxzL,2))*(gzzL*(2*gxyL*JacPDstandardNth1gxy - 
      gxxL*JacPDstandardNth1gyy) - 2*(gxyL*gxzL*JacPDstandardNth1gyz + 
      gyzL*(gxzL*JacPDstandardNth1gxy + gxyL*JacPDstandardNth1gxz - 
      gxxL*JacPDstandardNth1gyz)) - gyyL*(gzzL*JacPDstandardNth1gxx - 
      2*gxzL*JacPDstandardNth1gxz + gxxL*JacPDstandardNth1gzz) + 
      JacPDstandardNth1gzz*pow(gxyL,2) + JacPDstandardNth1gyy*pow(gxzL,2) + 
      JacPDstandardNth1gxx*pow(gyzL,2))) + (-(gzzL*JacPDstandardNth1gxx) + 
      2*gxzL*JacPDstandardNth1gxz - 
      gxxL*JacPDstandardNth1gzz)*(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + 
      pow(gyzL,2))))*pow(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + pow(gyzL,2)),-2);
    
    CCTK_REAL dsgu123 CCTK_ATTRIBUTE_UNUSED = ((gxyL*gxzL - 
      gxxL*gyzL)*(gzzL*(2*gxyL*JacPDstandardNth1gxy - 
      gxxL*JacPDstandardNth1gyy) - 2*(gxyL*gxzL*JacPDstandardNth1gyz + 
      gyzL*(gxzL*JacPDstandardNth1gxy + gxyL*JacPDstandardNth1gxz - 
      gxxL*JacPDstandardNth1gyz)) - gyyL*(gzzL*JacPDstandardNth1gxx - 
      2*gxzL*JacPDstandardNth1gxz + gxxL*JacPDstandardNth1gzz) + 
      JacPDstandardNth1gzz*pow(gxyL,2) + JacPDstandardNth1gyy*pow(gxzL,2) + 
      JacPDstandardNth1gxx*pow(gyzL,2)) - (-(gyzL*JacPDstandardNth1gxx) + 
      gxzL*JacPDstandardNth1gxy + gxyL*JacPDstandardNth1gxz - 
      gxxL*JacPDstandardNth1gyz)*(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + 
      pow(gyzL,2))))*pow(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + pow(gyzL,2)),-2);
    
    CCTK_REAL dsgu131 CCTK_ATTRIBUTE_UNUSED = (-((gxzL*gyyL - 
      gxyL*gyzL)*(gzzL*(2*gxyL*JacPDstandardNth1gxy - 
      gxxL*JacPDstandardNth1gyy) - 2*(gxyL*gxzL*JacPDstandardNth1gyz + 
      gyzL*(gxzL*JacPDstandardNth1gxy + gxyL*JacPDstandardNth1gxz - 
      gxxL*JacPDstandardNth1gyz)) - gyyL*(gzzL*JacPDstandardNth1gxx - 
      2*gxzL*JacPDstandardNth1gxz + gxxL*JacPDstandardNth1gzz) + 
      JacPDstandardNth1gzz*pow(gxyL,2) + JacPDstandardNth1gyy*pow(gxzL,2) + 
      JacPDstandardNth1gxx*pow(gyzL,2))) + (-(gyzL*JacPDstandardNth1gxy) + 
      gyyL*JacPDstandardNth1gxz + gxzL*JacPDstandardNth1gyy - 
      gxyL*JacPDstandardNth1gyz)*(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + 
      pow(gyzL,2))))*pow(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + pow(gyzL,2)),-2);
    
    CCTK_REAL dsgu132 CCTK_ATTRIBUTE_UNUSED = ((gxyL*gxzL - 
      gxxL*gyzL)*(gzzL*(2*gxyL*JacPDstandardNth1gxy - 
      gxxL*JacPDstandardNth1gyy) - 2*(gxyL*gxzL*JacPDstandardNth1gyz + 
      gyzL*(gxzL*JacPDstandardNth1gxy + gxyL*JacPDstandardNth1gxz - 
      gxxL*JacPDstandardNth1gyz)) - gyyL*(gzzL*JacPDstandardNth1gxx - 
      2*gxzL*JacPDstandardNth1gxz + gxxL*JacPDstandardNth1gzz) + 
      JacPDstandardNth1gzz*pow(gxyL,2) + JacPDstandardNth1gyy*pow(gxzL,2) + 
      JacPDstandardNth1gxx*pow(gyzL,2)) - (-(gyzL*JacPDstandardNth1gxx) + 
      gxzL*JacPDstandardNth1gxy + gxyL*JacPDstandardNth1gxz - 
      gxxL*JacPDstandardNth1gyz)*(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + 
      pow(gyzL,2))))*pow(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + pow(gyzL,2)),-2);
    
    CCTK_REAL dsgu133 CCTK_ATTRIBUTE_UNUSED = -(((-(gxxL*gyyL) + 
      pow(gxyL,2))*(gzzL*(2*gxyL*JacPDstandardNth1gxy - 
      gxxL*JacPDstandardNth1gyy) - 2*(gxyL*gxzL*JacPDstandardNth1gyz + 
      gyzL*(gxzL*JacPDstandardNth1gxy + gxyL*JacPDstandardNth1gxz - 
      gxxL*JacPDstandardNth1gyz)) - gyyL*(gzzL*JacPDstandardNth1gxx - 
      2*gxzL*JacPDstandardNth1gxz + gxxL*JacPDstandardNth1gzz) + 
      JacPDstandardNth1gzz*pow(gxyL,2) + JacPDstandardNth1gyy*pow(gxzL,2) + 
      JacPDstandardNth1gxx*pow(gyzL,2)) + (gyyL*JacPDstandardNth1gxx - 
      2*gxyL*JacPDstandardNth1gxy + 
      gxxL*JacPDstandardNth1gyy)*(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + 
      pow(gyzL,2))))*pow(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + pow(gyzL,2)),-2));
    
    CCTK_REAL dsgu211 CCTK_ATTRIBUTE_UNUSED = -(((-(gyyL*gzzL) + 
      pow(gyzL,2))*(gzzL*(2*gxyL*JacPDstandardNth2gxy - 
      gxxL*JacPDstandardNth2gyy) - 2*(gxyL*gxzL*JacPDstandardNth2gyz + 
      gyzL*(gxzL*JacPDstandardNth2gxy + gxyL*JacPDstandardNth2gxz - 
      gxxL*JacPDstandardNth2gyz)) - gyyL*(gzzL*JacPDstandardNth2gxx - 
      2*gxzL*JacPDstandardNth2gxz + gxxL*JacPDstandardNth2gzz) + 
      JacPDstandardNth2gzz*pow(gxyL,2) + JacPDstandardNth2gyy*pow(gxzL,2) + 
      JacPDstandardNth2gxx*pow(gyzL,2)) + (gzzL*JacPDstandardNth2gyy - 
      2*gyzL*JacPDstandardNth2gyz + 
      gyyL*JacPDstandardNth2gzz)*(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + 
      pow(gyzL,2))))*pow(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + pow(gyzL,2)),-2));
    
    CCTK_REAL dsgu212 CCTK_ATTRIBUTE_UNUSED = ((gxzL*gyzL - 
      gxyL*gzzL)*(gzzL*(2*gxyL*JacPDstandardNth2gxy - 
      gxxL*JacPDstandardNth2gyy) - 2*(gxyL*gxzL*JacPDstandardNth2gyz + 
      gyzL*(gxzL*JacPDstandardNth2gxy + gxyL*JacPDstandardNth2gxz - 
      gxxL*JacPDstandardNth2gyz)) - gyyL*(gzzL*JacPDstandardNth2gxx - 
      2*gxzL*JacPDstandardNth2gxz + gxxL*JacPDstandardNth2gzz) + 
      JacPDstandardNth2gzz*pow(gxyL,2) + JacPDstandardNth2gyy*pow(gxzL,2) + 
      JacPDstandardNth2gxx*pow(gyzL,2)) - (-(gzzL*JacPDstandardNth2gxy) + 
      gyzL*JacPDstandardNth2gxz + gxzL*JacPDstandardNth2gyz - 
      gxyL*JacPDstandardNth2gzz)*(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + 
      pow(gyzL,2))))*pow(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + pow(gyzL,2)),-2);
    
    CCTK_REAL dsgu213 CCTK_ATTRIBUTE_UNUSED = (-((gxzL*gyyL - 
      gxyL*gyzL)*(gzzL*(2*gxyL*JacPDstandardNth2gxy - 
      gxxL*JacPDstandardNth2gyy) - 2*(gxyL*gxzL*JacPDstandardNth2gyz + 
      gyzL*(gxzL*JacPDstandardNth2gxy + gxyL*JacPDstandardNth2gxz - 
      gxxL*JacPDstandardNth2gyz)) - gyyL*(gzzL*JacPDstandardNth2gxx - 
      2*gxzL*JacPDstandardNth2gxz + gxxL*JacPDstandardNth2gzz) + 
      JacPDstandardNth2gzz*pow(gxyL,2) + JacPDstandardNth2gyy*pow(gxzL,2) + 
      JacPDstandardNth2gxx*pow(gyzL,2))) + (-(gyzL*JacPDstandardNth2gxy) + 
      gyyL*JacPDstandardNth2gxz + gxzL*JacPDstandardNth2gyy - 
      gxyL*JacPDstandardNth2gyz)*(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + 
      pow(gyzL,2))))*pow(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + pow(gyzL,2)),-2);
    
    CCTK_REAL dsgu221 CCTK_ATTRIBUTE_UNUSED = ((gxzL*gyzL - 
      gxyL*gzzL)*(gzzL*(2*gxyL*JacPDstandardNth2gxy - 
      gxxL*JacPDstandardNth2gyy) - 2*(gxyL*gxzL*JacPDstandardNth2gyz + 
      gyzL*(gxzL*JacPDstandardNth2gxy + gxyL*JacPDstandardNth2gxz - 
      gxxL*JacPDstandardNth2gyz)) - gyyL*(gzzL*JacPDstandardNth2gxx - 
      2*gxzL*JacPDstandardNth2gxz + gxxL*JacPDstandardNth2gzz) + 
      JacPDstandardNth2gzz*pow(gxyL,2) + JacPDstandardNth2gyy*pow(gxzL,2) + 
      JacPDstandardNth2gxx*pow(gyzL,2)) - (-(gzzL*JacPDstandardNth2gxy) + 
      gyzL*JacPDstandardNth2gxz + gxzL*JacPDstandardNth2gyz - 
      gxyL*JacPDstandardNth2gzz)*(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + 
      pow(gyzL,2))))*pow(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + pow(gyzL,2)),-2);
    
    CCTK_REAL dsgu222 CCTK_ATTRIBUTE_UNUSED = (-((-(gxxL*gzzL) + 
      pow(gxzL,2))*(gzzL*(2*gxyL*JacPDstandardNth2gxy - 
      gxxL*JacPDstandardNth2gyy) - 2*(gxyL*gxzL*JacPDstandardNth2gyz + 
      gyzL*(gxzL*JacPDstandardNth2gxy + gxyL*JacPDstandardNth2gxz - 
      gxxL*JacPDstandardNth2gyz)) - gyyL*(gzzL*JacPDstandardNth2gxx - 
      2*gxzL*JacPDstandardNth2gxz + gxxL*JacPDstandardNth2gzz) + 
      JacPDstandardNth2gzz*pow(gxyL,2) + JacPDstandardNth2gyy*pow(gxzL,2) + 
      JacPDstandardNth2gxx*pow(gyzL,2))) + (-(gzzL*JacPDstandardNth2gxx) + 
      2*gxzL*JacPDstandardNth2gxz - 
      gxxL*JacPDstandardNth2gzz)*(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + 
      pow(gyzL,2))))*pow(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + pow(gyzL,2)),-2);
    
    CCTK_REAL dsgu223 CCTK_ATTRIBUTE_UNUSED = ((gxyL*gxzL - 
      gxxL*gyzL)*(gzzL*(2*gxyL*JacPDstandardNth2gxy - 
      gxxL*JacPDstandardNth2gyy) - 2*(gxyL*gxzL*JacPDstandardNth2gyz + 
      gyzL*(gxzL*JacPDstandardNth2gxy + gxyL*JacPDstandardNth2gxz - 
      gxxL*JacPDstandardNth2gyz)) - gyyL*(gzzL*JacPDstandardNth2gxx - 
      2*gxzL*JacPDstandardNth2gxz + gxxL*JacPDstandardNth2gzz) + 
      JacPDstandardNth2gzz*pow(gxyL,2) + JacPDstandardNth2gyy*pow(gxzL,2) + 
      JacPDstandardNth2gxx*pow(gyzL,2)) - (-(gyzL*JacPDstandardNth2gxx) + 
      gxzL*JacPDstandardNth2gxy + gxyL*JacPDstandardNth2gxz - 
      gxxL*JacPDstandardNth2gyz)*(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + 
      pow(gyzL,2))))*pow(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + pow(gyzL,2)),-2);
    
    CCTK_REAL dsgu231 CCTK_ATTRIBUTE_UNUSED = (-((gxzL*gyyL - 
      gxyL*gyzL)*(gzzL*(2*gxyL*JacPDstandardNth2gxy - 
      gxxL*JacPDstandardNth2gyy) - 2*(gxyL*gxzL*JacPDstandardNth2gyz + 
      gyzL*(gxzL*JacPDstandardNth2gxy + gxyL*JacPDstandardNth2gxz - 
      gxxL*JacPDstandardNth2gyz)) - gyyL*(gzzL*JacPDstandardNth2gxx - 
      2*gxzL*JacPDstandardNth2gxz + gxxL*JacPDstandardNth2gzz) + 
      JacPDstandardNth2gzz*pow(gxyL,2) + JacPDstandardNth2gyy*pow(gxzL,2) + 
      JacPDstandardNth2gxx*pow(gyzL,2))) + (-(gyzL*JacPDstandardNth2gxy) + 
      gyyL*JacPDstandardNth2gxz + gxzL*JacPDstandardNth2gyy - 
      gxyL*JacPDstandardNth2gyz)*(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + 
      pow(gyzL,2))))*pow(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + pow(gyzL,2)),-2);
    
    CCTK_REAL dsgu232 CCTK_ATTRIBUTE_UNUSED = ((gxyL*gxzL - 
      gxxL*gyzL)*(gzzL*(2*gxyL*JacPDstandardNth2gxy - 
      gxxL*JacPDstandardNth2gyy) - 2*(gxyL*gxzL*JacPDstandardNth2gyz + 
      gyzL*(gxzL*JacPDstandardNth2gxy + gxyL*JacPDstandardNth2gxz - 
      gxxL*JacPDstandardNth2gyz)) - gyyL*(gzzL*JacPDstandardNth2gxx - 
      2*gxzL*JacPDstandardNth2gxz + gxxL*JacPDstandardNth2gzz) + 
      JacPDstandardNth2gzz*pow(gxyL,2) + JacPDstandardNth2gyy*pow(gxzL,2) + 
      JacPDstandardNth2gxx*pow(gyzL,2)) - (-(gyzL*JacPDstandardNth2gxx) + 
      gxzL*JacPDstandardNth2gxy + gxyL*JacPDstandardNth2gxz - 
      gxxL*JacPDstandardNth2gyz)*(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + 
      pow(gyzL,2))))*pow(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + pow(gyzL,2)),-2);
    
    CCTK_REAL dsgu233 CCTK_ATTRIBUTE_UNUSED = -(((-(gxxL*gyyL) + 
      pow(gxyL,2))*(gzzL*(2*gxyL*JacPDstandardNth2gxy - 
      gxxL*JacPDstandardNth2gyy) - 2*(gxyL*gxzL*JacPDstandardNth2gyz + 
      gyzL*(gxzL*JacPDstandardNth2gxy + gxyL*JacPDstandardNth2gxz - 
      gxxL*JacPDstandardNth2gyz)) - gyyL*(gzzL*JacPDstandardNth2gxx - 
      2*gxzL*JacPDstandardNth2gxz + gxxL*JacPDstandardNth2gzz) + 
      JacPDstandardNth2gzz*pow(gxyL,2) + JacPDstandardNth2gyy*pow(gxzL,2) + 
      JacPDstandardNth2gxx*pow(gyzL,2)) + (gyyL*JacPDstandardNth2gxx - 
      2*gxyL*JacPDstandardNth2gxy + 
      gxxL*JacPDstandardNth2gyy)*(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + 
      pow(gyzL,2))))*pow(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + pow(gyzL,2)),-2));
    
    CCTK_REAL dsgu311 CCTK_ATTRIBUTE_UNUSED = -(((-(gyyL*gzzL) + 
      pow(gyzL,2))*(gzzL*(2*gxyL*JacPDstandardNth3gxy - 
      gxxL*JacPDstandardNth3gyy) - 2*(gxyL*gxzL*JacPDstandardNth3gyz + 
      gyzL*(gxzL*JacPDstandardNth3gxy + gxyL*JacPDstandardNth3gxz - 
      gxxL*JacPDstandardNth3gyz)) - gyyL*(gzzL*JacPDstandardNth3gxx - 
      2*gxzL*JacPDstandardNth3gxz + gxxL*JacPDstandardNth3gzz) + 
      JacPDstandardNth3gzz*pow(gxyL,2) + JacPDstandardNth3gyy*pow(gxzL,2) + 
      JacPDstandardNth3gxx*pow(gyzL,2)) + (gzzL*JacPDstandardNth3gyy - 
      2*gyzL*JacPDstandardNth3gyz + 
      gyyL*JacPDstandardNth3gzz)*(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + 
      pow(gyzL,2))))*pow(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + pow(gyzL,2)),-2));
    
    CCTK_REAL dsgu312 CCTK_ATTRIBUTE_UNUSED = ((gxzL*gyzL - 
      gxyL*gzzL)*(gzzL*(2*gxyL*JacPDstandardNth3gxy - 
      gxxL*JacPDstandardNth3gyy) - 2*(gxyL*gxzL*JacPDstandardNth3gyz + 
      gyzL*(gxzL*JacPDstandardNth3gxy + gxyL*JacPDstandardNth3gxz - 
      gxxL*JacPDstandardNth3gyz)) - gyyL*(gzzL*JacPDstandardNth3gxx - 
      2*gxzL*JacPDstandardNth3gxz + gxxL*JacPDstandardNth3gzz) + 
      JacPDstandardNth3gzz*pow(gxyL,2) + JacPDstandardNth3gyy*pow(gxzL,2) + 
      JacPDstandardNth3gxx*pow(gyzL,2)) - (-(gzzL*JacPDstandardNth3gxy) + 
      gyzL*JacPDstandardNth3gxz + gxzL*JacPDstandardNth3gyz - 
      gxyL*JacPDstandardNth3gzz)*(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + 
      pow(gyzL,2))))*pow(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + pow(gyzL,2)),-2);
    
    CCTK_REAL dsgu313 CCTK_ATTRIBUTE_UNUSED = (-((gxzL*gyyL - 
      gxyL*gyzL)*(gzzL*(2*gxyL*JacPDstandardNth3gxy - 
      gxxL*JacPDstandardNth3gyy) - 2*(gxyL*gxzL*JacPDstandardNth3gyz + 
      gyzL*(gxzL*JacPDstandardNth3gxy + gxyL*JacPDstandardNth3gxz - 
      gxxL*JacPDstandardNth3gyz)) - gyyL*(gzzL*JacPDstandardNth3gxx - 
      2*gxzL*JacPDstandardNth3gxz + gxxL*JacPDstandardNth3gzz) + 
      JacPDstandardNth3gzz*pow(gxyL,2) + JacPDstandardNth3gyy*pow(gxzL,2) + 
      JacPDstandardNth3gxx*pow(gyzL,2))) + (-(gyzL*JacPDstandardNth3gxy) + 
      gyyL*JacPDstandardNth3gxz + gxzL*JacPDstandardNth3gyy - 
      gxyL*JacPDstandardNth3gyz)*(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + 
      pow(gyzL,2))))*pow(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + pow(gyzL,2)),-2);
    
    CCTK_REAL dsgu321 CCTK_ATTRIBUTE_UNUSED = ((gxzL*gyzL - 
      gxyL*gzzL)*(gzzL*(2*gxyL*JacPDstandardNth3gxy - 
      gxxL*JacPDstandardNth3gyy) - 2*(gxyL*gxzL*JacPDstandardNth3gyz + 
      gyzL*(gxzL*JacPDstandardNth3gxy + gxyL*JacPDstandardNth3gxz - 
      gxxL*JacPDstandardNth3gyz)) - gyyL*(gzzL*JacPDstandardNth3gxx - 
      2*gxzL*JacPDstandardNth3gxz + gxxL*JacPDstandardNth3gzz) + 
      JacPDstandardNth3gzz*pow(gxyL,2) + JacPDstandardNth3gyy*pow(gxzL,2) + 
      JacPDstandardNth3gxx*pow(gyzL,2)) - (-(gzzL*JacPDstandardNth3gxy) + 
      gyzL*JacPDstandardNth3gxz + gxzL*JacPDstandardNth3gyz - 
      gxyL*JacPDstandardNth3gzz)*(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + 
      pow(gyzL,2))))*pow(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + pow(gyzL,2)),-2);
    
    CCTK_REAL dsgu322 CCTK_ATTRIBUTE_UNUSED = (-((-(gxxL*gzzL) + 
      pow(gxzL,2))*(gzzL*(2*gxyL*JacPDstandardNth3gxy - 
      gxxL*JacPDstandardNth3gyy) - 2*(gxyL*gxzL*JacPDstandardNth3gyz + 
      gyzL*(gxzL*JacPDstandardNth3gxy + gxyL*JacPDstandardNth3gxz - 
      gxxL*JacPDstandardNth3gyz)) - gyyL*(gzzL*JacPDstandardNth3gxx - 
      2*gxzL*JacPDstandardNth3gxz + gxxL*JacPDstandardNth3gzz) + 
      JacPDstandardNth3gzz*pow(gxyL,2) + JacPDstandardNth3gyy*pow(gxzL,2) + 
      JacPDstandardNth3gxx*pow(gyzL,2))) + (-(gzzL*JacPDstandardNth3gxx) + 
      2*gxzL*JacPDstandardNth3gxz - 
      gxxL*JacPDstandardNth3gzz)*(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + 
      pow(gyzL,2))))*pow(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + pow(gyzL,2)),-2);
    
    CCTK_REAL dsgu323 CCTK_ATTRIBUTE_UNUSED = ((gxyL*gxzL - 
      gxxL*gyzL)*(gzzL*(2*gxyL*JacPDstandardNth3gxy - 
      gxxL*JacPDstandardNth3gyy) - 2*(gxyL*gxzL*JacPDstandardNth3gyz + 
      gyzL*(gxzL*JacPDstandardNth3gxy + gxyL*JacPDstandardNth3gxz - 
      gxxL*JacPDstandardNth3gyz)) - gyyL*(gzzL*JacPDstandardNth3gxx - 
      2*gxzL*JacPDstandardNth3gxz + gxxL*JacPDstandardNth3gzz) + 
      JacPDstandardNth3gzz*pow(gxyL,2) + JacPDstandardNth3gyy*pow(gxzL,2) + 
      JacPDstandardNth3gxx*pow(gyzL,2)) - (-(gyzL*JacPDstandardNth3gxx) + 
      gxzL*JacPDstandardNth3gxy + gxyL*JacPDstandardNth3gxz - 
      gxxL*JacPDstandardNth3gyz)*(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + 
      pow(gyzL,2))))*pow(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + pow(gyzL,2)),-2);
    
    CCTK_REAL dsgu331 CCTK_ATTRIBUTE_UNUSED = (-((gxzL*gyyL - 
      gxyL*gyzL)*(gzzL*(2*gxyL*JacPDstandardNth3gxy - 
      gxxL*JacPDstandardNth3gyy) - 2*(gxyL*gxzL*JacPDstandardNth3gyz + 
      gyzL*(gxzL*JacPDstandardNth3gxy + gxyL*JacPDstandardNth3gxz - 
      gxxL*JacPDstandardNth3gyz)) - gyyL*(gzzL*JacPDstandardNth3gxx - 
      2*gxzL*JacPDstandardNth3gxz + gxxL*JacPDstandardNth3gzz) + 
      JacPDstandardNth3gzz*pow(gxyL,2) + JacPDstandardNth3gyy*pow(gxzL,2) + 
      JacPDstandardNth3gxx*pow(gyzL,2))) + (-(gyzL*JacPDstandardNth3gxy) + 
      gyyL*JacPDstandardNth3gxz + gxzL*JacPDstandardNth3gyy - 
      gxyL*JacPDstandardNth3gyz)*(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + 
      pow(gyzL,2))))*pow(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + pow(gyzL,2)),-2);
    
    CCTK_REAL dsgu332 CCTK_ATTRIBUTE_UNUSED = ((gxyL*gxzL - 
      gxxL*gyzL)*(gzzL*(2*gxyL*JacPDstandardNth3gxy - 
      gxxL*JacPDstandardNth3gyy) - 2*(gxyL*gxzL*JacPDstandardNth3gyz + 
      gyzL*(gxzL*JacPDstandardNth3gxy + gxyL*JacPDstandardNth3gxz - 
      gxxL*JacPDstandardNth3gyz)) - gyyL*(gzzL*JacPDstandardNth3gxx - 
      2*gxzL*JacPDstandardNth3gxz + gxxL*JacPDstandardNth3gzz) + 
      JacPDstandardNth3gzz*pow(gxyL,2) + JacPDstandardNth3gyy*pow(gxzL,2) + 
      JacPDstandardNth3gxx*pow(gyzL,2)) - (-(gyzL*JacPDstandardNth3gxx) + 
      gxzL*JacPDstandardNth3gxy + gxyL*JacPDstandardNth3gxz - 
      gxxL*JacPDstandardNth3gyz)*(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + 
      pow(gyzL,2))))*pow(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + pow(gyzL,2)),-2);
    
    CCTK_REAL dsgu333 CCTK_ATTRIBUTE_UNUSED = -(((-(gxxL*gyyL) + 
      pow(gxyL,2))*(gzzL*(2*gxyL*JacPDstandardNth3gxy - 
      gxxL*JacPDstandardNth3gyy) - 2*(gxyL*gxzL*JacPDstandardNth3gyz + 
      gyzL*(gxzL*JacPDstandardNth3gxy + gxyL*JacPDstandardNth3gxz - 
      gxxL*JacPDstandardNth3gyz)) - gyyL*(gzzL*JacPDstandardNth3gxx - 
      2*gxzL*JacPDstandardNth3gxz + gxxL*JacPDstandardNth3gzz) + 
      JacPDstandardNth3gzz*pow(gxyL,2) + JacPDstandardNth3gyy*pow(gxzL,2) + 
      JacPDstandardNth3gxx*pow(gyzL,2)) + (gyyL*JacPDstandardNth3gxx - 
      2*gxyL*JacPDstandardNth3gxy + 
      gxxL*JacPDstandardNth3gyy)*(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + 
      pow(gyzL,2))))*pow(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + pow(gyzL,2)),-2));
    
    CCTK_REAL dsbl11 CCTK_ATTRIBUTE_UNUSED = gxxL*dsb11 + gxyL*dsb12 + 
      gxzL*dsb13 + betaxL*dsg111 + betayL*dsg112 + betazL*dsg113;
    
    CCTK_REAL dsbl12 CCTK_ATTRIBUTE_UNUSED = gxyL*dsb11 + gyyL*dsb12 + 
      gyzL*dsb13 + betaxL*dsg121 + betayL*dsg122 + betazL*dsg123;
    
    CCTK_REAL dsbl13 CCTK_ATTRIBUTE_UNUSED = gxzL*dsb11 + gyzL*dsb12 + 
      gzzL*dsb13 + betaxL*dsg131 + betayL*dsg132 + betazL*dsg133;
    
    CCTK_REAL dsbl21 CCTK_ATTRIBUTE_UNUSED = gxxL*dsb21 + gxyL*dsb22 + 
      gxzL*dsb23 + betaxL*dsg211 + betayL*dsg212 + betazL*dsg213;
    
    CCTK_REAL dsbl22 CCTK_ATTRIBUTE_UNUSED = gxyL*dsb21 + gyyL*dsb22 + 
      gyzL*dsb23 + betaxL*dsg221 + betayL*dsg222 + betazL*dsg223;
    
    CCTK_REAL dsbl23 CCTK_ATTRIBUTE_UNUSED = gxzL*dsb21 + gyzL*dsb22 + 
      gzzL*dsb23 + betaxL*dsg231 + betayL*dsg232 + betazL*dsg233;
    
    CCTK_REAL dsbl31 CCTK_ATTRIBUTE_UNUSED = gxxL*dsb31 + gxyL*dsb32 + 
      gxzL*dsb33 + betaxL*dsg311 + betayL*dsg312 + betazL*dsg313;
    
    CCTK_REAL dsbl32 CCTK_ATTRIBUTE_UNUSED = gxyL*dsb31 + gyyL*dsb32 + 
      gyzL*dsb33 + betaxL*dsg321 + betayL*dsg322 + betazL*dsg323;
    
    CCTK_REAL dsbl33 CCTK_ATTRIBUTE_UNUSED = gxzL*dsb31 + gyzL*dsb32 + 
      gzzL*dsb33 + betaxL*dsg331 + betayL*dsg332 + betazL*dsg333;
    
    CCTK_REAL Liebg11 CCTK_ATTRIBUTE_UNUSED = 2*(gxxL*dsb11 + gxyL*dsb12 + 
      gxzL*dsb13) + betaxL*dsg111 + betayL*dsg211 + betazL*dsg311;
    
    CCTK_REAL Liebg12 CCTK_ATTRIBUTE_UNUSED = gyyL*dsb12 + gyzL*dsb13 + 
      gxxL*dsb21 + gxyL*(dsb11 + dsb22) + gxzL*dsb23 + betaxL*dsg112 + 
      betayL*dsg212 + betazL*dsg312;
    
    CCTK_REAL Liebg13 CCTK_ATTRIBUTE_UNUSED = gyzL*dsb12 + gzzL*dsb13 + 
      gxxL*dsb31 + gxyL*dsb32 + gxzL*(dsb11 + dsb33) + betaxL*dsg113 + 
      betayL*dsg213 + betazL*dsg313;
    
    CCTK_REAL Liebg21 CCTK_ATTRIBUTE_UNUSED = gyyL*dsb12 + gyzL*dsb13 + 
      gxxL*dsb21 + gxyL*(dsb11 + dsb22) + gxzL*dsb23 + betaxL*dsg121 + 
      betayL*dsg221 + betazL*dsg321;
    
    CCTK_REAL Liebg22 CCTK_ATTRIBUTE_UNUSED = 2*(gxyL*dsb21 + gyyL*dsb22 + 
      gyzL*dsb23) + betaxL*dsg122 + betayL*dsg222 + betazL*dsg322;
    
    CCTK_REAL Liebg23 CCTK_ATTRIBUTE_UNUSED = gxzL*dsb21 + gzzL*dsb23 + 
      gxyL*dsb31 + gyyL*dsb32 + gyzL*(dsb22 + dsb33) + betaxL*dsg123 + 
      betayL*dsg223 + betazL*dsg323;
    
    CCTK_REAL Liebg31 CCTK_ATTRIBUTE_UNUSED = gyzL*dsb12 + gzzL*dsb13 + 
      gxxL*dsb31 + gxyL*dsb32 + gxzL*(dsb11 + dsb33) + betaxL*dsg131 + 
      betayL*dsg231 + betazL*dsg331;
    
    CCTK_REAL Liebg32 CCTK_ATTRIBUTE_UNUSED = gxzL*dsb21 + gzzL*dsb23 + 
      gxyL*dsb31 + gyyL*dsb32 + gyzL*(dsb22 + dsb33) + betaxL*dsg132 + 
      betayL*dsg232 + betazL*dsg332;
    
    CCTK_REAL Liebg33 CCTK_ATTRIBUTE_UNUSED = 2*(gxzL*dsb31 + gyzL*dsb32 + 
      gzzL*dsb33) + betaxL*dsg133 + betayL*dsg233 + betazL*dsg333;
    
    CCTK_REAL Liebgu11 CCTK_ATTRIBUTE_UNUSED = betaxL*dsgu111 + 
      betayL*dsgu211 + betazL*dsgu311;
    
    CCTK_REAL Liebgu12 CCTK_ATTRIBUTE_UNUSED = betaxL*dsgu112 + 
      betayL*dsgu212 + betazL*dsgu312 + dsb12*gu11 + (-dsb11 + dsb22)*gu12 + 
      dsb32*gu13 - dsb21*gu22 - dsb31*gu23;
    
    CCTK_REAL Liebgu13 CCTK_ATTRIBUTE_UNUSED = betaxL*dsgu113 + 
      betayL*dsgu213 + betazL*dsgu313 + dsb13*gu11 + dsb23*gu12 + (-dsb11 + 
      dsb33)*gu13 - dsb21*gu23 - dsb31*gu33;
    
    CCTK_REAL Liebgu21 CCTK_ATTRIBUTE_UNUSED = betaxL*dsgu121 + 
      betayL*dsgu221 + betazL*dsgu321 - dsb12*gu11 + (dsb11 - dsb22)*gu12 - 
      dsb32*gu13 + dsb21*gu22 + dsb31*gu23;
    
    CCTK_REAL Liebgu22 CCTK_ATTRIBUTE_UNUSED = betaxL*dsgu122 + 
      betayL*dsgu222 + betazL*dsgu322;
    
    CCTK_REAL Liebgu23 CCTK_ATTRIBUTE_UNUSED = betaxL*dsgu123 + 
      betayL*dsgu223 + betazL*dsgu323 + dsb13*gu12 - dsb12*gu13 + dsb23*gu22 
      + (-dsb22 + dsb33)*gu23 - dsb32*gu33;
    
    CCTK_REAL Liebgu31 CCTK_ATTRIBUTE_UNUSED = betaxL*dsgu131 + 
      betayL*dsgu231 + betazL*dsgu331 - dsb13*gu11 - dsb23*gu12 + (dsb11 - 
      dsb33)*gu13 + dsb21*gu23 + dsb31*gu33;
    
    CCTK_REAL Liebgu32 CCTK_ATTRIBUTE_UNUSED = betaxL*dsgu132 + 
      betayL*dsgu232 + betazL*dsgu332 - dsb13*gu12 + dsb12*gu13 - dsb23*gu22 
      + (dsb22 - dsb33)*gu23 + dsb32*gu33;
    
    CCTK_REAL Liebgu33 CCTK_ATTRIBUTE_UNUSED = betaxL*dsgu133 + 
      betayL*dsgu233 + betazL*dsgu333;
    
    CCTK_REAL dtg11 CCTK_ATTRIBUTE_UNUSED = -2*alpL*kxxL + Liebg11;
    
    CCTK_REAL dtg12 CCTK_ATTRIBUTE_UNUSED = -2*alpL*kxyL + Liebg12;
    
    CCTK_REAL dtg13 CCTK_ATTRIBUTE_UNUSED = -2*alpL*kxzL + Liebg13;
    
    CCTK_REAL dtg21 CCTK_ATTRIBUTE_UNUSED = -2*alpL*kxyL + Liebg21;
    
    CCTK_REAL dtg22 CCTK_ATTRIBUTE_UNUSED = -2*alpL*kyyL + Liebg22;
    
    CCTK_REAL dtg23 CCTK_ATTRIBUTE_UNUSED = -2*alpL*kyzL + Liebg23;
    
    CCTK_REAL dtg31 CCTK_ATTRIBUTE_UNUSED = -2*alpL*kxzL + Liebg31;
    
    CCTK_REAL dtg32 CCTK_ATTRIBUTE_UNUSED = -2*alpL*kyzL + Liebg32;
    
    CCTK_REAL dtg33 CCTK_ATTRIBUTE_UNUSED = -2*alpL*kzzL + Liebg33;
    
    CCTK_REAL dtgu11 CCTK_ATTRIBUTE_UNUSED = -2*alpL*ku11 + Liebgu11;
    
    CCTK_REAL dtgu12 CCTK_ATTRIBUTE_UNUSED = -2*alpL*ku12 + Liebgu12;
    
    CCTK_REAL dtgu13 CCTK_ATTRIBUTE_UNUSED = -2*alpL*ku13 + Liebgu13;
    
    CCTK_REAL dtgu21 CCTK_ATTRIBUTE_UNUSED = -2*alpL*ku21 + Liebgu21;
    
    CCTK_REAL dtgu22 CCTK_ATTRIBUTE_UNUSED = -2*alpL*ku22 + Liebgu22;
    
    CCTK_REAL dtgu23 CCTK_ATTRIBUTE_UNUSED = -2*alpL*ku23 + Liebgu23;
    
    CCTK_REAL dtgu31 CCTK_ATTRIBUTE_UNUSED = -2*alpL*ku31 + Liebgu31;
    
    CCTK_REAL dtgu32 CCTK_ATTRIBUTE_UNUSED = -2*alpL*ku32 + Liebgu32;
    
    CCTK_REAL dtgu33 CCTK_ATTRIBUTE_UNUSED = -2*alpL*ku33 + Liebgu33;
    
    CCTK_REAL dtbl1 CCTK_ATTRIBUTE_UNUSED = dtbetaxL*gxxL + dtbetayL*gxyL 
      + dtbetazL*gxzL + betaxL*dtg11 + betayL*dtg12 + betazL*dtg13;
    
    CCTK_REAL dtbl2 CCTK_ATTRIBUTE_UNUSED = dtbetaxL*gxyL + dtbetayL*gyyL 
      + dtbetazL*gyzL + betaxL*dtg21 + betayL*dtg22 + betazL*dtg23;
    
    CCTK_REAL dtbl3 CCTK_ATTRIBUTE_UNUSED = dtbetaxL*gxzL + dtbetayL*gyzL 
      + dtbetazL*gzzL + betaxL*dtg31 + betayL*dtg32 + betazL*dtg33;
    
    CCTK_REAL Gdttt CCTK_ATTRIBUTE_UNUSED = 0.5*(-2*alpL*dtalpL + 
      dtbetaxL*bl1 + dtbetayL*bl2 + dtbetazL*bl3 + betaxL*dtbl1 + 
      betayL*dtbl2 + betazL*dtbl3);
    
    CCTK_REAL Gdtts1 CCTK_ATTRIBUTE_UNUSED = 0.5*(-2*alpL*dsa1 + bl1*dsb11 
      + bl2*dsb12 + bl3*dsb13 + betaxL*dsbl11 + betayL*dsbl12 + 
      betazL*dsbl13);
    
    CCTK_REAL Gdtts2 CCTK_ATTRIBUTE_UNUSED = 0.5*(-2*alpL*dsa2 + bl1*dsb21 
      + bl2*dsb22 + bl3*dsb23 + betaxL*dsbl21 + betayL*dsbl22 + 
      betazL*dsbl23);
    
    CCTK_REAL Gdtts3 CCTK_ATTRIBUTE_UNUSED = 0.5*(-2*alpL*dsa3 + bl1*dsb31 
      + bl2*dsb32 + bl3*dsb33 + betaxL*dsbl31 + betayL*dsbl32 + 
      betazL*dsbl33);
    
    CCTK_REAL Gdtss11 CCTK_ATTRIBUTE_UNUSED = dsbl11 - 0.5*dtg11;
    
    CCTK_REAL Gdtss12 CCTK_ATTRIBUTE_UNUSED = 0.5*(dsbl12 + dsbl21 - 
      dtg12);
    
    CCTK_REAL Gdtss13 CCTK_ATTRIBUTE_UNUSED = 0.5*(dsbl13 + dsbl31 - 
      dtg13);
    
    CCTK_REAL Gdtss21 CCTK_ATTRIBUTE_UNUSED = 0.5*(dsbl12 + dsbl21 - 
      dtg21);
    
    CCTK_REAL Gdtss22 CCTK_ATTRIBUTE_UNUSED = dsbl22 - 0.5*dtg22;
    
    CCTK_REAL Gdtss23 CCTK_ATTRIBUTE_UNUSED = 0.5*(dsbl23 + dsbl32 - 
      dtg23);
    
    CCTK_REAL Gdtss31 CCTK_ATTRIBUTE_UNUSED = 0.5*(dsbl13 + dsbl31 - 
      dtg31);
    
    CCTK_REAL Gdtss32 CCTK_ATTRIBUTE_UNUSED = 0.5*(dsbl23 + dsbl32 - 
      dtg32);
    
    CCTK_REAL Gdtss33 CCTK_ATTRIBUTE_UNUSED = dsbl33 - 0.5*dtg33;
    
    CCTK_REAL Gdstt1 CCTK_ATTRIBUTE_UNUSED = dtbl1 - Gdtts1;
    
    CCTK_REAL Gdstt2 CCTK_ATTRIBUTE_UNUSED = dtbl2 - Gdtts2;
    
    CCTK_REAL Gdstt3 CCTK_ATTRIBUTE_UNUSED = dtbl3 - Gdtts3;
    
    CCTK_REAL Gdsts11 CCTK_ATTRIBUTE_UNUSED = dsbl11 - Gdtss11;
    
    CCTK_REAL Gdsts12 CCTK_ATTRIBUTE_UNUSED = dsbl21 - Gdtss12;
    
    CCTK_REAL Gdsts13 CCTK_ATTRIBUTE_UNUSED = dsbl31 - Gdtss13;
    
    CCTK_REAL Gdsts21 CCTK_ATTRIBUTE_UNUSED = dsbl12 - Gdtss21;
    
    CCTK_REAL Gdsts22 CCTK_ATTRIBUTE_UNUSED = dsbl22 - Gdtss22;
    
    CCTK_REAL Gdsts23 CCTK_ATTRIBUTE_UNUSED = dsbl32 - Gdtss23;
    
    CCTK_REAL Gdsts31 CCTK_ATTRIBUTE_UNUSED = dsbl13 - Gdtss31;
    
    CCTK_REAL Gdsts32 CCTK_ATTRIBUTE_UNUSED = dsbl23 - Gdtss32;
    
    CCTK_REAL Gdsts33 CCTK_ATTRIBUTE_UNUSED = dsbl33 - Gdtss33;
    
    CCTK_REAL Gdsss111 CCTK_ATTRIBUTE_UNUSED = 0.5*dsg111;
    
    CCTK_REAL Gdsss112 CCTK_ATTRIBUTE_UNUSED = 0.5*dsg211;
    
    CCTK_REAL Gdsss113 CCTK_ATTRIBUTE_UNUSED = 0.5*dsg311;
    
    CCTK_REAL Gdsss121 CCTK_ATTRIBUTE_UNUSED = 0.5*dsg211;
    
    CCTK_REAL Gdsss122 CCTK_ATTRIBUTE_UNUSED = 0.5*(-dsg122 + dsg212 + 
      dsg221);
    
    CCTK_REAL Gdsss123 CCTK_ATTRIBUTE_UNUSED = 0.5*(-dsg123 + dsg213 + 
      dsg321);
    
    CCTK_REAL Gdsss131 CCTK_ATTRIBUTE_UNUSED = 0.5*dsg311;
    
    CCTK_REAL Gdsss132 CCTK_ATTRIBUTE_UNUSED = 0.5*(-dsg132 + dsg231 + 
      dsg312);
    
    CCTK_REAL Gdsss133 CCTK_ATTRIBUTE_UNUSED = 0.5*(-dsg133 + dsg313 + 
      dsg331);
    
    CCTK_REAL Gdsss211 CCTK_ATTRIBUTE_UNUSED = 0.5*(dsg112 + dsg121 - 
      dsg211);
    
    CCTK_REAL Gdsss212 CCTK_ATTRIBUTE_UNUSED = 0.5*dsg122;
    
    CCTK_REAL Gdsss213 CCTK_ATTRIBUTE_UNUSED = 0.5*(dsg123 - dsg213 + 
      dsg312);
    
    CCTK_REAL Gdsss221 CCTK_ATTRIBUTE_UNUSED = 0.5*dsg122;
    
    CCTK_REAL Gdsss222 CCTK_ATTRIBUTE_UNUSED = 0.5*dsg222;
    
    CCTK_REAL Gdsss223 CCTK_ATTRIBUTE_UNUSED = 0.5*dsg322;
    
    CCTK_REAL Gdsss231 CCTK_ATTRIBUTE_UNUSED = 0.5*(dsg132 - dsg231 + 
      dsg321);
    
    CCTK_REAL Gdsss232 CCTK_ATTRIBUTE_UNUSED = 0.5*dsg322;
    
    CCTK_REAL Gdsss233 CCTK_ATTRIBUTE_UNUSED = 0.5*(-dsg233 + dsg323 + 
      dsg332);
    
    CCTK_REAL Gdsss311 CCTK_ATTRIBUTE_UNUSED = 0.5*(dsg113 + dsg131 - 
      dsg311);
    
    CCTK_REAL Gdsss312 CCTK_ATTRIBUTE_UNUSED = 0.5*(dsg132 + dsg213 - 
      dsg312);
    
    CCTK_REAL Gdsss313 CCTK_ATTRIBUTE_UNUSED = 0.5*dsg133;
    
    CCTK_REAL Gdsss321 CCTK_ATTRIBUTE_UNUSED = 0.5*(dsg123 + dsg231 - 
      dsg321);
    
    CCTK_REAL Gdsss322 CCTK_ATTRIBUTE_UNUSED = 0.5*(dsg223 + dsg232 - 
      dsg322);
    
    CCTK_REAL Gdsss323 CCTK_ATTRIBUTE_UNUSED = 0.5*dsg233;
    
    CCTK_REAL Gdsss331 CCTK_ATTRIBUTE_UNUSED = 0.5*dsg133;
    
    CCTK_REAL Gdsss332 CCTK_ATTRIBUTE_UNUSED = 0.5*dsg233;
    
    CCTK_REAL Gdsss333 CCTK_ATTRIBUTE_UNUSED = 0.5*dsg333;
    
    CCTK_REAL Gttt CCTK_ATTRIBUTE_UNUSED = (betaxL*Gdstt1 + betayL*Gdstt2 
      + betazL*Gdstt3 - Gdttt)*pow(a2,-1);
    
    CCTK_REAL Gtts1 CCTK_ATTRIBUTE_UNUSED = (betaxL*Gdsts11 + 
      betayL*Gdsts21 + betazL*Gdsts31 - Gdtts1)*pow(a2,-1);
    
    CCTK_REAL Gtts2 CCTK_ATTRIBUTE_UNUSED = (betaxL*Gdsts12 + 
      betayL*Gdsts22 + betazL*Gdsts32 - Gdtts2)*pow(a2,-1);
    
    CCTK_REAL Gtts3 CCTK_ATTRIBUTE_UNUSED = (betaxL*Gdsts13 + 
      betayL*Gdsts23 + betazL*Gdsts33 - Gdtts3)*pow(a2,-1);
    
    CCTK_REAL Gtss11 CCTK_ATTRIBUTE_UNUSED = (betaxL*Gdsss111 + 
      betayL*Gdsss211 + betazL*Gdsss311 - Gdtss11)*pow(a2,-1);
    
    CCTK_REAL Gtss12 CCTK_ATTRIBUTE_UNUSED = (betaxL*Gdsss112 + 
      betayL*Gdsss212 + betazL*Gdsss312 - Gdtss12)*pow(a2,-1);
    
    CCTK_REAL Gtss13 CCTK_ATTRIBUTE_UNUSED = (betaxL*Gdsss113 + 
      betayL*Gdsss213 + betazL*Gdsss313 - Gdtss13)*pow(a2,-1);
    
    CCTK_REAL Gtss21 CCTK_ATTRIBUTE_UNUSED = (betaxL*Gdsss121 + 
      betayL*Gdsss221 + betazL*Gdsss321 - Gdtss21)*pow(a2,-1);
    
    CCTK_REAL Gtss22 CCTK_ATTRIBUTE_UNUSED = (betaxL*Gdsss122 + 
      betayL*Gdsss222 + betazL*Gdsss322 - Gdtss22)*pow(a2,-1);
    
    CCTK_REAL Gtss23 CCTK_ATTRIBUTE_UNUSED = (betaxL*Gdsss123 + 
      betayL*Gdsss223 + betazL*Gdsss323 - Gdtss23)*pow(a2,-1);
    
    CCTK_REAL Gtss31 CCTK_ATTRIBUTE_UNUSED = (betaxL*Gdsss131 + 
      betayL*Gdsss231 + betazL*Gdsss331 - Gdtss31)*pow(a2,-1);
    
    CCTK_REAL Gtss32 CCTK_ATTRIBUTE_UNUSED = (betaxL*Gdsss132 + 
      betayL*Gdsss232 + betazL*Gdsss332 - Gdtss32)*pow(a2,-1);
    
    CCTK_REAL Gtss33 CCTK_ATTRIBUTE_UNUSED = (betaxL*Gdsss133 + 
      betayL*Gdsss233 + betazL*Gdsss333 - Gdtss33)*pow(a2,-1);
    
    CCTK_REAL Gstt1 CCTK_ATTRIBUTE_UNUSED = (g4ss11*Gdstt1 + g4ss12*Gdstt2 
      + g4ss13*Gdstt3 + betaxL*Gdttt)*pow(a2,-1);
    
    CCTK_REAL Gstt2 CCTK_ATTRIBUTE_UNUSED = (g4ss21*Gdstt1 + g4ss22*Gdstt2 
      + g4ss23*Gdstt3 + betayL*Gdttt)*pow(a2,-1);
    
    CCTK_REAL Gstt3 CCTK_ATTRIBUTE_UNUSED = (g4ss31*Gdstt1 + g4ss32*Gdstt2 
      + g4ss33*Gdstt3 + betazL*Gdttt)*pow(a2,-1);
    
    CCTK_REAL Gsts11 CCTK_ATTRIBUTE_UNUSED = (g4ss11*Gdsts11 + 
      g4ss12*Gdsts21 + g4ss13*Gdsts31 + betaxL*Gdtts1)*pow(a2,-1);
    
    CCTK_REAL Gsts21 CCTK_ATTRIBUTE_UNUSED = (g4ss21*Gdsts11 + 
      g4ss22*Gdsts21 + g4ss23*Gdsts31 + betayL*Gdtts1)*pow(a2,-1);
    
    CCTK_REAL Gsts31 CCTK_ATTRIBUTE_UNUSED = (g4ss31*Gdsts11 + 
      g4ss32*Gdsts21 + g4ss33*Gdsts31 + betazL*Gdtts1)*pow(a2,-1);
    
    CCTK_REAL Gsts12 CCTK_ATTRIBUTE_UNUSED = (g4ss11*Gdsts12 + 
      g4ss12*Gdsts22 + g4ss13*Gdsts32 + betaxL*Gdtts2)*pow(a2,-1);
    
    CCTK_REAL Gsts22 CCTK_ATTRIBUTE_UNUSED = (g4ss21*Gdsts12 + 
      g4ss22*Gdsts22 + g4ss23*Gdsts32 + betayL*Gdtts2)*pow(a2,-1);
    
    CCTK_REAL Gsts32 CCTK_ATTRIBUTE_UNUSED = (g4ss31*Gdsts12 + 
      g4ss32*Gdsts22 + g4ss33*Gdsts32 + betazL*Gdtts2)*pow(a2,-1);
    
    CCTK_REAL Gsts13 CCTK_ATTRIBUTE_UNUSED = (g4ss11*Gdsts13 + 
      g4ss12*Gdsts23 + g4ss13*Gdsts33 + betaxL*Gdtts3)*pow(a2,-1);
    
    CCTK_REAL Gsts23 CCTK_ATTRIBUTE_UNUSED = (g4ss21*Gdsts13 + 
      g4ss22*Gdsts23 + g4ss23*Gdsts33 + betayL*Gdtts3)*pow(a2,-1);
    
    CCTK_REAL Gsts33 CCTK_ATTRIBUTE_UNUSED = (g4ss31*Gdsts13 + 
      g4ss32*Gdsts23 + g4ss33*Gdsts33 + betazL*Gdtts3)*pow(a2,-1);
    
    CCTK_REAL Gsss111 CCTK_ATTRIBUTE_UNUSED = (g4ss11*Gdsss111 + 
      g4ss12*Gdsss211 + g4ss13*Gdsss311 + betaxL*Gdtss11)*pow(a2,-1);
    
    CCTK_REAL Gsss211 CCTK_ATTRIBUTE_UNUSED = (g4ss21*Gdsss111 + 
      g4ss22*Gdsss211 + g4ss23*Gdsss311 + betayL*Gdtss11)*pow(a2,-1);
    
    CCTK_REAL Gsss311 CCTK_ATTRIBUTE_UNUSED = (g4ss31*Gdsss111 + 
      g4ss32*Gdsss211 + g4ss33*Gdsss311 + betazL*Gdtss11)*pow(a2,-1);
    
    CCTK_REAL Gsss112 CCTK_ATTRIBUTE_UNUSED = (g4ss11*Gdsss112 + 
      g4ss12*Gdsss212 + g4ss13*Gdsss312 + betaxL*Gdtss12)*pow(a2,-1);
    
    CCTK_REAL Gsss212 CCTK_ATTRIBUTE_UNUSED = (g4ss21*Gdsss112 + 
      g4ss22*Gdsss212 + g4ss23*Gdsss312 + betayL*Gdtss12)*pow(a2,-1);
    
    CCTK_REAL Gsss312 CCTK_ATTRIBUTE_UNUSED = (g4ss31*Gdsss112 + 
      g4ss32*Gdsss212 + g4ss33*Gdsss312 + betazL*Gdtss12)*pow(a2,-1);
    
    CCTK_REAL Gsss113 CCTK_ATTRIBUTE_UNUSED = (g4ss11*Gdsss113 + 
      g4ss12*Gdsss213 + g4ss13*Gdsss313 + betaxL*Gdtss13)*pow(a2,-1);
    
    CCTK_REAL Gsss213 CCTK_ATTRIBUTE_UNUSED = (g4ss21*Gdsss113 + 
      g4ss22*Gdsss213 + g4ss23*Gdsss313 + betayL*Gdtss13)*pow(a2,-1);
    
    CCTK_REAL Gsss313 CCTK_ATTRIBUTE_UNUSED = (g4ss31*Gdsss113 + 
      g4ss32*Gdsss213 + g4ss33*Gdsss313 + betazL*Gdtss13)*pow(a2,-1);
    
    CCTK_REAL Gsss121 CCTK_ATTRIBUTE_UNUSED = (g4ss11*Gdsss121 + 
      g4ss12*Gdsss221 + g4ss13*Gdsss321 + betaxL*Gdtss21)*pow(a2,-1);
    
    CCTK_REAL Gsss221 CCTK_ATTRIBUTE_UNUSED = (g4ss21*Gdsss121 + 
      g4ss22*Gdsss221 + g4ss23*Gdsss321 + betayL*Gdtss21)*pow(a2,-1);
    
    CCTK_REAL Gsss321 CCTK_ATTRIBUTE_UNUSED = (g4ss31*Gdsss121 + 
      g4ss32*Gdsss221 + g4ss33*Gdsss321 + betazL*Gdtss21)*pow(a2,-1);
    
    CCTK_REAL Gsss122 CCTK_ATTRIBUTE_UNUSED = (g4ss11*Gdsss122 + 
      g4ss12*Gdsss222 + g4ss13*Gdsss322 + betaxL*Gdtss22)*pow(a2,-1);
    
    CCTK_REAL Gsss222 CCTK_ATTRIBUTE_UNUSED = (g4ss21*Gdsss122 + 
      g4ss22*Gdsss222 + g4ss23*Gdsss322 + betayL*Gdtss22)*pow(a2,-1);
    
    CCTK_REAL Gsss322 CCTK_ATTRIBUTE_UNUSED = (g4ss31*Gdsss122 + 
      g4ss32*Gdsss222 + g4ss33*Gdsss322 + betazL*Gdtss22)*pow(a2,-1);
    
    CCTK_REAL Gsss123 CCTK_ATTRIBUTE_UNUSED = (g4ss11*Gdsss123 + 
      g4ss12*Gdsss223 + g4ss13*Gdsss323 + betaxL*Gdtss23)*pow(a2,-1);
    
    CCTK_REAL Gsss223 CCTK_ATTRIBUTE_UNUSED = (g4ss21*Gdsss123 + 
      g4ss22*Gdsss223 + g4ss23*Gdsss323 + betayL*Gdtss23)*pow(a2,-1);
    
    CCTK_REAL Gsss323 CCTK_ATTRIBUTE_UNUSED = (g4ss31*Gdsss123 + 
      g4ss32*Gdsss223 + g4ss33*Gdsss323 + betazL*Gdtss23)*pow(a2,-1);
    
    CCTK_REAL Gsss131 CCTK_ATTRIBUTE_UNUSED = (g4ss11*Gdsss131 + 
      g4ss12*Gdsss231 + g4ss13*Gdsss331 + betaxL*Gdtss31)*pow(a2,-1);
    
    CCTK_REAL Gsss231 CCTK_ATTRIBUTE_UNUSED = (g4ss21*Gdsss131 + 
      g4ss22*Gdsss231 + g4ss23*Gdsss331 + betayL*Gdtss31)*pow(a2,-1);
    
    CCTK_REAL Gsss331 CCTK_ATTRIBUTE_UNUSED = (g4ss31*Gdsss131 + 
      g4ss32*Gdsss231 + g4ss33*Gdsss331 + betazL*Gdtss31)*pow(a2,-1);
    
    CCTK_REAL Gsss132 CCTK_ATTRIBUTE_UNUSED = (g4ss11*Gdsss132 + 
      g4ss12*Gdsss232 + g4ss13*Gdsss332 + betaxL*Gdtss32)*pow(a2,-1);
    
    CCTK_REAL Gsss232 CCTK_ATTRIBUTE_UNUSED = (g4ss21*Gdsss132 + 
      g4ss22*Gdsss232 + g4ss23*Gdsss332 + betayL*Gdtss32)*pow(a2,-1);
    
    CCTK_REAL Gsss332 CCTK_ATTRIBUTE_UNUSED = (g4ss31*Gdsss132 + 
      g4ss32*Gdsss232 + g4ss33*Gdsss332 + betazL*Gdtss32)*pow(a2,-1);
    
    CCTK_REAL Gsss133 CCTK_ATTRIBUTE_UNUSED = (g4ss11*Gdsss133 + 
      g4ss12*Gdsss233 + g4ss13*Gdsss333 + betaxL*Gdtss33)*pow(a2,-1);
    
    CCTK_REAL Gsss233 CCTK_ATTRIBUTE_UNUSED = (g4ss21*Gdsss133 + 
      g4ss22*Gdsss233 + g4ss23*Gdsss333 + betayL*Gdtss33)*pow(a2,-1);
    
    CCTK_REAL Gsss333 CCTK_ATTRIBUTE_UNUSED = (g4ss31*Gdsss133 + 
      g4ss32*Gdsss233 + g4ss33*Gdsss333 + betazL*Gdtss33)*pow(a2,-1);
    
    CCTK_REAL u1 CCTK_ATTRIBUTE_UNUSED = (gxxL*velxL + gxyL*velyL + 
      gxzL*velzL)*w_lorentzL;
    
    CCTK_REAL u2 CCTK_ATTRIBUTE_UNUSED = (gxyL*velxL + gyyL*velyL + 
      gyzL*velzL)*w_lorentzL;
    
    CCTK_REAL u3 CCTK_ATTRIBUTE_UNUSED = (gxzL*velxL + gyzL*velyL + 
      gzzL*velzL)*w_lorentzL;
    
    CCTK_REAL u0u CCTK_ATTRIBUTE_UNUSED = w_lorentzL*pow(alpL,-1);
    
    CCTK_REAL u0l CCTK_ATTRIBUTE_UNUSED = -(a2*u0u) + betaxL*u1 + 
      betayL*u2 + betazL*u3;
    
    CCTK_REAL uup1 CCTK_ATTRIBUTE_UNUSED = (betaxL*u0l + g4ss11*u1 + 
      g4ss12*u2 + g4ss13*u3)*pow(a2,-1);
    
    CCTK_REAL uup2 CCTK_ATTRIBUTE_UNUSED = (betayL*u0l + g4ss21*u1 + 
      g4ss22*u2 + g4ss23*u3)*pow(a2,-1);
    
    CCTK_REAL uup3 CCTK_ATTRIBUTE_UNUSED = (betazL*u0l + g4ss31*u1 + 
      g4ss32*u2 + g4ss33*u3)*pow(a2,-1);
    
    CCTK_REAL dtu1 CCTK_ATTRIBUTE_UNUSED = pow(DDL + EEL,-1)*(-((DDrhsL + 
      EErhsL)*u1) + SS1rhsL*pow(1 + eosw,-1));
    
    CCTK_REAL dtu2 CCTK_ATTRIBUTE_UNUSED = pow(DDL + EEL,-1)*(-((DDrhsL + 
      EErhsL)*u2) + SS2rhsL*pow(1 + eosw,-1));
    
    CCTK_REAL dtu3 CCTK_ATTRIBUTE_UNUSED = pow(DDL + EEL,-1)*(-((DDrhsL + 
      EErhsL)*u3) + SS3rhsL*pow(1 + eosw,-1));
    
    CCTK_REAL dtu0u CCTK_ATTRIBUTE_UNUSED = -(dtalpL*u0u*pow(alpL,-1)) + 
      0.5*pow(a2,-1)*pow(u0u,-1)*((dtgu31*u1 + (dtgu23 + dtgu32)*u2)*u3 + 
      u1*((dtgu12 + dtgu21)*u2 + dtgu13*u3) + 2*(dtu1*(gu11*u1 + gu12*u2 + 
      gu13*u3) + dtu2*(gu12*u1 + gu22*u2 + gu23*u3) + dtu3*(gu13*u1 + gu23*u2 
      + gu33*u3)) + dtgu11*pow(u1,2) + dtgu22*pow(u2,2) + dtgu33*pow(u3,2));
    
    CCTK_REAL dtu0l CCTK_ATTRIBUTE_UNUSED = -(a2*dtu0u) + betaxL*dtu1 + 
      betayL*dtu2 + betazL*dtu3 - 2*alpL*dtalpL*u0u + dtbetaxL*u1 + 
      dtbetayL*u2 + dtbetazL*u3;
    
    CCTK_REAL dsu11 CCTK_ATTRIBUTE_UNUSED = 
      w_lorentzL*(velxL*JacPDstandardNth1gxx + velyL*JacPDstandardNth1gxy + 
      velzL*JacPDstandardNth1gxz + gxxL*JacPDstandardNth1velx + 
      gxyL*JacPDstandardNth1vely + gxzL*JacPDstandardNth1velz) + (gxxL*velxL 
      + gxyL*velyL + gxzL*velzL)*JacPDstandardNth1w_lorentz;
    
    CCTK_REAL dsu12 CCTK_ATTRIBUTE_UNUSED = 
      w_lorentzL*(velxL*JacPDstandardNth1gxy + velyL*JacPDstandardNth1gyy + 
      velzL*JacPDstandardNth1gyz + gxyL*JacPDstandardNth1velx + 
      gyyL*JacPDstandardNth1vely + gyzL*JacPDstandardNth1velz) + (gxyL*velxL 
      + gyyL*velyL + gyzL*velzL)*JacPDstandardNth1w_lorentz;
    
    CCTK_REAL dsu13 CCTK_ATTRIBUTE_UNUSED = 
      w_lorentzL*(velxL*JacPDstandardNth1gxz + velyL*JacPDstandardNth1gyz + 
      velzL*JacPDstandardNth1gzz + gxzL*JacPDstandardNth1velx + 
      gyzL*JacPDstandardNth1vely + gzzL*JacPDstandardNth1velz) + (gxzL*velxL 
      + gyzL*velyL + gzzL*velzL)*JacPDstandardNth1w_lorentz;
    
    CCTK_REAL dsu21 CCTK_ATTRIBUTE_UNUSED = 
      w_lorentzL*(velxL*JacPDstandardNth2gxx + velyL*JacPDstandardNth2gxy + 
      velzL*JacPDstandardNth2gxz + gxxL*JacPDstandardNth2velx + 
      gxyL*JacPDstandardNth2vely + gxzL*JacPDstandardNth2velz) + (gxxL*velxL 
      + gxyL*velyL + gxzL*velzL)*JacPDstandardNth2w_lorentz;
    
    CCTK_REAL dsu22 CCTK_ATTRIBUTE_UNUSED = 
      w_lorentzL*(velxL*JacPDstandardNth2gxy + velyL*JacPDstandardNth2gyy + 
      velzL*JacPDstandardNth2gyz + gxyL*JacPDstandardNth2velx + 
      gyyL*JacPDstandardNth2vely + gyzL*JacPDstandardNth2velz) + (gxyL*velxL 
      + gyyL*velyL + gyzL*velzL)*JacPDstandardNth2w_lorentz;
    
    CCTK_REAL dsu23 CCTK_ATTRIBUTE_UNUSED = 
      w_lorentzL*(velxL*JacPDstandardNth2gxz + velyL*JacPDstandardNth2gyz + 
      velzL*JacPDstandardNth2gzz + gxzL*JacPDstandardNth2velx + 
      gyzL*JacPDstandardNth2vely + gzzL*JacPDstandardNth2velz) + (gxzL*velxL 
      + gyzL*velyL + gzzL*velzL)*JacPDstandardNth2w_lorentz;
    
    CCTK_REAL dsu31 CCTK_ATTRIBUTE_UNUSED = 
      w_lorentzL*(velxL*JacPDstandardNth3gxx + velyL*JacPDstandardNth3gxy + 
      velzL*JacPDstandardNth3gxz + gxxL*JacPDstandardNth3velx + 
      gxyL*JacPDstandardNth3vely + gxzL*JacPDstandardNth3velz) + (gxxL*velxL 
      + gxyL*velyL + gxzL*velzL)*JacPDstandardNth3w_lorentz;
    
    CCTK_REAL dsu32 CCTK_ATTRIBUTE_UNUSED = 
      w_lorentzL*(velxL*JacPDstandardNth3gxy + velyL*JacPDstandardNth3gyy + 
      velzL*JacPDstandardNth3gyz + gxyL*JacPDstandardNth3velx + 
      gyyL*JacPDstandardNth3vely + gyzL*JacPDstandardNth3velz) + (gxyL*velxL 
      + gyyL*velyL + gyzL*velzL)*JacPDstandardNth3w_lorentz;
    
    CCTK_REAL dsu33 CCTK_ATTRIBUTE_UNUSED = 
      w_lorentzL*(velxL*JacPDstandardNth3gxz + velyL*JacPDstandardNth3gyz + 
      velzL*JacPDstandardNth3gzz + gxzL*JacPDstandardNth3velx + 
      gyzL*JacPDstandardNth3vely + gzzL*JacPDstandardNth3velz) + (gxzL*velxL 
      + gyzL*velyL + gzzL*velzL)*JacPDstandardNth3w_lorentz;
    
    CCTK_REAL dsu01 CCTK_ATTRIBUTE_UNUSED = -(a2*dtu0u) + betaxL*dtu1 + 
      betayL*dtu2 + betazL*dtu3 - 2*alpL*dtalpL*u0u + dtbetaxL*u1 + 
      dtbetayL*u2 + dtbetazL*u3;
    
    CCTK_REAL dsu02 CCTK_ATTRIBUTE_UNUSED = -(a2*dtu0u) + betaxL*dtu1 + 
      betayL*dtu2 + betazL*dtu3 - 2*alpL*dtalpL*u0u + dtbetaxL*u1 + 
      dtbetayL*u2 + dtbetazL*u3;
    
    CCTK_REAL dsu03 CCTK_ATTRIBUTE_UNUSED = -(a2*dtu0u) + betaxL*dtu1 + 
      betayL*dtu2 + betazL*dtu3 - 2*alpL*dtalpL*u0u + dtbetaxL*u1 + 
      dtbetayL*u2 + dtbetazL*u3;
    
    CCTK_REAL Dtut CCTK_ATTRIBUTE_UNUSED = dtu0l - Gttt*u0l - Gstt1*u1 - 
      Gstt2*u2 - Gstt3*u3;
    
    CCTK_REAL Gtsu1 CCTK_ATTRIBUTE_UNUSED = Gtts1*u0l + Gsts11*u1 + 
      Gsts21*u2 + Gsts31*u3;
    
    CCTK_REAL Gtsu2 CCTK_ATTRIBUTE_UNUSED = Gtts2*u0l + Gsts12*u1 + 
      Gsts22*u2 + Gsts32*u3;
    
    CCTK_REAL Gtsu3 CCTK_ATTRIBUTE_UNUSED = Gtts3*u0l + Gsts13*u1 + 
      Gsts23*u2 + Gsts33*u3;
    
    CCTK_REAL Dsut1 CCTK_ATTRIBUTE_UNUSED = dsu01 - Gtsu1;
    
    CCTK_REAL Dsut2 CCTK_ATTRIBUTE_UNUSED = dsu02 - Gtsu2;
    
    CCTK_REAL Dsut3 CCTK_ATTRIBUTE_UNUSED = dsu03 - Gtsu3;
    
    CCTK_REAL Dtus1 CCTK_ATTRIBUTE_UNUSED = dtu1 - Gtsu1;
    
    CCTK_REAL Dtus2 CCTK_ATTRIBUTE_UNUSED = dtu2 - Gtsu2;
    
    CCTK_REAL Dtus3 CCTK_ATTRIBUTE_UNUSED = dtu3 - Gtsu3;
    
    CCTK_REAL Dsus11 CCTK_ATTRIBUTE_UNUSED = dsu11 - Gtss11*u0l - 
      Gsss111*u1 - Gsss211*u2 - Gsss311*u3;
    
    CCTK_REAL Dsus12 CCTK_ATTRIBUTE_UNUSED = dsu12 - Gtss12*u0l - 
      Gsss112*u1 - Gsss212*u2 - Gsss312*u3;
    
    CCTK_REAL Dsus13 CCTK_ATTRIBUTE_UNUSED = dsu13 - Gtss13*u0l - 
      Gsss113*u1 - Gsss213*u2 - Gsss313*u3;
    
    CCTK_REAL Dsus21 CCTK_ATTRIBUTE_UNUSED = dsu21 - Gtss21*u0l - 
      Gsss121*u1 - Gsss221*u2 - Gsss321*u3;
    
    CCTK_REAL Dsus22 CCTK_ATTRIBUTE_UNUSED = dsu22 - Gtss22*u0l - 
      Gsss122*u1 - Gsss222*u2 - Gsss322*u3;
    
    CCTK_REAL Dsus23 CCTK_ATTRIBUTE_UNUSED = dsu23 - Gtss23*u0l - 
      Gsss123*u1 - Gsss223*u2 - Gsss323*u3;
    
    CCTK_REAL Dsus31 CCTK_ATTRIBUTE_UNUSED = dsu31 - Gtss31*u0l - 
      Gsss131*u1 - Gsss231*u2 - Gsss331*u3;
    
    CCTK_REAL Dsus32 CCTK_ATTRIBUTE_UNUSED = dsu32 - Gtss32*u0l - 
      Gsss132*u1 - Gsss232*u2 - Gsss332*u3;
    
    CCTK_REAL Dsus33 CCTK_ATTRIBUTE_UNUSED = dsu33 - Gtss33*u0l - 
      Gsss133*u1 - Gsss233*u2 - Gsss333*u3;
    
    CCTK_REAL htt CCTK_ATTRIBUTE_UNUSED = 2*(betaxL*(betayL*gxyL + 
      betazL*gxzL) + betayL*betazL*gyzL) - a2 + gxxL*pow(betaxL,2) + 
      gyyL*pow(betayL,2) + gzzL*pow(betazL,2) + pow(u0l,2);
    
    CCTK_REAL hts1 CCTK_ATTRIBUTE_UNUSED = bl1 + u0l*u1;
    
    CCTK_REAL hts2 CCTK_ATTRIBUTE_UNUSED = bl2 + u0l*u2;
    
    CCTK_REAL hts3 CCTK_ATTRIBUTE_UNUSED = bl3 + u0l*u3;
    
    CCTK_REAL hss11 CCTK_ATTRIBUTE_UNUSED = gxxL + pow(u1,2);
    
    CCTK_REAL hss12 CCTK_ATTRIBUTE_UNUSED = gxyL + u1*u2;
    
    CCTK_REAL hss13 CCTK_ATTRIBUTE_UNUSED = gxzL + u1*u3;
    
    CCTK_REAL hss22 CCTK_ATTRIBUTE_UNUSED = gyyL + pow(u2,2);
    
    CCTK_REAL hss23 CCTK_ATTRIBUTE_UNUSED = gyzL + u2*u3;
    
    CCTK_REAL hss33 CCTK_ATTRIBUTE_UNUSED = gzzL + pow(u3,2);
    
    CCTK_REAL hmtt CCTK_ATTRIBUTE_UNUSED = 1 + u0l*u0u;
    
    CCTK_REAL hmtsl1 CCTK_ATTRIBUTE_UNUSED = u0u*u1;
    
    CCTK_REAL hmtsl2 CCTK_ATTRIBUTE_UNUSED = u0u*u2;
    
    CCTK_REAL hmtsl3 CCTK_ATTRIBUTE_UNUSED = u0u*u3;
    
    CCTK_REAL hmtsu1 CCTK_ATTRIBUTE_UNUSED = u0l*uup1;
    
    CCTK_REAL hmtsu2 CCTK_ATTRIBUTE_UNUSED = u0l*uup2;
    
    CCTK_REAL hmtsu3 CCTK_ATTRIBUTE_UNUSED = u0l*uup3;
    
    CCTK_REAL hmss11 CCTK_ATTRIBUTE_UNUSED = 1 + u1*uup1;
    
    CCTK_REAL hmss21 CCTK_ATTRIBUTE_UNUSED = u1*uup2;
    
    CCTK_REAL hmss31 CCTK_ATTRIBUTE_UNUSED = u1*uup3;
    
    CCTK_REAL hmss12 CCTK_ATTRIBUTE_UNUSED = u2*uup1;
    
    CCTK_REAL hmss22 CCTK_ATTRIBUTE_UNUSED = 1 + u2*uup2;
    
    CCTK_REAL hmss32 CCTK_ATTRIBUTE_UNUSED = u2*uup3;
    
    CCTK_REAL hmss13 CCTK_ATTRIBUTE_UNUSED = u3*uup1;
    
    CCTK_REAL hmss23 CCTK_ATTRIBUTE_UNUSED = u3*uup2;
    
    CCTK_REAL hmss33 CCTK_ATTRIBUTE_UNUSED = 1 + u3*uup3;
    
    CCTK_REAL Dptut CCTK_ATTRIBUTE_UNUSED = hmtsu1*(Dsus11*hmtsu1 + 
      Dsus12*hmtsu2 + Dsus13*hmtsu3 + Dsut1*hmtt) + hmtsu2*(Dsus21*hmtsu1 + 
      Dsus22*hmtsu2 + Dsus23*hmtsu3 + Dsut2*hmtt) + hmtsu3*(Dsus31*hmtsu1 + 
      Dsus32*hmtsu2 + Dsus33*hmtsu3 + Dsut3*hmtt) + hmtt*(Dtus1*hmtsu1 + 
      Dtus2*hmtsu2 + Dtus3*hmtsu3 + Dtut*hmtt);
    
    CCTK_REAL Dptus1 CCTK_ATTRIBUTE_UNUSED = (Dsus11*hmss11 + 
      Dsus12*hmss21 + Dsus13*hmss31 + Dsut1*hmtsl1)*hmtsu1 + (Dsus21*hmss11 + 
      Dsus22*hmss21 + Dsus23*hmss31 + Dsut2*hmtsl1)*hmtsu2 + (Dsus31*hmss11 + 
      Dsus32*hmss21 + Dsus33*hmss31 + Dsut3*hmtsl1)*hmtsu3 + (Dtus1*hmss11 + 
      Dtus2*hmss21 + Dtus3*hmss31 + Dtut*hmtsl1)*hmtt;
    
    CCTK_REAL Dptus2 CCTK_ATTRIBUTE_UNUSED = (Dsus11*hmss12 + 
      Dsus12*hmss22 + Dsus13*hmss32 + Dsut1*hmtsl2)*hmtsu1 + (Dsus21*hmss12 + 
      Dsus22*hmss22 + Dsus23*hmss32 + Dsut2*hmtsl2)*hmtsu2 + (Dsus31*hmss12 + 
      Dsus32*hmss22 + Dsus33*hmss32 + Dsut3*hmtsl2)*hmtsu3 + (Dtus1*hmss12 + 
      Dtus2*hmss22 + Dtus3*hmss32 + Dtut*hmtsl2)*hmtt;
    
    CCTK_REAL Dptus3 CCTK_ATTRIBUTE_UNUSED = (Dsus11*hmss13 + 
      Dsus12*hmss23 + Dsus13*hmss33 + Dsut1*hmtsl3)*hmtsu1 + (Dsus21*hmss13 + 
      Dsus22*hmss23 + Dsus23*hmss33 + Dsut2*hmtsl3)*hmtsu2 + (Dsus31*hmss13 + 
      Dsus32*hmss23 + Dsus33*hmss33 + Dsut3*hmtsl3)*hmtsu3 + (Dtus1*hmss13 + 
      Dtus2*hmss23 + Dtus3*hmss33 + Dtut*hmtsl3)*hmtt;
    
    CCTK_REAL Dpsut1 CCTK_ATTRIBUTE_UNUSED = hmss11*(Dsus11*hmtsu1 + 
      Dsus12*hmtsu2 + Dsus13*hmtsu3 + Dsut1*hmtt) + hmss21*(Dsus21*hmtsu1 + 
      Dsus22*hmtsu2 + Dsus23*hmtsu3 + Dsut2*hmtt) + hmss31*(Dsus31*hmtsu1 + 
      Dsus32*hmtsu2 + Dsus33*hmtsu3 + Dsut3*hmtt) + hmtsl1*(Dtus1*hmtsu1 + 
      Dtus2*hmtsu2 + Dtus3*hmtsu3 + Dtut*hmtt);
    
    CCTK_REAL Dpsut2 CCTK_ATTRIBUTE_UNUSED = hmss12*(Dsus11*hmtsu1 + 
      Dsus12*hmtsu2 + Dsus13*hmtsu3 + Dsut1*hmtt) + hmss22*(Dsus21*hmtsu1 + 
      Dsus22*hmtsu2 + Dsus23*hmtsu3 + Dsut2*hmtt) + hmss32*(Dsus31*hmtsu1 + 
      Dsus32*hmtsu2 + Dsus33*hmtsu3 + Dsut3*hmtt) + hmtsl2*(Dtus1*hmtsu1 + 
      Dtus2*hmtsu2 + Dtus3*hmtsu3 + Dtut*hmtt);
    
    CCTK_REAL Dpsut3 CCTK_ATTRIBUTE_UNUSED = hmss13*(Dsus11*hmtsu1 + 
      Dsus12*hmtsu2 + Dsus13*hmtsu3 + Dsut1*hmtt) + hmss23*(Dsus21*hmtsu1 + 
      Dsus22*hmtsu2 + Dsus23*hmtsu3 + Dsut2*hmtt) + hmss33*(Dsus31*hmtsu1 + 
      Dsus32*hmtsu2 + Dsus33*hmtsu3 + Dsut3*hmtt) + hmtsl3*(Dtus1*hmtsu1 + 
      Dtus2*hmtsu2 + Dtus3*hmtsu3 + Dtut*hmtt);
    
    CCTK_REAL Dpsus11 CCTK_ATTRIBUTE_UNUSED = hmss11*(Dsus11*hmss11 + 
      Dsus12*hmss21 + Dsus13*hmss31 + Dsut1*hmtsl1) + hmss21*(Dsus21*hmss11 + 
      Dsus22*hmss21 + Dsus23*hmss31 + Dsut2*hmtsl1) + hmss31*(Dsus31*hmss11 + 
      Dsus32*hmss21 + Dsus33*hmss31 + Dsut3*hmtsl1) + hmtsl1*(Dtus1*hmss11 + 
      Dtus2*hmss21 + Dtus3*hmss31 + Dtut*hmtsl1);
    
    CCTK_REAL Dpsus12 CCTK_ATTRIBUTE_UNUSED = hmss11*(Dsus11*hmss12 + 
      Dsus12*hmss22 + Dsus13*hmss32 + Dsut1*hmtsl2) + hmss21*(Dsus21*hmss12 + 
      Dsus22*hmss22 + Dsus23*hmss32 + Dsut2*hmtsl2) + hmss31*(Dsus31*hmss12 + 
      Dsus32*hmss22 + Dsus33*hmss32 + Dsut3*hmtsl2) + hmtsl1*(Dtus1*hmss12 + 
      Dtus2*hmss22 + Dtus3*hmss32 + Dtut*hmtsl2);
    
    CCTK_REAL Dpsus13 CCTK_ATTRIBUTE_UNUSED = hmss11*(Dsus11*hmss13 + 
      Dsus12*hmss23 + Dsus13*hmss33 + Dsut1*hmtsl3) + hmss21*(Dsus21*hmss13 + 
      Dsus22*hmss23 + Dsus23*hmss33 + Dsut2*hmtsl3) + hmss31*(Dsus31*hmss13 + 
      Dsus32*hmss23 + Dsus33*hmss33 + Dsut3*hmtsl3) + hmtsl1*(Dtus1*hmss13 + 
      Dtus2*hmss23 + Dtus3*hmss33 + Dtut*hmtsl3);
    
    CCTK_REAL Dpsus21 CCTK_ATTRIBUTE_UNUSED = hmss12*(Dsus11*hmss11 + 
      Dsus12*hmss21 + Dsus13*hmss31 + Dsut1*hmtsl1) + hmss22*(Dsus21*hmss11 + 
      Dsus22*hmss21 + Dsus23*hmss31 + Dsut2*hmtsl1) + hmss32*(Dsus31*hmss11 + 
      Dsus32*hmss21 + Dsus33*hmss31 + Dsut3*hmtsl1) + (Dtus1*hmss11 + 
      Dtus2*hmss21 + Dtus3*hmss31 + Dtut*hmtsl1)*hmtsl2;
    
    CCTK_REAL Dpsus22 CCTK_ATTRIBUTE_UNUSED = hmss12*(Dsus11*hmss12 + 
      Dsus12*hmss22 + Dsus13*hmss32 + Dsut1*hmtsl2) + hmss22*(Dsus21*hmss12 + 
      Dsus22*hmss22 + Dsus23*hmss32 + Dsut2*hmtsl2) + hmss32*(Dsus31*hmss12 + 
      Dsus32*hmss22 + Dsus33*hmss32 + Dsut3*hmtsl2) + hmtsl2*(Dtus1*hmss12 + 
      Dtus2*hmss22 + Dtus3*hmss32 + Dtut*hmtsl2);
    
    CCTK_REAL Dpsus23 CCTK_ATTRIBUTE_UNUSED = hmss12*(Dsus11*hmss13 + 
      Dsus12*hmss23 + Dsus13*hmss33 + Dsut1*hmtsl3) + hmss22*(Dsus21*hmss13 + 
      Dsus22*hmss23 + Dsus23*hmss33 + Dsut2*hmtsl3) + hmss32*(Dsus31*hmss13 + 
      Dsus32*hmss23 + Dsus33*hmss33 + Dsut3*hmtsl3) + hmtsl2*(Dtus1*hmss13 + 
      Dtus2*hmss23 + Dtus3*hmss33 + Dtut*hmtsl3);
    
    CCTK_REAL Dpsus31 CCTK_ATTRIBUTE_UNUSED = hmss13*(Dsus11*hmss11 + 
      Dsus12*hmss21 + Dsus13*hmss31 + Dsut1*hmtsl1) + hmss23*(Dsus21*hmss11 + 
      Dsus22*hmss21 + Dsus23*hmss31 + Dsut2*hmtsl1) + hmss33*(Dsus31*hmss11 + 
      Dsus32*hmss21 + Dsus33*hmss31 + Dsut3*hmtsl1) + (Dtus1*hmss11 + 
      Dtus2*hmss21 + Dtus3*hmss31 + Dtut*hmtsl1)*hmtsl3;
    
    CCTK_REAL Dpsus32 CCTK_ATTRIBUTE_UNUSED = hmss13*(Dsus11*hmss12 + 
      Dsus12*hmss22 + Dsus13*hmss32 + Dsut1*hmtsl2) + hmss23*(Dsus21*hmss12 + 
      Dsus22*hmss22 + Dsus23*hmss32 + Dsut2*hmtsl2) + hmss33*(Dsus31*hmss12 + 
      Dsus32*hmss22 + Dsus33*hmss32 + Dsut3*hmtsl2) + (Dtus1*hmss12 + 
      Dtus2*hmss22 + Dtus3*hmss32 + Dtut*hmtsl2)*hmtsl3;
    
    CCTK_REAL Dpsus33 CCTK_ATTRIBUTE_UNUSED = hmss13*(Dsus11*hmss13 + 
      Dsus12*hmss23 + Dsus13*hmss33 + Dsut1*hmtsl3) + hmss23*(Dsus21*hmss13 + 
      Dsus22*hmss23 + Dsus23*hmss33 + Dsut2*hmtsl3) + hmss33*(Dsus31*hmss13 + 
      Dsus32*hmss23 + Dsus33*hmss33 + Dsut3*hmtsl3) + hmtsl3*(Dtus1*hmss13 + 
      Dtus2*hmss23 + Dtus3*hmss33 + Dtut*hmtsl3);
    
    CCTK_REAL Thetatt CCTK_ATTRIBUTE_UNUSED = Dptut;
    
    CCTK_REAL Thetats1 CCTK_ATTRIBUTE_UNUSED = 0.5*(Dpsut1 + Dptus1);
    
    CCTK_REAL Thetats2 CCTK_ATTRIBUTE_UNUSED = 0.5*(Dpsut2 + Dptus2);
    
    CCTK_REAL Thetats3 CCTK_ATTRIBUTE_UNUSED = 0.5*(Dpsut3 + Dptus3);
    
    CCTK_REAL Thetass11 CCTK_ATTRIBUTE_UNUSED = Dpsus11;
    
    CCTK_REAL Thetass12 CCTK_ATTRIBUTE_UNUSED = 0.5*(Dpsus12 + Dpsus21);
    
    CCTK_REAL Thetass13 CCTK_ATTRIBUTE_UNUSED = 0.5*(Dpsus13 + Dpsus31);
    
    CCTK_REAL Thetass22 CCTK_ATTRIBUTE_UNUSED = Dpsus22;
    
    CCTK_REAL Thetass23 CCTK_ATTRIBUTE_UNUSED = 0.5*(Dpsus23 + Dpsus32);
    
    CCTK_REAL Thetass33 CCTK_ATTRIBUTE_UNUSED = Dpsus33;
    
    ThetaL = 0.;
    
    CCTK_REAL sigmattL CCTK_ATTRIBUTE_UNUSED = 
      -0.333333333333333333333333333333*ThetaL*htt + Thetatt;
    
    CCTK_REAL sigmats1L CCTK_ATTRIBUTE_UNUSED = 
      -0.333333333333333333333333333333*ThetaL*hts1 + Thetats1;
    
    CCTK_REAL sigmats2L CCTK_ATTRIBUTE_UNUSED = 
      -0.333333333333333333333333333333*ThetaL*hts2 + Thetats2;
    
    CCTK_REAL sigmats3L CCTK_ATTRIBUTE_UNUSED = 
      -0.333333333333333333333333333333*ThetaL*hts3 + Thetats3;
    
    CCTK_REAL sigmass11L CCTK_ATTRIBUTE_UNUSED = 
      -0.333333333333333333333333333333*ThetaL*hss11 + Thetass11;
    
    CCTK_REAL sigmass12L CCTK_ATTRIBUTE_UNUSED = 
      -0.333333333333333333333333333333*ThetaL*hss12 + Thetass12;
    
    CCTK_REAL sigmass13L CCTK_ATTRIBUTE_UNUSED = 
      -0.333333333333333333333333333333*ThetaL*hss13 + Thetass13;
    
    CCTK_REAL sigmass22L CCTK_ATTRIBUTE_UNUSED = 
      -0.333333333333333333333333333333*ThetaL*hss22 + Thetass22;
    
    CCTK_REAL sigmass23L CCTK_ATTRIBUTE_UNUSED = 
      -0.333333333333333333333333333333*ThetaL*hss23 + Thetass23;
    
    CCTK_REAL sigmass33L CCTK_ATTRIBUTE_UNUSED = 
      -0.333333333333333333333333333333*ThetaL*hss33 + Thetass33;
    /* Copy local copies back to grid functions */
    sigmass11[index] = sigmass11L;
    sigmass12[index] = sigmass12L;
    sigmass13[index] = sigmass13L;
    sigmass22[index] = sigmass22L;
    sigmass23[index] = sigmass23L;
    sigmass33[index] = sigmass33L;
    sigmats1[index] = sigmats1L;
    sigmats2[index] = sigmats2L;
    sigmats3[index] = sigmats3L;
    sigmatt[index] = sigmattL;
    Theta[index] = ThetaL;
  }
  CCTK_ENDLOOP3(CT_Dust_Expansion);
}
extern "C" void CT_Dust_Expansion(CCTK_ARGUMENTS)
{
  #ifdef DECLARE_CCTK_ARGUMENTS_CT_Dust_Expansion
  DECLARE_CCTK_ARGUMENTS_CHECKED(CT_Dust_Expansion);
  #else
  DECLARE_CCTK_ARGUMENTS;
  #endif
  DECLARE_CCTK_PARAMETERS;
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering CT_Dust_Expansion_Body");
  }
  if (cctk_iteration % CT_Dust_Expansion_calc_every != CT_Dust_Expansion_calc_offset)
  {
    return;
  }
  
  const char* const groups[] = {
    "ADMBase::curv",
    "ADMBase::dtlapse",
    "ADMBase::dtshift",
    "ADMBase::lapse",
    "ADMBase::metric",
    "ADMBase::shift",
    "CT_Dust::CT_D",
    "CT_Dust::CT_Drhs",
    "CT_Dust::CT_E",
    "CT_Dust::CT_Erhs",
    "CT_Dust::CT_sigmass",
    "CT_Dust::CT_sigmats",
    "CT_Dust::CT_sigmatt",
    "CT_Dust::CT_Srhs",
    "CT_Dust::CT_Theta",
    "HydroBase::vel",
    "HydroBase::w_lorentz"};
  AssertGroupStorage(cctkGH, "CT_Dust_Expansion", 17, groups);
  
  switch (fdOrder)
  {
    case 2:
    {
      EnsureStencilFits(cctkGH, "CT_Dust_Expansion", 1, 1, 1);
      break;
    }
    
    case 4:
    {
      EnsureStencilFits(cctkGH, "CT_Dust_Expansion", 2, 2, 2);
      break;
    }
    
    case 6:
    {
      EnsureStencilFits(cctkGH, "CT_Dust_Expansion", 3, 3, 3);
      break;
    }
    
    case 8:
    {
      EnsureStencilFits(cctkGH, "CT_Dust_Expansion", 4, 4, 4);
      break;
    }
    default:
      CCTK_BUILTIN_UNREACHABLE();
  }
  
  LoopOverInterior(cctkGH, CT_Dust_Expansion_Body);
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Leaving CT_Dust_Expansion_Body");
  }
}

} // namespace CT_Dust
