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

namespace CosmoLapse {

extern "C" void CosmoLapse_dtLapseRhoPostRHS_SelectBCs(CCTK_ARGUMENTS)
{
  #ifdef DECLARE_CCTK_ARGUMENTS_CosmoLapse_dtLapseRhoPostRHS_SelectBCs
  DECLARE_CCTK_ARGUMENTS_CHECKED(CosmoLapse_dtLapseRhoPostRHS_SelectBCs);
  #else
  DECLARE_CCTK_ARGUMENTS;
  #endif
  DECLARE_CCTK_PARAMETERS;
  
  if (cctk_iteration % CosmoLapse_dtLapseRhoPostRHS_calc_every != CosmoLapse_dtLapseRhoPostRHS_calc_offset)
    return;
  CCTK_INT ierr CCTK_ATTRIBUTE_UNUSED = 0;
  ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "ADMBase::dtlapse","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for ADMBase::dtlapse.");
  ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN::ML_lapserhs","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for ML_BSSN::ML_lapserhs.");
  return;
}

static void CosmoLapse_dtLapseRhoPostRHS_Body(const cGH* restrict const cctkGH, const int dir, const int face, const CCTK_REAL normal[3], const CCTK_REAL tangentA[3], const CCTK_REAL tangentB[3], const int imin[3], const int imax[3], const int n_subblock_gfs, CCTK_REAL* restrict const subblock_gfs[])
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
  const CCTK_REAL p1o180dx2 CCTK_ATTRIBUTE_UNUSED = 0.00555555555555555555555555555556*pow(dx,-2);
  const CCTK_REAL p1o180dy2 CCTK_ATTRIBUTE_UNUSED = 0.00555555555555555555555555555556*pow(dy,-2);
  const CCTK_REAL p1o180dz2 CCTK_ATTRIBUTE_UNUSED = 0.00555555555555555555555555555556*pow(dz,-2);
  const CCTK_REAL p1o24dx CCTK_ATTRIBUTE_UNUSED = 0.0416666666666666666666666666667*pow(dx,-1);
  const CCTK_REAL p1o24dy CCTK_ATTRIBUTE_UNUSED = 0.0416666666666666666666666666667*pow(dy,-1);
  const CCTK_REAL p1o24dz CCTK_ATTRIBUTE_UNUSED = 0.0416666666666666666666666666667*pow(dz,-1);
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
  const CCTK_REAL p1odx2 CCTK_ATTRIBUTE_UNUSED = pow(dx,-2);
  const CCTK_REAL p1ody2 CCTK_ATTRIBUTE_UNUSED = pow(dy,-2);
  const CCTK_REAL p1odz2 CCTK_ATTRIBUTE_UNUSED = pow(dz,-2);
  const CCTK_REAL pm1o120dx CCTK_ATTRIBUTE_UNUSED = -0.00833333333333333333333333333333*pow(dx,-1);
  const CCTK_REAL pm1o120dy CCTK_ATTRIBUTE_UNUSED = -0.00833333333333333333333333333333*pow(dy,-1);
  const CCTK_REAL pm1o120dz CCTK_ATTRIBUTE_UNUSED = -0.00833333333333333333333333333333*pow(dz,-1);
  const CCTK_REAL pm1o12dx CCTK_ATTRIBUTE_UNUSED = -0.0833333333333333333333333333333*pow(dx,-1);
  const CCTK_REAL pm1o12dx2 CCTK_ATTRIBUTE_UNUSED = -0.0833333333333333333333333333333*pow(dx,-2);
  const CCTK_REAL pm1o12dy CCTK_ATTRIBUTE_UNUSED = -0.0833333333333333333333333333333*pow(dy,-1);
  const CCTK_REAL pm1o12dy2 CCTK_ATTRIBUTE_UNUSED = -0.0833333333333333333333333333333*pow(dy,-2);
  const CCTK_REAL pm1o12dz CCTK_ATTRIBUTE_UNUSED = -0.0833333333333333333333333333333*pow(dz,-1);
  const CCTK_REAL pm1o12dz2 CCTK_ATTRIBUTE_UNUSED = -0.0833333333333333333333333333333*pow(dz,-2);
  const CCTK_REAL pm1o16dx CCTK_ATTRIBUTE_UNUSED = -0.0625*pow(dx,-1);
  const CCTK_REAL pm1o16dy CCTK_ATTRIBUTE_UNUSED = -0.0625*pow(dy,-1);
  const CCTK_REAL pm1o16dz CCTK_ATTRIBUTE_UNUSED = -0.0625*pow(dz,-1);
  const CCTK_REAL pm1o256dx CCTK_ATTRIBUTE_UNUSED = -0.00390625*pow(dx,-1);
  const CCTK_REAL pm1o256dy CCTK_ATTRIBUTE_UNUSED = -0.00390625*pow(dy,-1);
  const CCTK_REAL pm1o256dz CCTK_ATTRIBUTE_UNUSED = -0.00390625*pow(dz,-1);
  const CCTK_REAL pm1o2dx CCTK_ATTRIBUTE_UNUSED = -0.5*pow(dx,-1);
  const CCTK_REAL pm1o2dy CCTK_ATTRIBUTE_UNUSED = -0.5*pow(dy,-1);
  const CCTK_REAL pm1o2dz CCTK_ATTRIBUTE_UNUSED = -0.5*pow(dz,-1);
  const CCTK_REAL pm1o4dx CCTK_ATTRIBUTE_UNUSED = -0.25*pow(dx,-1);
  const CCTK_REAL pm1o4dy CCTK_ATTRIBUTE_UNUSED = -0.25*pow(dy,-1);
  const CCTK_REAL pm1o4dz CCTK_ATTRIBUTE_UNUSED = -0.25*pow(dz,-1);
  const CCTK_REAL pm1o60dx CCTK_ATTRIBUTE_UNUSED = -0.0166666666666666666666666666667*pow(dx,-1);
  const CCTK_REAL pm1o60dy CCTK_ATTRIBUTE_UNUSED = -0.0166666666666666666666666666667*pow(dy,-1);
  const CCTK_REAL pm1o60dz CCTK_ATTRIBUTE_UNUSED = -0.0166666666666666666666666666667*pow(dz,-1);
  const CCTK_REAL pm1o6dx CCTK_ATTRIBUTE_UNUSED = -0.166666666666666666666666666667*pow(dx,-1);
  const CCTK_REAL pm1o6dy CCTK_ATTRIBUTE_UNUSED = -0.166666666666666666666666666667*pow(dy,-1);
  const CCTK_REAL pm1o6dz CCTK_ATTRIBUTE_UNUSED = -0.166666666666666666666666666667*pow(dz,-1);
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
  #pragma omp parallel
  CCTK_LOOP3(CosmoLapse_dtLapseRhoPostRHS,
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
    CCTK_REAL gxxL CCTK_ATTRIBUTE_UNUSED = gxx[index];
    CCTK_REAL gxyL CCTK_ATTRIBUTE_UNUSED = gxy[index];
    CCTK_REAL gxzL CCTK_ATTRIBUTE_UNUSED = gxz[index];
    CCTK_REAL gyyL CCTK_ATTRIBUTE_UNUSED = gyy[index];
    CCTK_REAL gyzL CCTK_ATTRIBUTE_UNUSED = gyz[index];
    CCTK_REAL gzzL CCTK_ATTRIBUTE_UNUSED = gzz[index];
    CCTK_REAL KtraceL CCTK_ATTRIBUTE_UNUSED = Ktrace[index];
    CCTK_REAL kxxL CCTK_ATTRIBUTE_UNUSED = kxx[index];
    CCTK_REAL kxyL CCTK_ATTRIBUTE_UNUSED = kxy[index];
    CCTK_REAL kxzL CCTK_ATTRIBUTE_UNUSED = kxz[index];
    CCTK_REAL kyyL CCTK_ATTRIBUTE_UNUSED = kyy[index];
    CCTK_REAL kyzL CCTK_ATTRIBUTE_UNUSED = kyz[index];
    CCTK_REAL kzzL CCTK_ATTRIBUTE_UNUSED = kzz[index];
    CCTK_REAL SS1L CCTK_ATTRIBUTE_UNUSED = SS1[index];
    CCTK_REAL SS2L CCTK_ATTRIBUTE_UNUSED = SS2[index];
    CCTK_REAL SS3L CCTK_ATTRIBUTE_UNUSED = SS3[index];
    CCTK_REAL tauL CCTK_ATTRIBUTE_UNUSED = tau[index];
    CCTK_REAL V1L CCTK_ATTRIBUTE_UNUSED = V1[index];
    CCTK_REAL V2L CCTK_ATTRIBUTE_UNUSED = V2[index];
    CCTK_REAL V3L CCTK_ATTRIBUTE_UNUSED = V3[index];
    
    
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
    CCTK_REAL PDstandardNth1DD CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2DD CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3DD CCTK_ATTRIBUTE_UNUSED;
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
    CCTK_REAL PDstandardNth1SS1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2SS1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3SS1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth1SS2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2SS2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3SS2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth1SS3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2SS3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3SS3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth1V1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2V1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3V1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth1V2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2V2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3V2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth1V3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2V3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3V3 CCTK_ATTRIBUTE_UNUSED;
    
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
        PDstandardNth1DD = PDstandardNthfdOrder21(&DD[index]);
        PDstandardNth2DD = PDstandardNthfdOrder22(&DD[index]);
        PDstandardNth3DD = PDstandardNthfdOrder23(&DD[index]);
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
        PDstandardNth1SS1 = PDstandardNthfdOrder21(&SS1[index]);
        PDstandardNth2SS1 = PDstandardNthfdOrder22(&SS1[index]);
        PDstandardNth3SS1 = PDstandardNthfdOrder23(&SS1[index]);
        PDstandardNth1SS2 = PDstandardNthfdOrder21(&SS2[index]);
        PDstandardNth2SS2 = PDstandardNthfdOrder22(&SS2[index]);
        PDstandardNth3SS2 = PDstandardNthfdOrder23(&SS2[index]);
        PDstandardNth1SS3 = PDstandardNthfdOrder21(&SS3[index]);
        PDstandardNth2SS3 = PDstandardNthfdOrder22(&SS3[index]);
        PDstandardNth3SS3 = PDstandardNthfdOrder23(&SS3[index]);
        PDstandardNth1V1 = PDstandardNthfdOrder21(&V1[index]);
        PDstandardNth2V1 = PDstandardNthfdOrder22(&V1[index]);
        PDstandardNth3V1 = PDstandardNthfdOrder23(&V1[index]);
        PDstandardNth1V2 = PDstandardNthfdOrder21(&V2[index]);
        PDstandardNth2V2 = PDstandardNthfdOrder22(&V2[index]);
        PDstandardNth3V2 = PDstandardNthfdOrder23(&V2[index]);
        PDstandardNth1V3 = PDstandardNthfdOrder21(&V3[index]);
        PDstandardNth2V3 = PDstandardNthfdOrder22(&V3[index]);
        PDstandardNth3V3 = PDstandardNthfdOrder23(&V3[index]);
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
        PDstandardNth1DD = PDstandardNthfdOrder41(&DD[index]);
        PDstandardNth2DD = PDstandardNthfdOrder42(&DD[index]);
        PDstandardNth3DD = PDstandardNthfdOrder43(&DD[index]);
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
        PDstandardNth1SS1 = PDstandardNthfdOrder41(&SS1[index]);
        PDstandardNth2SS1 = PDstandardNthfdOrder42(&SS1[index]);
        PDstandardNth3SS1 = PDstandardNthfdOrder43(&SS1[index]);
        PDstandardNth1SS2 = PDstandardNthfdOrder41(&SS2[index]);
        PDstandardNth2SS2 = PDstandardNthfdOrder42(&SS2[index]);
        PDstandardNth3SS2 = PDstandardNthfdOrder43(&SS2[index]);
        PDstandardNth1SS3 = PDstandardNthfdOrder41(&SS3[index]);
        PDstandardNth2SS3 = PDstandardNthfdOrder42(&SS3[index]);
        PDstandardNth3SS3 = PDstandardNthfdOrder43(&SS3[index]);
        PDstandardNth1V1 = PDstandardNthfdOrder41(&V1[index]);
        PDstandardNth2V1 = PDstandardNthfdOrder42(&V1[index]);
        PDstandardNth3V1 = PDstandardNthfdOrder43(&V1[index]);
        PDstandardNth1V2 = PDstandardNthfdOrder41(&V2[index]);
        PDstandardNth2V2 = PDstandardNthfdOrder42(&V2[index]);
        PDstandardNth3V2 = PDstandardNthfdOrder43(&V2[index]);
        PDstandardNth1V3 = PDstandardNthfdOrder41(&V3[index]);
        PDstandardNth2V3 = PDstandardNthfdOrder42(&V3[index]);
        PDstandardNth3V3 = PDstandardNthfdOrder43(&V3[index]);
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
        PDstandardNth1DD = PDstandardNthfdOrder61(&DD[index]);
        PDstandardNth2DD = PDstandardNthfdOrder62(&DD[index]);
        PDstandardNth3DD = PDstandardNthfdOrder63(&DD[index]);
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
        PDstandardNth1SS1 = PDstandardNthfdOrder61(&SS1[index]);
        PDstandardNth2SS1 = PDstandardNthfdOrder62(&SS1[index]);
        PDstandardNth3SS1 = PDstandardNthfdOrder63(&SS1[index]);
        PDstandardNth1SS2 = PDstandardNthfdOrder61(&SS2[index]);
        PDstandardNth2SS2 = PDstandardNthfdOrder62(&SS2[index]);
        PDstandardNth3SS2 = PDstandardNthfdOrder63(&SS2[index]);
        PDstandardNth1SS3 = PDstandardNthfdOrder61(&SS3[index]);
        PDstandardNth2SS3 = PDstandardNthfdOrder62(&SS3[index]);
        PDstandardNth3SS3 = PDstandardNthfdOrder63(&SS3[index]);
        PDstandardNth1V1 = PDstandardNthfdOrder61(&V1[index]);
        PDstandardNth2V1 = PDstandardNthfdOrder62(&V1[index]);
        PDstandardNth3V1 = PDstandardNthfdOrder63(&V1[index]);
        PDstandardNth1V2 = PDstandardNthfdOrder61(&V2[index]);
        PDstandardNth2V2 = PDstandardNthfdOrder62(&V2[index]);
        PDstandardNth3V2 = PDstandardNthfdOrder63(&V2[index]);
        PDstandardNth1V3 = PDstandardNthfdOrder61(&V3[index]);
        PDstandardNth2V3 = PDstandardNthfdOrder62(&V3[index]);
        PDstandardNth3V3 = PDstandardNthfdOrder63(&V3[index]);
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
        PDstandardNth1DD = PDstandardNthfdOrder81(&DD[index]);
        PDstandardNth2DD = PDstandardNthfdOrder82(&DD[index]);
        PDstandardNth3DD = PDstandardNthfdOrder83(&DD[index]);
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
        PDstandardNth1SS1 = PDstandardNthfdOrder81(&SS1[index]);
        PDstandardNth2SS1 = PDstandardNthfdOrder82(&SS1[index]);
        PDstandardNth3SS1 = PDstandardNthfdOrder83(&SS1[index]);
        PDstandardNth1SS2 = PDstandardNthfdOrder81(&SS2[index]);
        PDstandardNth2SS2 = PDstandardNthfdOrder82(&SS2[index]);
        PDstandardNth3SS2 = PDstandardNthfdOrder83(&SS2[index]);
        PDstandardNth1SS3 = PDstandardNthfdOrder81(&SS3[index]);
        PDstandardNth2SS3 = PDstandardNthfdOrder82(&SS3[index]);
        PDstandardNth3SS3 = PDstandardNthfdOrder83(&SS3[index]);
        PDstandardNth1V1 = PDstandardNthfdOrder81(&V1[index]);
        PDstandardNth2V1 = PDstandardNthfdOrder82(&V1[index]);
        PDstandardNth3V1 = PDstandardNthfdOrder83(&V1[index]);
        PDstandardNth1V2 = PDstandardNthfdOrder81(&V2[index]);
        PDstandardNth2V2 = PDstandardNthfdOrder82(&V2[index]);
        PDstandardNth3V2 = PDstandardNthfdOrder83(&V2[index]);
        PDstandardNth1V3 = PDstandardNthfdOrder81(&V3[index]);
        PDstandardNth2V3 = PDstandardNthfdOrder82(&V3[index]);
        PDstandardNth3V3 = PDstandardNthfdOrder83(&V3[index]);
        break;
      }
      default:
        CCTK_BUILTIN_UNREACHABLE();
    }
    /* Calculate temporaries and grid functions */
    CCTK_REAL JacPDstandardNth1alp CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1betax CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1betay CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1betaz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1DD CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1gxx CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1gxy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1gxz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1gyy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1gyz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1gzz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1SS1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1SS2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1SS3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth1V1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2alp CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2betax CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2betay CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2betaz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2DD CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2gxx CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2gxy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2gxz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2gyy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2gyz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2gzz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2SS1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2SS2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2SS3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth2V2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3alp CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3betax CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3betay CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3betaz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3DD CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3gxx CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3gxy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3gxz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3gyy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3gyz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3gzz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3SS1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3SS2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3SS3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDstandardNth3V3 CCTK_ATTRIBUTE_UNUSED;
    
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
      
      JacPDstandardNth1DD = J11L*PDstandardNth1DD + J21L*PDstandardNth2DD + 
        J31L*PDstandardNth3DD;
      
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
      
      JacPDstandardNth1SS1 = J11L*PDstandardNth1SS1 + J21L*PDstandardNth2SS1 
        + J31L*PDstandardNth3SS1;
      
      JacPDstandardNth1SS2 = J11L*PDstandardNth1SS2 + J21L*PDstandardNth2SS2 
        + J31L*PDstandardNth3SS2;
      
      JacPDstandardNth1SS3 = J11L*PDstandardNth1SS3 + J21L*PDstandardNth2SS3 
        + J31L*PDstandardNth3SS3;
      
      JacPDstandardNth1V1 = J11L*PDstandardNth1V1 + J21L*PDstandardNth2V1 + 
        J31L*PDstandardNth3V1;
      
      JacPDstandardNth2alp = J12L*PDstandardNth1alp + J22L*PDstandardNth2alp 
        + J32L*PDstandardNth3alp;
      
      JacPDstandardNth2betax = J12L*PDstandardNth1betax + 
        J22L*PDstandardNth2betax + J32L*PDstandardNth3betax;
      
      JacPDstandardNth2betay = J12L*PDstandardNth1betay + 
        J22L*PDstandardNth2betay + J32L*PDstandardNth3betay;
      
      JacPDstandardNth2betaz = J12L*PDstandardNth1betaz + 
        J22L*PDstandardNth2betaz + J32L*PDstandardNth3betaz;
      
      JacPDstandardNth2DD = J12L*PDstandardNth1DD + J22L*PDstandardNth2DD + 
        J32L*PDstandardNth3DD;
      
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
      
      JacPDstandardNth2SS1 = J12L*PDstandardNth1SS1 + J22L*PDstandardNth2SS1 
        + J32L*PDstandardNth3SS1;
      
      JacPDstandardNth2SS2 = J12L*PDstandardNth1SS2 + J22L*PDstandardNth2SS2 
        + J32L*PDstandardNth3SS2;
      
      JacPDstandardNth2SS3 = J12L*PDstandardNth1SS3 + J22L*PDstandardNth2SS3 
        + J32L*PDstandardNth3SS3;
      
      JacPDstandardNth2V2 = J12L*PDstandardNth1V2 + J22L*PDstandardNth2V2 + 
        J32L*PDstandardNth3V2;
      
      JacPDstandardNth3alp = J13L*PDstandardNth1alp + J23L*PDstandardNth2alp 
        + J33L*PDstandardNth3alp;
      
      JacPDstandardNth3betax = J13L*PDstandardNth1betax + 
        J23L*PDstandardNth2betax + J33L*PDstandardNth3betax;
      
      JacPDstandardNth3betay = J13L*PDstandardNth1betay + 
        J23L*PDstandardNth2betay + J33L*PDstandardNth3betay;
      
      JacPDstandardNth3betaz = J13L*PDstandardNth1betaz + 
        J23L*PDstandardNth2betaz + J33L*PDstandardNth3betaz;
      
      JacPDstandardNth3DD = J13L*PDstandardNth1DD + J23L*PDstandardNth2DD + 
        J33L*PDstandardNth3DD;
      
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
      
      JacPDstandardNth3SS1 = J13L*PDstandardNth1SS1 + J23L*PDstandardNth2SS1 
        + J33L*PDstandardNth3SS1;
      
      JacPDstandardNth3SS2 = J13L*PDstandardNth1SS2 + J23L*PDstandardNth2SS2 
        + J33L*PDstandardNth3SS2;
      
      JacPDstandardNth3SS3 = J13L*PDstandardNth1SS3 + J23L*PDstandardNth2SS3 
        + J33L*PDstandardNth3SS3;
      
      JacPDstandardNth3V3 = J13L*PDstandardNth1V3 + J23L*PDstandardNth2V3 + 
        J33L*PDstandardNth3V3;
    }
    else
    {
      JacPDstandardNth1alp = PDstandardNth1alp;
      
      JacPDstandardNth1betax = PDstandardNth1betax;
      
      JacPDstandardNth1betay = PDstandardNth1betay;
      
      JacPDstandardNth1betaz = PDstandardNth1betaz;
      
      JacPDstandardNth1DD = PDstandardNth1DD;
      
      JacPDstandardNth1gxx = PDstandardNth1gxx;
      
      JacPDstandardNth1gxy = PDstandardNth1gxy;
      
      JacPDstandardNth1gxz = PDstandardNth1gxz;
      
      JacPDstandardNth1gyy = PDstandardNth1gyy;
      
      JacPDstandardNth1gyz = PDstandardNth1gyz;
      
      JacPDstandardNth1gzz = PDstandardNth1gzz;
      
      JacPDstandardNth1SS1 = PDstandardNth1SS1;
      
      JacPDstandardNth1SS2 = PDstandardNth1SS2;
      
      JacPDstandardNth1SS3 = PDstandardNth1SS3;
      
      JacPDstandardNth1V1 = PDstandardNth1V1;
      
      JacPDstandardNth2alp = PDstandardNth2alp;
      
      JacPDstandardNth2betax = PDstandardNth2betax;
      
      JacPDstandardNth2betay = PDstandardNth2betay;
      
      JacPDstandardNth2betaz = PDstandardNth2betaz;
      
      JacPDstandardNth2DD = PDstandardNth2DD;
      
      JacPDstandardNth2gxx = PDstandardNth2gxx;
      
      JacPDstandardNth2gxy = PDstandardNth2gxy;
      
      JacPDstandardNth2gxz = PDstandardNth2gxz;
      
      JacPDstandardNth2gyy = PDstandardNth2gyy;
      
      JacPDstandardNth2gyz = PDstandardNth2gyz;
      
      JacPDstandardNth2gzz = PDstandardNth2gzz;
      
      JacPDstandardNth2SS1 = PDstandardNth2SS1;
      
      JacPDstandardNth2SS2 = PDstandardNth2SS2;
      
      JacPDstandardNth2SS3 = PDstandardNth2SS3;
      
      JacPDstandardNth2V2 = PDstandardNth2V2;
      
      JacPDstandardNth3alp = PDstandardNth3alp;
      
      JacPDstandardNth3betax = PDstandardNth3betax;
      
      JacPDstandardNth3betay = PDstandardNth3betay;
      
      JacPDstandardNth3betaz = PDstandardNth3betaz;
      
      JacPDstandardNth3DD = PDstandardNth3DD;
      
      JacPDstandardNth3gxx = PDstandardNth3gxx;
      
      JacPDstandardNth3gxy = PDstandardNth3gxy;
      
      JacPDstandardNth3gxz = PDstandardNth3gxz;
      
      JacPDstandardNth3gyy = PDstandardNth3gyy;
      
      JacPDstandardNth3gyz = PDstandardNth3gyz;
      
      JacPDstandardNth3gzz = PDstandardNth3gzz;
      
      JacPDstandardNth3SS1 = PDstandardNth3SS1;
      
      JacPDstandardNth3SS2 = PDstandardNth3SS2;
      
      JacPDstandardNth3SS3 = PDstandardNth3SS3;
      
      JacPDstandardNth3V3 = PDstandardNth3V3;
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
    
    CCTK_REAL dtgu11 CCTK_ATTRIBUTE_UNUSED = -2*alpL*(2*(kyzL*gu12*gu13 + 
      gu11*(kxyL*gu12 + kxzL*gu13)) + kxxL*pow(gu11,2) + kyyL*pow(gu12,2) + 
      kzzL*pow(gu13,2));
    
    CCTK_REAL dtgu12 CCTK_ATTRIBUTE_UNUSED = -2*alpL*(gu12*(kxxL*gu11 + 
      kxyL*gu12 + kxzL*gu13) + (kxyL*gu11 + kyyL*gu12 + kyzL*gu13)*gu22 + 
      (kxzL*gu11 + kyzL*gu12 + kzzL*gu13)*gu23);
    
    CCTK_REAL dtgu13 CCTK_ATTRIBUTE_UNUSED = -2*alpL*(gu11*(kxxL*gu13 + 
      kxyL*gu23 + kxzL*gu33) + gu12*(kxyL*gu13 + kyyL*gu23 + kyzL*gu33) + 
      gu13*(kxzL*gu13 + kyzL*gu23 + kzzL*gu33));
    
    CCTK_REAL dtgu22 CCTK_ATTRIBUTE_UNUSED = -2*alpL*(2*(kyzL*gu22*gu23 + 
      gu12*(kxyL*gu22 + kxzL*gu23)) + kxxL*pow(gu12,2) + kyyL*pow(gu22,2) + 
      kzzL*pow(gu23,2));
    
    CCTK_REAL dtgu23 CCTK_ATTRIBUTE_UNUSED = -2*alpL*(gu13*(kxxL*gu12 + 
      kxyL*gu22 + kxzL*gu23) + gu23*(kxyL*gu12 + kyyL*gu22 + kyzL*gu23) + 
      (kxzL*gu12 + kyzL*gu22 + kzzL*gu23)*gu33);
    
    CCTK_REAL dtgu33 CCTK_ATTRIBUTE_UNUSED = -2*alpL*(2*(kyzL*gu23*gu33 + 
      gu13*(kxyL*gu23 + kxzL*gu33)) + kxxL*pow(gu13,2) + kyyL*pow(gu23,2) + 
      kzzL*pow(gu33,2));
    
    CCTK_REAL dttau CCTK_ATTRIBUTE_UNUSED = pow(-2*betaxL*(betayL*gxyL + 
      betazL*gxzL) - 2*betayL*betazL*gyzL + pow(alpL,2) - gxxL*pow(betaxL,2) 
      - gyyL*pow(betayL,2) - gzzL*pow(betazL,2),0.5);
    
    CCTK_REAL u1 CCTK_ATTRIBUTE_UNUSED = SS1L*pow(DDL,-1);
    
    CCTK_REAL u2 CCTK_ATTRIBUTE_UNUSED = SS2L*pow(DDL,-1);
    
    CCTK_REAL u3 CCTK_ATTRIBUTE_UNUSED = SS3L*pow(DDL,-1);
    
    CCTK_REAL u0u CCTK_ATTRIBUTE_UNUSED = pow(alpL,-1)*pow(1 + 
      2*(gu23*u2*u3 + u1*(gu12*u2 + gu13*u3)) + gu11*pow(u1,2) + 
      gu22*pow(u2,2) + gu33*pow(u3,2),0.5);
    
    CCTK_REAL u0l CCTK_ATTRIBUTE_UNUSED = betaxL*u1 + betayL*u2 + 
      betazL*u3 - u0u*pow(alpL,2);
    
    CCTK_REAL dtD CCTK_ATTRIBUTE_UNUSED = -(V1L*JacPDstandardNth1DD) - 
      V2L*JacPDstandardNth2DD - V3L*JacPDstandardNth3DD - 
      DDL*(JacPDstandardNth1V1 + JacPDstandardNth2V2 + JacPDstandardNth3V3);
    
    CCTK_REAL SS0u CCTK_ATTRIBUTE_UNUSED = DDL*u0u;
    
    CCTK_REAL SS0l CCTK_ATTRIBUTE_UNUSED = DDL*u0l;
    
    CCTK_REAL SaSbdgab1 CCTK_ATTRIBUTE_UNUSED = 
      -(pow(alpL,-3)*pow(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + 
      pow(gyzL,2)),-2)*(pow(alpL,3)*(-2*SS1L*SS2L*((gzzL*JacPDstandardNth1gxy 
      - gyzL*JacPDstandardNth1gxz - gxzL*JacPDstandardNth1gyz + 
      gxyL*JacPDstandardNth1gzz)*(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + pow(gyzL,2))) + (gxzL*gyzL - 
      gxyL*gzzL)*(-2*(gyzL*(gxzL*JacPDstandardNth1gxy + 
      gxyL*JacPDstandardNth1gxz) + gxyL*gxzL*JacPDstandardNth1gyz) + 
      2*(gxyL*gzzL*JacPDstandardNth1gxy + gxzL*gyyL*JacPDstandardNth1gxz + 
      gxxL*gyzL*JacPDstandardNth1gyz) - gxxL*(gzzL*JacPDstandardNth1gyy + 
      gyyL*JacPDstandardNth1gzz) + JacPDstandardNth1gzz*pow(gxyL,2) + 
      JacPDstandardNth1gyy*pow(gxzL,2) + JacPDstandardNth1gxx*(-(gyyL*gzzL) + 
      pow(gyzL,2)))) + SS3L*(-2*SS2L*((gyzL*JacPDstandardNth1gxx - 
      gxzL*JacPDstandardNth1gxy - gxyL*JacPDstandardNth1gxz + 
      gxxL*JacPDstandardNth1gyz)*(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + pow(gyzL,2))) + (gxyL*gxzL - 
      gxxL*gyzL)*(-2*(gyzL*(gxzL*JacPDstandardNth1gxy + 
      gxyL*JacPDstandardNth1gxz) + gxyL*gxzL*JacPDstandardNth1gyz) + 
      2*(gxyL*gzzL*JacPDstandardNth1gxy + gxzL*gyyL*JacPDstandardNth1gxz + 
      gxxL*gyzL*JacPDstandardNth1gyz) - gxxL*(gzzL*JacPDstandardNth1gyy + 
      gyyL*JacPDstandardNth1gzz) + JacPDstandardNth1gzz*pow(gxyL,2) + 
      JacPDstandardNth1gyy*pow(gxzL,2) + JacPDstandardNth1gxx*(-(gyyL*gzzL) + 
      pow(gyzL,2)))) + 2*SS1L*((gyzL*JacPDstandardNth1gxy - 
      gyyL*JacPDstandardNth1gxz - gxzL*JacPDstandardNth1gyy + 
      gxyL*JacPDstandardNth1gyz)*(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + pow(gyzL,2))) + (gxzL*gyyL - 
      gxyL*gyzL)*(-2*(gyzL*(gxzL*JacPDstandardNth1gxy + 
      gxyL*JacPDstandardNth1gxz) + gxyL*gxzL*JacPDstandardNth1gyz) + 
      2*(gxyL*gzzL*JacPDstandardNth1gxy + gxzL*gyyL*JacPDstandardNth1gxz + 
      gxxL*gyzL*JacPDstandardNth1gyz) - gxxL*(gzzL*JacPDstandardNth1gyy + 
      gyyL*JacPDstandardNth1gzz) + JacPDstandardNth1gzz*pow(gxyL,2) + 
      JacPDstandardNth1gyy*pow(gxzL,2) + JacPDstandardNth1gxx*(-(gyyL*gzzL) + 
      pow(gyzL,2))))) + ((gzzL*JacPDstandardNth1gyy - 
      2*gyzL*JacPDstandardNth1gyz + 
      gyyL*JacPDstandardNth1gzz)*(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + pow(gyzL,2))) + (-(gyyL*gzzL) + 
      pow(gyzL,2))*(-2*(gyzL*(gxzL*JacPDstandardNth1gxy + 
      gxyL*JacPDstandardNth1gxz) + gxyL*gxzL*JacPDstandardNth1gyz) + 
      2*(gxyL*gzzL*JacPDstandardNth1gxy + gxzL*gyyL*JacPDstandardNth1gxz + 
      gxxL*gyzL*JacPDstandardNth1gyz) - gxxL*(gzzL*JacPDstandardNth1gyy + 
      gyyL*JacPDstandardNth1gzz) + JacPDstandardNth1gzz*pow(gxyL,2) + 
      JacPDstandardNth1gyy*pow(gxzL,2) + JacPDstandardNth1gxx*(-(gyyL*gzzL) + 
      pow(gyzL,2))))*pow(SS1L,2) + ((gzzL*JacPDstandardNth1gxx - 
      2*gxzL*JacPDstandardNth1gxz + 
      gxxL*JacPDstandardNth1gzz)*(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + pow(gyzL,2))) + (-(gxxL*gzzL) + 
      pow(gxzL,2))*(-2*(gyzL*(gxzL*JacPDstandardNth1gxy + 
      gxyL*JacPDstandardNth1gxz) + gxyL*gxzL*JacPDstandardNth1gyz) + 
      2*(gxyL*gzzL*JacPDstandardNth1gxy + gxzL*gyyL*JacPDstandardNth1gxz + 
      gxxL*gyzL*JacPDstandardNth1gyz) - gxxL*(gzzL*JacPDstandardNth1gyy + 
      gyyL*JacPDstandardNth1gzz) + JacPDstandardNth1gzz*pow(gxyL,2) + 
      JacPDstandardNth1gyy*pow(gxzL,2) + JacPDstandardNth1gxx*(-(gyyL*gzzL) + 
      pow(gyzL,2))))*pow(SS2L,2) + ((gyyL*JacPDstandardNth1gxx - 
      2*gxyL*JacPDstandardNth1gxy + 
      gxxL*JacPDstandardNth1gyy)*(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + pow(gyzL,2))) + (-(gxxL*gyyL) + 
      pow(gxyL,2))*(-2*(gyzL*(gxzL*JacPDstandardNth1gxy + 
      gxyL*JacPDstandardNth1gxz) + gxyL*gxzL*JacPDstandardNth1gyz) + 
      2*(gxyL*gzzL*JacPDstandardNth1gxy + gxzL*gyyL*JacPDstandardNth1gxz + 
      gxxL*gyzL*JacPDstandardNth1gyz) - gxxL*(gzzL*JacPDstandardNth1gyy + 
      gyyL*JacPDstandardNth1gzz) + JacPDstandardNth1gzz*pow(gxyL,2) + 
      JacPDstandardNth1gyy*pow(gxzL,2) + JacPDstandardNth1gxx*(-(gyyL*gzzL) + 
      pow(gyzL,2))))*pow(SS3L,2)) + (2*(SS1L*(2*betaxL*JacPDstandardNth1alp - 
      alpL*JacPDstandardNth1betax) + SS2L*(2*betayL*JacPDstandardNth1alp - 
      alpL*JacPDstandardNth1betay) + SS3L*(2*betazL*JacPDstandardNth1alp - 
      alpL*JacPDstandardNth1betaz))*SS0l - 
      2*JacPDstandardNth1alp*pow(SS0l,2))*pow(-2*gxyL*gxzL*gyzL + 
      gzzL*pow(gxyL,2) + gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + 
      pow(gyzL,2)),2)));
    
    CCTK_REAL SaSbdgab2 CCTK_ATTRIBUTE_UNUSED = 
      -(pow(alpL,-3)*pow(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + 
      pow(gyzL,2)),-2)*(pow(alpL,3)*(-2*SS1L*SS2L*((gzzL*JacPDstandardNth2gxy 
      - gyzL*JacPDstandardNth2gxz - gxzL*JacPDstandardNth2gyz + 
      gxyL*JacPDstandardNth2gzz)*(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + pow(gyzL,2))) + (gxzL*gyzL - 
      gxyL*gzzL)*(-2*(gyzL*(gxzL*JacPDstandardNth2gxy + 
      gxyL*JacPDstandardNth2gxz) + gxyL*gxzL*JacPDstandardNth2gyz) + 
      2*(gxyL*gzzL*JacPDstandardNth2gxy + gxzL*gyyL*JacPDstandardNth2gxz + 
      gxxL*gyzL*JacPDstandardNth2gyz) - gxxL*(gzzL*JacPDstandardNth2gyy + 
      gyyL*JacPDstandardNth2gzz) + JacPDstandardNth2gzz*pow(gxyL,2) + 
      JacPDstandardNth2gyy*pow(gxzL,2) + JacPDstandardNth2gxx*(-(gyyL*gzzL) + 
      pow(gyzL,2)))) + SS3L*(-2*SS2L*((gyzL*JacPDstandardNth2gxx - 
      gxzL*JacPDstandardNth2gxy - gxyL*JacPDstandardNth2gxz + 
      gxxL*JacPDstandardNth2gyz)*(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + pow(gyzL,2))) + (gxyL*gxzL - 
      gxxL*gyzL)*(-2*(gyzL*(gxzL*JacPDstandardNth2gxy + 
      gxyL*JacPDstandardNth2gxz) + gxyL*gxzL*JacPDstandardNth2gyz) + 
      2*(gxyL*gzzL*JacPDstandardNth2gxy + gxzL*gyyL*JacPDstandardNth2gxz + 
      gxxL*gyzL*JacPDstandardNth2gyz) - gxxL*(gzzL*JacPDstandardNth2gyy + 
      gyyL*JacPDstandardNth2gzz) + JacPDstandardNth2gzz*pow(gxyL,2) + 
      JacPDstandardNth2gyy*pow(gxzL,2) + JacPDstandardNth2gxx*(-(gyyL*gzzL) + 
      pow(gyzL,2)))) + 2*SS1L*((gyzL*JacPDstandardNth2gxy - 
      gyyL*JacPDstandardNth2gxz - gxzL*JacPDstandardNth2gyy + 
      gxyL*JacPDstandardNth2gyz)*(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + pow(gyzL,2))) + (gxzL*gyyL - 
      gxyL*gyzL)*(-2*(gyzL*(gxzL*JacPDstandardNth2gxy + 
      gxyL*JacPDstandardNth2gxz) + gxyL*gxzL*JacPDstandardNth2gyz) + 
      2*(gxyL*gzzL*JacPDstandardNth2gxy + gxzL*gyyL*JacPDstandardNth2gxz + 
      gxxL*gyzL*JacPDstandardNth2gyz) - gxxL*(gzzL*JacPDstandardNth2gyy + 
      gyyL*JacPDstandardNth2gzz) + JacPDstandardNth2gzz*pow(gxyL,2) + 
      JacPDstandardNth2gyy*pow(gxzL,2) + JacPDstandardNth2gxx*(-(gyyL*gzzL) + 
      pow(gyzL,2))))) + ((gzzL*JacPDstandardNth2gyy - 
      2*gyzL*JacPDstandardNth2gyz + 
      gyyL*JacPDstandardNth2gzz)*(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + pow(gyzL,2))) + (-(gyyL*gzzL) + 
      pow(gyzL,2))*(-2*(gyzL*(gxzL*JacPDstandardNth2gxy + 
      gxyL*JacPDstandardNth2gxz) + gxyL*gxzL*JacPDstandardNth2gyz) + 
      2*(gxyL*gzzL*JacPDstandardNth2gxy + gxzL*gyyL*JacPDstandardNth2gxz + 
      gxxL*gyzL*JacPDstandardNth2gyz) - gxxL*(gzzL*JacPDstandardNth2gyy + 
      gyyL*JacPDstandardNth2gzz) + JacPDstandardNth2gzz*pow(gxyL,2) + 
      JacPDstandardNth2gyy*pow(gxzL,2) + JacPDstandardNth2gxx*(-(gyyL*gzzL) + 
      pow(gyzL,2))))*pow(SS1L,2) + ((gzzL*JacPDstandardNth2gxx - 
      2*gxzL*JacPDstandardNth2gxz + 
      gxxL*JacPDstandardNth2gzz)*(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + pow(gyzL,2))) + (-(gxxL*gzzL) + 
      pow(gxzL,2))*(-2*(gyzL*(gxzL*JacPDstandardNth2gxy + 
      gxyL*JacPDstandardNth2gxz) + gxyL*gxzL*JacPDstandardNth2gyz) + 
      2*(gxyL*gzzL*JacPDstandardNth2gxy + gxzL*gyyL*JacPDstandardNth2gxz + 
      gxxL*gyzL*JacPDstandardNth2gyz) - gxxL*(gzzL*JacPDstandardNth2gyy + 
      gyyL*JacPDstandardNth2gzz) + JacPDstandardNth2gzz*pow(gxyL,2) + 
      JacPDstandardNth2gyy*pow(gxzL,2) + JacPDstandardNth2gxx*(-(gyyL*gzzL) + 
      pow(gyzL,2))))*pow(SS2L,2) + ((gyyL*JacPDstandardNth2gxx - 
      2*gxyL*JacPDstandardNth2gxy + 
      gxxL*JacPDstandardNth2gyy)*(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + pow(gyzL,2))) + (-(gxxL*gyyL) + 
      pow(gxyL,2))*(-2*(gyzL*(gxzL*JacPDstandardNth2gxy + 
      gxyL*JacPDstandardNth2gxz) + gxyL*gxzL*JacPDstandardNth2gyz) + 
      2*(gxyL*gzzL*JacPDstandardNth2gxy + gxzL*gyyL*JacPDstandardNth2gxz + 
      gxxL*gyzL*JacPDstandardNth2gyz) - gxxL*(gzzL*JacPDstandardNth2gyy + 
      gyyL*JacPDstandardNth2gzz) + JacPDstandardNth2gzz*pow(gxyL,2) + 
      JacPDstandardNth2gyy*pow(gxzL,2) + JacPDstandardNth2gxx*(-(gyyL*gzzL) + 
      pow(gyzL,2))))*pow(SS3L,2)) + (2*(SS1L*(2*betaxL*JacPDstandardNth2alp - 
      alpL*JacPDstandardNth2betax) + SS2L*(2*betayL*JacPDstandardNth2alp - 
      alpL*JacPDstandardNth2betay) + SS3L*(2*betazL*JacPDstandardNth2alp - 
      alpL*JacPDstandardNth2betaz))*SS0l - 
      2*JacPDstandardNth2alp*pow(SS0l,2))*pow(-2*gxyL*gxzL*gyzL + 
      gzzL*pow(gxyL,2) + gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + 
      pow(gyzL,2)),2)));
    
    CCTK_REAL SaSbdgab3 CCTK_ATTRIBUTE_UNUSED = 
      -(pow(alpL,-3)*pow(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + 
      pow(gyzL,2)),-2)*(pow(alpL,3)*(-2*SS1L*SS2L*((gzzL*JacPDstandardNth3gxy 
      - gyzL*JacPDstandardNth3gxz - gxzL*JacPDstandardNth3gyz + 
      gxyL*JacPDstandardNth3gzz)*(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + pow(gyzL,2))) + (gxzL*gyzL - 
      gxyL*gzzL)*(-2*(gyzL*(gxzL*JacPDstandardNth3gxy + 
      gxyL*JacPDstandardNth3gxz) + gxyL*gxzL*JacPDstandardNth3gyz) + 
      2*(gxyL*gzzL*JacPDstandardNth3gxy + gxzL*gyyL*JacPDstandardNth3gxz + 
      gxxL*gyzL*JacPDstandardNth3gyz) - gxxL*(gzzL*JacPDstandardNth3gyy + 
      gyyL*JacPDstandardNth3gzz) + JacPDstandardNth3gzz*pow(gxyL,2) + 
      JacPDstandardNth3gyy*pow(gxzL,2) + JacPDstandardNth3gxx*(-(gyyL*gzzL) + 
      pow(gyzL,2)))) + SS3L*(-2*SS2L*((gyzL*JacPDstandardNth3gxx - 
      gxzL*JacPDstandardNth3gxy - gxyL*JacPDstandardNth3gxz + 
      gxxL*JacPDstandardNth3gyz)*(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + pow(gyzL,2))) + (gxyL*gxzL - 
      gxxL*gyzL)*(-2*(gyzL*(gxzL*JacPDstandardNth3gxy + 
      gxyL*JacPDstandardNth3gxz) + gxyL*gxzL*JacPDstandardNth3gyz) + 
      2*(gxyL*gzzL*JacPDstandardNth3gxy + gxzL*gyyL*JacPDstandardNth3gxz + 
      gxxL*gyzL*JacPDstandardNth3gyz) - gxxL*(gzzL*JacPDstandardNth3gyy + 
      gyyL*JacPDstandardNth3gzz) + JacPDstandardNth3gzz*pow(gxyL,2) + 
      JacPDstandardNth3gyy*pow(gxzL,2) + JacPDstandardNth3gxx*(-(gyyL*gzzL) + 
      pow(gyzL,2)))) + 2*SS1L*((gyzL*JacPDstandardNth3gxy - 
      gyyL*JacPDstandardNth3gxz - gxzL*JacPDstandardNth3gyy + 
      gxyL*JacPDstandardNth3gyz)*(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + pow(gyzL,2))) + (gxzL*gyyL - 
      gxyL*gyzL)*(-2*(gyzL*(gxzL*JacPDstandardNth3gxy + 
      gxyL*JacPDstandardNth3gxz) + gxyL*gxzL*JacPDstandardNth3gyz) + 
      2*(gxyL*gzzL*JacPDstandardNth3gxy + gxzL*gyyL*JacPDstandardNth3gxz + 
      gxxL*gyzL*JacPDstandardNth3gyz) - gxxL*(gzzL*JacPDstandardNth3gyy + 
      gyyL*JacPDstandardNth3gzz) + JacPDstandardNth3gzz*pow(gxyL,2) + 
      JacPDstandardNth3gyy*pow(gxzL,2) + JacPDstandardNth3gxx*(-(gyyL*gzzL) + 
      pow(gyzL,2))))) + ((gzzL*JacPDstandardNth3gyy - 
      2*gyzL*JacPDstandardNth3gyz + 
      gyyL*JacPDstandardNth3gzz)*(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + pow(gyzL,2))) + (-(gyyL*gzzL) + 
      pow(gyzL,2))*(-2*(gyzL*(gxzL*JacPDstandardNth3gxy + 
      gxyL*JacPDstandardNth3gxz) + gxyL*gxzL*JacPDstandardNth3gyz) + 
      2*(gxyL*gzzL*JacPDstandardNth3gxy + gxzL*gyyL*JacPDstandardNth3gxz + 
      gxxL*gyzL*JacPDstandardNth3gyz) - gxxL*(gzzL*JacPDstandardNth3gyy + 
      gyyL*JacPDstandardNth3gzz) + JacPDstandardNth3gzz*pow(gxyL,2) + 
      JacPDstandardNth3gyy*pow(gxzL,2) + JacPDstandardNth3gxx*(-(gyyL*gzzL) + 
      pow(gyzL,2))))*pow(SS1L,2) + ((gzzL*JacPDstandardNth3gxx - 
      2*gxzL*JacPDstandardNth3gxz + 
      gxxL*JacPDstandardNth3gzz)*(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + pow(gyzL,2))) + (-(gxxL*gzzL) + 
      pow(gxzL,2))*(-2*(gyzL*(gxzL*JacPDstandardNth3gxy + 
      gxyL*JacPDstandardNth3gxz) + gxyL*gxzL*JacPDstandardNth3gyz) + 
      2*(gxyL*gzzL*JacPDstandardNth3gxy + gxzL*gyyL*JacPDstandardNth3gxz + 
      gxxL*gyzL*JacPDstandardNth3gyz) - gxxL*(gzzL*JacPDstandardNth3gyy + 
      gyyL*JacPDstandardNth3gzz) + JacPDstandardNth3gzz*pow(gxyL,2) + 
      JacPDstandardNth3gyy*pow(gxzL,2) + JacPDstandardNth3gxx*(-(gyyL*gzzL) + 
      pow(gyzL,2))))*pow(SS2L,2) + ((gyyL*JacPDstandardNth3gxx - 
      2*gxyL*JacPDstandardNth3gxy + 
      gxxL*JacPDstandardNth3gyy)*(-2*gxyL*gxzL*gyzL + gzzL*pow(gxyL,2) + 
      gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + pow(gyzL,2))) + (-(gxxL*gyyL) + 
      pow(gxyL,2))*(-2*(gyzL*(gxzL*JacPDstandardNth3gxy + 
      gxyL*JacPDstandardNth3gxz) + gxyL*gxzL*JacPDstandardNth3gyz) + 
      2*(gxyL*gzzL*JacPDstandardNth3gxy + gxzL*gyyL*JacPDstandardNth3gxz + 
      gxxL*gyzL*JacPDstandardNth3gyz) - gxxL*(gzzL*JacPDstandardNth3gyy + 
      gyyL*JacPDstandardNth3gzz) + JacPDstandardNth3gzz*pow(gxyL,2) + 
      JacPDstandardNth3gyy*pow(gxzL,2) + JacPDstandardNth3gxx*(-(gyyL*gzzL) + 
      pow(gyzL,2))))*pow(SS3L,2)) + (2*(SS1L*(2*betaxL*JacPDstandardNth3alp - 
      alpL*JacPDstandardNth3betax) + SS2L*(2*betayL*JacPDstandardNth3alp - 
      alpL*JacPDstandardNth3betay) + SS3L*(2*betazL*JacPDstandardNth3alp - 
      alpL*JacPDstandardNth3betaz))*SS0l - 
      2*JacPDstandardNth3alp*pow(SS0l,2))*pow(-2*gxyL*gxzL*gyzL + 
      gzzL*pow(gxyL,2) + gyyL*pow(gxzL,2) + gxxL*(-(gyyL*gzzL) + 
      pow(gyzL,2)),2)));
    
    CCTK_REAL dtS1 CCTK_ATTRIBUTE_UNUSED = -(V1L*JacPDstandardNth1SS1) - 
      V2L*JacPDstandardNth2SS1 - V3L*JacPDstandardNth3SS1 - 
      SS1L*(JacPDstandardNth1V1 + JacPDstandardNth2V2 + JacPDstandardNth3V3) 
      - 0.5*SaSbdgab1*pow(SS0u,-1);
    
    CCTK_REAL dtS2 CCTK_ATTRIBUTE_UNUSED = -(V1L*JacPDstandardNth1SS2) - 
      V2L*JacPDstandardNth2SS2 - V3L*JacPDstandardNth3SS2 - 
      SS2L*(JacPDstandardNth1V1 + JacPDstandardNth2V2 + JacPDstandardNth3V3) 
      - 0.5*SaSbdgab2*pow(SS0u,-1);
    
    CCTK_REAL dtS3 CCTK_ATTRIBUTE_UNUSED = -(V1L*JacPDstandardNth1SS3) - 
      V2L*JacPDstandardNth2SS3 - V3L*JacPDstandardNth3SS3 - 
      SS3L*(JacPDstandardNth1V1 + JacPDstandardNth2V2 + JacPDstandardNth3V3) 
      - 0.5*SaSbdgab3*pow(SS0u,-1);
    
    CCTK_REAL Snorm CCTK_ATTRIBUTE_UNUSED = 2*(SS1L*(SS2L*gu12 + 
      SS3L*gu13) + SS2L*SS3L*gu23) + gu11*pow(SS1L,2) + gu22*pow(SS2L,2) + 
      gu33*pow(SS3L,2);
    
    CCTK_REAL dtSnorm CCTK_ATTRIBUTE_UNUSED = 2*(dtS1*(SS1L*gu11 + 
      SS2L*gu12 + SS3L*gu13) + SS1L*(SS2L*dtgu12 + SS3L*dtgu13 + dtS3*gu13) + 
      dtS2*(SS1L*gu12 + SS2L*gu22 + SS3L*gu23) + SS2L*(SS3L*dtgu23 + 
      dtS3*gu23) + SS3L*dtS3*gu33) + dtgu11*pow(SS1L,2) + dtgu22*pow(SS2L,2) 
      + dtgu33*pow(SS3L,2);
    
    CCTK_REAL dtF CCTK_ATTRIBUTE_UNUSED = 2*DDL*dtD + dtSnorm;
    
    CCTK_REAL F CCTK_ATTRIBUTE_UNUSED = Snorm + pow(DDL,2);
    
    CCTK_REAL partialtalpha CCTK_ATTRIBUTE_UNUSED = 
      -(alpL*pow(alphaN,-1)*(alpL*KtraceL + 2*(dtD*pow(DDL,-1) + 
      dttau*pow(tauL,-1)) - 0.5*dtF*pow(F,-1)));
    
    CCTK_REAL dtalpL CCTK_ATTRIBUTE_UNUSED = partialtalpha;
    
    CCTK_REAL alpharhsL CCTK_ATTRIBUTE_UNUSED = partialtalpha;
    /* Copy local copies back to grid functions */
    alpharhs[index] = alpharhsL;
    dtalp[index] = dtalpL;
  }
  CCTK_ENDLOOP3(CosmoLapse_dtLapseRhoPostRHS);
}
extern "C" void CosmoLapse_dtLapseRhoPostRHS(CCTK_ARGUMENTS)
{
  #ifdef DECLARE_CCTK_ARGUMENTS_CosmoLapse_dtLapseRhoPostRHS
  DECLARE_CCTK_ARGUMENTS_CHECKED(CosmoLapse_dtLapseRhoPostRHS);
  #else
  DECLARE_CCTK_ARGUMENTS;
  #endif
  DECLARE_CCTK_PARAMETERS;
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering CosmoLapse_dtLapseRhoPostRHS_Body");
  }
  if (cctk_iteration % CosmoLapse_dtLapseRhoPostRHS_calc_every != CosmoLapse_dtLapseRhoPostRHS_calc_offset)
  {
    return;
  }
  
  const char* const groups[] = {
    "ADMBase::curv",
    "ADMBase::dtlapse",
    "ADMBase::lapse",
    "ADMBase::metric",
    "ADMBase::shift",
    "CosmoLapse::CL_Ktrace",
    "CT_Dust::CT_D",
    "CT_Dust::CT_S",
    "CT_Dust::CT_V",
    "ML_BSSN::ML_lapserhs",
    "CosmoLapse::propertime"};
  AssertGroupStorage(cctkGH, "CosmoLapse_dtLapseRhoPostRHS", 11, groups);
  
  switch (fdOrder)
  {
    case 2:
    {
      EnsureStencilFits(cctkGH, "CosmoLapse_dtLapseRhoPostRHS", 1, 1, 1);
      break;
    }
    
    case 4:
    {
      EnsureStencilFits(cctkGH, "CosmoLapse_dtLapseRhoPostRHS", 2, 2, 2);
      break;
    }
    
    case 6:
    {
      EnsureStencilFits(cctkGH, "CosmoLapse_dtLapseRhoPostRHS", 3, 3, 3);
      break;
    }
    
    case 8:
    {
      EnsureStencilFits(cctkGH, "CosmoLapse_dtLapseRhoPostRHS", 4, 4, 4);
      break;
    }
    default:
      CCTK_BUILTIN_UNREACHABLE();
  }
  
  LoopOverInterior(cctkGH, CosmoLapse_dtLapseRhoPostRHS_Body);
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Leaving CosmoLapse_dtLapseRhoPostRHS_Body");
  }
}

} // namespace CosmoLapse
