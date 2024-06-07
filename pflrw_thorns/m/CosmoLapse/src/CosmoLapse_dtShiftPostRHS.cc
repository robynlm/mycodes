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

extern "C" void CosmoLapse_dtShiftPostRHS_SelectBCs(CCTK_ARGUMENTS)
{
  #ifdef DECLARE_CCTK_ARGUMENTS_CosmoLapse_dtShiftPostRHS_SelectBCs
  DECLARE_CCTK_ARGUMENTS_CHECKED(CosmoLapse_dtShiftPostRHS_SelectBCs);
  #else
  DECLARE_CCTK_ARGUMENTS;
  #endif
  DECLARE_CCTK_PARAMETERS;
  
  if (cctk_iteration % CosmoLapse_dtShiftPostRHS_calc_every != CosmoLapse_dtShiftPostRHS_calc_offset)
    return;
  CCTK_INT ierr CCTK_ATTRIBUTE_UNUSED = 0;
  ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "ADMBase::dtshift","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for ADMBase::dtshift.");
  ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN::ML_dtshiftrhs","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for ML_BSSN::ML_dtshiftrhs.");
  ierr = KrancBdy_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN::ML_shiftrhs","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for ML_BSSN::ML_shiftrhs.");
  return;
}

static void CosmoLapse_dtShiftPostRHS_Body(const cGH* restrict const cctkGH, const int dir, const int face, const CCTK_REAL normal[3], const CCTK_REAL tangentA[3], const CCTK_REAL tangentB[3], const int imin[3], const int imax[3], const int n_subblock_gfs, CCTK_REAL* restrict const subblock_gfs[])
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
  CCTK_LOOP3(CosmoLapse_dtShiftPostRHS,
    i,j,k, imin0,imin1,imin2, imax0,imax1,imax2,
    cctk_ash[0],cctk_ash[1],cctk_ash[2])
  {
    const ptrdiff_t index CCTK_ATTRIBUTE_UNUSED = di*i + dj*j + dk*k;
    /* Assign local copies of grid functions */
    
    CCTK_REAL alpL CCTK_ATTRIBUTE_UNUSED = alp[index];
    CCTK_REAL B1L CCTK_ATTRIBUTE_UNUSED = B1[index];
    CCTK_REAL B2L CCTK_ATTRIBUTE_UNUSED = B2[index];
    CCTK_REAL B3L CCTK_ATTRIBUTE_UNUSED = B3[index];
    CCTK_REAL betaxL CCTK_ATTRIBUTE_UNUSED = betax[index];
    CCTK_REAL betayL CCTK_ATTRIBUTE_UNUSED = betay[index];
    CCTK_REAL betazL CCTK_ATTRIBUTE_UNUSED = betaz[index];
    CCTK_REAL Xt1L CCTK_ATTRIBUTE_UNUSED = Xt1[index];
    CCTK_REAL Xt1rhsL CCTK_ATTRIBUTE_UNUSED = Xt1rhs[index];
    CCTK_REAL Xt2L CCTK_ATTRIBUTE_UNUSED = Xt2[index];
    CCTK_REAL Xt2rhsL CCTK_ATTRIBUTE_UNUSED = Xt2rhs[index];
    CCTK_REAL Xt3L CCTK_ATTRIBUTE_UNUSED = Xt3[index];
    CCTK_REAL Xt3rhsL CCTK_ATTRIBUTE_UNUSED = Xt3rhs[index];
    
    
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
    CCTK_REAL PDupwindNthSymm1B1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDupwindNthSymm2B1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDupwindNthSymm3B1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDupwindNthAnti1B1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDupwindNthAnti2B1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDupwindNthAnti3B1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDupwindNthSymm1B2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDupwindNthSymm2B2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDupwindNthSymm3B2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDupwindNthAnti1B2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDupwindNthAnti2B2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDupwindNthAnti3B2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDupwindNthSymm1B3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDupwindNthSymm2B3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDupwindNthSymm3B3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDupwindNthAnti1B3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDupwindNthAnti2B3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDupwindNthAnti3B3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDupwindNthSymm1betax CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDupwindNthSymm2betax CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDupwindNthSymm3betax CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDupwindNthAnti1betax CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDupwindNthAnti2betax CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDupwindNthAnti3betax CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDupwindNthSymm1betay CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDupwindNthSymm2betay CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDupwindNthSymm3betay CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDupwindNthAnti1betay CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDupwindNthAnti2betay CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDupwindNthAnti3betay CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDupwindNthSymm1betaz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDupwindNthSymm2betaz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDupwindNthSymm3betaz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDupwindNthAnti1betaz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDupwindNthAnti2betaz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDupwindNthAnti3betaz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDupwindNthSymm1Xt1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDupwindNthSymm2Xt1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDupwindNthSymm3Xt1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDupwindNthAnti1Xt1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDupwindNthAnti2Xt1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDupwindNthAnti3Xt1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDupwindNthSymm1Xt2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDupwindNthSymm2Xt2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDupwindNthSymm3Xt2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDupwindNthAnti1Xt2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDupwindNthAnti2Xt2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDupwindNthAnti3Xt2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDupwindNthSymm1Xt3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDupwindNthSymm2Xt3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDupwindNthSymm3Xt3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDupwindNthAnti1Xt3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDupwindNthAnti2Xt3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDupwindNthAnti3Xt3 CCTK_ATTRIBUTE_UNUSED;
    
    switch (fdOrder)
    {
      case 2:
      {
        PDupwindNthSymm1B1 = PDupwindNthSymmfdOrder21(&B1[index]);
        PDupwindNthSymm2B1 = PDupwindNthSymmfdOrder22(&B1[index]);
        PDupwindNthSymm3B1 = PDupwindNthSymmfdOrder23(&B1[index]);
        PDupwindNthAnti1B1 = PDupwindNthAntifdOrder21(&B1[index]);
        PDupwindNthAnti2B1 = PDupwindNthAntifdOrder22(&B1[index]);
        PDupwindNthAnti3B1 = PDupwindNthAntifdOrder23(&B1[index]);
        PDupwindNthSymm1B2 = PDupwindNthSymmfdOrder21(&B2[index]);
        PDupwindNthSymm2B2 = PDupwindNthSymmfdOrder22(&B2[index]);
        PDupwindNthSymm3B2 = PDupwindNthSymmfdOrder23(&B2[index]);
        PDupwindNthAnti1B2 = PDupwindNthAntifdOrder21(&B2[index]);
        PDupwindNthAnti2B2 = PDupwindNthAntifdOrder22(&B2[index]);
        PDupwindNthAnti3B2 = PDupwindNthAntifdOrder23(&B2[index]);
        PDupwindNthSymm1B3 = PDupwindNthSymmfdOrder21(&B3[index]);
        PDupwindNthSymm2B3 = PDupwindNthSymmfdOrder22(&B3[index]);
        PDupwindNthSymm3B3 = PDupwindNthSymmfdOrder23(&B3[index]);
        PDupwindNthAnti1B3 = PDupwindNthAntifdOrder21(&B3[index]);
        PDupwindNthAnti2B3 = PDupwindNthAntifdOrder22(&B3[index]);
        PDupwindNthAnti3B3 = PDupwindNthAntifdOrder23(&B3[index]);
        PDupwindNthSymm1betax = PDupwindNthSymmfdOrder21(&betax[index]);
        PDupwindNthSymm2betax = PDupwindNthSymmfdOrder22(&betax[index]);
        PDupwindNthSymm3betax = PDupwindNthSymmfdOrder23(&betax[index]);
        PDupwindNthAnti1betax = PDupwindNthAntifdOrder21(&betax[index]);
        PDupwindNthAnti2betax = PDupwindNthAntifdOrder22(&betax[index]);
        PDupwindNthAnti3betax = PDupwindNthAntifdOrder23(&betax[index]);
        PDupwindNthSymm1betay = PDupwindNthSymmfdOrder21(&betay[index]);
        PDupwindNthSymm2betay = PDupwindNthSymmfdOrder22(&betay[index]);
        PDupwindNthSymm3betay = PDupwindNthSymmfdOrder23(&betay[index]);
        PDupwindNthAnti1betay = PDupwindNthAntifdOrder21(&betay[index]);
        PDupwindNthAnti2betay = PDupwindNthAntifdOrder22(&betay[index]);
        PDupwindNthAnti3betay = PDupwindNthAntifdOrder23(&betay[index]);
        PDupwindNthSymm1betaz = PDupwindNthSymmfdOrder21(&betaz[index]);
        PDupwindNthSymm2betaz = PDupwindNthSymmfdOrder22(&betaz[index]);
        PDupwindNthSymm3betaz = PDupwindNthSymmfdOrder23(&betaz[index]);
        PDupwindNthAnti1betaz = PDupwindNthAntifdOrder21(&betaz[index]);
        PDupwindNthAnti2betaz = PDupwindNthAntifdOrder22(&betaz[index]);
        PDupwindNthAnti3betaz = PDupwindNthAntifdOrder23(&betaz[index]);
        PDupwindNthSymm1Xt1 = PDupwindNthSymmfdOrder21(&Xt1[index]);
        PDupwindNthSymm2Xt1 = PDupwindNthSymmfdOrder22(&Xt1[index]);
        PDupwindNthSymm3Xt1 = PDupwindNthSymmfdOrder23(&Xt1[index]);
        PDupwindNthAnti1Xt1 = PDupwindNthAntifdOrder21(&Xt1[index]);
        PDupwindNthAnti2Xt1 = PDupwindNthAntifdOrder22(&Xt1[index]);
        PDupwindNthAnti3Xt1 = PDupwindNthAntifdOrder23(&Xt1[index]);
        PDupwindNthSymm1Xt2 = PDupwindNthSymmfdOrder21(&Xt2[index]);
        PDupwindNthSymm2Xt2 = PDupwindNthSymmfdOrder22(&Xt2[index]);
        PDupwindNthSymm3Xt2 = PDupwindNthSymmfdOrder23(&Xt2[index]);
        PDupwindNthAnti1Xt2 = PDupwindNthAntifdOrder21(&Xt2[index]);
        PDupwindNthAnti2Xt2 = PDupwindNthAntifdOrder22(&Xt2[index]);
        PDupwindNthAnti3Xt2 = PDupwindNthAntifdOrder23(&Xt2[index]);
        PDupwindNthSymm1Xt3 = PDupwindNthSymmfdOrder21(&Xt3[index]);
        PDupwindNthSymm2Xt3 = PDupwindNthSymmfdOrder22(&Xt3[index]);
        PDupwindNthSymm3Xt3 = PDupwindNthSymmfdOrder23(&Xt3[index]);
        PDupwindNthAnti1Xt3 = PDupwindNthAntifdOrder21(&Xt3[index]);
        PDupwindNthAnti2Xt3 = PDupwindNthAntifdOrder22(&Xt3[index]);
        PDupwindNthAnti3Xt3 = PDupwindNthAntifdOrder23(&Xt3[index]);
        break;
      }
      
      case 4:
      {
        PDupwindNthSymm1B1 = PDupwindNthSymmfdOrder41(&B1[index]);
        PDupwindNthSymm2B1 = PDupwindNthSymmfdOrder42(&B1[index]);
        PDupwindNthSymm3B1 = PDupwindNthSymmfdOrder43(&B1[index]);
        PDupwindNthAnti1B1 = PDupwindNthAntifdOrder41(&B1[index]);
        PDupwindNthAnti2B1 = PDupwindNthAntifdOrder42(&B1[index]);
        PDupwindNthAnti3B1 = PDupwindNthAntifdOrder43(&B1[index]);
        PDupwindNthSymm1B2 = PDupwindNthSymmfdOrder41(&B2[index]);
        PDupwindNthSymm2B2 = PDupwindNthSymmfdOrder42(&B2[index]);
        PDupwindNthSymm3B2 = PDupwindNthSymmfdOrder43(&B2[index]);
        PDupwindNthAnti1B2 = PDupwindNthAntifdOrder41(&B2[index]);
        PDupwindNthAnti2B2 = PDupwindNthAntifdOrder42(&B2[index]);
        PDupwindNthAnti3B2 = PDupwindNthAntifdOrder43(&B2[index]);
        PDupwindNthSymm1B3 = PDupwindNthSymmfdOrder41(&B3[index]);
        PDupwindNthSymm2B3 = PDupwindNthSymmfdOrder42(&B3[index]);
        PDupwindNthSymm3B3 = PDupwindNthSymmfdOrder43(&B3[index]);
        PDupwindNthAnti1B3 = PDupwindNthAntifdOrder41(&B3[index]);
        PDupwindNthAnti2B3 = PDupwindNthAntifdOrder42(&B3[index]);
        PDupwindNthAnti3B3 = PDupwindNthAntifdOrder43(&B3[index]);
        PDupwindNthSymm1betax = PDupwindNthSymmfdOrder41(&betax[index]);
        PDupwindNthSymm2betax = PDupwindNthSymmfdOrder42(&betax[index]);
        PDupwindNthSymm3betax = PDupwindNthSymmfdOrder43(&betax[index]);
        PDupwindNthAnti1betax = PDupwindNthAntifdOrder41(&betax[index]);
        PDupwindNthAnti2betax = PDupwindNthAntifdOrder42(&betax[index]);
        PDupwindNthAnti3betax = PDupwindNthAntifdOrder43(&betax[index]);
        PDupwindNthSymm1betay = PDupwindNthSymmfdOrder41(&betay[index]);
        PDupwindNthSymm2betay = PDupwindNthSymmfdOrder42(&betay[index]);
        PDupwindNthSymm3betay = PDupwindNthSymmfdOrder43(&betay[index]);
        PDupwindNthAnti1betay = PDupwindNthAntifdOrder41(&betay[index]);
        PDupwindNthAnti2betay = PDupwindNthAntifdOrder42(&betay[index]);
        PDupwindNthAnti3betay = PDupwindNthAntifdOrder43(&betay[index]);
        PDupwindNthSymm1betaz = PDupwindNthSymmfdOrder41(&betaz[index]);
        PDupwindNthSymm2betaz = PDupwindNthSymmfdOrder42(&betaz[index]);
        PDupwindNthSymm3betaz = PDupwindNthSymmfdOrder43(&betaz[index]);
        PDupwindNthAnti1betaz = PDupwindNthAntifdOrder41(&betaz[index]);
        PDupwindNthAnti2betaz = PDupwindNthAntifdOrder42(&betaz[index]);
        PDupwindNthAnti3betaz = PDupwindNthAntifdOrder43(&betaz[index]);
        PDupwindNthSymm1Xt1 = PDupwindNthSymmfdOrder41(&Xt1[index]);
        PDupwindNthSymm2Xt1 = PDupwindNthSymmfdOrder42(&Xt1[index]);
        PDupwindNthSymm3Xt1 = PDupwindNthSymmfdOrder43(&Xt1[index]);
        PDupwindNthAnti1Xt1 = PDupwindNthAntifdOrder41(&Xt1[index]);
        PDupwindNthAnti2Xt1 = PDupwindNthAntifdOrder42(&Xt1[index]);
        PDupwindNthAnti3Xt1 = PDupwindNthAntifdOrder43(&Xt1[index]);
        PDupwindNthSymm1Xt2 = PDupwindNthSymmfdOrder41(&Xt2[index]);
        PDupwindNthSymm2Xt2 = PDupwindNthSymmfdOrder42(&Xt2[index]);
        PDupwindNthSymm3Xt2 = PDupwindNthSymmfdOrder43(&Xt2[index]);
        PDupwindNthAnti1Xt2 = PDupwindNthAntifdOrder41(&Xt2[index]);
        PDupwindNthAnti2Xt2 = PDupwindNthAntifdOrder42(&Xt2[index]);
        PDupwindNthAnti3Xt2 = PDupwindNthAntifdOrder43(&Xt2[index]);
        PDupwindNthSymm1Xt3 = PDupwindNthSymmfdOrder41(&Xt3[index]);
        PDupwindNthSymm2Xt3 = PDupwindNthSymmfdOrder42(&Xt3[index]);
        PDupwindNthSymm3Xt3 = PDupwindNthSymmfdOrder43(&Xt3[index]);
        PDupwindNthAnti1Xt3 = PDupwindNthAntifdOrder41(&Xt3[index]);
        PDupwindNthAnti2Xt3 = PDupwindNthAntifdOrder42(&Xt3[index]);
        PDupwindNthAnti3Xt3 = PDupwindNthAntifdOrder43(&Xt3[index]);
        break;
      }
      
      case 6:
      {
        PDupwindNthSymm1B1 = PDupwindNthSymmfdOrder61(&B1[index]);
        PDupwindNthSymm2B1 = PDupwindNthSymmfdOrder62(&B1[index]);
        PDupwindNthSymm3B1 = PDupwindNthSymmfdOrder63(&B1[index]);
        PDupwindNthAnti1B1 = PDupwindNthAntifdOrder61(&B1[index]);
        PDupwindNthAnti2B1 = PDupwindNthAntifdOrder62(&B1[index]);
        PDupwindNthAnti3B1 = PDupwindNthAntifdOrder63(&B1[index]);
        PDupwindNthSymm1B2 = PDupwindNthSymmfdOrder61(&B2[index]);
        PDupwindNthSymm2B2 = PDupwindNthSymmfdOrder62(&B2[index]);
        PDupwindNthSymm3B2 = PDupwindNthSymmfdOrder63(&B2[index]);
        PDupwindNthAnti1B2 = PDupwindNthAntifdOrder61(&B2[index]);
        PDupwindNthAnti2B2 = PDupwindNthAntifdOrder62(&B2[index]);
        PDupwindNthAnti3B2 = PDupwindNthAntifdOrder63(&B2[index]);
        PDupwindNthSymm1B3 = PDupwindNthSymmfdOrder61(&B3[index]);
        PDupwindNthSymm2B3 = PDupwindNthSymmfdOrder62(&B3[index]);
        PDupwindNthSymm3B3 = PDupwindNthSymmfdOrder63(&B3[index]);
        PDupwindNthAnti1B3 = PDupwindNthAntifdOrder61(&B3[index]);
        PDupwindNthAnti2B3 = PDupwindNthAntifdOrder62(&B3[index]);
        PDupwindNthAnti3B3 = PDupwindNthAntifdOrder63(&B3[index]);
        PDupwindNthSymm1betax = PDupwindNthSymmfdOrder61(&betax[index]);
        PDupwindNthSymm2betax = PDupwindNthSymmfdOrder62(&betax[index]);
        PDupwindNthSymm3betax = PDupwindNthSymmfdOrder63(&betax[index]);
        PDupwindNthAnti1betax = PDupwindNthAntifdOrder61(&betax[index]);
        PDupwindNthAnti2betax = PDupwindNthAntifdOrder62(&betax[index]);
        PDupwindNthAnti3betax = PDupwindNthAntifdOrder63(&betax[index]);
        PDupwindNthSymm1betay = PDupwindNthSymmfdOrder61(&betay[index]);
        PDupwindNthSymm2betay = PDupwindNthSymmfdOrder62(&betay[index]);
        PDupwindNthSymm3betay = PDupwindNthSymmfdOrder63(&betay[index]);
        PDupwindNthAnti1betay = PDupwindNthAntifdOrder61(&betay[index]);
        PDupwindNthAnti2betay = PDupwindNthAntifdOrder62(&betay[index]);
        PDupwindNthAnti3betay = PDupwindNthAntifdOrder63(&betay[index]);
        PDupwindNthSymm1betaz = PDupwindNthSymmfdOrder61(&betaz[index]);
        PDupwindNthSymm2betaz = PDupwindNthSymmfdOrder62(&betaz[index]);
        PDupwindNthSymm3betaz = PDupwindNthSymmfdOrder63(&betaz[index]);
        PDupwindNthAnti1betaz = PDupwindNthAntifdOrder61(&betaz[index]);
        PDupwindNthAnti2betaz = PDupwindNthAntifdOrder62(&betaz[index]);
        PDupwindNthAnti3betaz = PDupwindNthAntifdOrder63(&betaz[index]);
        PDupwindNthSymm1Xt1 = PDupwindNthSymmfdOrder61(&Xt1[index]);
        PDupwindNthSymm2Xt1 = PDupwindNthSymmfdOrder62(&Xt1[index]);
        PDupwindNthSymm3Xt1 = PDupwindNthSymmfdOrder63(&Xt1[index]);
        PDupwindNthAnti1Xt1 = PDupwindNthAntifdOrder61(&Xt1[index]);
        PDupwindNthAnti2Xt1 = PDupwindNthAntifdOrder62(&Xt1[index]);
        PDupwindNthAnti3Xt1 = PDupwindNthAntifdOrder63(&Xt1[index]);
        PDupwindNthSymm1Xt2 = PDupwindNthSymmfdOrder61(&Xt2[index]);
        PDupwindNthSymm2Xt2 = PDupwindNthSymmfdOrder62(&Xt2[index]);
        PDupwindNthSymm3Xt2 = PDupwindNthSymmfdOrder63(&Xt2[index]);
        PDupwindNthAnti1Xt2 = PDupwindNthAntifdOrder61(&Xt2[index]);
        PDupwindNthAnti2Xt2 = PDupwindNthAntifdOrder62(&Xt2[index]);
        PDupwindNthAnti3Xt2 = PDupwindNthAntifdOrder63(&Xt2[index]);
        PDupwindNthSymm1Xt3 = PDupwindNthSymmfdOrder61(&Xt3[index]);
        PDupwindNthSymm2Xt3 = PDupwindNthSymmfdOrder62(&Xt3[index]);
        PDupwindNthSymm3Xt3 = PDupwindNthSymmfdOrder63(&Xt3[index]);
        PDupwindNthAnti1Xt3 = PDupwindNthAntifdOrder61(&Xt3[index]);
        PDupwindNthAnti2Xt3 = PDupwindNthAntifdOrder62(&Xt3[index]);
        PDupwindNthAnti3Xt3 = PDupwindNthAntifdOrder63(&Xt3[index]);
        break;
      }
      
      case 8:
      {
        PDupwindNthSymm1B1 = PDupwindNthSymmfdOrder81(&B1[index]);
        PDupwindNthSymm2B1 = PDupwindNthSymmfdOrder82(&B1[index]);
        PDupwindNthSymm3B1 = PDupwindNthSymmfdOrder83(&B1[index]);
        PDupwindNthAnti1B1 = PDupwindNthAntifdOrder81(&B1[index]);
        PDupwindNthAnti2B1 = PDupwindNthAntifdOrder82(&B1[index]);
        PDupwindNthAnti3B1 = PDupwindNthAntifdOrder83(&B1[index]);
        PDupwindNthSymm1B2 = PDupwindNthSymmfdOrder81(&B2[index]);
        PDupwindNthSymm2B2 = PDupwindNthSymmfdOrder82(&B2[index]);
        PDupwindNthSymm3B2 = PDupwindNthSymmfdOrder83(&B2[index]);
        PDupwindNthAnti1B2 = PDupwindNthAntifdOrder81(&B2[index]);
        PDupwindNthAnti2B2 = PDupwindNthAntifdOrder82(&B2[index]);
        PDupwindNthAnti3B2 = PDupwindNthAntifdOrder83(&B2[index]);
        PDupwindNthSymm1B3 = PDupwindNthSymmfdOrder81(&B3[index]);
        PDupwindNthSymm2B3 = PDupwindNthSymmfdOrder82(&B3[index]);
        PDupwindNthSymm3B3 = PDupwindNthSymmfdOrder83(&B3[index]);
        PDupwindNthAnti1B3 = PDupwindNthAntifdOrder81(&B3[index]);
        PDupwindNthAnti2B3 = PDupwindNthAntifdOrder82(&B3[index]);
        PDupwindNthAnti3B3 = PDupwindNthAntifdOrder83(&B3[index]);
        PDupwindNthSymm1betax = PDupwindNthSymmfdOrder81(&betax[index]);
        PDupwindNthSymm2betax = PDupwindNthSymmfdOrder82(&betax[index]);
        PDupwindNthSymm3betax = PDupwindNthSymmfdOrder83(&betax[index]);
        PDupwindNthAnti1betax = PDupwindNthAntifdOrder81(&betax[index]);
        PDupwindNthAnti2betax = PDupwindNthAntifdOrder82(&betax[index]);
        PDupwindNthAnti3betax = PDupwindNthAntifdOrder83(&betax[index]);
        PDupwindNthSymm1betay = PDupwindNthSymmfdOrder81(&betay[index]);
        PDupwindNthSymm2betay = PDupwindNthSymmfdOrder82(&betay[index]);
        PDupwindNthSymm3betay = PDupwindNthSymmfdOrder83(&betay[index]);
        PDupwindNthAnti1betay = PDupwindNthAntifdOrder81(&betay[index]);
        PDupwindNthAnti2betay = PDupwindNthAntifdOrder82(&betay[index]);
        PDupwindNthAnti3betay = PDupwindNthAntifdOrder83(&betay[index]);
        PDupwindNthSymm1betaz = PDupwindNthSymmfdOrder81(&betaz[index]);
        PDupwindNthSymm2betaz = PDupwindNthSymmfdOrder82(&betaz[index]);
        PDupwindNthSymm3betaz = PDupwindNthSymmfdOrder83(&betaz[index]);
        PDupwindNthAnti1betaz = PDupwindNthAntifdOrder81(&betaz[index]);
        PDupwindNthAnti2betaz = PDupwindNthAntifdOrder82(&betaz[index]);
        PDupwindNthAnti3betaz = PDupwindNthAntifdOrder83(&betaz[index]);
        PDupwindNthSymm1Xt1 = PDupwindNthSymmfdOrder81(&Xt1[index]);
        PDupwindNthSymm2Xt1 = PDupwindNthSymmfdOrder82(&Xt1[index]);
        PDupwindNthSymm3Xt1 = PDupwindNthSymmfdOrder83(&Xt1[index]);
        PDupwindNthAnti1Xt1 = PDupwindNthAntifdOrder81(&Xt1[index]);
        PDupwindNthAnti2Xt1 = PDupwindNthAntifdOrder82(&Xt1[index]);
        PDupwindNthAnti3Xt1 = PDupwindNthAntifdOrder83(&Xt1[index]);
        PDupwindNthSymm1Xt2 = PDupwindNthSymmfdOrder81(&Xt2[index]);
        PDupwindNthSymm2Xt2 = PDupwindNthSymmfdOrder82(&Xt2[index]);
        PDupwindNthSymm3Xt2 = PDupwindNthSymmfdOrder83(&Xt2[index]);
        PDupwindNthAnti1Xt2 = PDupwindNthAntifdOrder81(&Xt2[index]);
        PDupwindNthAnti2Xt2 = PDupwindNthAntifdOrder82(&Xt2[index]);
        PDupwindNthAnti3Xt2 = PDupwindNthAntifdOrder83(&Xt2[index]);
        PDupwindNthSymm1Xt3 = PDupwindNthSymmfdOrder81(&Xt3[index]);
        PDupwindNthSymm2Xt3 = PDupwindNthSymmfdOrder82(&Xt3[index]);
        PDupwindNthSymm3Xt3 = PDupwindNthSymmfdOrder83(&Xt3[index]);
        PDupwindNthAnti1Xt3 = PDupwindNthAntifdOrder81(&Xt3[index]);
        PDupwindNthAnti2Xt3 = PDupwindNthAntifdOrder82(&Xt3[index]);
        PDupwindNthAnti3Xt3 = PDupwindNthAntifdOrder83(&Xt3[index]);
        break;
      }
      default:
        CCTK_BUILTIN_UNREACHABLE();
    }
    /* Calculate temporaries and grid functions */
    CCTK_REAL JacPDupwindNthAnti1B1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDupwindNthAnti1B2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDupwindNthAnti1B3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDupwindNthAnti1betax CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDupwindNthAnti1betay CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDupwindNthAnti1betaz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDupwindNthAnti1Xt1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDupwindNthAnti1Xt2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDupwindNthAnti1Xt3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDupwindNthAnti2B1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDupwindNthAnti2B2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDupwindNthAnti2B3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDupwindNthAnti2betax CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDupwindNthAnti2betay CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDupwindNthAnti2betaz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDupwindNthAnti2Xt1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDupwindNthAnti2Xt2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDupwindNthAnti2Xt3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDupwindNthAnti3B1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDupwindNthAnti3B2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDupwindNthAnti3B3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDupwindNthAnti3betax CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDupwindNthAnti3betay CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDupwindNthAnti3betaz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDupwindNthAnti3Xt1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDupwindNthAnti3Xt2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDupwindNthAnti3Xt3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDupwindNthSymm1B1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDupwindNthSymm1B2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDupwindNthSymm1B3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDupwindNthSymm1betax CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDupwindNthSymm1betay CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDupwindNthSymm1betaz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDupwindNthSymm1Xt1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDupwindNthSymm1Xt2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDupwindNthSymm1Xt3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDupwindNthSymm2B1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDupwindNthSymm2B2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDupwindNthSymm2B3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDupwindNthSymm2betax CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDupwindNthSymm2betay CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDupwindNthSymm2betaz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDupwindNthSymm2Xt1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDupwindNthSymm2Xt2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDupwindNthSymm2Xt3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDupwindNthSymm3B1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDupwindNthSymm3B2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDupwindNthSymm3B3 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDupwindNthSymm3betax CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDupwindNthSymm3betay CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDupwindNthSymm3betaz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDupwindNthSymm3Xt1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDupwindNthSymm3Xt2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL JacPDupwindNthSymm3Xt3 CCTK_ATTRIBUTE_UNUSED;
    
    if (usejacobian)
    {
      JacPDupwindNthSymm1B1 = J11L*PDupwindNthSymm1B1 + 
        J21L*PDupwindNthSymm2B1 + J31L*PDupwindNthSymm3B1;
      
      JacPDupwindNthSymm1B2 = J11L*PDupwindNthSymm1B2 + 
        J21L*PDupwindNthSymm2B2 + J31L*PDupwindNthSymm3B2;
      
      JacPDupwindNthSymm1B3 = J11L*PDupwindNthSymm1B3 + 
        J21L*PDupwindNthSymm2B3 + J31L*PDupwindNthSymm3B3;
      
      JacPDupwindNthSymm1betax = J11L*PDupwindNthSymm1betax + 
        J21L*PDupwindNthSymm2betax + J31L*PDupwindNthSymm3betax;
      
      JacPDupwindNthSymm1betay = J11L*PDupwindNthSymm1betay + 
        J21L*PDupwindNthSymm2betay + J31L*PDupwindNthSymm3betay;
      
      JacPDupwindNthSymm1betaz = J11L*PDupwindNthSymm1betaz + 
        J21L*PDupwindNthSymm2betaz + J31L*PDupwindNthSymm3betaz;
      
      JacPDupwindNthSymm1Xt1 = J11L*PDupwindNthSymm1Xt1 + 
        J21L*PDupwindNthSymm2Xt1 + J31L*PDupwindNthSymm3Xt1;
      
      JacPDupwindNthSymm1Xt2 = J11L*PDupwindNthSymm1Xt2 + 
        J21L*PDupwindNthSymm2Xt2 + J31L*PDupwindNthSymm3Xt2;
      
      JacPDupwindNthSymm1Xt3 = J11L*PDupwindNthSymm1Xt3 + 
        J21L*PDupwindNthSymm2Xt3 + J31L*PDupwindNthSymm3Xt3;
      
      JacPDupwindNthSymm2B1 = J12L*PDupwindNthSymm1B1 + 
        J22L*PDupwindNthSymm2B1 + J32L*PDupwindNthSymm3B1;
      
      JacPDupwindNthSymm2B2 = J12L*PDupwindNthSymm1B2 + 
        J22L*PDupwindNthSymm2B2 + J32L*PDupwindNthSymm3B2;
      
      JacPDupwindNthSymm2B3 = J12L*PDupwindNthSymm1B3 + 
        J22L*PDupwindNthSymm2B3 + J32L*PDupwindNthSymm3B3;
      
      JacPDupwindNthSymm2betax = J12L*PDupwindNthSymm1betax + 
        J22L*PDupwindNthSymm2betax + J32L*PDupwindNthSymm3betax;
      
      JacPDupwindNthSymm2betay = J12L*PDupwindNthSymm1betay + 
        J22L*PDupwindNthSymm2betay + J32L*PDupwindNthSymm3betay;
      
      JacPDupwindNthSymm2betaz = J12L*PDupwindNthSymm1betaz + 
        J22L*PDupwindNthSymm2betaz + J32L*PDupwindNthSymm3betaz;
      
      JacPDupwindNthSymm2Xt1 = J12L*PDupwindNthSymm1Xt1 + 
        J22L*PDupwindNthSymm2Xt1 + J32L*PDupwindNthSymm3Xt1;
      
      JacPDupwindNthSymm2Xt2 = J12L*PDupwindNthSymm1Xt2 + 
        J22L*PDupwindNthSymm2Xt2 + J32L*PDupwindNthSymm3Xt2;
      
      JacPDupwindNthSymm2Xt3 = J12L*PDupwindNthSymm1Xt3 + 
        J22L*PDupwindNthSymm2Xt3 + J32L*PDupwindNthSymm3Xt3;
      
      JacPDupwindNthSymm3B1 = J13L*PDupwindNthSymm1B1 + 
        J23L*PDupwindNthSymm2B1 + J33L*PDupwindNthSymm3B1;
      
      JacPDupwindNthSymm3B2 = J13L*PDupwindNthSymm1B2 + 
        J23L*PDupwindNthSymm2B2 + J33L*PDupwindNthSymm3B2;
      
      JacPDupwindNthSymm3B3 = J13L*PDupwindNthSymm1B3 + 
        J23L*PDupwindNthSymm2B3 + J33L*PDupwindNthSymm3B3;
      
      JacPDupwindNthSymm3betax = J13L*PDupwindNthSymm1betax + 
        J23L*PDupwindNthSymm2betax + J33L*PDupwindNthSymm3betax;
      
      JacPDupwindNthSymm3betay = J13L*PDupwindNthSymm1betay + 
        J23L*PDupwindNthSymm2betay + J33L*PDupwindNthSymm3betay;
      
      JacPDupwindNthSymm3betaz = J13L*PDupwindNthSymm1betaz + 
        J23L*PDupwindNthSymm2betaz + J33L*PDupwindNthSymm3betaz;
      
      JacPDupwindNthSymm3Xt1 = J13L*PDupwindNthSymm1Xt1 + 
        J23L*PDupwindNthSymm2Xt1 + J33L*PDupwindNthSymm3Xt1;
      
      JacPDupwindNthSymm3Xt2 = J13L*PDupwindNthSymm1Xt2 + 
        J23L*PDupwindNthSymm2Xt2 + J33L*PDupwindNthSymm3Xt2;
      
      JacPDupwindNthSymm3Xt3 = J13L*PDupwindNthSymm1Xt3 + 
        J23L*PDupwindNthSymm2Xt3 + J33L*PDupwindNthSymm3Xt3;
      
      JacPDupwindNthAnti1B1 = J11L*PDupwindNthAnti1B1 + 
        J21L*PDupwindNthAnti2B1 + J31L*PDupwindNthAnti3B1;
      
      JacPDupwindNthAnti1B2 = J11L*PDupwindNthAnti1B2 + 
        J21L*PDupwindNthAnti2B2 + J31L*PDupwindNthAnti3B2;
      
      JacPDupwindNthAnti1B3 = J11L*PDupwindNthAnti1B3 + 
        J21L*PDupwindNthAnti2B3 + J31L*PDupwindNthAnti3B3;
      
      JacPDupwindNthAnti1betax = J11L*PDupwindNthAnti1betax + 
        J21L*PDupwindNthAnti2betax + J31L*PDupwindNthAnti3betax;
      
      JacPDupwindNthAnti1betay = J11L*PDupwindNthAnti1betay + 
        J21L*PDupwindNthAnti2betay + J31L*PDupwindNthAnti3betay;
      
      JacPDupwindNthAnti1betaz = J11L*PDupwindNthAnti1betaz + 
        J21L*PDupwindNthAnti2betaz + J31L*PDupwindNthAnti3betaz;
      
      JacPDupwindNthAnti1Xt1 = J11L*PDupwindNthAnti1Xt1 + 
        J21L*PDupwindNthAnti2Xt1 + J31L*PDupwindNthAnti3Xt1;
      
      JacPDupwindNthAnti1Xt2 = J11L*PDupwindNthAnti1Xt2 + 
        J21L*PDupwindNthAnti2Xt2 + J31L*PDupwindNthAnti3Xt2;
      
      JacPDupwindNthAnti1Xt3 = J11L*PDupwindNthAnti1Xt3 + 
        J21L*PDupwindNthAnti2Xt3 + J31L*PDupwindNthAnti3Xt3;
      
      JacPDupwindNthAnti2B1 = J12L*PDupwindNthAnti1B1 + 
        J22L*PDupwindNthAnti2B1 + J32L*PDupwindNthAnti3B1;
      
      JacPDupwindNthAnti2B2 = J12L*PDupwindNthAnti1B2 + 
        J22L*PDupwindNthAnti2B2 + J32L*PDupwindNthAnti3B2;
      
      JacPDupwindNthAnti2B3 = J12L*PDupwindNthAnti1B3 + 
        J22L*PDupwindNthAnti2B3 + J32L*PDupwindNthAnti3B3;
      
      JacPDupwindNthAnti2betax = J12L*PDupwindNthAnti1betax + 
        J22L*PDupwindNthAnti2betax + J32L*PDupwindNthAnti3betax;
      
      JacPDupwindNthAnti2betay = J12L*PDupwindNthAnti1betay + 
        J22L*PDupwindNthAnti2betay + J32L*PDupwindNthAnti3betay;
      
      JacPDupwindNthAnti2betaz = J12L*PDupwindNthAnti1betaz + 
        J22L*PDupwindNthAnti2betaz + J32L*PDupwindNthAnti3betaz;
      
      JacPDupwindNthAnti2Xt1 = J12L*PDupwindNthAnti1Xt1 + 
        J22L*PDupwindNthAnti2Xt1 + J32L*PDupwindNthAnti3Xt1;
      
      JacPDupwindNthAnti2Xt2 = J12L*PDupwindNthAnti1Xt2 + 
        J22L*PDupwindNthAnti2Xt2 + J32L*PDupwindNthAnti3Xt2;
      
      JacPDupwindNthAnti2Xt3 = J12L*PDupwindNthAnti1Xt3 + 
        J22L*PDupwindNthAnti2Xt3 + J32L*PDupwindNthAnti3Xt3;
      
      JacPDupwindNthAnti3B1 = J13L*PDupwindNthAnti1B1 + 
        J23L*PDupwindNthAnti2B1 + J33L*PDupwindNthAnti3B1;
      
      JacPDupwindNthAnti3B2 = J13L*PDupwindNthAnti1B2 + 
        J23L*PDupwindNthAnti2B2 + J33L*PDupwindNthAnti3B2;
      
      JacPDupwindNthAnti3B3 = J13L*PDupwindNthAnti1B3 + 
        J23L*PDupwindNthAnti2B3 + J33L*PDupwindNthAnti3B3;
      
      JacPDupwindNthAnti3betax = J13L*PDupwindNthAnti1betax + 
        J23L*PDupwindNthAnti2betax + J33L*PDupwindNthAnti3betax;
      
      JacPDupwindNthAnti3betay = J13L*PDupwindNthAnti1betay + 
        J23L*PDupwindNthAnti2betay + J33L*PDupwindNthAnti3betay;
      
      JacPDupwindNthAnti3betaz = J13L*PDupwindNthAnti1betaz + 
        J23L*PDupwindNthAnti2betaz + J33L*PDupwindNthAnti3betaz;
      
      JacPDupwindNthAnti3Xt1 = J13L*PDupwindNthAnti1Xt1 + 
        J23L*PDupwindNthAnti2Xt1 + J33L*PDupwindNthAnti3Xt1;
      
      JacPDupwindNthAnti3Xt2 = J13L*PDupwindNthAnti1Xt2 + 
        J23L*PDupwindNthAnti2Xt2 + J33L*PDupwindNthAnti3Xt2;
      
      JacPDupwindNthAnti3Xt3 = J13L*PDupwindNthAnti1Xt3 + 
        J23L*PDupwindNthAnti2Xt3 + J33L*PDupwindNthAnti3Xt3;
    }
    else
    {
      JacPDupwindNthSymm1B1 = PDupwindNthSymm1B1;
      
      JacPDupwindNthSymm1B2 = PDupwindNthSymm1B2;
      
      JacPDupwindNthSymm1B3 = PDupwindNthSymm1B3;
      
      JacPDupwindNthSymm1betax = PDupwindNthSymm1betax;
      
      JacPDupwindNthSymm1betay = PDupwindNthSymm1betay;
      
      JacPDupwindNthSymm1betaz = PDupwindNthSymm1betaz;
      
      JacPDupwindNthSymm1Xt1 = PDupwindNthSymm1Xt1;
      
      JacPDupwindNthSymm1Xt2 = PDupwindNthSymm1Xt2;
      
      JacPDupwindNthSymm1Xt3 = PDupwindNthSymm1Xt3;
      
      JacPDupwindNthSymm2B1 = PDupwindNthSymm2B1;
      
      JacPDupwindNthSymm2B2 = PDupwindNthSymm2B2;
      
      JacPDupwindNthSymm2B3 = PDupwindNthSymm2B3;
      
      JacPDupwindNthSymm2betax = PDupwindNthSymm2betax;
      
      JacPDupwindNthSymm2betay = PDupwindNthSymm2betay;
      
      JacPDupwindNthSymm2betaz = PDupwindNthSymm2betaz;
      
      JacPDupwindNthSymm2Xt1 = PDupwindNthSymm2Xt1;
      
      JacPDupwindNthSymm2Xt2 = PDupwindNthSymm2Xt2;
      
      JacPDupwindNthSymm2Xt3 = PDupwindNthSymm2Xt3;
      
      JacPDupwindNthSymm3B1 = PDupwindNthSymm3B1;
      
      JacPDupwindNthSymm3B2 = PDupwindNthSymm3B2;
      
      JacPDupwindNthSymm3B3 = PDupwindNthSymm3B3;
      
      JacPDupwindNthSymm3betax = PDupwindNthSymm3betax;
      
      JacPDupwindNthSymm3betay = PDupwindNthSymm3betay;
      
      JacPDupwindNthSymm3betaz = PDupwindNthSymm3betaz;
      
      JacPDupwindNthSymm3Xt1 = PDupwindNthSymm3Xt1;
      
      JacPDupwindNthSymm3Xt2 = PDupwindNthSymm3Xt2;
      
      JacPDupwindNthSymm3Xt3 = PDupwindNthSymm3Xt3;
      
      JacPDupwindNthAnti1B1 = PDupwindNthAnti1B1;
      
      JacPDupwindNthAnti1B2 = PDupwindNthAnti1B2;
      
      JacPDupwindNthAnti1B3 = PDupwindNthAnti1B3;
      
      JacPDupwindNthAnti1betax = PDupwindNthAnti1betax;
      
      JacPDupwindNthAnti1betay = PDupwindNthAnti1betay;
      
      JacPDupwindNthAnti1betaz = PDupwindNthAnti1betaz;
      
      JacPDupwindNthAnti1Xt1 = PDupwindNthAnti1Xt1;
      
      JacPDupwindNthAnti1Xt2 = PDupwindNthAnti1Xt2;
      
      JacPDupwindNthAnti1Xt3 = PDupwindNthAnti1Xt3;
      
      JacPDupwindNthAnti2B1 = PDupwindNthAnti2B1;
      
      JacPDupwindNthAnti2B2 = PDupwindNthAnti2B2;
      
      JacPDupwindNthAnti2B3 = PDupwindNthAnti2B3;
      
      JacPDupwindNthAnti2betax = PDupwindNthAnti2betax;
      
      JacPDupwindNthAnti2betay = PDupwindNthAnti2betay;
      
      JacPDupwindNthAnti2betaz = PDupwindNthAnti2betaz;
      
      JacPDupwindNthAnti2Xt1 = PDupwindNthAnti2Xt1;
      
      JacPDupwindNthAnti2Xt2 = PDupwindNthAnti2Xt2;
      
      JacPDupwindNthAnti2Xt3 = PDupwindNthAnti2Xt3;
      
      JacPDupwindNthAnti3B1 = PDupwindNthAnti3B1;
      
      JacPDupwindNthAnti3B2 = PDupwindNthAnti3B2;
      
      JacPDupwindNthAnti3B3 = PDupwindNthAnti3B3;
      
      JacPDupwindNthAnti3betax = PDupwindNthAnti3betax;
      
      JacPDupwindNthAnti3betay = PDupwindNthAnti3betay;
      
      JacPDupwindNthAnti3betaz = PDupwindNthAnti3betaz;
      
      JacPDupwindNthAnti3Xt1 = PDupwindNthAnti3Xt1;
      
      JacPDupwindNthAnti3Xt2 = PDupwindNthAnti3Xt2;
      
      JacPDupwindNthAnti3Xt3 = PDupwindNthAnti3Xt3;
    }
    
    CCTK_REAL partialtbeta1 CCTK_ATTRIBUTE_UNUSED = B1L + 
      IfThen(betaFullLieDeriv != 0,betaxL*JacPDupwindNthAnti1betax + 
      betayL*JacPDupwindNthAnti2betax + betazL*JacPDupwindNthAnti3betax + 
      JacPDupwindNthSymm1betax*fabs(betaxL) + 
      JacPDupwindNthSymm2betax*fabs(betayL) + 
      JacPDupwindNthSymm3betax*fabs(betazL),0);
    
    CCTK_REAL partialtbeta2 CCTK_ATTRIBUTE_UNUSED = B2L + 
      IfThen(betaFullLieDeriv != 0,betaxL*JacPDupwindNthAnti1betay + 
      betayL*JacPDupwindNthAnti2betay + betazL*JacPDupwindNthAnti3betay + 
      JacPDupwindNthSymm1betay*fabs(betaxL) + 
      JacPDupwindNthSymm2betay*fabs(betayL) + 
      JacPDupwindNthSymm3betay*fabs(betazL),0);
    
    CCTK_REAL partialtbeta3 CCTK_ATTRIBUTE_UNUSED = B3L + 
      IfThen(betaFullLieDeriv != 0,betaxL*JacPDupwindNthAnti1betaz + 
      betayL*JacPDupwindNthAnti2betaz + betazL*JacPDupwindNthAnti3betaz + 
      JacPDupwindNthSymm1betaz*fabs(betaxL) + 
      JacPDupwindNthSymm2betaz*fabs(betayL) + 
      JacPDupwindNthSymm3betaz*fabs(betazL),0);
    
    CCTK_REAL partialtB1 CCTK_ATTRIBUTE_UNUSED = -(B1L*betaEta) + 
      IfThen(betaFullLieDeriv != 0,betaxL*JacPDupwindNthAnti1B1 + 
      betayL*JacPDupwindNthAnti2B1 + betazL*JacPDupwindNthAnti3B1 + 
      JacPDupwindNthSymm1B1*fabs(betaxL) + JacPDupwindNthSymm2B1*fabs(betayL) 
      + JacPDupwindNthSymm3B1*fabs(betazL),0) + betaXi*(Xt1rhsL - 
      IfThen(betaFullLieDeriv != 0,betaxL*JacPDupwindNthAnti1Xt1 + 
      betayL*JacPDupwindNthAnti2Xt1 + betazL*JacPDupwindNthAnti3Xt1 + 
      JacPDupwindNthSymm1Xt1*fabs(betaxL) + 
      JacPDupwindNthSymm2Xt1*fabs(betayL) + 
      JacPDupwindNthSymm3Xt1*fabs(betazL),0))*pow(alpL,betaP);
    
    CCTK_REAL partialtB2 CCTK_ATTRIBUTE_UNUSED = -(B2L*betaEta) + 
      IfThen(betaFullLieDeriv != 0,betaxL*JacPDupwindNthAnti1B2 + 
      betayL*JacPDupwindNthAnti2B2 + betazL*JacPDupwindNthAnti3B2 + 
      JacPDupwindNthSymm1B2*fabs(betaxL) + JacPDupwindNthSymm2B2*fabs(betayL) 
      + JacPDupwindNthSymm3B2*fabs(betazL),0) + betaXi*(Xt2rhsL - 
      IfThen(betaFullLieDeriv != 0,betaxL*JacPDupwindNthAnti1Xt2 + 
      betayL*JacPDupwindNthAnti2Xt2 + betazL*JacPDupwindNthAnti3Xt2 + 
      JacPDupwindNthSymm1Xt2*fabs(betaxL) + 
      JacPDupwindNthSymm2Xt2*fabs(betayL) + 
      JacPDupwindNthSymm3Xt2*fabs(betazL),0))*pow(alpL,betaP);
    
    CCTK_REAL partialtB3 CCTK_ATTRIBUTE_UNUSED = -(B3L*betaEta) + 
      IfThen(betaFullLieDeriv != 0,betaxL*JacPDupwindNthAnti1B3 + 
      betayL*JacPDupwindNthAnti2B3 + betazL*JacPDupwindNthAnti3B3 + 
      JacPDupwindNthSymm1B3*fabs(betaxL) + JacPDupwindNthSymm2B3*fabs(betayL) 
      + JacPDupwindNthSymm3B3*fabs(betazL),0) + betaXi*(Xt3rhsL - 
      IfThen(betaFullLieDeriv != 0,betaxL*JacPDupwindNthAnti1Xt3 + 
      betayL*JacPDupwindNthAnti2Xt3 + betazL*JacPDupwindNthAnti3Xt3 + 
      JacPDupwindNthSymm1Xt3*fabs(betaxL) + 
      JacPDupwindNthSymm2Xt3*fabs(betayL) + 
      JacPDupwindNthSymm3Xt3*fabs(betazL),0))*pow(alpL,betaP);
    
    CCTK_REAL dtbetaxL CCTK_ATTRIBUTE_UNUSED = partialtbeta1;
    
    CCTK_REAL dtbetayL CCTK_ATTRIBUTE_UNUSED = partialtbeta2;
    
    CCTK_REAL dtbetazL CCTK_ATTRIBUTE_UNUSED = partialtbeta3;
    
    CCTK_REAL beta1rhsL CCTK_ATTRIBUTE_UNUSED = partialtbeta1;
    
    CCTK_REAL beta2rhsL CCTK_ATTRIBUTE_UNUSED = partialtbeta2;
    
    CCTK_REAL beta3rhsL CCTK_ATTRIBUTE_UNUSED = partialtbeta3;
    
    CCTK_REAL B1rhsL CCTK_ATTRIBUTE_UNUSED = partialtB1;
    
    CCTK_REAL B2rhsL CCTK_ATTRIBUTE_UNUSED = partialtB2;
    
    CCTK_REAL B3rhsL CCTK_ATTRIBUTE_UNUSED = partialtB3;
    /* Copy local copies back to grid functions */
    B1rhs[index] = B1rhsL;
    B2rhs[index] = B2rhsL;
    B3rhs[index] = B3rhsL;
    beta1rhs[index] = beta1rhsL;
    beta2rhs[index] = beta2rhsL;
    beta3rhs[index] = beta3rhsL;
    dtbetax[index] = dtbetaxL;
    dtbetay[index] = dtbetayL;
    dtbetaz[index] = dtbetazL;
  }
  CCTK_ENDLOOP3(CosmoLapse_dtShiftPostRHS);
}
extern "C" void CosmoLapse_dtShiftPostRHS(CCTK_ARGUMENTS)
{
  #ifdef DECLARE_CCTK_ARGUMENTS_CosmoLapse_dtShiftPostRHS
  DECLARE_CCTK_ARGUMENTS_CHECKED(CosmoLapse_dtShiftPostRHS);
  #else
  DECLARE_CCTK_ARGUMENTS;
  #endif
  DECLARE_CCTK_PARAMETERS;
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering CosmoLapse_dtShiftPostRHS_Body");
  }
  if (cctk_iteration % CosmoLapse_dtShiftPostRHS_calc_every != CosmoLapse_dtShiftPostRHS_calc_offset)
  {
    return;
  }
  
  const char* const groups[] = {
    "ADMBase::dtshift",
    "ADMBase::lapse",
    "ADMBase::shift",
    "ML_BSSN::ML_dtshift",
    "ML_BSSN::ML_dtshiftrhs",
    "ML_BSSN::ML_Gamma",
    "ML_BSSN::ML_Gammarhs",
    "ML_BSSN::ML_shiftrhs"};
  AssertGroupStorage(cctkGH, "CosmoLapse_dtShiftPostRHS", 8, groups);
  
  switch (fdOrder)
  {
    case 2:
    {
      EnsureStencilFits(cctkGH, "CosmoLapse_dtShiftPostRHS", 2, 2, 2);
      break;
    }
    
    case 4:
    {
      EnsureStencilFits(cctkGH, "CosmoLapse_dtShiftPostRHS", 3, 3, 3);
      break;
    }
    
    case 6:
    {
      EnsureStencilFits(cctkGH, "CosmoLapse_dtShiftPostRHS", 4, 4, 4);
      break;
    }
    
    case 8:
    {
      EnsureStencilFits(cctkGH, "CosmoLapse_dtShiftPostRHS", 5, 5, 5);
      break;
    }
    default:
      CCTK_BUILTIN_UNREACHABLE();
  }
  
  LoopOverInterior(cctkGH, CosmoLapse_dtShiftPostRHS_Body);
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Leaving CosmoLapse_dtShiftPostRHS_Body");
  }
}

} // namespace CosmoLapse
