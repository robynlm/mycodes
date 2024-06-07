SetEnhancedTimes[False];

(******************************************************************************)
(******************************************************************************)
(******************************************************************************)
(* Options *)
(******************************************************************************)
(******************************************************************************)
(******************************************************************************)

derivOrder = 4;
maxTimelevels = 4;

prefix = "CT_";
suffix = "" <> If [derivOrder!=4, "_O" <> ToString[derivOrder], ""];

CTThorn = prefix <> "Dust" <> suffix;

(******************************************************************************)
(******************************************************************************)
(******************************************************************************)
(* Derivatives *)
(******************************************************************************)
(******************************************************************************)
(******************************************************************************)

KD = KroneckerDelta;

derivatives =
{
  PDstandardNth[i_]    -> StandardCenteredDifferenceOperator[1,fdOrder/2,i],
  PDstandardNth[i_,i_] -> StandardCenteredDifferenceOperator[2,fdOrder/2,i],
  PDstandardNth[i_,j_] -> StandardCenteredDifferenceOperator[1,fdOrder/2,i] *
                          StandardCenteredDifferenceOperator[1,fdOrder/2,j],
  PDdissipationNth[i_] ->
    spacing[i]^(fdOrder+1) / 2^(fdOrder+2) *
    StandardCenteredDifferenceOperator[fdOrder+2,fdOrder/2+1,i],
  
(* PD: These come from my mathematica notebook
   "Upwind-Kranc-Convert.nb" that converts upwinding finite
   differencing operators generated by
   StandardUpwindDifferenceOperator into this form *)

  Sequence@@Flatten[Table[
   {PDupwindNth[i] -> Switch[fdOrder,
      2, (dir[i]*(-3 + 4*shift[i]^dir[i] - shift[i]^(2*dir[i])))/(2*spacing[i]),
      4, (dir[i]*(-10 - 3/shift[i]^dir[i] + 18*shift[i]^dir[i] -
          6*shift[i]^(2*dir[i]) + shift[i]^(3*dir[i])))/(12*spacing[i]),
      6, (dir[i]*(-35 + 2/shift[i]^(2*dir[i]) - 24/shift[i]^dir[i] + 80*shift[i]^dir[i] -
          30*shift[i]^(2*dir[i]) + 8*shift[i]^(3*dir[i]) - shift[i]^(4*dir[i])))/(60*spacing[i]),
      8, (dir[i]*(-378 - 5/shift[i]^(3*dir[i]) + 60/shift[i]^(2*dir[i]) - 420/shift[i]^dir[i] +
          1050*shift[i]^dir[i] - 420*shift[i]^(2*dir[i]) + 140*shift[i]^(3*dir[i]) - 30*shift[i]^(4*dir[i]) +
          3*shift[i]^(5*dir[i])))/(840*spacing[i])],

    PDupwindNthAnti[i] -> Switch[fdOrder,
      2, (+1 shift[i]^(-2) -4 shift[i]^(-1) +0 shift[i]^( 0) +4 shift[i]^(+1) -1 shift[i]^(+2)) / (4 spacing[i]),
      4, (-1 shift[i]^(-3) +6 shift[i]^(-2) -21 shift[i]^(-1 )+0 shift[i]^( 0) +21 shift[i]^(+1)
          -6 shift[i]^(+2) +1 shift[i]^(+3)) / (24 spacing[i]),
      6, (+1 shift[i]^(-4) -8 shift[i]^(-3) +32 shift[i]^(-2) -104 shift[i]^(-1) +0 shift[i]^( 0)
          +104 shift[i]^(+1) -32 shift[i]^(+2) +8 shift[i]^(+3) -1 shift[i]^(+4)) / (120 spacing[i]),
      8, (-3 shift[i]^(-5) +30 shift[i]^(-4) -145 shift[i]^(-3) +480 shift[i]^(-2) -1470 shift[i]^(-1)
          +0 shift[i]^( 0) +1470 shift[i]^(+1) -480 shift[i]^(+2) +145 shift[i]^(+3) -30 shift[i]^(+4)
          +3 shift[i]^(+5)) / (1680 spacing[i])],

    PDupwindNthSymm[i] -> Switch[fdOrder,
     2, (-1 shift[i]^(-2) +4 shift[i]^(-1) -6 shift[i]^( 0) +4 shift[i]^(+1) -1 shift[i]^(+2)) / (4 spacing[i]),
     4, (+1 shift[i]^(-3) -6 shift[i]^(-2) +15 shift[i]^(-1) -20 shift[i]^( 0) +15 shift[i]^(+1)
         -6 shift[i]^(+2) +1 shift[i]^(+3)) / (24 spacing[i]),
     6, (-1 shift[i]^(-4) +8 shift[i]^(-3) - 28 shift[i]^(-2)+56 shift[i]^(-1)-70 shift[i]^( 0)
         +56 shift[i]^(+1) -28 shift[i]^(+2) +8 shift[i]^(+3) -1 shift[i]^(+4)) / (120 spacing[i]),
     8, (+3 shift[i]^(-5) -30 shift[i]^(-4) +135 shift[i]^(-3) -360 shift[i]^(-2) +630 shift[i]^(-1)
         -756 shift[i]^( 0) +630 shift[i]^(+1) -360 shift[i]^(+2) +135 shift[i]^(+3) -30 shift[i]^(+4)
         +3 shift[i]^(+5)) / (1680 spacing[i])],

    (* TODO: make these higher order stencils *)
    PDonesided[i] -> dir[i] (-1 + shift[i]^dir[i]) / spacing[i]} /. i->j, {j,1,3}],1]
};

PD     = PDstandardNth;

If [splitUpwindDerivs,
    Upwind[dir_, var_, idx_] := dir PDua[var,idx] + Abs[dir] PDus[var,idx],
    Upwind[dir_, var_, idx_] := dir PDu[var,idx]];

(******************************************************************************)
(******************************************************************************)
(******************************************************************************)
(* Tensors *)
(******************************************************************************)
(******************************************************************************)
(******************************************************************************)

(* Register the tensor quantities with the TensorTools package *)
Map [DefineTensor,
     {DD, EE, SS, V,
      rho, eps, prs, W, u,
      T, T0, t, t0, rhodp,
      g, gu, detg, k, a, b, dir,
      dtS, SaSbdgab, SSloc, Sstress,
      bl, g4ss, ku, dsa, dsb, dsg, dsgu, dsbl, Liebg, Liebgu, dta, dtb, dtg, dtgu, dtbl,
      Gdtts, Gdtss, Gdstt, Gdsts, Gdsss, Gtts, Gtss, Gstt, Gsts, Gsss,
      uup, dtu, dsu, dsu0, hts, hss, hmtsl, hmtsu, hmss, huts, huss,
      Gtsu, Dsut, Dtus, Dsus, Dptus, Dpsut, Dpsus,
      Thetats, Thetass, Theta, sigmatt, sigmats, sigmass
}];

Map [AssertSymmetricIncreasing,
     {g[la,lb], T[la,lb], k[la,lb], sigmass[la,lb]}];
Map [AssertSymmetricIncreasing,
     {gu[ua,ub]}];

(* Use the CartGrid3D variable names *)
x1=x; x2=y; x3=z;

(* Use the ADMBase variable names *)
g11=gxx; g12=gxy; g22=gyy; g13=gxz; g23=gyz; g33=gzz;
a=alp;
b1=betax; b2=betay; b3=betaz;
dta=dtalp;
dtb1=dtbetax; dtb2=dtbetay; dtb3=dtbetaz;
k11=kxx; k12=kxy; k22=kyy; k13=kxz; k23=kyz; k33=kzz;
bssntrK=trK;

(* Use the TmunuBase variable names *)
T00=eTtt;
T01=eTtx; T02=eTty; T03=eTtz;
T11=eTxx; T12=eTxy; T22=eTyy; T13=eTxz; T23=eTyz; T33=eTzz;

(******************************************************************************)
(******************************************************************************)
(******************************************************************************)
(* Expressions *)
(******************************************************************************)
(******************************************************************************)
(******************************************************************************)

pi = N[Pi,40];

detgExpr  = Det [MatrixOfComponents [g [la,lb]]];

(******************************************************************************)
(******************************************************************************)
(******************************************************************************)
(* Groups *)
(******************************************************************************)
(******************************************************************************)
(******************************************************************************)

(* Define a group and declare its name *)
DefineGroup1[name_, tensor_, timelevels_:1] :=
  Module[{group},
         group = CreateGroupFromTensor[tensor];
         group = AddGroupExtra[group, Timelevels -> timelevels];
         group = SetGroupName[group, name];
         group];
       
evolvedGroups =
   {DefineGroup1[prefix <> "D", DD    , maxTimelevels],
    DefineGroup1[prefix <> "E", EE    , maxTimelevels],
    DefineGroup1[prefix <> "S", SS[la], maxTimelevels]};

evaluatedGroups =
   {SetGroupName [CreateGroupFromTensor [rho  ], prefix <> "rho"],
    SetGroupName [CreateGroupFromTensor [eps  ], prefix <> "eps"],
    SetGroupName [CreateGroupFromTensor [prs  ], prefix <> "prs"],
    SetGroupName [CreateGroupFromTensor [u[la]], prefix <> "u"],
    SetGroupName [CreateGroupFromTensor [V[ua]], prefix <> "V"],
    SetGroupName [CreateGroupFromTensor [W    ], prefix <> "W"],
    SetGroupName [CreateGroupFromTensor [rhodp], prefix <> "rhodp"],
    SetGroupName [CreateGroupFromTensor [Theta         ], prefix <> "Theta"],
    SetGroupName [CreateGroupFromTensor [sigmatt       ], prefix <> "sigmatt"],
    SetGroupName [CreateGroupFromTensor [sigmats[la]   ], prefix <> "sigmats"],
    SetGroupName [CreateGroupFromTensor [sigmass[la,lb]], prefix <> "sigmass"]};
(*
evaluatedGroups =
   {DefineGroup1[prefix <> "rho"    , rho   , maxTimelevels],
    DefineGroup1[prefix <> "eps"    , eps   , maxTimelevels],
    DefineGroup1[prefix <> "prs"    , prs   , maxTimelevels],
    DefineGroup1[prefix <> "u"      , u[la] , maxTimelevels],
    DefineGroup1[prefix <> "V"      , V[ua] , maxTimelevels],
    DefineGroup1[prefix <> "W"      , W     , maxTimelevels],
    DefineGroup1[prefix <> "rhodp"  , rhodp , maxTimelevels],
    DefineGroup1[prefix <> "Theta"  , Theta          , maxTimelevels],
    DefineGroup1[prefix <> "sigmatt", sigmatt        , maxTimelevels],
    DefineGroup1[prefix <> "sigmats", sigmats[la]    , maxTimelevels],
    DefineGroup1[prefix <> "sigmass", sigmass[la,lb] , maxTimelevels]};
*)

declaredGroups = Join [evolvedGroups, evaluatedGroups];
declaredGroupNames = Map [First, declaredGroups];

inheritedImplementations =
  {"ADMBase", "TmunuBase", "ML_BSSN"};

extraGroups =
  {{"Grid::coordinates", {x, y, z, r}},
   {"ADMBase::metric", {gxx, gxy, gxz, gyy, gyz, gzz}},
   {"ADMBase::lapse", {alp}},
   {"ADMBase::dtlapse", {dtalp}},
   {"ADMBase::shift", {betax, betay, betaz}},
   {"ADMBase::dtshift", {dtbetax, dtbetay, dtbetaz}},
   {"ADMBase::curv", {kxx, kxy, kxz, kyy, kyz, kzz}},
   {"ML_BSSN::ML_trace_curv", {trK}},
   {"TmunuBase::stress_energy_scalar", {eTtt}},
   {"TmunuBase::stress_energy_vector", {eTtx, eTty, eTtz}},
   {"TmunuBase::stress_energy_tensor", {eTxx, eTxy, eTxz, eTyy, eTyz, eTzz}}
};

groups = Join [declaredGroups, extraGroups];

(******************************************************************************)
(******************************************************************************)
(******************************************************************************)
(* Initial data *)
(******************************************************************************)
(******************************************************************************)
(******************************************************************************)

initialMinkowskiCalc =
{
  Name -> CTThorn <> "_Minkowski",
  Schedule -> {"IN CCTK_INITIAL after ADMBase_Initial before " <> CTThorn <> "_setCTTrhs"},
  ConditionalOnKeyword -> {"my_initial_data", "Minkowski"},
  Equations -> 
  {
    rho        -> 0,
    rhodp      -> 0,
    eps        -> 0,
    prs        -> eosw rho,
    u[la]      -> 0
  }
};

initialFLRWPertCalc =
{
  Name -> CTThorn <> "_FLRW_Pert",
  Schedule -> {"IN CCTK_INITIAL after ADMBase_Initial before " <> CTThorn <> "_setCTTrhs"},
  ConditionalOnKeyword -> {"my_initial_data", "pFLRW"},
  Shorthands -> {rhocrit, a2, aa, HH},
  Equations -> 
  {
    HH         -> "pflrw_H0" "pflrw_t0" / t,
    aa         -> "pflrw_a0" ( t / "pflrw_t0" )^(2/3),

    rhocrit    -> 3 HH^2 / (8 Pi),
 
    rho        -> "pflrw_omegaM" rhocrit ( 1 + Sum[ "(pflrw_ax["<>ToString[i]<>"])" Sin["pflrw_kx["<>ToString[i]<>"]" x], {i,0,19}]
			        	     + Sum[ "(pflrw_ay["<>ToString[i]<>"])" Sin["pflrw_ky["<>ToString[i]<>"]" y], {i,0,19}]
                                             + Sum[ "(pflrw_az["<>ToString[i]<>"])" Sin["pflrw_kz["<>ToString[i]<>"]" z], {i,0,19}] ),
    rhodp      -> rho aa^3,
    eps        -> 0,
    prs        -> eosw rho,
    u[la]      -> 0,

    a2         -> aa^2,

    g11        -> a2,
    g12        -> 0,
    g13        -> 0,
    g22        -> a2,
    g23        -> 0,
    g33        -> a2,

    k11        -> - a2 HH,
    k12        -> 0,
    k13        -> 0,
    k22        -> - a2 HH,
    k23        -> 0,
    k33        -> - a2 HH
  }
};

(******************************************************************************)
(******************************************************************************)
(******************************************************************************)
(* Fluid evolution *)
(******************************************************************************)
(******************************************************************************)
(******************************************************************************)

(******************************************************************************)
(* Convert from primitives *)
(******************************************************************************)
convertFromPrimitivesCalc =
{
  Name -> CTThorn <> "_convertFromPrimitives",
  Schedule -> {"IN CCTK_INITIAL after ADMBase_Initial after " <> CTThorn <> "_FLRW_Pert after " <> CTThorn <> "_Minkowski after CT_MultiLevel"},
  ConditionalOnKeyword -> {"formalism", "Wilson"},
  Shorthands -> {detg, gu[ua,ub], u0, h, rho0},
  Equations ->
  {
    detg      -> detgExpr,
    gu[ua,ub] -> 1/detg detgExpr MatrixInverse [g[ua,ub]],
    u0        -> Sqrt[gu[ua,ub] u[la] u[lb] + 1] / a,
    h         -> (1 + eps) (1 + eosw),
    rho0      -> rho / ( 1 + eps ),

    W  -> a Sqrt[detg] u0,
    DD -> W rho0,
    EE -> W rho0 eps,
    SS[la] -> W rho0 h u[la],
    V[ua]  -> gu[ua,ub] u[lb] / u0 - b[ua]
  }
};
       
convertFromPrimitivesCalcVal =
{
  Name -> CTThorn <> "_convertFromPrimitivesVal",
  Schedule -> {"IN CCTK_INITIAL after ADMBase_Initial after " <> CTThorn <> "_FLRW_Pert after " <> CTThorn <> "_Minkowski after CT_MultiLevel"},
  ConditionalOnKeyword -> {"formalism", "Valencia"},
  Shorthands -> {detg, gu[ua,ub], u0, W2},
  Equations ->
  {
    detg      -> detgExpr,
    gu[ua,ub] -> 1/detg detgExpr MatrixInverse [g[ua,ub]],
      
    (* fluid velocity *)
    u0 -> Sqrt[gu[ua,ub] u[la] u[lb] + 1] / a,
    W  -> a u0,
    W2 -> W^2,
    
    (* conserved quantities *)
    DD -> 0,  (* neglect *)
    EE -> Sqrt[detg] ((rho + prs) W2 - prs),
    SS[la] -> Sqrt[detg] ((rho + prs) W u[la]),
    V[ua] -> (gu[ua,ub] u[lb]) / W
  }
};

(******************************************************************************)
(* Convert to primitives *)
(******************************************************************************)

convertToPrimitivesCalc =
{
  Name -> CTThorn <> "_convertToPrimitives",
  Schedule -> {"in MoL_PostStep before SetTmunu after ADMBase_SetADMVars"},
  ConditionalOnKeyword -> {"formalism", "Wilson"},
  Where -> Everywhere,
  Shorthands -> {detg, gu[ua,ub], u0, h, rho0},
  Equations ->
  {
    detg      -> detgExpr,
    gu[ua,ub] -> 1/detg detgExpr MatrixInverse [g[ua,ub]],
    eps       -> EE / DD,
    h         -> (1 + eps) (1 + eosw),
    u[la]     -> SS[la] / (DD h),
    u0        -> Sqrt[gu[ua,ub] u[la] u[lb] + 1] / a,
    V[ua]     -> gu[ua, ub] u[lb] / u0 - b[ua],
    W         -> a u0 Sqrt[detg],
    rho0      -> DD / W,
    rho       -> rho0 (1 + eps),
    prs       -> eosw rho,
    rhodp     -> rho Sqrt[detg]
  }
};
       
convertToPrimitivesCalcVal =
{
  Name -> CTThorn <> "_convertToPrimitivesVal",
  Schedule -> {"in MoL_PostStep before SetTmunu after ADMBase_SetADMVars"},
  ConditionalOnKeyword -> {"formalism", "Valencia"},
  Where -> Everywhere,
  Shorthands -> {detg, gu[ua,ub], EEloc, SSloc[la], S2},
  Equations ->
  {
    detg      -> detgExpr,
    gu[ua,ub] -> 1/detg detgExpr MatrixInverse [g[ua,ub]],
      
    (* conservated variable without gdet *)
    EEloc -> EE / Sqrt[detg],
    SSloc[la] -> SS[la] / Sqrt[detg],
    S2 -> gu[ua,ub] SSloc[la] SSloc[lb],
      
    (* pressure *)
    prs -> (EEloc ( eosw - 1 ) + Sqrt[EEloc^2 ( eosw + 1 )^2 - 4 S2 eosw])/2,
    rho -> prs / eosw,
    rhodp -> rho Sqrt[detg],
    eps -> 0.0, (* neglect *)
      
    (* fluid velocity *)
    W -> Sqrt[(EEloc + prs) / (rho + prs)],
    u[la] -> SSloc[la] / ((rho + prs) W),
    V[ua] -> gu[ua,ub] u[lb] / W
  }
};

(******************************************************************************)
(* Evolution equations *)
(******************************************************************************)

evolCalc =
{
  Name -> CTThorn <> "_RHS",
  Schedule -> {"IN MoL_CalcRHS"},
  ConditionalOnKeyword -> {"formalism", "Wilson"},
  Where -> InteriorNoSync,
  Shorthands -> {dir[ua], detg, gu[ua,ub], u0u, u0l, dtD, dtW, dtS[la], dtL, SS0u, SS0l, h, rho0, SaSbdgab[la], ShiS0},
  Equations ->
  {
    dir[ua] -> Sign[b[ua]],

    detg      -> detgExpr,
    gu[ua,ub] -> 1/detg detgExpr MatrixInverse [g[ua,ub]],
    u0u       -> Sqrt[gu[ua,ub] u[la] u[lb] + 1] / a,
    u0l       -> -a^2 u0u + b[ua] u[la],
    h         -> (1 + eps) (1 + eosw),
    rho0      -> rho / (1 + eps),
    SS0u      -> W rho0 h u0u,
    SS0l      -> W rho0 h u0l,

    dtD -> - PD[DD V[ua], la],
    dot[DD] -> dtD,
      
    SaSbdgab[la] -> SS0l SS0l PD[-1/a^2, la]
                    + 2 SS0l SS[lc] PD[b[uc]/a^2, la]
                    + SS[lb] SS[lc] PD[MatrixInverse[g[ub,uc]], la],
    dtS[la] -> - PD[SS[la] V[ub], lb]
               - SaSbdgab[la] / (2 SS0u)
               - a Sqrt[detg] PD[prs, la],
    dot[SS[la]] -> dtS[la],
    
    dtL -> (-2 dtD gu[ua,ub] SS[la] SS[lb] / DD
            -2 a gu[ua,uc] gu[ub,ud] k[lc,ld] SS[la] SS[lb]
            + gu[ua,ub] (dtS[la] SS[lb]
                         + SS[la] dtS[lb])) / DD^2,
    dtW -> (sqrt[detg] dtL / (2 u0u a)) - a^2 u0u bssntrK sqrt[detg],
    dot[EE] -> - PD[EE V[ua], la] - prs dtW - prs PD[W V[ua], la]
  }
};
       
evolCalcVal =
{
  Name -> CTThorn <> "_RHSVal",
  Schedule -> {"IN MoL_CalcRHS"},
  ConditionalOnKeyword -> {"formalism", "Valencia"},
  Where -> InteriorNoSync,
  Shorthands -> {dir[ua], detg, gu[ua,ub], Sstress[ua,ub]},
  Equations ->
  {
    dir[ua] -> Sign[b[ua]],

    detg      -> detgExpr,
    gu[ua,ub] -> 1/detg detgExpr MatrixInverse [g[ua,ub]],
      
    Sstress[ua,ub] -> (rho + prs) W^2 V[ua] V[ub] + prs gu[ua,ub],

    dot[DD] -> 0.0, (* neglect *)
    dot[SS[la]] -> (- PD[a detgExpr g[la,lc] ((rho + prs) W^2 V[ub] V[uc]
                                              + prs MatrixInverse[g[ub,uc]]), lb]
                    + PD[b[ub] SS[la], lb]
                    + a Sqrt[detg] Sstress[ub,uc] PD[g[lb,lc], la] / 2
                    + SS[lb] PD[b[ub], la]
                    - EE PD[a, la]),
    dot[EE] -> (- PD[a MatrixInverse[g[ub,uc]] SS[lc], lb]
                + PD[b[ub] EE, lb]
                + a Sqrt[detg] Sstress[ub,uc] k[lb,lc]
                - gu[ub,uc] SS[lc] PD[a, lb])
  }
};

RHSStaticBoundaryCalc =
{
  Name -> CTThorn <> "_RHSStaticBoundary",
  Schedule -> {"IN MoL_CalcRHS"},
  ConditionalOnKeyword -> {"my_rhs_boundary_condition", "static"},
  Where -> Boundary,
  Equations -> 
  {
    dot[DD]       -> 0,
    dot[EE]       -> 0,
    dot[SS[la]]   -> 0
  }
};

(* Initialise the RHS variables in analysis in case they are going to
   be output - the noninterior points cannot be filled, so we define
   them to be zero *)
initRHSCalc =
{
  Name -> CTThorn <> "_InitRHS",
  Schedule -> {"AT analysis BEFORE " <> CTThorn <> "_RHS"},
  Where -> Everywhere,
  Equations -> 
  {
    dot[DD]       -> 0,
    dot[EE]       -> 0,
    dot[SS[la]]   -> 0
  }
};

(******************************************************************************)
(* Populate the energy-momentum tensor *)
(******************************************************************************)

addToTmunuCalc =
{
  Name -> CTThorn <> "_addToTmunu",
  Schedule -> {"IN AddToTmunu"},
  ConditionalOnKeyword -> {"coupling", "yes"},
  Where -> Everywhere,
  Shorthands -> {detg, gu[ua,ub], u0u, u0l, bsq, rhocrit, rhof},
  Equations ->
  {
    detg      -> detgExpr,
    gu[ua,ub] -> 1/detg detgExpr MatrixInverse [g[ua,ub]],
    u0u       -> Sqrt[gu[ua,ub] u[la] u[lb] + 1] / a,
    bsq       -> g[la,lb] b[ua] b[ub],
    u0l       -> -a^2 u0u + b[ua] u[la],

    T00       -> (rho + prs) u0l^2 + (prs - Lambda/(8 Pi)) (- a^2 + bsq),
    T0[la]    -> (rho + prs) u0l u[la] + (prs - Lambda/(8 Pi)) g[la,lb] b[ub],
    T[la,lb]  -> (rho + prs) u[la] u[lb] + (prs - Lambda/(8 Pi)) g[la,lb]
  }
};

setCTTrhs =
{
  Name -> CTThorn <> "_setCTTrhs",
  Schedule -> {"AT INITIAL before CT_MultiLevel"},
  Where -> Interior,
  Shorthands -> {detg, gu[ua,ub], u0u, u0l, bsq, t00, t0[la], t[la,lb]},
  Equations -> 
  {
    detg      -> detgExpr,
    gu[ua,ub] -> 1/detg detgExpr MatrixInverse [g[ua,ub]],
    u0u       -> Sqrt[gu[ua,ub] u[la] u[lb] + 1] / a,
    bsq       -> g[la,lb] b[ua] b[ub],
    u0l       -> -a^2 u0u + b[ua] u[la],

    t00       -> (rho + prs) u0l^2 + (prs - Lambda/(8 Pi)) (- a^2 + bsq),
    t0[la]    -> (rho + prs) u0l u[la] + (prs - Lambda/(8 Pi)) g[la,lb] b[ub],
    t[la,lb]  -> (rho + prs) u[la] u[lb] + (prs - Lambda/(8 Pi)) g[la,lb]
  }
};

(******************************************************************************)
(* Boundary conditions *)
(******************************************************************************)

boundaryCalc =
{
  Name -> CTThorn <> "_boundary",
  Schedule -> {"IN MoL_PostStep"},
  ConditionalOnKeyword -> {"my_boundary_condition", "Minkowski"},
  Where -> Boundary,
  Equations -> 
  {
    DD         -> 0,
    EE         -> 0,
    SS[la]     -> 0
  }
};

(******************************************************************************)
(******************************************************************************)
(******************************************************************************)
(* Calcultae expansion *)
(******************************************************************************)
(******************************************************************************)
(******************************************************************************)

calcExpansion =
{
  Name -> CTThorn <> "_Expansion",
  Schedule -> {"AT Poststep"},
  ConditionalOnKeyword -> {"calc_expansion", "yes"},
  Where -> InteriorNoSync,
  Shorthands -> {dir[ua], detg, gu[ua,ub], a2, bl[la], g4ss[ua,ub], ku[ua,ub],
                 dsa[la], dsb[la,ub], dsg[lc,la,lb], dsgu[lc,ua,ub], dsbl[la,lb],
                 Liebg[la,lb], Liebgu[ua,ub], dtg[la,lb], dtgu[ua,ub], dtbl[la],
                 Gdttt, Gdtts[la], Gdtss[la,lb], Gdstt[la], Gdsts[la,lb], Gdsss[la,lb,lc],
                 Gttt, Gtts[la], Gtss[la,lb], Gstt[ua], Gsts[ua,lc], Gsss[ua,lc,ld],
                 u0u, u0l, uup[ua], dtu[la], dtu0u, dtu0l, dsu[la,lb], dsu0[la],
                 Dtut, Gtsu[la], Dsut[la], Dtus[la], Dsus[la,lb],
                 htt, hts[la], hss[la,lb], hmtt, hmtsl[la], hmtsu[ua], hmss[ua,lb], hutt, huts[ua], huss[ua,ub],
                 Dptut, Dptus[lb], Dpsut[la], Dpsus[la,lb],
                 Thetatt, Thetats[la], Thetass[la,lb]},
  Equations ->
  {
      (******************** metric ********************)
      dir[ua] -> Sign[b[ua]],
      detg -> detgExpr,
      gu[ua,ub] -> 1/detg detgExpr MatrixInverse [g[ua,ub]],
      a2 -> a^2,
      bl[la] -> b[ub] g[la,lb],
      g4ss[ua,ub] -> a2 gu[ua,ub] - b[ua] b[ub],
      ku[ua,ub] -> gu[ua,uc] gu[ub,ud] k[lc,ld],
      
      (* space derivatives *)
      dsa[la] -> PD[a, la],
      dsb[la,ub] -> PD[b[ub], la],
      dsg[lc,la,lb] -> PD[g[la,lb], lc],
      dsgu[lc,ua,ub] -> PD[MatrixInverse[g[ua,ub]], lc],
      dsbl[la,lb] -> dsb[la,uc] g[lb,lc] + b[uc] dsg[la,lb,lc],
      
      (* time derivatives *)
      Liebg[la,lb] -> b[uc] dsg[lc,la,lb] + dsb[la,uc] g[lc,lb] + dsb[lb,uc] g[la,lc],
      Liebgu[ua,ub] -> b[uc] dsgu[lc,ua,ub] - dsb[lc,ua] gu[uc,ub] + dsb[lc,ub] gu[ua,uc],
      dtg[la,lb] -> Liebg[la,lb] - 2 a k[la,lb],
      dtgu[ua,ub] -> Liebgu[ua,ub] - 2 a ku[ua,ub],
      dtbl[la] -> dtb[ub] g[la,lb] + b[ub] dtg[la,lb],
      
      (* 4D Christoffel symbols ddd *)
      Gdttt           -> - a dta + (dtb[ua] bl[la] + b[ua] dtbl[la]) / 2,
      Gdtts[la]       -> - a dsa[la] + (dsb[la,ub] bl[lb] + b[ub] dsbl[la,lb]) / 2,
      Gdtss[la,lb]    -> (dsbl[la,lb] + dsbl[lb,la] - dtg[la,lb]) / 2,
      Gdstt[la]       -> dtbl[la] - Gdtts[la],
      Gdsts[la,lb]    -> dsbl[lb,la] - Gdtss[la,lb],
      Gdsss[la,lb,lc] -> (dsg[lb,la,lc] + dsg[lc,lb,la] - dsg[la,lb,lc]) / 2,
      
      (* 4D Christoffel symbols udd *)
      Gttt           -> (-Gdttt + b[uc] Gdstt[lc]) / a2,
      Gtts[la]       -> (-Gdtts[la] + b[uc] Gdsts[lc,la]) / a2,
      Gtss[la,lb]    -> (-Gdtss[la,lb] + b[uc] Gdsss[lc,la,lb]) / a2,
      Gstt[ua]       -> (b[ua] Gdttt + g4ss[ua,ub] Gdstt[lb]) / a2,
      Gsts[ua,lc]    -> (b[ua] Gdtts[lc] + g4ss[ua,ub] Gdsts[lb,lc]) / a2,
      Gsss[ua,lc,ld] -> (b[ua] Gdtss[lc,ld] + g4ss[ua,ub] Gdsss[lb,lc,ld]) / a2,
      
      (******************** fluid velocity ********************)
      u0u -> Sqrt[gu[ua,ub] u[la] u[lb] + 1] / a,
      u0l -> -a2 u0u + b[ua] u[la],
      uup[ua] -> (b[ua] u0l + g4ss[ua,ub] u[lb]) / a2,
      
      (* time derivative *)
      dtu[la] -> ( (dot[SS[la]] / (1 + eosw))  - u[la] (dot[DD] + dot[EE]) ) / (DD + EE),
      dtu0u -> - (dta u0u / a) + ((dtgu[ua,ub] u[la] u[lb]
                                   + gu[ua,ub] dtu[la] u[lb]
                                   + gu[ua,ub] u[la] dtu[lb]) / (2 a2 u0u)),
      dtu0l -> -2 a dta u0u - a2 dtu0u + dtb[ua] u[la] + b[ua] dtu[la],
      
      (* space derivative *)
      dsu[la,lb] -> PD[u[lb], la],
      dsu0[la] -> - 2 a dta u0u - a2 dtu0u + dtb[ua] u[la] + b[ua] dtu[la],
      
      (* Covariant derivative *)
      Dtut -> dtu0l - Gttt u0l - Gstt[ua] u[la],
      Gtsu[la] -> Gtts[la] u0l + Gsts[uc,la] u[lc],
      Dsut[la] -> dsu0[la] - Gtsu[la],
      Dtus[la] -> dtu[la] - Gtsu[la],
      Dsus[la,lb] -> dsu[la,lb] - Gtss[la,lb] u0l - Gsss[uc,la,lb] u[lc],
      
      (******************** fluid projection ********************)
      (* indices down *)
      htt -> (-a2 + g[la,lb] b[ua] b[ub]) + u0l u0l,
      hts[la] -> bl[la] + u0l u[la],
      hss[la,lb] -> g[la,lb] + u[la] u[lb],
      
      (* mixed indices *)
      hmtt -> 1 + u0u u0l,
      hmtsl[la] -> u0u u[la], (* t up *)
      hmtsu[ua] -> u0l uup[ua], (* t down *)
      hmss[ua,lb] -> KroneckerDelta[ua,lb] + uup[ua] u[lb],
      
      (* indices up *)
      hutt -> - (1/a2) + u0u u0u,
      huts[ua] -> (b[ua]/a2) + u0u uup[ua],
      huss[ua,ub] -> (g4ss[ua,ub]/a2) + uup[ua] uup[ub],
      
      (* projected covariant derivative *)
      (* Dpu[la,lb] -> (hmt[la] (hmt[lb] Dtut + hms[ud,lb] Dtus[ld])
                        + hms[uc,la] (hmt[lb] Dsut[lc] + hms[ud,lb] Dsus[lc,ld])), *)
      Dptut -> (hmtt (hmtt Dtut + hmtsu[ud] Dtus[ld])
                + hmtsu[uc] (hmtt Dsut[lc] + hmtsu[ud] Dsus[lc,ld])),
      Dptus[lb] -> (hmtt (hmtsl[lb] Dtut + hmss[ud,lb] Dtus[ld])
                    + hmtsu[uc] (hmtsl[lb] Dsut[lc] + hmss[ud,lb] Dsus[lc,ld])),
      Dpsut[la] -> (hmtsl[la] (hmtt Dtut + hmtsu[ud] Dtus[ld])
                    + hmss[uc,la] (hmtt Dsut[lc] + hmtsu[ud] Dsus[lc,ld])),
      Dpsus[la,lb] -> (hmtsl[la] (hmtsl[lb] Dtut + hmss[ud,lb] Dtus[ld])
                       + hmss[uc,la] (hmtsl[lb] Dsut[lc] + hmss[ud,lb] Dsus[lc,ld])),
      
      (******************** Expansion ********************)
      Thetatt -> Dptut,
      Thetats[la] -> (Dptus[la] + Dpsut[la]) / 2,
      Thetass[la,lb] -> (Dpsus[la,lb] + Dpsus[lb,la]) / 2,
      
      (* expansion scalar *)
      Theta -> hutt Thetatt + 2 huts[ua] Thetats[la] + huss[ua,ub] Thetass[la,lb],
      
      (* shear *)
      sigmatt -> Thetatt - htt Theta / 3,
      sigmats[la] -> Thetats[la] - hts[la] Theta / 3,
      sigmass[la,lb] -> Thetass[la,lb] - hss[la,lb] Theta / 3
  }
};

calcExpansionBoundary =
{
  Name -> CTThorn <> "_ExpansionBoundary",
  Schedule -> {"AT Poststep"},
  ConditionalOnKeyword -> {"calc_expansion", "yes"},
  Where -> Boundary,
  Equations ->
  {
      Theta -> 0,
      sigmatt -> 0,
      sigmats[la] -> 0,
      sigmass[la,lb] -> 0
  }
};

(******************************************************************************)
(******************************************************************************)
(******************************************************************************)
(* Parameters *)
(******************************************************************************)
(******************************************************************************)
(******************************************************************************)


extendedKeywordParameters =
{
  {
    Name -> "ADMBase::initial_data",
    AllowedValues -> {CTThorn}
  },
  {
    Name -> "ADMBase::initial_lapse",
    AllowedValues -> {CTThorn}
  },
  {
    Name -> "ADMBase::initial_shift",
    AllowedValues -> {CTThorn} 
  }
};

keywordParameters =
{
  {
    Name -> "my_initial_data",
    Visibility -> "restricted",
    AllowedValues -> {"Minkowski", "pFLRW"},
    Default -> "Minkowski"
  },
    {
      Name -> "formalism",
      Visibility -> "restricted",
      Description -> "Wilson for dust and Valencia for radiation",
      AllowedValues -> {"Wilson", "Valencia"},
      Default -> "Wilson"
    },
  {
    Name -> "my_rhs_boundary_condition",
    Visibility -> "restricted",
    AllowedValues -> {"none", "static", "radiative"},
    Default -> "none"
  },
  {
    Name -> "my_boundary_condition",
    AllowedValues -> {"none", "Minkowski"},
    Default -> "none"
  },
  {
    Name -> "coupling",
    AllowedValues -> {"yes", "no"},
    Default -> "no"
  },
  {
    Name -> "calc_expansion",
    AllowedValues -> {"yes", "no"},
    Default -> "no"
  }
};

intParameters =
{
  {
    Name -> nmodes,
    Description -> "Number of Fourier modes in each direction", 
    Default -> 0
  },
  {
    Name -> fdOrder,
    Default -> derivOrder,
    AllowedValues -> {2,4,6,8}
  }
};

realParameters =
{
  {
    Name -> eosw,
    Description -> "Ratio of pressure to density",
    Default -> 0
  },
  {
    Name -> Lambda,
    Description -> "Cosmological constant",
    Default -> 0
  },
  {
    Name -> "pflrw_t0",
    Description -> "Initial time in RW models",
    Default -> 1
  },
  {
    Name -> "pflrw_a0",
    Description -> "Initial value of the scale factor in RW models",
    Default -> 1
  },
  {
    Name -> "pflrw_H0",
    Description -> "Initial value of the Hubble rate in RW models",
    Default -> 1
  },
  {
    Name -> "pflrw_omegaM",
    Description -> "Initial value of the matter density parameter",
    Default -> 1(*,
    AllowedValues -> {1}*)
  },
  {
    Name -> "pflrw_kx[20]",
    Description -> "Wave numbers of the x-modes",
    Default -> 0
  },
  {
    Name -> "pflrw_ky[20]",
    Description -> "Wave numbers of the y-modes",
    Default -> 0
  },
  {
    Name -> "pflrw_kz[20]",
    Description -> "Wave numbers of the z-modes",
    Default -> 0
  },
  {
    Name -> "pflrw_ax[20]",
    Description -> "Amplitudes of the x-modes",
    Default -> 0
  },
  {
    Name -> "pflrw_ay[20]",
    Description -> "Amplitudes of the y-modes",
    Default -> 0
  },
  {
    Name -> "pflrw_az[20]",
    Description -> "Amplitudes of the z-modes",
    Default -> 0
  }
};

(******************************************************************************)
(******************************************************************************)
(******************************************************************************)
(* Construct the thorns *)
(******************************************************************************)
(******************************************************************************)
(******************************************************************************)

calculations =
{
  initialMinkowskiCalc,
  initialFLRWPertCalc,
  convertFromPrimitivesCalc,
  convertFromPrimitivesCalcVal,
  evolCalc,
  evolCalcVal,
  initRHSCalc,
  RHSStaticBoundaryCalc,
  boundaryCalc,
  convertToPrimitivesCalc,
  convertToPrimitivesCalcVal,
  addToTmunuCalc,
  setCTTrhs,
  calcExpansion,
  calcExpansionBoundary
};

CreateKrancThornTT [groups, ".", CTThorn,
  Calculations -> calculations,
  DeclaredGroups -> declaredGroupNames,
  PartialDerivatives -> derivatives,
  EvolutionTimelevels -> maxTimelevels,
  DefaultEvolutionTimelevels -> Min[3,maxTimelevels],
  UseJacobian -> True,
  UseLoopControl -> True,
  UseVectors -> False,
  InheritedImplementations -> inheritedImplementations,
  ExtendedKeywordParameters -> extendedKeywordParameters,
  KeywordParameters -> keywordParameters,
  IntParameters -> intParameters,
  RealParameters -> realParameters
];