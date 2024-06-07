SetEnhancedTimes[False];

thorn = "CosmoLapse";
prefix = "CL_";

(******************************************************************************)
(******************************************************************************)
(******************************************************************************)
(* Tensors *)
(******************************************************************************)
(******************************************************************************)
(******************************************************************************)

DefineTensor1[name_, indices_:{}, syms_:{}] :=
  Module[{tensor},
         DefineTensor[name];
         tensor = If[indices=={}, name, name[Sequence@@indices]];
         AssertSymmetricIncreasing[tensor, Sequence@@#]& /@ syms;
         tensor];

DefineTensor1[dir, {ua}];

(* Local *)
DefineTensor1[tau];
DefineTensor1[tauini];
DefineTensor1[Ktrace];
DefineTensor1[Kini];
DefineTensor1[Ktraceaverage];
DefineTensor1[Ktraceaverageini];
DefineTensor1[arhoini];

(* Spatial metric and extrinsic curvature *)
DefineTensor1[g, {la,lb}, {{la,lb}}];
DefineTensor1[gu, {ua,ub}, {{ua,ub}}];
DefineTensor1[dtgu, {ua,ub}, {{ua,ub}}];
DefineTensor1[k, {la,lb}, {{la,lb}}];
Map [AssertSymmetricIncreasing,
     {g[la,lb], k[la,lb]}];
Map [AssertSymmetricIncreasing,
     {gu[ua,ub]}];

(* Lapse *)
DefineTensor1[a];
DefineTensor1[admdtalpha];
DefineTensor1[mlbalpharhs];
DefineTensor1[partialtalpha];

(* Shift *)
DefineTensor1[b, {ua}];
DefineTensor1[admdtbeta, {ua}];
DefineTensor1[mlbbetarhs, {ua}];
DefineTensor1[partialtbeta, {ua}];
DefineTensor1[mlbB, {ua}];
DefineTensor1[mlbBrhs, {ua}];
DefineTensor1[partialtB, {ua}];

(* Gamma *)
DefineTensor1[mlbXt, {ua}];
DefineTensor1[mlbXtrhs, {ua}];
SetTensorAttribute[mlbXt,  TensorWeight, +2/3];

(* ADMBase variables *)
g11=gxx; g12=gxy; g22=gyy; g13=gxz; g23=gyz; g33=gzz;
k11=kxx; k12=kxy; k22=kyy; k13=kxz; k23=kyz; k33=kzz;
a=alp;
admdtalpha=dtalp;
b1=betax; b2=betay; b3=betaz;
admdtbeta1=dtbetax; admdtbeta2=dtbetay; admdtbeta3=dtbetaz;
 
(* ML_BSSN variables *)
mlbalpha=alpha;
mlbalpharhs=alpharhs;
mlbbetarhs1=beta1rhs; mlbbetarhs2=beta2rhs; mlbbetarhs3=beta3rhs;

mlbB1=B1; mlbB2=B2; mlbB3=B3;
mlbBrhs1=B1rhs; mlbBrhs2=B2rhs; mlbBrhs3=B3rhs;

mlbXt1=Xt1; mlbXt2=Xt2; mlbXt3=Xt3;
mlbXtrhs1=Xt1rhs; mlbXtrhs2=Xt2rhs; mlbXtrhs3=Xt3rhs;

(* Fluid *)
DefineTensor1[rhoL];
rhoL = rho;
DefineTensor1[DDL];
DDL = DD;
DefineTensor1[SSL, {la}];
SSL1 = SS1; SSL2 = SS2; SSL3 = SS3;
DefineTensor1[u, {la}];
DefineTensor1[SaSbdgab, {la}];
DefineTensor1[dtS, {la}];
DefineTensor1[VL, {ua}];
VL1 = V1; VL2 = V2; VL3 = V3;

(******************************************************************************)
(******************************************************************************)
(******************************************************************************)
(* Groups *)
(******************************************************************************)
(******************************************************************************)
(******************************************************************************)

inheritedImplementations =
   {"ADMBase", "ML_BSSN", "Cactus", "ICPertFLRW", "LocalReduce"};

ThornGroups = {SetGroupName [CreateGroupFromTensor [tau], "propertime"],
               SetGroupName [CreateGroupFromTensor [tauini], "propertimeinitial"],
               SetGroupName [CreateGroupFromTensor [Ktrace], "CL_Ktrace"],
               SetGroupName [CreateGroupFromTensor [Kini], "Kinitial"],
               SetGroupName [CreateGroupFromTensor [Ktraceaverage], "CL_Ktraceaverage"],
               SetGroupName [CreateGroupFromTensor [Ktraceaverageini], "CL_Ktraceaverageinitial"],
               SetGroupName [CreateGroupFromTensor [arhoini], "arhoinitialfactor"]};

extraGroups =
{
    {"ADMBase::metric",        {gxx, gxy, gxz, gyy, gyz, gzz}},
    {"ADMBase::curv",          {kxx, kxy, kxz, kyy, kyz, kzz}},
    {"ADMBase::lapse",         {alp}},
    {"ADMBase::dtlapse",       {dtalp}},
    {"ADMBase::shift",         {betax, betay, betaz}},
    {"ADMBase::dtshift",       {dtbetax, dtbetay, dtbetaz}},
    {"ML_BSSN::ML_lapse",      {alpha}},
    {"ML_BSSN::ML_lapserhs",   {alpharhs}},
    {"ML_BSSN::ML_shiftrhs",   {beta1rhs, beta2rhs, beta3rhs}},
    {"ML_BSSN::ML_dtshift",    {B1, B2, B3}},
    {"ML_BSSN::ML_dtshiftrhs", {B1rhs, B2rhs, B3rhs}},
    {"ML_BSSN::ML_Gamma",      {Xt1, Xt2, Xt3}},
    {"ML_BSSN::ML_Gammarhs",   {Xt1rhs, Xt2rhs, Xt3rhs}},
    {"CT_Dust::CT_rho",        {rho}},
    {"CT_Dust::CT_D",          {DD}},
    {"CT_Dust::CT_S",          {SS1, SS2, SS3}},
    {"CT_Dust::CT_V",          {V1, V2, V3}}
};

declaredGroupNames = Map [First, ThornGroups];
groups = Join [ThornGroups, extraGroups];

(******************************************************************************)
(******************************************************************************)
(******************************************************************************)
(* Definitions *)
(******************************************************************************)
(******************************************************************************)
(******************************************************************************)

(* Functions & Definitions *)

detgExpr = Det[MatrixOfComponents[g[la,lb]]];

(* Derivatives *)

SCDO = StandardCenteredDifferenceOperator;
SUDO = StandardUpwindDifferenceOperator;

SUDO1[p_, m1_Integer, m2_Integer, i_Integer] :=
  Module[{dir},
         dir[j_Integer] := Symbol["dir"<>ToString[j]];
         dir[i] SUDO[p, m1, m2, i] /. shift[j_] -> shift[j]^dir[j]];

SUDOsymm[p_, m1_Integer, m2_Integer, i_Integer] :=
  1/2 (SUDO[p, m1, m2, i] - SUDO[p, m2, m1, i])
SUDOanti[p_, m1_Integer, m2_Integer, i_Integer] :=
  1/2 (SUDO[p, m1, m2, i] + SUDO[p, m2, m1, i])

derivatives = {
    PDstandardNth[i_]    -> SCDO[1, fdOrder/2, i],
    PDstandardNth[i_,i_] -> SCDO[2, fdOrder/2, i],
    PDstandardNth[i_,j_] -> SCDO[1, fdOrder/2, i] SCDO[1, fdOrder/2, j],

    PDupwindNth[i_]     -> SUDO1[1, fdOrder/2-1, fdOrder/2+1, i],
    PDupwindNthSymm[i_] -> SUDOsymm[1, fdOrder/2-1, fdOrder/2+1, i],
    PDupwindNthAnti[i_] -> SUDOanti[1, fdOrder/2-1, fdOrder/2+1, i],

    (* Note: for stability (and to reduce the stencil radius), we lower
       the order of accuracy here *)
    PDonesided[i_] -> SUDO1[1, 0, fdOrder/2+1, i],

    PDdissipationNth[i_] -> ((-1)^(fdOrder/2)
                             spacing[i]^(fdOrder+1) / 2^(fdOrder+2)
                             SCDO[fdOrder+2, fdOrder/2+1, i])};

PD  = PDstandardNth;
PDu = PDupwindNth;
PDua   = PDupwindNthAnti;
PDus   = PDupwindNthSymm;
PDo    = PDonesided;
PDdiss = PDdissipationNth;

splitUpwindDerivsKranc = True;
If[splitUpwindDerivsKranc,
   Upwind[dir_, var_, idx_] := dir PDua[var,idx] + Abs[dir] PDus[var,idx],
   Upwind[dir_, var_, idx_] := dir PDu[var,idx]];
 
 (* Split a calculation *)
 PartialCalculation[calc_, suffix_, updates_, vars_] :=
   Module[
     {name, calc1, replacements, calc2, vars1, patterns, eqs, calc3},
     (* Add suffix to name *)
     name  = lookupDefault[calc, Name, ""] <> suffix;
     calc1 = mapReplaceAdd[calc, Name, name];
     (* Replace some entries in the calculation *)
     replacements = updates //. (lhs_ -> rhs_) -> (mapReplaceAdd[#, lhs, rhs]&);
     calc2        = calc1 // Composition@@replacements;
     (* Remove unnecessary equations *)
     vars1    = Join[vars, lookupDefault[calc2, Shorthands, {}]];
     patterns = Replace[vars1, {Tensor[n_,__]      ->     Tensor[n,__] ,
                                dot[Tensor[n_,__]] -> dot[Tensor[n,__]]}, 1];
     eqs      = FilterRules[lookup[calc, Equations], patterns];
     calc3    = mapReplace[calc2, Equations, eqs];
     calc3];

(******************************************************************************)
(******************************************************************************)
(******************************************************************************)
(* Master calculations *)
(******************************************************************************)
(******************************************************************************)
(******************************************************************************)

(***** Gauge evolution equations *****)

MasterCalc = {
    Shorthands -> {dir[ua], detg, gu[ua,ub],
                   partialtalpha, Kexp, Kn,
                   partialtbeta[ua], partialtB[ua]},
    Equations -> {
        dir[ua]   -> Sign[b[ua]],
      
        (***** alpha *****)
        Kexp -> IfThen[Kexpression==1, - 2.0 / tauini,
                       IfThen[Kexpression==2, - 2.0 / tau,
                              IfThen[Kexpression==3, Kini tauini / tau,
                                     IfThen[Kexpression==4, Ktraceaverage,
                                            IfThen[Kexpression==5, Kini Ktraceaverage / Ktraceaverageini,
                                                   0.0]]]]],
        Kn -> IfThen[Knorm==1, 2.0 / tau,
                     IfThen[Knorm==2, tauini / tau,
                            IfThen[Knorm==3, Ktraceaverage / Ktraceaverageini,
                                   1.0]]],
        partialtalpha  -> (- alphaF a^alphaN (Ktrace - Kexp) / Kn
                      + IfThen[alphaFullLieDeriv!=0, Upwind[b[ub], a, lb], 0]),
        
        (***** beta *****)
        partialtbeta[ua] -> mlbB[ua] + IfThen[betaFullLieDeriv!=0, Upwind[b[ub], b[ua], lb], 0],
        partialtB[ua] -> (betaXi a^betaP
                          (mlbXtrhs[ua] - IfThen[betaFullLieDeriv!=0, Upwind[b[ub], mlbXt[ua], lb], 0])
                          - betaEta mlbB[ua]
                          + IfThen[betaFullLieDeriv!=0, Upwind[b[ub], mlbB[ua], lb], 0]),

        (***** Update ADMBase and ML_BSSN *****)
        admdtalpha  -> partialtalpha,
        mlbalpharhs -> partialtalpha,
        admdtbeta[ua] -> partialtbeta[ua],
        mlbbetarhs[ua] -> partialtbeta[ua],
        mlbBrhs[ua] -> partialtB[ua]
    }};

MasterRhoCalc = {
    Shorthands -> {dir[ua], detg, gu[ua,ub], dtgu[ua,ub],
                   dttau, u[la], u0u, u0l, dtD,
                   SS0u, SS0l, SaSbdgab[la], dtS[la],
                   Snorm, dtSnorm, dtF, F,
                   partialtalpha},
    Equations -> {
        dir[ua]   -> Sign[b[ua]],
        detg      -> detgExpr,
        gu[ua,ub] -> detgExpr/detg MatrixInverse[g[ua,ub]],
        dtgu[ua, ub] -> - 2 a gu[ua, uc] gu[ub, ud] k[lc, ld], (* [TO DO: add beta] *)
        
        dttau -> (a^2 - b[ua] b[ub] g[la,lb])^(1/2),
        u[la] -> SSL[la] / DDL,  (* [TO DO: add h] *)
        u0u -> Sqrt[gu[ua,ub] u[la] u[lb] + 1] / a,
        u0l -> -a^2 u0u + b[ua] u[la],
        dtD -> - PD[DDL VL[ua], la],
        
        (* velocity *)
        SS0u -> DDL u0u,  (* [TO DO: add h] *)
        SS0l -> DDL u0l,  (* [TO DO: add h] *)
        SaSbdgab[la] -> SS0l SS0l PD[-1/a^2, la]
                     + 2 SS0l SSL[lc] PD[b[uc]/a^2, la]
                     + SSL[lb] SSL[lc] PD[MatrixInverse[g[ub,uc]], la],
        dtS[la] -> - PD[SSL[la] VL[ub], lb] - SaSbdgab[la] / (2 SS0u),  (* [TO DO: add pressure] *)
        
        Snorm -> gu[ua, ub] SSL[la] SSL[lb],
        dtSnorm -> dtgu[ua, ub] SSL[la] SSL[lb] + gu[ua, ub] dtS[la] SSL[lb] + gu[ua, ub] SSL[la] dtS[lb],
        
        dtF -> dtSnorm + 2 dtD DDL,
        F -> Snorm + DDL DDL,
        
        partialtalpha -> - a ( 2 (dttau/tau) + 2 (dtD/DDL) + a Ktrace - (1/2) (dtF/F) ) / alphaN,
        (* [TO DO: add beta] *)

        (***** Update ADMBase and ML_BSSN *****)
        admdtalpha  -> partialtalpha,
        mlbalpharhs -> partialtalpha
    }};

(********************* Lapse *********************)
 
dtLapsePostStep = PartialCalculation[
    MasterCalc, "",
     {
         Name -> thorn <> "_dtLapsePostStep",
         Schedule -> {"IN MoL_PostStep AFTER ML_BSSN_ADMBaseBoundaryScalar BEFORE ADMBase_SetADMVars"},
         ConditionalOnKeyword -> {"dtalpha", "yes"},
         Where -> InteriorNoSync
     },
     {admdtalpha, mlbalpharhs}];
 
dtLapsePostRHS = PartialCalculation[
    MasterCalc, "",
     {
         Name -> thorn <> "_dtLapsePostRHS",
         Schedule -> {"IN MoL_PostRHS"},
         ConditionalOnKeyword -> {"dtalpha", "yes"},
         Where -> InteriorNoSync
     },
     {admdtalpha, mlbalpharhs}];

dtLapsePostStepBoundary = 
    {
         Name -> thorn <> "_dtLapsePostStepBoundary",
         Schedule -> {"IN MoL_PostStep AFTER " <> thorn <> "_dtLapsePostStep"},
         Where -> Boundary,
         Equations ->
         {
            admdtalpha -> 0.0,
            mlbalpharhs -> 0.0
         }
     };
 
dtLapsePostRHSBoundary =
    {
         Name -> thorn <> "_dtLapsePostRHSBoundary",
         Schedule -> {"IN MoL_PostRHS AFTER " <> thorn <> "_dtLapsePostRHS"},
         Where -> Boundary,
         Equations ->
         {
            admdtalpha -> 0.0,
            mlbalpharhs -> 0.0
         }
     };

AlphaFuncRho =
{
  Name -> thorn <> "_LapseFuncRho",
  Schedule -> {"IN MoL_PostStep BEFORE ML_BSSN_ADMBaseEverywhere"},
  ConditionalOnKeyword -> {"alpha_func_of_rho", "yes"},
  Where -> Everywhere,
  Equations ->
  {
    a -> ( arhoini / (tau tau rhoL) ) ^ ( 1 / alphaN ),
    mlbalpha -> a
  }
};

dtLapseRhoPostStep = PartialCalculation[
    MasterRhoCalc, "",
     {
         Name -> thorn <> "_dtLapseRhoPostStep",
         Schedule -> {"IN MoL_PostStep AFTER ML_BSSN_ADMBaseBoundaryScalar BEFORE ADMBase_SetADMVars"},
         ConditionalOnKeyword -> {"alpha_func_of_rho", "yes"},
         Where -> InteriorNoSync
     },
     {admdtalpha, mlbalpharhs}];
 
dtLapseRhoPostRHS = PartialCalculation[
    MasterRhoCalc, "",
     {
         Name -> thorn <> "_dtLapseRhoPostRHS",
         Schedule -> {"IN MoL_PostRHS"},
         ConditionalOnKeyword -> {"alpha_func_of_rho", "yes"},
         Where -> InteriorNoSync
     },
     {admdtalpha, mlbalpharhs}];

(********************* Shift *********************)

dtShiftPostStep = PartialCalculation[
    MasterCalc, "",
     {
         Name -> thorn <> "_dtShiftPostStep",
         Schedule -> {"IN MoL_PostStep AFTER ML_BSSN_ADMBaseBoundaryScalar BEFORE ADMBase_SetADMVars"},
         ConditionalOnKeyword -> {"dtshift", "yes"},
         Where -> InteriorNoSync
     },
     {admdtbeta[ua], mlbbetarhs[ua], mlbBrhs[ua]}];
 
dtShiftPostRHS = PartialCalculation[
    MasterCalc, "",
     {
         Name -> thorn <> "_dtShiftPostRHS",
         Schedule -> {"IN MoL_PostRHS"},
         ConditionalOnKeyword -> {"dtshift", "yes"},
         Where -> InteriorNoSync
     },
     {admdtbeta[ua], mlbbetarhs[ua], mlbBrhs[ua]}];

dtShiftPostStepBoundary =
    {
         Name -> thorn <> "_dtShiftPostStepBoundary",
         Schedule -> {"IN MoL_PostStep AFTER " <> thorn <> "_dtShiftPostStep"},
         ConditionalOnKeyword -> {"dtshift", "yes"},
         Where -> Boundary,
         Equations ->
         {
            admdtbeta[ua] -> 0.0,
            mlbbetarhs[ua] -> 0.0,
            mlbBrhs[ua] -> 0.0
         }
     };
 
dtShiftPostRHSBoundary =
    {
         Name -> thorn <> "_dtShiftPostRHSBoundary",
         Schedule -> {"IN MoL_PostRHS AFTER " <> thorn <> "_dtShiftPostRHS"},
         ConditionalOnKeyword -> {"dtshift", "yes"},
         Where -> Boundary,
         Equations ->
         {
            admdtbeta[ua] -> 0.0,
            mlbbetarhs[ua] -> 0.0,
            mlbBrhs[ua] -> 0.0
         }
     };

(******************* Proper time *******************)

initialise =
{
  Name -> thorn <> "_InitialTau",
  Schedule -> {"IN CCTK_INITIAL"},
  Where -> Everywhere,
  Shorthands -> {detg, gu[ua,ub], sum, count},
  Equations ->
  {
    tauini -> t,
    tau -> tauini,
      
    detg      -> detgExpr,
    gu[ua,ub] -> detgExpr/detg MatrixInverse[g[ua,ub]],
    Kini -> gu[ua,ub] k[la,lb],
    arhoini -> tauini tauini rhoL,
    sum -> Kini,
    count -> 1.0,
    Ktraceaverageini -> sum / count
  }
};

evolCalc =
{
  Name -> thorn <> "_RHSTau",
  Schedule -> {"IN MoL_CalcRHS"},
  Where -> Everywhere,
  Equations ->
  {
    dot[tau] -> (a^2 - b[ua] b[ub] g[la,lb])^(1/2)
  }
};
       
initRHSCalc =
{
  Name -> thorn <> "_InitRHSTau",
  Schedule -> {"AT analysis BEFORE " <> thorn <> "_RHSTau"},
  Where -> Everywhere,
  Equations ->
  {
    dot[tau] -> 0
  }
};

(******************* Proper time *******************)
traceofK =
{
  Name -> thorn <> "_TraceK",
  Schedule -> {"IN MoL_PostStep BEFORE ML_BSSN_ADMBaseBoundaryScalar"},
  Where -> Everywhere,
  Shorthands -> {detg, gu[ua,ub]},
  Equations ->
  {
    detg      -> detgExpr,
    gu[ua,ub] -> detgExpr/detg MatrixInverse[g[ua,ub]],
    Ktrace -> gu[ua,ub] k[la,lb]
  }
};

averageofK =
{
  Name -> thorn <> "_AverageK",
  Schedule -> {"IN MoL_PostStep BEFORE ML_BSSN_ADMBaseBoundaryScalar AFTER " <> thorn <> "_TraceK"},
  Where -> Everywhere,
  Shorthands -> {sum, count},
  Equations ->
  {
      sum -> Ktrace,
      count -> 1.0,
      Ktraceaverage -> sum / count
  }
};
 
calculations = {
    dtLapsePostStep,
    dtLapsePostRHS,
    dtLapsePostStepBoundary,
    dtLapsePostRHSBoundary,
    dtLapseRhoPostStep,
    dtLapseRhoPostRHS,
    dtShiftPostStep,
    dtShiftPostRHS,
    dtShiftPostStepBoundary,
    dtShiftPostRHSBoundary,
    traceofK,
    averageofK,
    initialise,
    evolCalc,
    initRHSCalc};

(******************************************************************************)
(******************************************************************************)
(******************************************************************************)
(* Parameters *)
(******************************************************************************)
(******************************************************************************)
(******************************************************************************)

keywordParameters =
{
  {
    Name -> "dtalpha",
    Visibility -> "restricted",
    AllowedValues -> {"yes", "no"},
    Default -> "no"
  },
  {
    Name -> "dtshift",
    Visibility -> "restricted",
    AllowedValues -> {"yes", "no"},
    Default -> "no"
  },
  {
    Name -> "alpha_func_of_rho",
    AllowedValues -> {"yes", "no"},
    Default -> "no"
  }
};

realParameters = {
    {
        Name -> alphaF,
        Description -> "d/dt alpha = - F alpha^N (K - Kexpression)",
        Default -> 0.0
    },
    {
        Name -> betaXi,
        Description -> "d/dt B = alpha^P Xi d/dt Xt^i - Eta B^i",
        Default -> 0.0
    },
    {
        Name -> betaEta,
        Description -> "d/dt B = alpha^P Xi d/dt Xt^i - Eta B^i",
        Default -> 0.0
    }};
  
intParameters = {
    {
        Name -> alphaN,
        Description -> "d/dt alpha = - F alpha^N (K - Kexpression) / Knorm  OR  alpha = ( tauIN^2 rhoIN / tau^2 rho )^( 1/N ) where it can not be 1",
        Default -> 1
    },
    {
       Name -> Kexpression,
       Description -> "d/dt alpha = - F alpha^N (K - Kexpression) / Knorm",
       AllowedValues -> {{Value -> 0, Description -> "zero"},
                         {Value -> 1, Description -> "initial_background_K"},
                         {Value -> 2, Description -> "background_K"},
                         {Value -> 3, Description -> "initial_K * initial_tau / tau"},
                         {Value -> 4, Description -> "mean_K"},
                         {Value -> 5, Description -> "initial_K * initial_mean_K / mean_K"}},
       Default -> 0
    },
    {
       Name -> Knorm,
       Description -> "d/dt alpha = - F alpha^N (K - Kexpression) / Knorm",
       AllowedValues -> {{Value -> 0, Description -> "Knorm = 1"},
                         {Value -> 1, Description -> "Knorm = - background_K"},
                         {Value -> 2, Description -> "Knorm = initial_tau / tau"},
                         {Value -> 3, Description -> "Knorm = initial_mean_K / mean_K"}},
       Default -> 0
    },
    {
        Name -> alphaFullLieDeriv,
        Description -> "Full Lie derivative, include beta^k partial_k term in alpha evo eq",
        AllowedValues -> {{Value -> 0, Description -> "no, off"},
                          {Value -> 1, Description -> "yes, on"}},
        Default -> 0
     },
     {
        Name -> betaP,
        Description -> "d/dt B = alpha^P Xi d/dt Xt^i - Eta B^i",
        Default -> 1
     },
     {
        Name -> betaFullLieDeriv,
        Description -> "Full Lie derivative, include beta^k partial_k term in beta evo eq",
        AllowedValues -> {{Value -> 0, Description -> "no off"},
                          {Value -> 1, Description -> "yes on"}},
        Default -> 0
     },
     {
        Name -> fdOrder,
        Description -> "Finite differencing order",
        AllowedValues -> {2, 4, 6, 8},
        Default -> 4
     }};

(******************************************************************************)
(******************************************************************************)
(******************************************************************************)
(* Construct the thorn *)
(******************************************************************************)
(******************************************************************************)
(******************************************************************************)

CreateKrancThornTT[
  groups, ".", thorn,
  Calculations -> calculations,
  PartialDerivatives -> derivatives,
  DeclaredGroups -> declaredGroupNames,
  EvolutionTimelevels -> 4,
  DefaultEvolutionTimelevels -> 3,
  UseJacobian -> True,
  UseLoopControl -> True,
  UseVectors -> False,
  InheritedImplementations -> inheritedImplementations,
  KeywordParameters -> keywordParameters,
  IntParameters -> intParameters,
  RealParameters -> realParameters
];
