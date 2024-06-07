SetEnhancedTimes[False];

(******************************************************************************)
(******************************************************************************)
(******************************************************************************)
(* Options *)
(******************************************************************************)
(******************************************************************************)
(******************************************************************************)

createCode[derivOrder_, useJacobian_, evolutionTimelevels_] :=
Module[{},
ThornName = "CosmoLapse";
prefix = "CL_";

(******************************************************************************)
(******************************************************************************)
(******************************************************************************)
(* Tensors *)
(******************************************************************************)
(******************************************************************************)
(******************************************************************************)

Map [DefineTensor, {g, gu, k, a, b, tau}];
Map [AssertSymmetricIncreasing, {g[la,lb], k[la,lb], gu[ua,ub]}];

(* Use the ADMBase variable names *)
g11=gxx; g12=gxy; g22=gyy; g13=gxz; g23=gyz; g33=gzz;
k11=kxx; k12=kxy; k22=kyy; k13=kxz; k23=kyz; k33=kzz;
a=alp;
b1=betax; b2=betay; b3=betaz;

(******************************************************************************)
(******************************************************************************)
(******************************************************************************)
(* Groups *)
(******************************************************************************)
(******************************************************************************)
(******************************************************************************)

ThornGroups = {SetGroupName [CreateGroupFromTensor [tau], "propertime"]};

extraGroups =
  {{"ADMBase::metric",  {gxx, gxy, gxz, gyy, gyz, gzz}},
   {"ADMBase::curv",    {kxx, kxy, kxz, kyy, kyz, kzz}},
   {"ADMBase::lapse",   {alp}},
   {"ADMBase::shift",   {betax, betay, betaz}}
};

declaredGroupNames = Map [First, ThornGroups];
groups = Join [ThornGroups, extraGroups];
       
(******************************************************************************)
(******************************************************************************)
(******************************************************************************)
(* Initial data *)
(******************************************************************************)
(******************************************************************************)
(******************************************************************************)

initialise =
{
  Name -> prefix <> "Initial",
  Schedule -> {"IN CCTK_INITIAL"},
  Where -> Everywhere,
  Equations ->
  {
    tau -> t
  }
};

(******************************************************************************)
(******************************************************************************)
(******************************************************************************)
(* Evolution equations *)
(******************************************************************************)
(******************************************************************************)
(******************************************************************************)

evolCalc =
{
  Name -> prefix <> "RHS",
  Schedule -> {"IN MoL_CalcRHS"},
  Where -> Everywhere,
  Equations ->
  {
    dot[tau] -> (a^2 - b[ua] b[ub] g[la,lb])^(1/2)
  }
};
       
initRHSCalc =
{
  Name -> prefix <> "InitRHS",
  Schedule -> {"AT analysis BEFORE " <> prefix <> "RHS"},
  Where -> Everywhere,
  Equations -> 
  {
    dot[tau] -> 0
  }
};
       
(******************************************************************************)
(******************************************************************************)
(******************************************************************************)
(* Parameters *)
(******************************************************************************)
(******************************************************************************)
(******************************************************************************)

extendedKeywordParameters = {
    {
        Name -> "ADMBase::lapse_evolution_method",
        AllowedValues -> {thorn}
    },
    {
        Name -> "ADMBase::dtlapse_evolution_method",
        AllowedValues -> {thorn}
    }};

keywordParameters = {
    {
        Name -> "CosmoLapse_Kexpression",
        Description -> "d/dt alpha = - F alpha^N (K - Kexpression)",
        AllowedValues -> {{Value -> "zero", Description -> ""},
                          {Value -> "initial_background_K", Description -> ""},
                          {Value -> "background_K", Description -> ""}},
        Default -> "zero"
    }};

realParameters = {
    {
        Name -> alphaF,
        Description -> "d/dt alpha = - F alpha^N (K - Kexpression)",
        Default -> 0.0
    }};
         
intParameters = {
    {
        Name -> alphaN,
        Description -> "d/dt alpha = - F alpha^N (K - Kexpression)",
        Default -> 1
    },
    {
        Name -> FullLieDeriv,
        Description -> "Full Lie derivative, include beta^k partial_k term",
        AllowedValues -> {{Value -> 1, Description -> "yes on"},
                          {Value -> 0, Description -> "no off"}},
        Default -> 0
    }};

(******************************************************************************)
(******************************************************************************)
(******************************************************************************)
(* Construct the thorn *)
(******************************************************************************)
(******************************************************************************)
(******************************************************************************)

inheritedImplementations = {"ADMBase", "ML_BSSN", "Cactus", "ICPertFLRW"};

calculations =
{
  initialise,
  evolCalc,
  initRHSCalc
};

CreateKrancThornTT [groups, ".", ThornName,
  Calculations -> calculations,
  DeclaredGroups -> declaredGroupNames,
  EvolutionTimelevels -> evolutionTimelevels,
  DefaultEvolutionTimelevels -> 3,
  UseJacobian -> useJacobian,
  UseLoopControl -> True,
  UseVectors -> False,
  InheritedImplementations -> inheritedImplementations,
  ExtendedKeywordParameters -> extendedKeywordParameters,
  KeywordParameters -> keywordParameters,
  RealParameters -> realParameters,
  IntParameters -> intParameters
];

];

createCode[4, True, 3];
