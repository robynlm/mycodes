Get["KrancThorn`"];

TName = "MyIC";

groups =
  {{"Grid::coordinates", {x, y, z, r}},
   {"ADMBase::metric",  {gxx, gxy, gxz, gyy, gyz, gzz}},
   {"ADMBase::lapse",   {alp}},
   {"ADMBase::shift",   {betax, betay, betaz}},
   {"ADMBase::curv",  {kxx, kxy, kxz, kyy, kyz, kzz}},
   {"HydroBase::rho",  {rho}},
   {"HydroBase::press",  {press}},
   {"HydroBase::eps",  {eps}},
   {"HydroBase::vel",  {vel}}
};

ai2 = "ai"^2;
Hprop = 2 / ( 3 "ti" );
kappa = 8 Pi "G";
rhoflrw = 3 Hprop^2 "Omega_mi" / kappa;
kflrw = -ai2 Hprop;

initialFLRW = 
{
  Name -> TName <> "_FLRW",
  Schedule -> {"IN ADMBase_InitialData"},
  ConditionalOnKeyword -> {"my_initial_data", "FLRW"},
  Equations -> 
  {
    gxx -> ai2,
    gxy -> 0,
    gxz -> 0,
    gyy -> ai2,
    gyz -> 0,
    gzz -> ai2,
    
    kxx -> kflrw,
    kxy -> 0,
    kxz -> 0,
    kyy -> kflrw,
    kyz -> 0,
    kzz -> kflrw,

    rho -> rhoflrw,
    press -> 0,
    eps -> 0,
    vel -> 0
  }
};

initialsynch = 
{
  Name -> TName <> "_Synch",
  Schedule -> {"IN ADMBase_InitialData"},
  ConditionalOnKeyword -> {"Alpha_and_Beta", "Synchronous"},
  Equations -> 
  {
    alp -> 1,
    betax -> 0,
    betay -> 0,
    betaz -> 0
  }
};
(* !! Change their evolution with ML *)

realparam = 
{
  {
    Name -> "ti",
    Description -> "Initial Time",
    Default -> 1
  },
  {
    Name -> "Omega_mi",
    Description -> "Initial Matter Density Parameter",
    Default -> 1
  },
  {
    Name -> "ai",
    Description -> "Initial Scale Factor",
    Default -> 1
  },
  {
    Name -> "G",
    Description -> "Gravitational Constant",
    Default -> 1
  }
};

extendedKeywordParameters =
{
  {
    Name -> "ADMBase::initial_data",
    AllowedValues -> {TName}
  },
  {
    Name -> "ADMBase::initial_lapse",
    AllowedValues -> {TName}
  },
  {
    Name -> "ADMBase::initial_shift",
    AllowedValues -> {TName} 
  }};(*,
  {
    Name -> "HydroBase::initial_hydro",
    AllowedValues -> {TName}
  }
};*)

keywordParameters =
{
  {
    Name -> "Alpha_and_Beta",
    AllowedValues -> {"Synchronous"},
    Default -> "Synchronous"
  },
  {
    Name -> "my_initial_data",
    AllowedValues -> {"FLRW"},
    Default -> "FLRW"
  }
};

calculations = {initialFLRW, initialsynch};

CreateKrancThornTT [groups, ".", TName,
  Calculations -> calculations,
  DefaultEvolutionTimelevels -> 3,
  UseLoopControl -> True,
  UseVectors -> False,
  InheritedImplementations -> {"ADMBase", "HydroBase"},
  ExtendedKeywordParameters -> extendedKeywordParameters,
  KeywordParameters -> keywordParameters,
  RealParameters -> realparam
];




