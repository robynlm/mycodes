Get["KrancThorn`"];

SetEnhancedTimes[False];

TName = "IC_pFLRW"

(***********************************************************************)
(* Groups *)
(***********************************************************************)

groups =
  {{"Grid::coordinates", {x, y, z, r}},
   {"ADMBase::metric",  {gxx, gxy, gxz, gyy, gyz, gzz}},
   {"ADMBase::curv",  {kxx, kxy, kxz, kyy, kyz, kzz}},
   {"HydroBase::rho",  {rho}},
   {"HydroBase::press",  {press}},
   {"HydroBase::eps",  {eps}},
   {"HydroBase::vel",  {vel}}
};

(***********************************************************************)
(* Initial data *)
(***********************************************************************)

ai2 = "ai"^2;
Hprop = 2 / ( 3 "ti" );
Hconf = "ai" Hprop;
f     = "Omega_mi"^(6/11);
F     = f + 3 "Omega_mi" / 2;
iFH   = 1 / ( F Hconf^2 );
iFHp  = f / ( F Hconf );

kx := "nbrWx" 2 Pi /"L";
ky := "nbrWy" 2 Pi /"L";
kz := "nbrWz" 2 Pi /"L";

Rc[x_, y_, z_] = "Ax" Sin[x kx] + "Ay" Sin[y ky] + "Az" Sin[z kz];
ddRc[x_, y_, z_] = -( "Ax" kx^2 Sin[x kx] + "Ay" ky^2 Sin[y ky] + "Az" kz^2 Sin[y kz] );


psi1[x_, y_, z_] := iFH ddRc[x, y, z] / 3 + Rc[x, y, z];
psi1p[x_, y_, z_] := iFHp ddRc[x, y, z] / 3;

chiRterm11[x_, y_, z_] := "Ay" ky^2 Sin[y ky] + "Az" kz^2 Sin[y kz];
chiRterm22[x_, y_, z_] := "Ax" kx^2 Sin[x kx] + "Az" kz^2 Sin[y kz];
chiRterm33[x_, y_, z_] := "Ax" kx^2 Sin[x kx] + "Ay" ky^2 Sin[y ky];


kappa = 8 Pi "G";
rhoflrw = 3 Hprop^2 "Omega_mi" / kappa;
rhothetaterm[x_, y_, z_] = (-psi1p[x, y, z] - iFHp chiRterm11[x, y, z])^2 + (-psi1p[x, y, z] - iFHp chiRterm22[x, y, z])^2 + (-psi1p[x, y, z] - iFHp chiRterm33[x, y, z])^2;

delta[x_, y_, z_] := ( 9 psi1p[x, y, z]^2 - rhothetaterm[x, y, z] - 12 Hconf psi1p[x, y, z] + 4 ddRc[x, y, z]) / (6 Hconf^2 "Omega_mi");

initialCalc = 
{
  Name -> "ICCalc",
  Schedule -> {"IN HydroBase_Initial"},
  Equations -> 
  {
    gxx -> ai2 ( 1 - 2 psi1[x, y, z] - 2 iFH chiRterm11[x, y, z] ),
    gxy -> 0,
    gxz -> 0,
    gyy -> ai2 ( 1 - 2 psi1[x, y, z] - 2 iFH chiRterm22[x, y, z] ),
    gyz -> 0,
    gzz -> ai2 ( 1 - 2 psi1[x, y, z] - 2 iFH chiRterm33[x, y, z] ),
    
    kxx -> - Hconf - psi1p[x, y, z] - iFHp chiRterm11[x, y, z],
    kxy -> 0,
    kxz -> 0,
    kyy -> - Hconf - psi1p[x, y, z] - iFHp chiRterm22[x, y, z],
    kyz -> 0,
    kzz -> - Hconf - psi1p[x, y, z] - iFHp chiRterm33[x, y, z],

    rho -> rhoflrw ( 1 + delta[x, y, z] ),
    press -> 0,
    eps -> 0,
    vel -> 0
  }
};

(******************************************************************************)
(* Parameters *)
(******************************************************************************)

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
    Default -> 6
  },
  {
    Name -> "L",
    Description -> "Box size",
    Default -> 6
  },
  {
    Name -> "G",
    Description -> "Gravitational Constant",
    Default -> 1
  },
  {
    Name -> "Ax",
    Description -> "Amplitude of x-mode",
    Default -> 0
  },
  {
    Name -> "Ay",
    Description -> "Amplitude of y-mode",
    Default -> 0
  },
  {
    Name -> "Az",
    Description -> "Amplitude of z-mode",
    Default -> 0
  },
  {
    Name -> "nbrWx",
    Description -> "Number of waves along x",
    Default -> 1
  },
  {
    Name -> "nbrWy",
    Description -> "Number of waves along y",
    Default -> 1
  },
  {
    Name -> "nbrWz",
    Description -> "Number of waves along z",
    Default -> 1
  }
};

extendedKeywordParameters =
{
  { Name -> "ADMBase::initial_data", AllowedValues -> {TName} }};(*,
  { Name -> "HydroBase::initial_hydro", AllowedValues -> {TName} }
};*)

(******************************************************************************)
(* Construct the thorns *)
(******************************************************************************)

CreateKrancThornTT [groups, ".", TName,
  Calculations -> {initialCalc},
  DefaultEvolutionTimelevels -> 3,
  UseLoopControl -> True,
  UseVectors -> False,
  InheritedImplementations -> {"ADMBase", "HydroBase"},
  ExtendedKeywordParameters -> extendedKeywordParameters,
  RealParameters -> realparam
];

