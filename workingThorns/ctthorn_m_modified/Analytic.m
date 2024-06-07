(**************************************************************************************)
(*										      *)
(* Copyright 2013-2016 Eloisa Bentivegna					      *)
(*										      *)
(* Analytic.m is a simple Kranc script used to create grid functions which will be    *)
(* used as coefficients (or initial guesses, or exact solutions) by CT_MultiLevel.    *)
(* It is distributed under the General Public License, version 3 or higher.           *)
(*										      *)
(* The runmath.sh and copy-if-changed.sh have been adapted from the same-name         *)
(* scripts in McLachlan; please see McLachlan's licensing and copyright notices       *)
(* regarding this material.                                                           *)
(*										      *)
(**************************************************************************************)

Get["KrancThorn`"];

SetEnhancedTimes[False];
SetSourceLanguage["C"];

(******************************************************************************)
(* Options *)
(******************************************************************************)

createCode[derivOrder_] :=
Module[{},

prefix = "CT_";
suffix =
  ""
  <> If [derivOrder!=4, "_O" <> ToString[derivOrder], ""]
  ;

ThornName = prefix <> "Analytic" <> suffix;

(******************************************************************************)
(* Derivatives *)
(******************************************************************************)

KD = KroneckerDelta;

derivatives =
{
  PDstandardNth[i_]    -> StandardCenteredDifferenceOperator[1,derivOrder/2,i],
  PDstandardNth[i_,i_] -> StandardCenteredDifferenceOperator[2,derivOrder/2,i],
  PDstandardNth[i_,j_] -> StandardCenteredDifferenceOperator[1,derivOrder/2,i] *
                          StandardCenteredDifferenceOperator[1,derivOrder/2,j]
};

FD  = PDstandardNth;


(******************************************************************************)
(* Tensors *)
(******************************************************************************)

(* Register the tensor quantities with the TensorTools package *)
Map [DefineTensor,
     {testinipsi, testinixx, testinixy, testinixz, epsi, elaplacian,
      testcxx, testcxy, testcxz, testcyy, testcyz, testczz,
      testcx, testcy, testcz,
      testc0, testc1, testc2, testc3, testc4, 
      testa0, testa1, testa2, testa3, testa4, 
      testW, testK, testXx, testXy, testXz, testZ, testdxK, testdyK, testdzK,
      testAxx, testAxy, testAxz, testAyy, testAyz, testAzz,
      x, y, z, r}];

(* Use the CartGrid3D variable names *)
x1=x; x2=y; x3=z;

(******************************************************************************)
(* Expressions *)
(******************************************************************************)

pi = N[Pi,40]; 

(******************************************************************************)
(* Groups *)
(******************************************************************************)

evolvedGroups = {};
evaluatedGroups =
  {SetGroupName [CreateGroupFromTensor [testinipsi  ], prefix <> "testinipsi"],
   SetGroupName [CreateGroupFromTensor [testinixx   ], prefix <> "testinixx"],
   SetGroupName [CreateGroupFromTensor [testinixy   ], prefix <> "testinixy"],
   SetGroupName [CreateGroupFromTensor [testinixz   ], prefix <> "testinixz"],
   SetGroupName [CreateGroupFromTensor [epsi     ], prefix <> "epsi"   ],
   SetGroupName [CreateGroupFromTensor [elaplacian ], prefix <> "elaplacian"   ],
   SetGroupName [CreateGroupFromTensor [testcxx  ], prefix <> "testcxx"],
   SetGroupName [CreateGroupFromTensor [testcxy  ], prefix <> "testcxy"],
   SetGroupName [CreateGroupFromTensor [testcxz  ], prefix <> "testcxz"],
   SetGroupName [CreateGroupFromTensor [testcyy  ], prefix <> "testcyy"],
   SetGroupName [CreateGroupFromTensor [testcyz  ], prefix <> "testcyz"],
   SetGroupName [CreateGroupFromTensor [testczz  ], prefix <> "testczz"],
   SetGroupName [CreateGroupFromTensor [testcx   ], prefix <> "testcx"],
   SetGroupName [CreateGroupFromTensor [testcy   ], prefix <> "testcy"],
   SetGroupName [CreateGroupFromTensor [testcz   ], prefix <> "testcz"],
   SetGroupName [CreateGroupFromTensor [testc0   ], prefix <> "testc0"],
   SetGroupName [CreateGroupFromTensor [testc1   ], prefix <> "testc1"],
   SetGroupName [CreateGroupFromTensor [testc2   ], prefix <> "testc2"],
   SetGroupName [CreateGroupFromTensor [testc3   ], prefix <> "testc3"],
   SetGroupName [CreateGroupFromTensor [testc4   ], prefix <> "testc4"],
   SetGroupName [CreateGroupFromTensor [testa0   ], prefix <> "testa0"],
   SetGroupName [CreateGroupFromTensor [testa1   ], prefix <> "testa1"],
   SetGroupName [CreateGroupFromTensor [testa2   ], prefix <> "testa2"],
   SetGroupName [CreateGroupFromTensor [testa3   ], prefix <> "testa3"],
   SetGroupName [CreateGroupFromTensor [testa4   ], prefix <> "testa4"],
   SetGroupName [CreateGroupFromTensor [testW    ], prefix <> "testW"],
   SetGroupName [CreateGroupFromTensor [testK    ], prefix <> "testK"],
   SetGroupName [CreateGroupFromTensor [testdxK  ], prefix <> "testdxK"],
   SetGroupName [CreateGroupFromTensor [testdyK  ], prefix <> "testdyK"],
   SetGroupName [CreateGroupFromTensor [testdzK  ], prefix <> "testdzK"],
   SetGroupName [CreateGroupFromTensor [testXx   ], prefix <> "testXx"],
   SetGroupName [CreateGroupFromTensor [testXy   ], prefix <> "testXy"],
   SetGroupName [CreateGroupFromTensor [testXz   ], prefix <> "testXz"],
   SetGroupName [CreateGroupFromTensor [testZ    ], prefix <> "testZ"],
   SetGroupName [CreateGroupFromTensor [testAxx  ], prefix <> "testAxx"],
   SetGroupName [CreateGroupFromTensor [testAxy  ], prefix <> "testAxy"],
   SetGroupName [CreateGroupFromTensor [testAxz  ], prefix <> "testAxz"],
   SetGroupName [CreateGroupFromTensor [testAyy  ], prefix <> "testAyy"],
   SetGroupName [CreateGroupFromTensor [testAyz  ], prefix <> "testAyz"],
   SetGroupName [CreateGroupFromTensor [testAzz  ], prefix <> "testAzz"]
  };

declaredGroups = Join [evolvedGroups, evaluatedGroups];

declaredGroupNames = Map [First, declaredGroups];



extraGroups =
  {{"Grid::coordinates", {x, y, z, r}}
};



groups = Join [declaredGroups, extraGroups];

(******************************************************************************)
(* Constraint equations *)
(******************************************************************************)

HTheta[r_] := (Sign[r] + 1) / 2;
rhoBall[r_] := HTheta[rBall-r];
psiIn[r_] := rhoBall[0] (r^2 - 3 rBall^2) / 6;
psiOut[r_] := - Sign[x-rBall] rhoBall[0] rBall^3 / (3 r);
psiBall[r_] := HTheta[rBall-r] rhoBall[0] (r^2 - 3 rBall^2) / 6 - HTheta[r-rBall] rhoBall[0] rBall^3 / (3 (r+0.1));

gaussian[x_, y_, z_] := Sum[ "(amp["<>ToString[i]<>"])" Exp[-1/2 ((x-"(x0["<>ToString[i]<>"])")/"(sigmax["<>ToString[i]<>"])")^2] 
                                                        Exp[-1/2 ((y-"(y0["<>ToString[i]<>"])")/"(sigmay["<>ToString[i]<>"])")^2] 
                                                        Exp[-1/2 ((z-"(z0["<>ToString[i]<>"])")/"(sigmaz["<>ToString[i]<>"])")^2], {i,0,4}];
sine[x_, y_, z_] := Sin[2 pi (kx x + phasex)] Sin[2 pi (ky y + phasey)] Sin[2 pi (kz z + phasez)];
sing[x_, y_, z_] := 0.5 ( massa / (Sqrt[(x-xa)^2+(y-ya)^2+(z-za)^2+eps^2]) + massb / (Sqrt[(x-xb)^2+(y-yb)^2+(z-zb)^2+eps^2]));
soln[x_, y_, z_] := ampG gaussian[x,y,z] + ampS sine[x,y,z] + ampC + ampSg sing[x,y,z];
dxxsoln[x_, y_, z_] := D[soln[x,y,z], {x,2}];
dyysoln[x_, y_, z_] := D[soln[x,y,z], {y,2}];
dzzsoln[x_, y_, z_] := D[soln[x,y,z], {z,2}];

rad[x_, y_, z_] := Sqrt[x^2+y^2+z^2+eps^2];
W[r_] := pi/2 + ArcTan[10(r-5)];
Axx[x_, y_, z_] := 0;
Axy[x_, y_, z_] := 0;
Axz[x_, y_, z_] := 0;
Ayy[x_, y_, z_] := 0;
Ayz[x_, y_, z_] := 0;
Azz[x_, y_, z_] := 0;

W[r_] := HTheta[r-l] HTheta[l+sigma-r-eps] ((r-l-sigma)^6/sigma^6 - 1)^6 + HTheta[r-l-sigma+eps];
sfact[r_] := sing[r/Sqrt[3],r/Sqrt[3],r/Sqrt[3]]; 
mW[r_] := sfact[r] W[r];
dW[r_] := (HTheta[rr-l] HTheta[l+sigma-rr] D[((rr-l-sigma)^6/sigma^6 - 1)^6,rr])/.rr->r;
dxW[x_, y_, z_] := x dW[rad[x,y,z]] / rad[x,y,z];
dyW[x_, y_, z_] := y dW[rad[x,y,z]] / rad[x,y,z];
dzW[x_, y_, z_] := z dW[rad[x,y,z]] / rad[x,y,z];
d2W[r_] := HTheta[r-l] HTheta[l+sigma-r] (D[rr^2 D[sfact[rrr] ((rrr-l-sigma)^6/sigma^6 - 1)^6,rrr]/.rrr->rr,rr]/(rr^2+eps^2))/.rr->r;
sLmW[r_] := d2W[r];
lapW[r_] := HTheta[r-l] HTheta[l+sigma-r] (D[rr^2 D[((rrr-l-sigma)^6/sigma^6 - 1)^6,rrr]/.rrr->rr,rr]/(rr^2+eps^2))/.rr->r;

AnalyticPoissonCalc =
{
  Name -> ThornName <> "_Poisson_Calc",
  Schedule -> {"AT CCTK_INITIAL before CT_MultiLevel"},
  ConditionalOnKeyword -> {"free_data", "Poisson"},
  Where -> Everywhere,
  Equations -> 
  {
    testinipsi -> ampC,

    epsi -> psiBall[r],

    testcxx -> 1,
    testcxy -> 0,
    testcxz -> 0,
    testcyy -> 1,
    testcyz -> 0,
    testczz -> 1,

    testcx -> 0,
    testcy -> 0,
    testcz -> 0,

    testc1 -> 0,
    testc2 -> 0,
    testc3 -> 0,
    testc4 -> 0,

    testc0 -> - rhoBall[r],

    testW -> 0,
    testK -> 0,
    testXx -> 0,
    testXy -> 0,
    testXz -> 0,
    testZ  -> 0
  }
};

AnalyticExactCalc =
{
  Name -> ThornName <> "_Exact_Calc",
  Schedule -> {"AT CCTK_INITIAL before CT_MultiLevel"},
  ConditionalOnKeyword -> {"free_data", "exact"},
  Where -> Everywhere,
  NoSimplify -> True,
  Shorthands -> {d2soln},
  Equations -> 
  {
    d2soln -> dxxsoln[x,y,z] + dyysoln[x,y,z] + dzzsoln[x,y,z],

    epsi -> soln[x,y,z],
    elaplacian -> d2soln,

    testinipsi -> ampI soln[x,y,z],

    testcxx -> 1,
    testcxy -> 0,
    testcxz -> 0,
    testcyy -> 1,
    testcyz -> 0,
    testczz -> 1,

    testcx -> 0,
    testcy -> 0,
    testcz -> 0,

    testc1 -> ampC1,
    testc2 -> 0,
    testc3 -> 0,
    testc4 -> 0,

    testc0 -> -d2soln - testc1 soln[x,y,z] - testc2 (soln[x,y,z])^2,
    testc1 -> ampC1
  }
};

AnalyticExpandingLatticeCalc =
{
  Name -> ThornName <> "_ExpandingLattice_Calc",
  Schedule -> {"AT CCTK_INITIAL before CT_MultiLevel"},
  ConditionalOnKeyword -> {"free_data", "Expanding BH lattice"},
  Where -> Everywhere,
  Shorthands -> {d2mW},
  NoSimplify -> True,
  Equations -> 
  {
    d2mW -> sLmW[r],
    testinipsi -> ampI,
    epsi -> 1,

    testcxx -> 1,
    testcyy -> 1,
    testczz -> 1,

    testc0 -> -d2mW,
    testc1 -> - (Kc W[r])^2 / 12,
    testc2 -> Kc,
    testc3 -> lapW[r],

    testa0 -> sfact[r] (1-W[r]),

    testW -> W[r],
    testK -> Kc W[r],
    testdxK -> Kc dxW[x,y,z],
    testdyK -> Kc dyW[x,y,z],
    testdzK -> Kc dzW[x,y,z],
    testXx -> 0,
    testXy -> 0,
    testXz -> 0
  }
};

AnalyticBYCalc =
{
  Name -> ThornName <> "_BY_Calc",
  Schedule -> {"AT CCTK_INITIAL before CT_MultiLevel"},
  ConditionalOnKeyword -> {"free_data", "Bowen-York"},
  Where -> Everywhere,
  Shorthands -> {d2soln, ra, nax, nay, naz, rb, nbx, nby, nbz, nPa, nPb},
  Equations -> 
  {
    d2soln -> dxxsoln[x,y,z] + dyysoln[x,y,z] + dzzsoln[x,y,z],

    ra -> Sqrt[(x-xa)^2+(y-ya)^2+(z-za)^2+eps^2],
    nax -> (x-xa) / ra,
    nay -> (y-ya) / ra,
    naz -> (z-za) / ra,

    rb -> Sqrt[(x-xb)^2+(y-yb)^2+(z-zb)^2+eps^2],
    nbx -> (x-xb) / rb,
    nby -> (y-yb) / rb,
    nbz -> (z-zb) / rb,

    nPa -> nax Pax + nay Pay + naz Paz,
    nPb -> nbx Pbx + nby Pby + nbz Pbz,

    epsi -> 1,

    testinipsi -> ampI, 

    testcxx -> 1,
    testcyy -> 1,
    testczz -> 1,

    testa0 -> sing[x,y,z],
    testa1 -> gaussian[x,y,z],

    testXx -> -0.25 ( 7 Pax + nax nPa ) / ra - 0.25 ( 7 Pbx + nbx nPb ) / rb,
    testXy -> -0.25 ( 7 Pay + nay nPa ) / ra - 0.25 ( 7 Pby + nby nPb ) / rb,
    testXz -> -0.25 ( 7 Paz + naz nPa ) / ra - 0.25 ( 7 Pbz + nbz nPb ) / rb,
    testAxx -> 1.5 ( 2 Pax nax - (1 - nax nax) nPa ) / ra^2 + 1.5 ( 2 Pbx nbx - (1 - nbx nbx) nPb ) / rb^2,
    testAxy -> 1.5 ( Pax nay + Pay nax + nax nay nPa ) / ra^2 + 1.5 ( Pbx nby + Pby nbx + nbx nby nPb ) / rb^2, 
    testAxz -> 1.5 ( Pax naz + Paz nax + nax naz nPa ) / ra^2 + 1.5 ( Pbx nbz + Pbz nbx + nbx nbz nPb ) / rb^2, 
    testAyy -> 1.5 ( 2 Pay nay - (1- nay nay) nPa ) / ra^2 + 1.5 ( 2 Pby nby - (1- nby nby) nPb ) / rb^2,
    testAyz -> 1.5 ( Pay naz + Paz nay + nay naz nPa ) / ra^2 + 1.5 ( Pby nbz + Pbz nby + nby nbz nPb ) / rb^2, 
    testAzz -> 1.5 ( 2 Paz naz - (1- naz naz) nPa ) / ra^2 + 1.5 ( 2 Pbz nbz - (1- nbz nbz) nPb ) / rb^2,
    testZ  -> 0,
    testK  -> 0,
    testdxK -> 0,
    testdyK -> 0,
    testdzK -> 0
  }
};

AnalyticExactBoundary =
{
  Name -> ThornName <> "_ExactBoundary",
  Schedule -> {"AT CCTK_INITIAL before CT_MultiLevel after CT_Analytic_Calc"},
  ConditionalOnKeyword -> {"free_data", "exact"},
  Where -> Boundary,
  NoSimplify -> True,
  Equations -> 
  {
    testinipsi -> soln[x,y,z]
  }
};

(* psi *)
soln[x_, y_, z_] := ampC + ampS Sin[2 pi (kx x + phasex)] Sin[2 pi (ky y + phasey)] Sin[2 pi (kz z + phasez)] + ampG gaussian[x,y,z];
  (* Second derivatives *)
dxxsoln[x_, y_, z_] := D[soln[x,y,z], {x,2}];
dyysoln[x_, y_, z_] := D[soln[x,y,z], {y,2}];
dzzsoln[x_, y_, z_] := D[soln[x,y,z], {z,2}];

(* X^i *)
solnXx[x_, y_, z_] := ampV Sin[2 pi kx x] Sin[2 pi ky y] Sin[2 pi kz z] + ampVG (Exp[-((x-vecA)^2+y^2+z^2)/(2 sigma^2)]-Exp[-((x+vecA)^2+y^2+z^2)/(2 sigma^2)]);
solnXy[x_, y_, z_] := ampV Sin[2 pi kx x] Sin[2 pi ky y] Sin[2 pi kz z] + ampVG (Exp[-(x^2+(y-vecA)^2+z^2)/(2 sigma^2)]-Exp[-(x^2+(y+vecA)^2+z^2)/(2 sigma^2)]);
solnXz[x_, y_, z_] := ampV Sin[2 pi kx x] Sin[2 pi ky y] Sin[2 pi kz z] + ampVG (Exp[-(x^2+y^2+(z-vecA)^2)/(2 sigma^2)]-Exp[-(x^2+y^2+(z+vecA)^2)/(2 sigma^2)]);
  (* First derivatives *)
dxsolnXx[x_, y_, z_] := D[solnXx[x,y,z], x];
dysolnXx[x_, y_, z_] := D[solnXx[x,y,z], y];
dzsolnXx[x_, y_, z_] := D[solnXx[x,y,z], z];
dxsolnXy[x_, y_, z_] := D[solnXy[x,y,z], x];
dysolnXy[x_, y_, z_] := D[solnXy[x,y,z], y];
dzsolnXy[x_, y_, z_] := D[solnXy[x,y,z], z];
dxsolnXz[x_, y_, z_] := D[solnXz[x,y,z], x];
dysolnXz[x_, y_, z_] := D[solnXz[x,y,z], y];
dzsolnXz[x_, y_, z_] := D[solnXz[x,y,z], z];
  (* Second derivatives *)
dxxsolnXx[x_, y_, z_] := D[solnXx[x,y,z], {x,2}];
dyysolnXx[x_, y_, z_] := D[solnXx[x,y,z], {y,2}];
dzzsolnXx[x_, y_, z_] := D[solnXx[x,y,z], {z,2}];
dxxsolnXy[x_, y_, z_] := D[solnXy[x,y,z], {x,2}];
dyysolnXy[x_, y_, z_] := D[solnXy[x,y,z], {y,2}];
dzzsolnXy[x_, y_, z_] := D[solnXy[x,y,z], {z,2}];
dxxsolnXz[x_, y_, z_] := D[solnXz[x,y,z], {x,2}];
dyysolnXz[x_, y_, z_] := D[solnXz[x,y,z], {y,2}];
dzzsolnXz[x_, y_, z_] := D[solnXz[x,y,z], {z,2}];
  (* Divergence and its derivatives *)
divsolnX[x_, y_, z_] := D[solnXx[x,y,z], x] + D[solnXy[x,y,z], y] + D[solnXz[x,y,z], z];
dxdivsolnX[x_, y_, z_] := D[divsolnX[x,y,z], x];
dydivsolnX[x_, y_, z_] := D[divsolnX[x,y,z], y];
dzdivsolnX[x_, y_, z_] := D[divsolnX[x,y,z], z];
  (* Square of Aij *)
Axx[x_, y_, z_] := 2 dxsolnXx[x, y, z] - 2/3 divsolnX[x, y, z];
Ayy[x_, y_, z_] := 2 dysolnXy[x, y, z] - 2/3 divsolnX[x, y, z];
Azz[x_, y_, z_] := 2 dzsolnXz[x, y, z] - 2/3 divsolnX[x, y, z];
Axy[x_, y_, z_] := dxsolnXy[x, y, z] + dysolnXx[x, y, z];
Axz[x_, y_, z_] := dxsolnXz[x, y, z] + dzsolnXx[x, y, z];
Ayz[x_, y_, z_] := dzsolnXy[x, y, z] + dysolnXz[x, y, z];
L2[x_, y_, z_] := Axx[x, y, z] Axx[x, y, z] + Ayy[x, y, z] Ayy[x, y, z] + Azz[x, y, z] Azz[x, y, z] + 2 Axy[x, y, z] Axy[x, y, z] + 2 Axz[x, y, z] Axz[x, y, z] + 2 Ayz[x, y, z] Ayz[x, y, z]; 
  (* Trace of Kij and its derivatives *)
(*Kfunc[x_, y_, z_] := Kc gaussian[x,y,z];*)
Kfunc[x_, y_, z_] := Kc;
Kfunce[x_, y_, z_] := Ke;
dxKfunce[x_, y_, z_] := D[Kfunce[x,y,z],x];
dyKfunce[x_, y_, z_] := D[Kfunce[x,y,z],y];
dzKfunce[x_, y_, z_] := D[Kfunce[x,y,z],z];
AnalyticLumpCalc =
{
  Name -> ThornName <> "_Lump_Calc",
  Schedule -> {"AT CCTK_INITIAL before CT_MultiLevel"},
  ConditionalOnKeyword -> {"free_data", "Lump"},
  Where -> Everywhere,
  Shorthands -> {d2soln, d2solnXx, d2solnXy, d2solnXz},
  NoSimplify -> True,
  Equations -> 
  {
    d2soln -> dxxsoln[x,y,z] + dyysoln[x,y,z] + dzzsoln[x,y,z],
    d2solnXx -> dxxsolnXx[x,y,z] + dyysolnXx[x,y,z] + dzzsolnXx[x,y,z],
    d2solnXy -> dxxsolnXy[x,y,z] + dyysolnXy[x,y,z] + dzzsolnXy[x,y,z],
    d2solnXz -> dxxsolnXz[x,y,z] + dyysolnXz[x,y,z] + dzzsolnXz[x,y,z],

    epsi -> soln[x,y,z],
    elaplacian -> d2soln,

    testinipsi -> ampI (*soln[x,y,z]*),
    testinixx  -> 0,
    testinixy  -> 0,
    testinixz  -> 0,

    testa0 -> 0,

    testcxx -> 1,
    testcxy -> 0,
    testcxz -> 0,
    testcyy -> 1,
    testcyz -> 0,
    testczz -> 1,

    testcx -> 0,
    testcy -> 0,
    testcz -> 0,

    testc0 -> - Kfunc[x,y,z]^2 / 12,
    testc1 ->   Kfunce[x,y,z]^2 / 12 (soln[x,y,z])^5 - L2[x,y,z] / 8 (soln[x,y,z])^(-7) - d2soln,
    testc2 -> - d2solnXx - dxdivsolnX[x,y,z] / 3 + 2 (soln[x,y,z])^6 dxKfunce[x,y,z] / 3,
    testc3 -> - d2solnXy - dydivsolnX[x,y,z] / 3 + 2 (soln[x,y,z])^6 dyKfunce[x,y,z] / 3,
    testc4 -> - d2solnXz - dzdivsolnX[x,y,z] / 3 + 2 (soln[x,y,z])^6 dzKfunce[x,y,z] / 3,

    testXx -> solnXx[x,y,z], 
    testXy -> solnXy[x,y,z], 
    testXz -> solnXz[x,y,z],
    testK  -> Kfunc[x,y,z],
    testdxK -> dxKfunce[x,y,z],
    testdyK -> dyKfunce[x,y,z],
    testdzK -> dzKfunce[x,y,z]
  }
};

AnalyticLumpBoundary =
{
  Name -> ThornName <> "_LumpBoundary",
  Schedule -> {"AT CCTK_INITIAL before CT_MultiLevel after CT_Analytic_Calc"},
  ConditionalOnKeyword -> {"free_data", "Lump"},
  Where -> Boundary,
  NoSimplify -> True,
  Equations -> 
  {
    testinipsi -> soln[x,y,z],
    testinixx  -> solnXx[x,y,z],
    testinixy  -> solnXy[x,y,z],
    testinixz  -> solnXz[x,y,z]
  }
};

(******************************************************************************)
(* Implementations *)
(******************************************************************************)

(******************************************************************************)
(* Parameters *)
(******************************************************************************)

intParameters =
{
  {
    Name -> imaxF,
    Description -> "Max number of Fourier modes to include in x direction",
    Default -> 1
  },
  {
    Name -> jmaxF,
    Description -> "Max number of Fourier modes to include in y direction",
    Default -> 1
  },
  {
    Name -> kmaxF,
    Description -> "Max number of Fourier modes to include in z direction",
    Default -> 1
  }
};

realParameters =
{
  {
    Name -> kx,
    Description -> "Wavelength parameter along x",
    Default -> 0
  },
  {
    Name -> ky,
    Description -> "Wavelength parameter along y",
    Default -> 0
  },
  {
    Name -> kz,
    Description -> "Wavelength parameter along z",
    Default -> 0
  },
  {
    Name -> ampG,
    Description -> "Coefficient of the gaussian term in the exact solution",
    Default -> 0
  },
  {
    Name -> ampS,
    Description -> "Coefficient of the sine term in the exact solution",
    Default -> 0
  },
  {
    Name -> ampC,
    Description -> "Constant coefficient in the exact solution",
    Default -> 0
  },
  {
    Name -> ampI,
    Description -> "Multiplication factor between initial guess and exact solution",
    Default -> 0
  },
  {
    Name -> ampC1,
    Description -> "Initial value for testc1",
    Default -> 0
  },
  {
    Name -> ampSg,
    Description -> "Coefficient of the 1/r term in the exact solution",
    Default -> 0
  },
  {
    Name -> ampV,
    Description -> "Coefficient of the vector part in the exact solution",
    Default -> 0
  },
  {
    Name -> ampVG,
    Description -> "Coefficient of the vector part in the exact solution (gaussian term)",
    Default -> 0
  },
  {
    Name -> sigma,
    Description -> "Width of transition function in extrinsic curvature",
    Default -> 1
  },
  {
    Name -> l,
    Description -> "Location of transition function in extrinsic curvature",
    Default -> 0
  },
  {
    Name -> phasex,
    Description -> "Phase in the initial data for psi along x",
    Default -> 0
  },
  {
    Name -> phasey,
    Description -> "Phase in the initial data for psi along y",
    Default -> 0
  },
  {
    Name -> phasez,
    Description -> "Phase in the initial data for psi along z",
    Default -> 0
  },
  {
    Name -> Kc,
    Description -> "Coefficient of extrinsic curvature",
    Default -> 0
  },
  {
    Name -> Ke,
    Description -> "Coefficient of extrinsic curvature in exact solution",
    Default -> 0
  },
  {
    Name -> massa,
    Description -> "mass of first black hole",
    Default -> 0
  },
  {
    Name -> massb,
    Description -> "mass of second black hole",
    Default -> 0
  },
  {
    Name -> xa,
    Description -> "x-coordinate of first black hole for BY initial data",
    Default -> 0
  },
  {
    Name -> ya,
    Description -> "y-coordinate of first black hole for BY initial data",
    Default -> 0
  },
  {
    Name -> za,
    Description -> "z-coordinate of first black hole for BY initial data",
    Default -> 0
  },
  {
    Name -> xb,
    Description -> "x-coordinate of second black hole for BY initial data",
    Default -> 0
  },
  {
    Name -> yb,
    Description -> "y-coordinate of second black hole for BY initial data",
    Default -> 0
  },
  {
    Name -> zb,
    Description -> "z-coordinate of second black hole for BY initial data",
    Default -> 0
  },
  {
    Name -> Pax,
    Description -> "x-component of linear momentum of first black hole for BY initial data",
    Default -> 0
  },
  {
    Name -> Pay,
    Description -> "y-component of linear momentum of first black hole for BY initial data",
    Default -> 0
  },
  {
    Name -> Paz,
    Description -> "z-component of linear momentum of first black hole for BY initial data",
    Default -> 0
  },
  {
    Name -> Pbx,
    Description -> "x-component of linear momentum of second black hole for BY initial data",
    Default -> 0
  },
  {
    Name -> Pby,
    Description -> "y-component of linear momentum of second black hole for BY initial data",
    Default -> 0
  },
  {
    Name -> Pbz,
    Description -> "z-component of linear momentum of second black hole for BY initial data",
    Default -> 0
  },
  {
    Name -> eps,
    Description -> "Smoothing factor",
    Default -> 1`*^-6
  },
  {
    Name -> edgeL,
    Description -> "Coordinate length of cell edge",
    Default -> 10
  },
  {
    Name -> rBall,
    Description -> "Coordinate radius of ball of density for Poisson's equation",
    Default -> 1
  },
  {
    Name -> vecA,
    Description -> "Coordinate center of gaussian representing the X^i vector in the CTT decomposition of the constraints",
    Default -> 1
  },
  {
    Name -> "amp[5]",
    Description -> "Initial amplitude of peaks",
    Default -> 0
  },
  {
    Name -> "x0[5]",
    Description -> "Initial x-locations of peaks",
    Default -> 0
  },
  {
    Name -> "y0[5]",
    Description -> "Initial y-locations of peaks",
    Default -> 0
  },
  {
    Name -> "z0[5]",
    Description -> "Initial z-locations of peaks",
    Default -> 0
  },
  {
    Name -> "sigmax[5]",
    Description -> "x-spreads of initial gaussians",
    Default -> 1
  },
  {
    Name -> "sigmay[5]",
    Description -> "y-spreads of initial gaussians",
    Default -> 1
  },
  {
    Name -> "sigmaz[5]",
    Description -> "z-spreads of initial gaussians",
    Default -> 1
  }
};

keywordParameters =
{
  {
    Name -> "free_data",
    Description -> "How to set the free data for the extrinsic curvature?",
    AllowedValues -> {"exact", "Expanding BH lattice", "Bowen-York", "Poisson", "Lump"},
    Default -> "exact"
  }
};

(******************************************************************************)
(* Construct the thorns *)
(******************************************************************************)

calculations =
{
  AnalyticPoissonCalc,
  AnalyticExactCalc,
  AnalyticExpandingLatticeCalc,
  AnalyticBYCalc,
  AnalyticLumpCalc,
  AnalyticExactBoundary,
  AnalyticLumpBoundary
};

CreateKrancThornTT [groups, ".", ThornName,
  Calculations -> calculations,
  DeclaredGroups -> declaredGroupNames,
  PartialDerivatives -> derivatives,
  UseLoopControl -> True,
  UseVectors -> False,
  IntParameters -> intParameters,
  RealParameters -> realParameters,
  KeywordParameters -> keywordParameters
];

];



(******************************************************************************)
(* Options *)
(******************************************************************************)

(* derivative order: 2, 4, 6, 8, ... *)
(* useGlobalDerivs: False or True *)
(* timelevels: 2 or 3
   (keep this at 3; this is better chosen with a run-time parameter) *)
(* matter: 0 or 1
   (matter seems cheap; it should be always enabled) *)

createCode[4];
