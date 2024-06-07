#! /bin/bash

#rm -rf CT_Dust
#/Users/robynmunoz/EinsteinToolkit/Cactus/repos/Kranc/Bin/kranc Dust.m
#/Users/rlm36AA/codes/ET/ET_2023_05/Cactus/repos/Kranc/Bin/kranc Dust.m

cd CT_Dust
sed -i'.bak' 's/wlorentz/w_lorentz/g' *
sed -i'.bak' 's/shares: ADMBase/shares: HydroBase\n\nEXTENDS CCTK_KEYWORD initial_hydro "initial_hydro"\n{\n  "CT_Dust" :: ""\n}\n\nEXTENDS CCTK_KEYWORD evolution_method "evolution_method"\n{\n  "CT_Dust" :: ""\n}\n\n\nshares: ADMBase/g' param.ccl
sed -i'.bak' 's/HydroBase::velx(Everywhere)/HydroBase::vel(Everywhere)/g' schedule.ccl
sed -i'.bak' 's/WRITES: HydroBase::vely(Everywhere)//g' schedule.ccl
sed -i'.bak' 's/WRITES: HydroBase::velz(Everywhere)//g' schedule.ccl
sed -i'.bak' 's/READS: HydroBase::vely(Everywhere)//g' schedule.ccl
sed -i'.bak' 's/READS: HydroBase::velz(Everywhere)//g' schedule.ccl
rm *.bak

cd src
sed -i'.bak' 's/wlorentz/w_lorentz/g' *
sed -i'.bak' 's/const int imax2=imax\[2\]/const int imax2=imax\[2\];\n  const int N = cctk_ash\[0\]\*cctk_ash\[1\]\*cctk_ash\[2\]/g' *
sed -i'.bak' 's/velx\[index\]/vel\[index\]/g' *
sed -i'.bak' 's/vely\[index\]/vel\[index+N\]/g' *
sed -i'.bak' 's/velz\[index\]/vel\[index+2*N\]/g' *
rm *.bak
