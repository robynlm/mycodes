#!/bin/bash

FILE=CosmoLapse/src/CosmoLapse_AverageK.cc

export OLD="imax\[2\]\;"
export NEW="imax\[2\]\;\\
  CCTK_REAL sum CCTK_ATTRIBUTE_UNUSED \= 0.\;\\
  CCTK_REAL count CCTK_ATTRIBUTE_UNUSED \= 0.\;"
sed -i "" "s/${OLD}/${NEW}/" $FILE

export OLD="CCTK_REAL sum CCTK_ATTRIBUTE_UNUSED \= KtraceL\;"
export NEW="sum \= sum \+ KtraceL\;"
sed -i "" "s/${OLD}/${NEW}/" $FILE

export OLD="CCTK_REAL count CCTK_ATTRIBUTE_UNUSED \= 1.\;"
export NEW="count \= count \+ 1.\;"
sed -i "" "s/${OLD}/${NEW}/" $FILE

export OLD="  CCTK_REAL KtraceaverageL CCTK_ATTRIBUTE_UNUSED \= sum\*pow(count\,\-1)\;"
export NEW="\}\\
  CCTK_ENDLOOP3\(CosmoLapse_AverageK\)\;\\
  \\
  CCTK_REAL KtraceaverageL CCTK_ATTRIBUTE_UNUSED \= sum\*pow\(count\,\-1\)\;\\
  \\
  \#pragma omp parallel\\
  CCTK_LOOP3\(CosmoLapse_AverageK\,\\
    i\,j\,k\,imin0\,imin1\,imin2\,imax0\,imax1\,imax2\,\\
    cctk_ash\[0\]\,cctk_ash\[1\]\,cctk_ash\[2\]\)\\
  \{\\
    const ptrdiff_t index CCTK_ATTRIBUTE_UNUSED \= di\*i \+ dj\*j \+ dk\*k\;"
sed -i "" "s/${OLD}/${NEW}/" $FILE

##########################################################################################
##########################################################################################
##########################################################################################

FILE=CosmoLapse/src/CosmoLapse_InitialTau.cc

export OLD="imax\[2\]\;"
export NEW="imax\[2\]\;\\
  CCTK_REAL sum CCTK_ATTRIBUTE_UNUSED \= 0.\;\\
  CCTK_REAL count CCTK_ATTRIBUTE_UNUSED \= 0.\;"
sed -i "" "s/${OLD}/${NEW}/" $FILE

export OLD="CCTK_REAL sum CCTK_ATTRIBUTE_UNUSED \= KtraceL\;"
export NEW="sum \= sum \+ KiniL\;"
sed -i "" "s/${OLD}/${NEW}/" $FILE

export OLD="CCTK_REAL count CCTK_ATTRIBUTE_UNUSED \= 1.\;"
export NEW="count \= count \+ 1.\;"
sed -i "" "s/${OLD}/${NEW}/" $FILE

export OLD="CCTK_REAL KtraceaverageiniL CCTK_ATTRIBUTE_UNUSED \= sum\*pow(count\,\-1)\;"
export NEW=" "
sed -i "" "s/${OLD}/${NEW}/" $FILE

export OLD="Ktraceaverageini\[index\] \= KtraceaverageiniL\;"
export NEW=" "
sed -i "" "s/${OLD}/${NEW}/" $FILE

export OLD="CCTK_ENDLOOP3(CosmoLapse_InitialTau)\;"
export NEW="CCTK_ENDLOOP3\(CosmoLapse_InitialTau\)\;\\
  \\
  CCTK_REAL KtraceaverageiniL CCTK_ATTRIBUTE_UNUSED \= sum\*pow\(count\,\-1\)\;\\
  \\
  \#pragma omp parallel\\
  CCTK_LOOP3\(CosmoLapse_InitialTau\,\\
    i\,j\,k\,imin0\,imin1\,imin2\,imax0\,imax1\,imax2\,\\
    cctk_ash\[0\]\,cctk_ash\[1\]\,cctk_ash\[2\]\)\\
  \{\\
    const ptrdiff_t index CCTK_ATTRIBUTE_UNUSED \= di\*i \+ dj\*j \+ dk\*k\;\\
    Ktraceaverageini\[index\] \= KtraceaverageiniL\;\\
  \}\\
  CCTK_ENDLOOP3\(CosmoLapse_InitialTau\)\;"
sed -i "" "s/${OLD}/${NEW}/" $FILE

