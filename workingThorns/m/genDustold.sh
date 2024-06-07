#! /bin/bash

\rm -rf CT_Dust
./runmath.sh Dust.m
#./dust.patchmac
./copy-if-changed.sh CT_Dust ../CT_Dust
