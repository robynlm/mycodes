#! /bin/bash

script=ProperTime.m

\rm -rf ProperTime
set -e
error=$(basename $script .m).err
output=$(basename $script .m).out
rm -f $output
/Users/robynmunoz/EinsteinToolkit/Cactus/repos/Kranc/Bin/kranc -v $script | tee $error

