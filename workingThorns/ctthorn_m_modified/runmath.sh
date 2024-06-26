#! /bin/bash

# Abort on errors
set -e

script=$1

if test -z "$script"; then
    echo "Usage:"
    echo "$0 <script.m>"
    exit 2
fi

error=$(basename $script .m).err
output=$(basename $script .m).out

rm -f $output

# Run Kranc to regenerate the code
/home/robynm/ET/code/2020_05/Cactus/repos/Kranc/Bin/kranc -v $script | tee $error
[ $PIPESTATUS -eq 0 ] || exit $PIPESTATUS

mv $error $output
