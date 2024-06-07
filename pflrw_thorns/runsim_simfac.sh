#!/bin/bash
cd ~/code/2020_05/Cactus

RMPI=12
RNUMTHREADS=15
RPROCS=$(( $RNUMTHREADS * $RMPI))
PARDIR=~/mycodes/pflrw_thorns/par/$1.par

./simfactory/bin/sim submit $1 --machine=sciama --parfile=$PARDIR --procs=$RPROCS  --num-threads=$RNUMTHREADS
