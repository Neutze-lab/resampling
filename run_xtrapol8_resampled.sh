#!/usr/bin/env bash

###############################################################
### run_xtrapol8_resampled                 version 240724    ##
### (C) 2023 Adams Vallejos                                  ##
###############################################################

BASEDIR="/PATH/TO/resampling"
METHOD='B'

LIGHT="${@: -1}"
DARK="dark"

INDIR="$BASEDIR/output/${LIGHT}/{METHOD}"


X8(){
        phenix.python /PATH/TO/Xtrapol8/Fextr.py $@
}

for MTZIN in `find ${INDIR} -name '*_trunc.mtz'|sort`
do
    PREFIX=${MTZIN:0:-4}
    DARKMTZ="${MTZIN//$LIGHT/$DARK}"
    DARKMODEL=$(ls ${DARKMTZ:0:-4}*_001.pdb)
    echo Working on ${MTZIN}
    X8 \
    x8_params.phil \
    input.reference_mtz=${DARKMTZ} \
    input.triggered_mtz=$MTZIN \
    input.reference_pdb=$DARKMODEL \
    output.outdir=$(dirname $MTZIN)/xtra
done

echo $i Finished
