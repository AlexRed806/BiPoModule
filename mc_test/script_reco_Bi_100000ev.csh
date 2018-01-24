#!/bin/csh

cd $HOME
source ConfigNemoLyon.csh
cd analysis/mc_test/

flreconstruct -p reco_config_Bi -i Bi_100000ev.brio -o Bi_100000ev_reco.brio
