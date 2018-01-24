#!/bin/csh

cd $HOME
source ConfigNemoLyon.csh
cd analysis/mc_test/

flreconstruct -p reco_config_BiPo -i BiPo_100000ev.brio -o BiPo_100000ev_reco.brio

