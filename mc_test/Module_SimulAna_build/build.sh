#!/bin/sh

make -j4

cp -r libSimulAna.* $FALAISE_PATH/lib64/Falaise/modules/
