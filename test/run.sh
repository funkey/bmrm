#!/bin/bash

TRAIN=../build/linear-bmrm/linear-bmrm-train
CONF=test.conf

if [ "$1" == "gdb" ]
then gdb --args $TRAIN $CONF
else $TRAIN $CONF
fi
