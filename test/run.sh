#!/bin/bash

TRAIN=../build/linear-bmrm/linear-bmrm-train
CONF=test.conf

if [ "$1" == "gdb" ]
then gdb --args $TRAIN $CONF
else if [ "$1" == "valgrind" ]
then valgrind --leak-check=full $TRAIN $CONF
else
$TRAIN $CONF
fi
fi
