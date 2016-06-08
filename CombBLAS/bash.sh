#!/bin/bash

./Applications/ppcl TESTDATA/SCALE16BTW-TRANSBOOL  
mpirun -np -4  ./Applications/ppcl TESTDATA/SCALE16BTW-TRANSBOOL  
mpirun -np -8  ./Applications/ppcl TESTDATA/SCALE16BTW-TRANSBOOL  
