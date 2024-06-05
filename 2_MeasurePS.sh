#!/bin/bash
 
 
gfortran -I/usr/local/include -L/usr/local/lib /Users/fei/WSP/Scie/Proj12/Prog/PSrand.f90 /Users/fei/WSP/Scie/Proj12/Prog/PSestFun.f90  -O3 -o Proj -lm -lgsl -lgslcblas -lfftw3
./Proj
gfortran -I/usr/local/include -L/usr/local/lib /Users/fei/WSP/Scie/Proj12/Prog/PSdata.f90 /Users/fei/WSP/Scie/Proj12/Prog/PSestFun.f90  -O3 -o Proj -lm -lgsl -lgslcblas -lfftw3
./Proj
