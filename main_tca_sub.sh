#!/bin/bash 

datafile=input.dat
for U1 in 8.0 10.0 12.0 14.0 16.0  
do
mkdir -p U_$U1
cd U_$U1

for tp in 0.0 #0.2 0.4 0.6 0.8 1.0 
do 
mkdir -p tp_$tp
cd tp_$tp

for strnth in  0.5 1.0
do
mkdir -p strnth_$strnth
cd strnth_$strnth

for seed in 87959 #98352 79869 69857 96387
do
 
mkdir -p seed_$seed

cd seed_$seed


cp ../../../../ionic_hub_tca_mf_mc.f90 .
cp ../../../../mtfort90.f .


#!******************************************************************************

echo   "8                          !d "             > $datafile
echo   "4                          !ds "            >> $datafile
echo   "-1.0d0                     !t1     "        >> $datafile
echo   $tp"0d0                      !t2     "        >> $datafile
echo   "14                         !temp_max "      >> $datafile
echo   "4000                        !MCSW "          >> $datafile
echo   "5                          !intrvl   "      >> $datafile
echo   $U1"d0                      !U1 "             >> $datafile  
echo   "1.0d0                      !filling "       >> $datafile
echo   "0.050d0                       !gama "          >> $datafile 
echo   "0.020d0                       !gama_m "        >> $datafile 
echo   "$seed                      !seed "          >> $datafile
echo   $strnth"0d0                        !strnth  "       >> $datafile
#!******************************************************************************

gfortran -o run.x  mtfort90.f ionic_hub_tca_mf_mc.f90 -llapack -lblas
nohup ./run.x  &

cd ..
done

cd ..
done 


cd ..
done

cd ..
done
