#!/bin/bash

##usage 
## ./run.sh example_atom_position.dat

##example of converting file format
##atomic.dat format example 2 atoms , 2 sets of data
#step: 50   ##molecular dynamics steps, actually meaningless in msd program
#0.3 0.3 0.3
#0.4 0.4 0.4 
#step: 100  
#0.31 0.31 0.31
#0.43 0.43 0.43 

 
num_atom=3
#cat $1 | awk 'BEGIN{a=0;b="'$num_atom'"+1} {if(b=="'$num_atom'"+1){b=0 ;a=a+50; print $2 " " a }else{print  $1 " " $2 " " $3 " "};b=b+1}' >atomic.dat

if [ ! -n "$1" ]; then
   echo "atome filename IS NULL"
   echo "usage:"
   echo "      ./run.sh example_atom_position.dat"
   exit 232
else
    echo "atome file name is: $1 "
fi



echo 'recompile?(y/n):'
read ichose
if [ $ichose = 'y' ]; then
 echo 'method ?(1/2):'
 read ichose
 if [ $ichose = '1' ]; then
  ifort mean_square_displacement.f90 -o msd.exec
 elif [ $ichose = '2' ] ; then
  ifort difussion_coefficient.f90 -o msd.exec
 else
  echo "bad choice! exit"
  exit 222
 fi
fi

echo 'convert file format?(y/n):'
read ichose
if [ $ichose = 'y' ]; then
cat $1 | awk 'BEGIN{a=0;b="'$num_atom'"+1} {if(b=="'$num_atom'"+1){b=0 ;a=a+50; print a " " $1 }else{print  $1 " " $2 " " $3 " "};b=b+1}' >atomic.dat
else
cp $1 atomic.dat
fi


./msd.exec

cat  output_info.txt 
