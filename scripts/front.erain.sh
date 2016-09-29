#!/bin/csh -f 
#Check global variables in prog and prog1 !!!!!
#eg  mlon=145 mlat=73 for NCEP            !!!!!

set d = erain
set lev = $3

 set prog=front.$d
 gfortran -o ${prog} ${prog}.f

 set prog1=seg2.$d
 gfortran -o $prog1 $prog1.f /work/irudeva/include/libeda.a

set Din = ../output/$d/$lev
set Dout = ../output/$d/$lev

set x = all # modified crit I - actually I and VII
set cdv = 2 #critical dv (_C option in front_seg1)

#set H = ( "0S" "0N" )   #Hemisphere
#set nH = 1

#while ($nH <= 1 )
set H = $2

#@ y = 2007
set y=$1

#while ($y <= 2007)
set farea = $Dout/farea$lev.$d.$y.${H}.dat

 echo 'start ' $prog
 ./${prog} -C $cdv -dv $Din/dv${lev}.$d.$y.${H}.cmp  -a $farea 
#
 ./${prog1} -i $farea -o $Dout/front${lev}.$d.$y.$H.dat -h $H
#  \rm -f $farea
#
#  @ y ++
#end

# @ nH ++
#end
  \rm -f ${prog} ${prog1}
exit
