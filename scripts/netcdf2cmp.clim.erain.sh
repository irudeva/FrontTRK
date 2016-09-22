:
d=erain
prog=netcdf2cmp.clim.$d
#gfortran -o $prog $prog.f /work/irudeva/include/libnetcdf.a
gfortran -o $prog $prog.f  -lnetcdff -lnetcdf
:
Dout=../../$d
HS=0S
var=( "u" "v" )
nvar=( "U500" "V500" )

#while [ $n -le 1 ]
#do 
#y=1980
fcmp=$Dout/uv500.$d.clim.$HS.cmp
echo $fcmp
:
#py=$[$y - 1]
#ny=$[$y + 1]
sy=1979
ey=2011
start=${sy}01
end=${ey}12
fin=/work/irudeva/DATA/ERAint/500/$d.uv_H500.$sy-${ey}mon.nc
:
./$prog <<mark
$fin
$HS
$start
$end
${var[0]}
${var[1]}
$fcmp
${nvar[0]}
${nvar[1]}
mark
:
#n=$[$n+1]
#done
:
rm -f $prog

