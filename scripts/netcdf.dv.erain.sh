:
d=erain
prog=netcdf.dv.$d
#gfortran -o $prog $prog.f /work/irudeva/include/libnetcdf.a
gfortran -o $prog $prog.f  -lnetcdff -lnetcdf
:
#lev="10m"
lev=$3
Dout=../../$d/${lev}
#HS=0N
HS=$2
if [ $lev -eq "10" ]; then
Dnc="/work/irudeva/DATA/ERAint/10m/$d.uv_10m."
else
Dnc="/work/irudeva/DATA/ERAint/$lev/$d.uv_$lev."
fi
var=( "u"${lev} "v"${lev} )
miss=-99.99
y=$1
sy=$4
ey=$5
lmon=$6
#y=2007
y=$1
#while [ $y -le $ey ]
#do 
fcmp=$Dout/dv${lev}.$d.$y.$HS.cmp
#echo $fcmp
:
py=$[$y - 1]
ny=$[$y + 1]
start=${py}1101
if [ $y -eq $sy ];then
start=${y}0101;fi
end=${ny}0131
if [ $y -gt $ey ];then
end=${y}0430;fi
:
./$prog <<mark
$HS
$lev
$start
$end
$ey
$lmon
${var[0]}
${var[1]}
$Dnc
$fcmp
DV$lev
$miss
mark
:
#y=$[$y+1]
#done
:
rm -f $prog

