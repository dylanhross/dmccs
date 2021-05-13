#!/bin/zsh

n=$1

rm ./*_meas_pred_ccs.png

echo "CIPR"                     
./plot_ccs.py $n cipr 165

echo "ENOX"
./plot_ccs.py $n enox 160

echo "ENRO"
./plot_ccs.py $n enro 155

echo "LEVO"
./plot_ccs.py $n levo 150

echo "LOME"
./plot_ccs.py $n lome 170

echo "NORF"
./plot_ccs.py $n norf 165

echo "ORBI"
./plot_ccs.py $n orbi 180

echo "PEFL"
./plot_ccs.py $n pefl 150

