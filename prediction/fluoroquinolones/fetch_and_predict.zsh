#!/bin/zsh

n=$1

rm ./*.xyzmq

echo "CIPR"
scp 'xucube:~/gen3d/cipr_p*_*.xyzmq' .
./cipr1.py "$n"                         
./predict.py cipr 

echo "ENOX"
scp 'xucube:~/gen3d/enox_p*_*.xyzmq' .
./enox1.py "$n"                        
./predict.py enox 

echo "ENRO"
scp 'xucube:~/gen3d/enro_p*_*.xyzmq' .
./enro1.py "$n"                       
./predict.py enro

echo "LEVO"
scp 'xucube:~/gen3d/levo_p*_*.xyzmq' .
./levo1.py "$n"                         
./predict.py levo

echo "LOME"
scp 'xucube:~/gen3d/lome_p*_*.xyzmq' .
./lome1.py "$n"                         
./predict.py lome

echo "NORF"
scp 'xucube:~/gen3d/norf_p*_*.xyzmq' .
./norf1.py "$n"                         
./predict.py norf 

echo "ORBI"
scp 'xucube:~/gen3d/orbi_p*_*.xyzmq' .
./orbi1.py "$n"                         
./predict.py orbi

echo "PEFL"
scp 'xucube:~/gen3d/pefl_p*_*.xyzmq' .
./pefl1.py "$n"                        
./predict.py pefl 

