#!/bin/bash

dir=$(pwd)

for t in 500 600;do
	cd $dir
	mkdir $dir/$t
	for i in $(seq 0.00 0.02 1.00); do
	   	cd $dir/$t
		QVALUE=$(echo $i)
		#echo $QVALUE
		mkdir "$dir"/$t/"$i"
		cd $dir/$t/"$i"
		cp -- "$dir"'/'* "$dir"'/'"$t"'/'"$i"
		randomseed=$(shuf -i 2000-65000 -n 1)
		sed -i "s/TEMPERATURE/$t/g" foxp1.in
		sed -i "s/QVALUE/$i/g" fix_qbias_coeff.data
		sed -i "s/KHARMONIC/1300/g" fix_qbias_coeff.data
		sed -i "s/RANDOMSEED/$randomseed/g" foxp1.in
		sbatch sbatch_mazinger.sh
	done	
done
