#!/bin/bash

# codigos de funcoes do VNS
# FUNCTIONS_CODE=("0" "1" "2" "3" "4" "5" "6" "7" "8" "9" "11" "12" "13" "14" "15" "16" "17" "18" "19" "20" "21" "22" "23" "25" "26" "27" "28" "30" "31" "32" "33" "34" "35" "36" "37" "38" "39" "40" "41" "42" "43")
# FUNCTIONS=( "RASTRIGIN" "SCHAFFER_F7" "GRIEWANK" "ACKLEY" "ROSENBROCK" "SPHERE" "MPE" "SCHAFFER_F6" "SCHWEFELS226" "STEP" "PENALIZED1" "LEVY" "ZAKHAROV" "EGG_HOLDER" "HOLZMAN" "MICHALEWITZ" "PENALIZED2" "POWEL" "RANA" "SHUBERT" "STRETCHEDV" "MULTIMOD"  "SCHWEFEL222" "SHIFTED_SPHERE" "SFHITED_SCHWEFEL221" "SHIFTED_ROSENBROCK" "SHIFTED_RASTRIGIN" "SHIFTED_ACKLEY" "SHIFTED_SCHWEFEL222" "SHIFTED_SCHWEFEL12" "SHIFTED_EXTENDED_f10" "SHIFTED_BOHACHEVSKY" "SHIFTED_SCHAFFER" "HYBRID_1"  "HYBRID_2" "HYBRID_3" "HYBRID_4" "HYBRID_5" "HYBRID_6" "HYBRID_7" "HYBRID_8")

FUNCTIONS_CODE=("0" "1" "2" "3" "4" "5" "7" "12" "13" "23")
FUNCTIONS=( "RASTRIGIN" "SCHAFFER_F7" "GRIEWANK" "ACKLEY" "ROSENBROCK" "SPHERE"  "SCHAFFER_F6"  "LEVY" "ZAKHAROV"  "SCHWEFEL222" )
DIMENSIONS=("250")

METHODS_CODE=( "4" )
METHODS=( "RVNS" )

P_CROSS=("0.9" )

	echo "Creating dirs ..."
	echo "Generating input file ..."

	k=0;
	for p in "${P_CROSS[@]}"
	do
		j=0;
		for m in "${METHODS[@]}"
		do
			i=0;
			for f in "${FUNCTIONS[@]}"
			do
				for d in "${DIMENSIONS[@]}"
				do
					mkdir -p $m/$f/$d
					cp input.in $m/$f/$d/
					cd $m/$f/$d/
					sed -i.bak s/FUNCTIONS/${FUNCTIONS_CODE[$i]}/g input.in
					rm input.in.bak
					sed -i.bak s/METHODS/${METHODS_CODE[$j]}/g input.in
					rm input.in.bak
					sed -i.bak s/P_CROSS/${P_CROSS[$k]}/g input.in
					rm input.in.bak
					cd ../../../
				done
				let i++;
			done
			let j++;
		done
		let k++;
	done
	
	echo "Generation complete!"
	echo "Runing all ..."
	for m in "${METHODS[@]}"
	do
		for f in "${FUNCTIONS[@]}"
		do
			# mkdir Convergencia
			#wait finish before run more testes
			while pgrep "algorithm" > /dev/null; do sleep 10; done
				gnome-terminal -x sh -c "{ time ./algorithm $m/$f/250/input.in > $m/$f/250/output.out ; } 2>> $m/$f/250/output.out ;"
				sleep 1
			while pgrep "algorithm" > /dev/null; do sleep 10; done

			# media=`ls Convergencia/*Media.data`
			# melhor=`ls Convergencia/*Melhor.data`
			# saida_img=$m-$f-250.png
			# ./gnuplot_grafico.sh $saida_img $media $melhor
			# mv $saida_img Convergencia
			# mv Convergencia $m/$f/250/
		done
	done
