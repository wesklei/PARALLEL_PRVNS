#!/bin/bash

display_usage() { 
	echo "This script must be run with this arguments." 
	echo "   $1 saida.png"
	echo "   $2 arquivo_media"
	echo "   $3 arquivo melhor"
}

if [  $# -le 2 ] 
then 
	display_usage
	exit 1
fi 

echo "
set bmargin 7
unset colorbox
set terminal png enhanced font 'Verdana,10'
set output '$1'
set style line 2  lc rgb '#0404E9' lt 1 pt 7 
set style line 3  lc rgb '#B22C2C' lt 1 pt 7
set xlabel 'Geracoes'
set ylabel 'FO'
set title 'PRVNS'
plot '$2' using 1:2 title 'Media' ls 2 with lines, '$3' using 1:2 title 'Melhor' ls 3 with lines" | gnuplot


