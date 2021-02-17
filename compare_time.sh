#!/usr/bin/env bash

# Script to compare two galsims with respect to time.
# To compare two binaries, first make with "make" and then do
# changes to the source. Then run "make galsim_2" to get the second binary.
# Then run this script to see how they compare.

# Command to time and print seconds:
# ( /usr/bin/time -f '%e' -p ./galsim 3000 input_data/ellipse_N_10000.gal 200 0.00005 0 ) 2>&1 | awk '/real/{print $2}'
if [[ -e res_N$1_P$2 ]]; then
    rm res_N$1_P$2;
fi
echo N, gal1, gal2
for (( i = 0; i <= $1; i=i+$2 )); do
    no1=$((/usr/bin/time -f '%e' -p ./galsim $i input_data/ellipse_N_10000.gal 200 0.00005 0 ) 2>&1 | awk '/real/{print $2}');
    no2=$((/usr/bin/time -f '%e' -p ./galsim_2 $i input_data/ellipse_N_10000.gal 200 0.00005 0 ) 2>&1 | awk '/real/{print $2}');
    echo $i, $no1, $no2 | tee -a time_comp/res_N$1_$2.csv;
done

play ~/Ding-sound-effect.mp3;
