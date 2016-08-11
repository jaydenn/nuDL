#!/bin/bash

rm -rf file.list
let i=1

for source in "TAMU_reactor  1" "Sun  1"
do

for det in "#GERMANIUM  20000" "#SILICON  10000"
do

for N in 0 1
do

let E=1

if [ ${N} -ne 0 ]; then
let E=0
fi

for BSM in 1 2 3 4
do

filename=config${i}.dat

cat > ${filename} << EOF
//Conudl configuration options
1                  // running mode: 1=rate calculation, 2=disc. limit evolution, 3=exclusion limits
./results/CN_      // root of output file names, can include directories
${N}               // do nuclear scattering? 1=yes, 0=no
${E}               // do electron scattering? 1=yes, 0=no
${BSM}             // which BSM to consider? (0=SM, 1=scalar, 2=pseudoscalar, 3=vector, 4=axialvector)
1e-10	//initial coupling
1e-3    // mediator mass in GeV
//Source
${source}   1
// Detector  |  Exposure (kg*days) (or starting exposure in evolution mode)
${det}
//Asimov data set/random Monte-Carlo (1/0)
1
//log bins
1
EOF
echo "conf/${filename}" >> file.list

let i=i+1

done

done

done

done
