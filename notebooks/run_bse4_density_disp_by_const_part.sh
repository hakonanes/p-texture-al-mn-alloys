#!/bin/bash

declare -a SAMPLE=("0s"\ "175c"\ "300c"\ "325c")
declare -a DSET=("1"\ "2"\ "3")
declare -a THICKNESS_PFZ=(0.5\ 1.0\ 1.5\ 2.0) # um
declare -a CONSTITUENT_THRESHOLD=(0.7\ 0.8\ 0.9\ 1.0) # um

for i in ${SAMPLE}; do
    for j in ${DSET}; do
        for CT in ${CONSTITUENT_THRESHOLD}; do
            for TPFZ in ${THICKNESS_PFZ}; do
                echo ${i} ${j} ${CT} ${TPFZ}
                python bse4_density_disp_by_const_part.py ${i} ${j} ${CT} ${TPFZ}
            done
        done
    done
done
