#!/bin/sh
cd Simulation/v3

# Run DIFseq to generate the parameter values at each iteration in the folder "MCESM_iter_K5", and the posterior inference in the folder "Inference_K5" 

../../src/DIFseq_MCESM -d./  -r./RawCountData/ -p simulation -v 3 -K 3 -O2000 -M4000 -Q10 -E0.01 -s 4139 -c 16

../../src/DIFseq_MCESM -d./  -r./RawCountData/ -p simulation -v 3 -K 4 -O2000 -M4000 -Q10 -E0.01 -s 7697 -c 16

../../src/DIFseq_MCESM -d./  -r./RawCountData/ -p simulation -v 3 -K 5 -O2000 -M4000 -Q10 -E0.01 -s 9934 -c 16

../../src/DIFseq_MCESM -d./  -r./RawCountData/ -p simulation -v 3 -K 6 -O2000 -M4000 -Q10 -E0.01 -s 1113 -c 16

../../src/DIFseq_MCESM -d./  -r./RawCountData/ -p simulation -v 3 -K 7 -O2000 -M4000 -Q10 -E0.01 -s 8545 -c 16

../../src/DIFseq_MCESM -d./  -r./RawCountData/ -p simulation -v 3 -K 8 -O2000 -M4000 -Q10 -E0.01 -s 7521 -c 16


