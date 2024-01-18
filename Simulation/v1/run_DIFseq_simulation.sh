#!/bin/sh
cd Simulation/v1

# Run DIFseq to generate the parameter values at each iteration in the folder "MCESM_iter_K5", and the posterior inference in the folder "Inference_K5" 

../../src/DIFseq_MCESM -d./  -r./RawCountData/ -p simulation -v 1 -K 3 -O2000 -M1000 -Q10 -E0.01 -s 1660 -c 16

../../src/DIFseq_MCESM -d./  -r./RawCountData/ -p simulation -v 1 -K 4 -O2000 -M1000 -Q10 -E0.01 -s 8186 -c 16

../../src/DIFseq_MCESM -d./  -r./RawCountData/ -p simulation -v 1 -K 5 -O2000 -M1000 -Q10 -E0.01 -s 6230 -c 16

../../src/DIFseq_MCESM -d./  -r./RawCountData/ -p simulation -v 1 -K 6 -O2000 -M1000 -Q10 -E0.01 -s 5998 -c 16

../../src/DIFseq_MCESM -d./  -r./RawCountData/ -p simulation -v 1 -K 7 -O2000 -M1000 -Q10 -E0.01 -s 8641 -c 16

../../src/DIFseq_MCESM -d./  -r./RawCountData/ -p simulation -v 1 -K 8 -O2000 -M1000 -Q10 -E0.01 -s 1021 -c 16
