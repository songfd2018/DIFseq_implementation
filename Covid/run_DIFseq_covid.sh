#!/bin/sh
cd Covid/

# Run DIFseq to generate the parameter values at each iteration in the folder "MCESM_iter_K12", and the posterior inference in the folder "Inference_K12" 

../../src/DIFseq_MCESM -d./  -r../RawCountData/ -p covid -v 1 -K 8 -O2000 -M4000 -Q10 -s 5946 -c 16

../../src/DIFseq_MCESM -d./  -r../RawCountData/ -p covid -v 1 -K 9 -O2000 -M4000 -Q10 -s 5138 -c 16

../../src/DIFseq_MCESM -d./  -r../RawCountData/ -p covid -v 1 -K 10 -O2000 -M4000 -Q10 -s 6828 -c 16

../../src/DIFseq_MCESM -d./  -r../RawCountData/ -p covid -v 1 -K 11 -O2000 -M4000 -Q10 -s 1572 -c 16

../../src/DIFseq_MCESM -d./  -r../RawCountData/ -p covid -v 1 -K 12 -O2000 -M4000 -Q10 -s 8367 -c 16

../../src/DIFseq_MCESM -d./  -r../RawCountData/ -p covid -v 1 -K 13 -O2000 -M4000 -Q10 -s 1761 -c 16

../../src/DIFseq_MCESM -d./  -r../RawCountData/ -p covid -v 1 -K 14 -O2000 -M4000 -Q10 -s 3499 -c 16

../../src/DIFseq_MCESM -d./  -r../RawCountData/ -p covid -v 1 -K 15 -O2000 -M4000 -Q10 -s 4374 -c 16


