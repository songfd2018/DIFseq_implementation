#!/bin/sh
cd Pancreas/

# Run DIFseq to generate the parameter values at each iteration in the folder "MCESM_iter_K12", and the posterior inference in the folder "Inference_K12" 

../../src/DIFseq_MCESM -d./  -r./RawCountData/ -p pancreas -v 1 -K 3 -O2000 -M2000 -Q10 -s 7118 -c 16

../../src/DIFseq_MCESM -d./  -r./RawCountData/ -p pancreas -v 1 -K 4 -O2000 -M2000 -Q10 -s 5959 -c 16

../../src/DIFseq_MCESM -d./  -r./RawCountData/ -p pancreas -v 1 -K 5 -O2000 -M2000 -Q10 -s 6151 -c 16

../../src/DIFseq_MCESM -d./  -r./RawCountData/ -p pancreas -v 1 -K 6 -O2000 -M2000 -Q10 -s 8567 -c 16

../../src/DIFseq_MCESM -d./  -r./RawCountData/ -p pancreas -v 1 -K 7 -O2000 -M2000 -Q10 -s 6807 -c 16

../../src/DIFseq_MCESM -d./  -r./RawCountData/ -p pancreas -v 1 -K 8 -O2000 -M2000 -Q10 -s 1448 -c 16

../../src/DIFseq_MCESM -d./  -r./RawCountData/ -p pancreas -v 1 -K 9 -O2000 -M2000 -Q10 -s 2885 -c 16

../../src/DIFseq_MCESM -d./  -r./RawCountData/ -p pancreas -v 1 -K 10 -O2000 -M2000 -Q10 -s 9029 -c 16

../../src/DIFseq_MCESM -d./  -r./RawCountData/ -p pancreas -v 1 -K 11 -O2000 -M2000 -Q10 -s 1171 -c 16

../../src/DIFseq_MCESM -d./  -r./RawCountData/ -p pancreas -v 1 -K 12 -O2000 -M2000 -Q10 -s 1533 -c 16



