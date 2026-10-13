# test-analysis

How to build:
1. `source /cvmfs/sft.cern.ch/lcg/views/LCG_107/x86_64-el9-gcc11-opt/setup.sh`
2. Go to test-analysis directory
3. `mkdir build`
4. `cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DUSE_CONAN=OFF`
5. `cd build/src`
6. `make`

Make sure the enviroment is set by command #1 before `make`, otherwise correctionlib does not work.
# hmm-analyzer
