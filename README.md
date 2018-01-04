## Dijet Analysis Package

Common repo for dijet analysis code, starting with my di-jet imbalance measurement.

### Dependencies
1) ROOT (tested with ROOT6, https://root.cern.ch)
2) FastJet (http://fastjet.fr)
3) Pythia8 (http://home.thep.lu.se/~torbjorn/Pythia.html)
4) Boost C++ (http://www.boost.org)
5) eventStructuredAu (https://github.com/kkauder/eventStructuredAu)

### To Install
```mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=~/dijet_analysis (or wherever you choose to install it)
make
make test (optional but recommended)
make install
```

Default cmake install location is /usr/local, but since this isn't a system utility, I override the default
and make the install location ${CMAKE_BINARY_DIR}/install
