# GBP: An open source library for generalized belief propagation decoding of quantum codes

A C++ library implementing a Generalized Belief Propagation [YFM-GBP](https://ieeexplore.ieee.org/abstract/document/14590449). It is currently restricted to 
- two-layer Regiongraphs with one check and its adjacent qubits in top- and intersection qubits in bottom-regions
and provides special features
- hard-decision based on beliefs over top-regions

## Build Library

### Prerequsites
Make sure to have a C++17 compatible compiler and CMAKE at least Version 3.16.
Libraries used are
- [JSON for Modern C++](https://github.com/nlohmann/json) for I/O handling
- [NTL](https://github.com/libntl/ntl) for computation of GF(2) ranks
- [gmp](https://gmplib.org) in conjunction with NTL
- [xtensor](https://github.com/xtensor-stack/xtensor), [xtensor-blas](https://github.com/xtensor-stack/xtensor-blas) and their dependencies as data structure
- [lemon graph library](https://lemon.cs.elte.hu/trac/lemon) as graph structure
- [libDAI](https://staff.fnwi.uva.nl/j.m.mooij/libDAI/) for factors

The json header only library included in the `include` directory should work right away. In order to work with the given CMakeLists files, best install NTL with gmp, xtensor et al. and lemon into the root directory, such that header files are in `include`.


### Build
Then simply run from the repository root

```
mkdir build
cd build
cmake ..
make
```

## Install and Run Example Simulations

We provide examples illustrating the usage of the library. They can be found within the `sim` directory and are compiled with the library. In order to get the binaries to the simulation folder just run

```
make install
```

The examples are
- `sim.cpp` for binary GBP decoding

### Run examples
The example scripts use input files to set the properties of decoding. Exemplary input files are given in `sim/input/input_files`.
The binaries use two command line parameters to specify the input file and the output directory, an exemplary run command is

```
./sim input/input_files/input_test.json output
```

## Quantum Codes

Some codes are provided in `sim/input/codes`, in the [alist](http://www.inference.org.uk/mackay/codes/alist.html) format and numpy `.npy`.
- `random_irregular`,`random_regular`: codes constructed from random classical base codes with the hypergraph product construction ([arXiv:0903.0566](https://arxiv.org/abs/0903.0566))
- `surface`: topological surface codes (e.g. [arXiv:quant-ph/0110143](https://arxiv.org/abs/quant-ph/0110143))
- `xzzx`: twisted surface codes for biased channels ([arXiv:2009.07851](https://arxiv.org/abs/2009.07851))

