This directory contains the code used to produce and analyse the results in our manuscript entitled: *"A framework for evaluating biosignature potential against the abiotic baseline on ocean worlds"*. All the relevant code for that manuscript is in the IsotopeCycling and Chirality directories.

This code can replicate the isotope cycling models described in the manuscript. Some datasets or the bulk chemistry are from Higgins et al., (2024), JGR:Planets. The datasets for that work have been archived in figshare (doi: 10.6084/m9.figshare.22557706) and can be extracted and placed into the `data` directory.

This code is designed to be used in the directory layout specified in the [main NutMEG-Implementations README](https://github.com/pmhiggins/NutMEG-Implementations/README.md). For best compatability, use the nm_v2 branch of NutMEG. The only additional dependency of this code compared to NutMEG's dependencies is `sympy`, which can be installed via `pip`.

- `IsotopeSpeciation/ChemicalIsotopeSpeciation.py` contains a class and methods for computing the carbonate isotope speciation, partitioning isotopes between CO2, HCO3 and CO3.

- `OceanTransport/OceanAdvectiveZone.py` contains a class and methods for computing the compositional profile (in isotopes and bulk chemistry) in the fast, advective region of Encelaus's ocean (between x1 and x2 in the manuscript). It also contains some implementation code to compute ranges of net isotope difference across our expected Enceladus ocean parameter space.

- `OceanTransport/OceanDiffusiveZone.py` contains the implementation of the dimensionless compositional profile for a stratified layer in Enceladus' ocean where diffusion is the only mechanism of vertical transport (between x2 and x3 in the manuscript).

- `FracFactors/MethanogenInterpolation/GroppFits.py` includes the interpolation procedure for biological methanogenesis isotope enrichment factors.

- `SteadyStateMassBalance/EnceladusSSMBdf.py` contains a general class which can represent the modules of the framework for CO2 or CH4, but is generalisable to any species and any arbitrary number of processes happening inside a module.

- `SteadyStateMassBalance/endmember_fractionation.py` contains a simple example deployment of the modules, used to create Fig 3a-d in the manuscript.

- The `Chirality` directory contains an example of computing the timescales to complete racemization of biotic and abiotic mixture in Enceladus' seafloor and oceans. The files `Chiral_summary.py` and `ChiralPlotter` contain code to replicate Fig 3e-h and Extended Data Fig 8 respectively.
