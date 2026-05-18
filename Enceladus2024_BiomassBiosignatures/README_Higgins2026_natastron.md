This directory contains the code used to produce and analyse the results in our manuscript entitled: *"A framework for evaluating biosignature potential against the abiotic baseline on ocean worlds"*. All the relevant code for that manuscript is in the IsotopeCycling and Chirality directories. **This implementation uses an early build of NutMEG v2 which is not compatible with the final release NutMEG v2. Please use the specific version of NutMEG [here](https://github.com/pmhiggins/NutMEG/tree/2742434ceb64f7bd47d7eb32a0a7f9a217b25bef) and ensure you install reaktoro v2 (not v1) to replicate.** Alternatively, a working build of both this code and NutMEG has been [archived on zenodo with doi 10.5281/zenodo.19820889 ](https://doi.org/10.5281/zenodo.19820889).

This code can replicate the isotope cycling models described in the manuscript. Some datasets for the bulk chemistry are from Higgins et al., (2024), JGR:Planets. The datasets for that work have been archived in figshare (doi: 10.6084/m9.figshare.22557706) and can be extracted and placed into the `data` directory. You can also regenerate them (Follow the Higgins2024 readme). None of the files below should take more than a few minutes to run on typical modern laptop hardware.

This code is designed to be used in the directory layout specified in the [main NutMEG-Implementations README](https://github.com/pmhiggins/NutMEG-Implementations/README.md) (ensure you follow instructions for NutMEG v1). For best compatability, use the version of NutMEG as of [this commit](https://github.com/pmhiggins/NutMEG/tree/2742434ceb64f7bd47d7eb32a0a7f9a217b25bef). The only additional dependency of this code compared to NutMEG's dependencies is `sympy`, which can be installed via PyPI.

- `IsotopeSpeciation/ChemicalIsotopeSpeciation.py` contains a class and methods for computing the carbonate isotope speciation, partitioning isotopes between CO2, HCO3 and CO3. It replicates Extended Data Fig 2.

- `OceanTransport/OceanAdvectiveZone.py` contains a class and methods for computing the compositional profile (in isotopes and bulk chemistry) in the fast, advective region of Encelaus's ocean (between x1 and x2 in the manuscript). It also contains some implementation code to compute ranges of net isotope difference across our expected Enceladus ocean parameter space. It replicates Extended Data Fig 4.

- `OceanTransport/OceanDiffusiveZone.py` contains the implementation of the dimensionless compositional profile for a stratified layer in Enceladus' ocean where diffusion is the only mechanism of vertical transport (between x2 and x3 in the manuscript). It replicates Extended Data Figs 5 and 6.

- `FracFactors/MethanogenInterpolation/GroppFits.py` includes the interpolation procedure for biological methanogenesis isotope enrichment factors. It replicates Extended Data Fig 7.

- `SteadyStateMassBalance/EnceladusSSMBdf.py` contains a general class which can represent the modules of the framework for CO2 or CH4, but is generalisable to any species and any arbitrary number of processes happening inside a module.

- `SteadyStateMassBalance/endmember_fractionation.py` contains a simple example deployment of the modules, used to create Fig 3 in the manuscript.

- The `Chirality` directory contains an example of computing the timescales to complete racemization of biotic and abiotic mixture in Enceladus' seafloor and oceans. The files `Chiral_summary.py` and `ChiralPlotter` contain code to replicate Fig 4 and Extended Data Fig 8 respectively.
