Welcome! This directory contains the code used to produce and analyse the results in our manuscript entitled: *"A bioenergetic model to predict habitability, biomass and biosignatures in astrobiology and extreme conditions"*. There are a few helper classes, and some then some for the plotting and analysis. The code is designed to be used in the director layout specified in the [main NutMEG-Implementations README](https://github.com/pmhiggins/NutMEG-Implementations/README.md).

`methanogen_implementer.py` contains the `efficiencies` class, which is a helper class to use the `NutMEmatcher` module of NutMEG to estimate the maintenance requirements of the empirical methanogens. It also sets up and performs the growth prediction simulations specified by the plotting files below.

`unique_efficiencies.py` contains methods to generate the parameters of the typical optimal methanogen [TOM] and some other specific empirical methanogens. These can then be used to generate `horde` objects for input to `NutMEG.ecosystem`.

`methanogen_extractor.py` contains methods for reading from and writing to the data directory for easier extraction of methanogen properties than NutMEG's native database. The native database is used for the growth simulation results.

`theory_emp_match.py` is for generating Figure 2 in the manuscript. It predicts the maintenance requirements for the empirical methanogens and the TOM.

`EnergyNutrientLimitation.py` is for generating Figures 3 and 4 in the manuscript and performing the necessary simulations.

`maintenance_Esynthmesh.py` is for generating Figure 5 in the manuscript.

`supplement.py` is for generating the supplementary figures.

`GCanimation.py` is for generating the growth curve animation (also a supplementary figure).

**Please note:** The data directory contains some of the data used in the manuscript but excludes the growth simulations due to the large NutMEG database this would create (several GB). When run, specifically `EnergyNutrientLimitation.py` this code will provide identical results but may take some time.
