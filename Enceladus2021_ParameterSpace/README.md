Welcome! This directory contains the code used to produce and analyse the results in our manuscript entitled: *"Instantaneous habitable windows in the parameter space of Enceladus' ocean"*. The code is designed to be used in the director layout specified in the [main NutMEG-Implementations README](https://github.com/pmhiggins/NutMEG-Implementations/README.md).

There are a number of helper classes to gather the data into the formats presented in the manuscript, but the core use of NutMEG here was to compute the Gibbs free energy of methanogenesis and methanogenic power supplies and demands for a wide number of physicochemical and organismic parameters.

This code makes use of the in-built NutMEG saved_organism `TypicalOptimalMethanogen` and the saved_system `Enceladus`.

# Plotting tools

`EnergyPowerPlots.py` contains calls to relevant functions to plot Figures 1, 2, 3 and 6, and supplemental figures S5, S6, S7, S10 and S11. The keywords of these can be changed (or more can be added) to change the style of plot (contour vs mesh) or geochemical calculations (with/without speciation, salt level etc.)

`SampleHabitabilityPlots.py` similarly contains calls to functions for generating figures 4, 5 and Supplemental figures S8 and S9. This file can also be used to perform variance calculations with different sample sizes.

`EnceladusPlotStyles.py` contains methods for generating the plots in `EnergyPowerPlots.py`

`SamplingPlotStyles.py` contains methods for generating the plots in `SampleHabitabilityPlots.py`, and the supplemental variance plots.

`colormapping.py` contains the `cmapper` class for custom colormaps.

# Code for energy calculations

`EncergyCalculations.py` contains methods for building and Enceladus-loving `TypicalOptimalMethanogen`, extracting its bioenergetic parameters, and computing its microbially accessible power supply.

`QuotientUncertainties.py` contains methods for computing the quotient of methanogenesis when activities of CO2, H2 and CH4 have such large errors that the `uncertainties` package cannot adequately deal with them. Using these methods upper, lower and nominal bounds for the quotient, power supply, free energy available and rate constant can be calculated for various defined uncertainties in salt content.

`SpecitationGrids.py` is for extracting information from the carbonate speciation data files.

`EnceladusGrids.py` is for building grids of bioenergetic output parameters across the parameter space.

# Code for sampling

`Sampler.py` contains methods to sample across the parameter space.

`GetSamples.py` contains methods to plot on multiple axes samples across the parameter space.
