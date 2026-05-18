# NutMEG-Implementations
Repository for example or published codes using [NutMEG](http://github.com/pmhiggins/NutMEG).

**Note:** Not all Implementations have been updated to work with NutMEG v2.0, which is currently in development. We encourage users to use the appropriate NutMEG version as identified in the published manuscripts until v2.0 releases.

These are designed such that one can clone NutMEG-Implementations and position it in the same directory as their NutMEG clone. All code in NutMEG-Implementations appends the NutMEG directory to sys.path by default, so if you structure your projects differently you'll need to change it accordingly.

Example structure:

    projects/
    |--NutMEG/
    |   |--NutMEG/
    |   `--docs/
    `--NutMEG-Implementations/
        |--Venus/
            `--EVenus.py

Then to run some code:

    cd projects/NutMEG-Implementations/Venus
    python3 EVenus.py

## Implementations Quick Look

Below is a brief overview of each Implementation. Each directory has its own readME file breaking down its contents, please refer to those or the papers below for more information.

### Competition_Example
A guided example for using NutMEG to simulate competition between two fictional organisms: a methanogen and a sulfate reducer. Full guide available in [NutMEG's documentation](https://nutmeg-astrobiology.readthedocs.io).

### Venus
 Code for assessing the energetic availability in Venus' clouds.  Outputs of this were used in the publication: 'Biologically available Chemical Energy in the Temperate but Uninhabitable Venusian Cloud Layer: What do we Want to Know?' Cockell et al. 2021 *Astrobiology* DOI: [10.1089/ast.2020.2280](https://doi.org/10.1089/ast.2020.2280)

### TOM
Code for examining the growth and biosignature behaviour of methanogens.  Outputs of this were used in the publication: Higgins P.M. and Cockell C.S. (2020) 'A bioenergetic model to predict habitability, biomass and biosignatures in extreme conditions and astrobiology' *J. R. Soc. Interface* 17 (171) pp. 20200588 DOI: [10.1098/rsif.2020.0588](https://doi.org/10.1098/rsif.2020.0588)


### Enceladus2021_ParameterSpace
Code for mapping the habitability of the parameter space of Enceladus' ocean. Outputs of this were used in the publication: Higgins P.M.,  Glein C.R., and Cockell C.S. (2021) 'Instantaneous Habitable Windows in the Parameter Space of Enceladus' Ocean' *JGR: Planets* 126 (11) pp. e2021JE006951 DOI: [10.1029/2021JE006951](https://doi.org/10.1029/2021JE006951)

### Enceladus2024_BiomassBiosignatures
- Code for analysing the sustainable steady state biomass and turnover of biotic carbon in Enceladus' ocean used in the publication Higgins et al., (2024) 'Quantifying uncertainty in biomass and production of biotic carbon in Enceladus' notional methanogenic biosphere' *JGR: Planets* 129 (3) pp. e2023JE008166 DOI: [10.1029/2021JE006951](https://doi.org/10.1029/2023JE008166)
- Code for modelling isotope and enantiomeric cycling that owes to abiotic processes on Enceladus (the abiotic baseline) and possible biotic contributions above that baseline. Used in the manuscript Higgins et al., (submitted) 'A framework for evaluating biosignature potential
against the abiotic baseline on ocean worlds'. **This implementation uses an early build of NutMEG v2. Please use the version of NutMEG [here](https://github.com/pmhiggins/NutMEG/tree/2742434ceb64f7bd47d7eb32a0a7f9a217b25bef)**. Please refer to the [specific README](https://github.com/pmhiggins/NutMEG-Implementations/blob/IsotopeCycling/Enceladus2024_BiomassBiosignatures/README_Higgins2026_natastron.md) for more information.


### Thesis
Some useful code to accompany my PhD thesis. Thesis published at the University of Edinburgh: Higgins, P.M. (2022) 'Modelling extraterrestrial habitability, biomass and biosignatures through the bioenergetic lens' *The University of Edinburgh* DOI: [10.7488/era/2078](https://doi.org/10.7488/era/2078 )
