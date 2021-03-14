# NutMEG-Implementations
Repository for example or published codes using NutMEG.

These should be designed such that one can clone NutMEG-Implementations and position it in the same directory as their NutMEG clone. All code in NutMEG-Implementations appends the NutMEG directory to sys.path by default, so if you structure your projects differently you'll need to change it accordingly.

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

**Competition_Example** A guided example for using NutMEG to simulate competition between two fictional organisms: a methanogen and a sulfate reducer. Full guide available in [NutMEG's documentation](https://nutmeg-astrobiology.readthedocs.io).

**Venus** Code for assessing the energetic availability in Venus' clouds.  Outputs of this were used in the publication: 'Biologically available Chemical Energy in the Temperate but Uninhabitable Venusian Cloud Layer: What do we Want to Know?' Cockell et al. 2021 Astrobiology [DOI](https://doi.org/10.1089/ast.2020.2280)

**TOM** Code for examining the growth and biosignature behaviour of methanogens.  Outputs of this were used in the publication: 'A bioenergetic model to predict habitability, biomass and biosignatures in extreme conditions and astrobiolgy' Higgins & Cockell (2020) J. R. Soc. Interface [DOI](https://doi.org/10.1098/rsif.2020.0588)
