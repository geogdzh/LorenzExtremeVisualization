# LorenzExtremeVisualization

Repository containing all code to generate and analyze the data used in "The Evolving Butterfly: Statistics in a Changing Attractor" (Geogdzhayev, Souza, and Ferrari, 2024). The repository is structured as follows:

### Data Generation
generate_data.jl generates all of the data used in the paper. The individual sections of the script can be called separately to only generate some sets of data. All generated data will be stored in a data/ folder.

### Analysis
analyze_data.jl runs three separate scripts, each of which generate a subset of the plots shown in the paper. The different scipts can be run independently. All generated figures are saved into a figs/ folder.
