# Source code and data for Bayer, Lautenbach, Arneth *Benefits and trade-offs of optimizing global land use for food, water, and carbon*

data 
contains the data used for the analysis

The code consists of three stages: 
the first one is written in R, the other two in python 2.7.
Python 2.7 can be e.g.used conveniently based on anaconda, see [documentation](https://docs.anaconda.com/free/anaconda/configurations/switch-environment/).

Python modules used: 

  - numpy
  - scipy
  - matplotlib
  - seaborn
  - fnmatch
  - itertools
  - re
  - pandas
  - sys
  - os

R-packages used:

- emoa
- foreign
- geometry
- lattice
- sp
- maptools
- hexbin
- raster
- sf
- tidyverse

## Running the code

Each stage returns results and stores them under results. For the users convenience, some of the results as used in the publication are provided.

For python you should be able to use *python* if you activate tthe conda environment with the python 2.7 interpreter but this might differ for other settings.

### Biome level

tradeOffsBiomeLevel.R

Performs a brute force approach by running all biome level combinations and when calculating the solutions that span the Pareto-frontier. Results are used by the fpu level optimization.

### FPU level

Start optimization based on results stored in results/biomes

python lpj_pareto_fpu_forage.py -r <runid> --mutation_rate= --corFac= --tournament_size= --max_generations= --time_period=<2042|2099> --rcp=<26|60> --pop_extra=' 

- the run-id is used to create a subdirectory to store the results in
- time_period and pop_extra specify which scenario to use
- pop_extra specifies the number of random seeds to be added to the population read from the biome level optimization. Intended to increase the search space a bit
- mutation_rate specifies the mutation rate of the GA (default: 0.0125)
- tournament-size specifies how many candidates are drawn for the selection of the fittest (defaults to 2)
- correction factor: a factor multiplied with the constraints (the ecosystem service provisioning for the reference sitaution (landuse 2017)). As it is hard at the fpu level to be better than the current conditions it might be good to lower the constraints a bit . Defults to 0.98 (i.e. .98*current food prvisdion/carbon storag/water supply/forage provision are used as constraints)
- max_generation: for how long shall the optimization run - defualts to 500 which takes quite a while

Result: 
  - results/fpu/run_<run-ID>/
  - config4run<run-ID>.txt - configuration if the optimization run
  - reportArchiveNSGAII.txt - reports on the GA archive (i.e. the solutions at the Pareto-frontier) across the different generations
  - reportNSGAII.txt - reports on the GA population across the different generations
  - nsGA_1.pkl the genomes of the population. These will be used by the postprocessing scripts to generate input for the cell level optimization (note: the pickle files have been compressed)
  
postprocessing:

run createGenomeTableFromFront_fpu.py for the post processing of a single run or prepareSeedsFromSeveralRuns to combine the pareto fronts across several runs

the script does not use commandline parameters so far but has to be adjusted in the code to run for the 4 different scenarios

### Cell level

#### Helper code

code/escpy-0.4 the evolutionary algorithm module used as the basis for the optimization code. A couple of files are replace in the FPU and cell level code directiors (*my_archivers* used instead of *archivers*)

``ecspy`` -- A framework for creating evolutionary computations in Python.


Original code from   

* Homepage: http://ecspy.googlecode.com
* Email: aaron.lee.garrett@gmail.com
