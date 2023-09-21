# Test functionality
# load seed candidates at FPU level

import numpy, numpy.random
import pandas as pd
import sys, os,  fnmatch
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as mplcol


#codePath = "d:/Nutzer/development/python/GA/ecspy-0.4/"
#sys.path.append(codePath)
basePath = "d:/Nutzer/projekte_neu/ess_tradeoffs/LPJ_pareto/"
#sys.path.append(basePath + "code/version_fpu_updated_july_2019/")
inPath = basePath + "data/update_july_2019/"
sys.path.append(basePath + "code/version_fpu_updated_july_2019/")

import lpj_pareto_functions as lpf

rcp= 26
time_period_long= "2033-2042"
inPath = basePath + "data/update_july_2019/" # + ESquant_final_v2 ??
inStartGenomesFolderBase = basePath + "results/version_july_2019/biomes/fpuGenomes/"
inStartGenomesFolder = inStartGenomesFolderBase + "fpuGenomes_" + str(rcp) + "_" + str(time_period_long) + "/"
nSeeds = len(fnmatch.filter(os.listdir(inStartGenomesFolder), 'fpuGenome_*.csv'))

fpuGrid = pd.read_csv(basePath + "data/update_july_2019/fpu2biomeGrid.csv", sep=";")

fnamePart1= inPath + "/ESquant_final_v2/ESquant_rcp" + str(rcp) + "_" + str(time_period_long) + "_INT_9426_v8b_4gcms_"
fnamePart3 = "_cor4ProtArea.txt"
Cstorage = pd.read_csv( fnamePart1 + "Cstor" + fnamePart3, delim_whitespace=True)
FoodProv = pd.read_csv(fnamePart1 + "FoPro" + fnamePart3, delim_whitespace=True)
WaSu = pd.read_csv( fnamePart1 + "WaSu" + fnamePart3, delim_whitespace=True)
Forage = pd.read_csv(fnamePart1 + "Forage" + fnamePart3, delim_whitespace=True)

# convert to numpy for evaluation of objective functions
Cstorage_np = Cstorage.values
WaSu_np = WaSu.values
FoodProv_np = FoodProv.values
Forage_np = Forage.values

# list of objectives used in the objective function
listOfObjectiveLookupTables = [Cstorage_np, FoodProv_np, WaSu_np, Forage_np]

# maximum values for each objective
# corrected for additional column (pasture), endPos=12 includes this colum
# before it was endPos=11
Cmax = getMaximumObjectiveValue(Cstorage, startPos=2, endPos=12)
Fmax = getMaximumObjectiveValue(FoodProv, startPos=2, endPos=12)
Wmax = getMaximumObjectiveValue(WaSu, startPos=2, endPos=12)
Foragemax = getMaximumObjectiveValue(Forage, startPos=2, endPos=12)

# minimum values for each objective
Cmin = getMinimumObjectiveValue(Cstorage, startPos=2, endPos=12)
Fmin = getMinimumObjectiveValue(FoodProv, startPos=2, endPos=12)
Wmin = getMinimumObjectiveValue(WaSu, startPos=2, endPos=12)
Foragemin = getMinimumObjectiveValue(Forage, startPos=2, endPos=12)

if Wmin < 0:
    Wmin = 0

maximumObjectives = (Cmax, Fmax, Wmax, Foragemax)
minimumObjectives = (Cmin, Fmin, Wmin, Foragemin)

foodMinThreshold = 0  # has already been done in the corrrected inpuit files, 13.08.2018

# fpuOrder is the reference with respect to the order in the geome
# should be the same but better careful than sad

fn = inStartGenomesFolder + "fpuGenome_1.csv"
aFpu=  pd.read_csv(fn, sep=";")
fpuOrder = aFpu.newFPU
selectedGenomesList =[]


for i in range(0, nSeeds):

    fn = inStartGenomesFolder + "fpuGenome_" + str(i +1) + ".csv"
    try:
        aFpuGenome = pd.read_csv(fn, sep=";")
        #genom_array = fpuGenome.PFT  # if the field is named differently this needs to be changed...needs a better solution
    except (OSError, IOError) as e:
        print "generator_predef: \n" + e.message
        raise e

    #aFpuGridWithGenome = pd.merge(fpuGrid, aFpuGenome, on=('newFPU'))
    ## need to order accordingly to the objectives
    #aFpuGridWithGenome.sort_values(by='ID', inplace=True)
    #aFpuGridWithGenome.reset_index(drop=True, inplace=True)

    #selectedGenomesList.append(aFpuGridWithGenome['PFT'])
    selectedGenomesList.append(aFpuGrid['PFT'])
    if(i==2):
        break

# join first genome with FPUgrid to allow plotting in space
aFPUGenome = selectedGenomesList[0]
aFpuGridWithGenome = pd.merge(fpuGrid, aFpuGenome, on=('newFPU'))

# plot map
cm = plt.cm.get_cmap('Set1')
#colors = ['xkcd:dark olive green', 'xkcd:sea green', 'xkcd:sea blue', 'xkcd:orange', 'xkcd:light plum',
#          'xkcd:dark plum']
#cm = mplcol.ListedColormap(cm)

titleStr = "test"

labels = ('PNV', 'Wheat', 'Rice', 'Corn', 'Soybean', 'Tropical roots')

fig, ax = plt.subplots()
fig.set_size_inches(14, 8)

s = ax.scatter(aFpuGridWithGenome['Lon'], aFpuGridWithGenome['Lat'], c=aFpuGridWithGenome['PFT'], s=1.5, cmap=cm, label=labels)
ax.set_title(titleStr)
cbar = fig.colorbar(mappable=s, ax=ax)
cbar.set_ticks(range(0, 6))
cbar.set_ticklabels(labels)
plt.show()
fig.savefig('c:/sven/test.png',    dpi=100)

plt.close(fig)

# test objective function
corFac= .7
tupleOfConstraints = (1.325e+6 * corFac, 0.226601e+11 * corFac, 3.32955027e+7 * corFac, 3.4772044e+6)

fpuGenome = {'newFPU': fpuOrder, 'PFT': aFPUGenome}

fpuGenome = pd.DataFrame(fpuGenome)
fpuGridWithGenome = pd.merge(fpuGrid, fpuGenome, on=('newFPU'))
# need to order accordingly to the objectives
# need to order accordingly to the objectives
fpuGridWithGenome.sort_values(by='ID', inplace=True)
fpuGridWithGenome.reset_index(drop=True, inplace=True)

luList = fpuGridWithGenome.PFT

fodPos = 1  # position of food provisioning in the lookup tables list

indicators, constraintViolations = lpf.returnValuesSumAndConstraints(
            luList=luList,
            listOfObjectiveLookupTables=listOfObjectiveLookupTables,
            foodMinThreshold=foodMinThreshold,
            foodPos=fodPos,
            tupleOfConstraints=tupleOfConstraints,
            maximumObjectives=maximumObjectives,
            minimumObjectives=minimumObjectives)