####################
# Create one big genome table from front archive
# here for the FPU level archive
###################

import numpy as np
import pandas as pd
import sys, os
#import h5py
import pickle

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns

codePath = "c:/sven/Nutzer/development/python/GA/ecspy-0.4/"
sys.path.append(codePath)
basePath = "c:/sven/Nutzer/projekte_neu/ess_tradeoffs/LPJ_pareto/"
sys.path.append(basePath + "code/version_cell_updated_july_2019/")
outPath= basePath + 'results/version_july_2019/fpu/genome4cellLevel/'
inPathData = basePath + "data/update_july_2019/"
path2Archive = basePath + 'results/version_july_2019/fpu/run_1_26_2042/'

import ec_arrayGenome_4ConstraintsAndPredefinedGenomes as ec_arrayGenome #ec_arrayGenome

#%% read from picklefile

picklefile = open(path2Archive + "nsGA_1.pkl")
fpuLevelPop = pickle.load(picklefile)
picklefile.close()

# fpuGrid is used to bring the fpu level genome to the cell level biome
fpuGrid = pd.read_csv(basePath + "data/update_july_2019/fpu2biomeGrid.csv", sep=";")
fn= basePath + "results/version_july_2019/biomes/fpuGenomes/fpuGenomes_26_2033-2042/fpuGenome_1.csv"
aFpuGenome = pd.read_csv(fn, sep=";")

i= 0
for ind in fpuLevelPop:
    theGenome = ind.candidate
    if i==0:
        fpuTable = pd.DataFrame({'FPU': aFpuGenome['newFPU']})
    colName = 'sol_' + str(i)
    # print(colName)
    fpuTable = fpuTable.assign(**{colName: ind.candidate})
    i = i + 1

fpuTable.shape

# merge with cell level information

fpuTableGrid = pd.merge(fpuGrid, fpuTable, left_on=('newFPU'), right_on= ('FPU'))

# need to order accordingly to the objectives
fpuTableGrid.sort_values(by='ID', inplace=True)
fpuTableGrid.reset_index(drop=True, inplace=True)

# join with biom info
fn= inPathData + "biomes_lai_tslice2008-2017_rcp26_INT_4gcms_0p4_FINAL_ATbiomes_9426.txt"
biome= pd.read_csv(fn, delim_whitespace= True)

fpuTableGridBiome= pd.merge(fpuTableGrid, biome, left_on=['Lon', 'Lat'], right_on= ['Lon', 'Lat'])

fpuTableGridBiome.to_csv(path_or_buf= outPath + "fpuFrontGenomeTable.csv", sep=",", index= False)

# serialise genome for import in cell level optimization
#list(fpuTableGrid.columns)
theGenomesFPU= fpuTableGrid.drop(labels=['Lon', 'Lat', 'ATbiome', 'FPU_code', 'newFPU', 'ID', 'LONLAT', 'FPU'], axis=1)

theGenomesFPU.head()

outFn = outPath + "fpuGenomes4CellLevel_run_1_26_2042.pkl"
theGenomesFPU.to_pickle(path=outFn )

# count occurence of different pfts/cfts for each cell
assignedLU = pd.DataFrame({'Lat': fpuTableGrid['Lat'], 'Lon': fpuTableGrid['Lon'], 'newFPU': fpuTableGrid['newFPU']})

assignedLU['nPnv']= theGenomesFPU.isin(['0']).sum(axis=1)
assignedLU['nC3cer']= theGenomesFPU.isin(['1']).sum(axis=1)
assignedLU['nC3cerI']= theGenomesFPU.isin(['2']).sum(axis=1)
assignedLU['nC3crps']= theGenomesFPU.isin(['3']).sum(axis=1)
assignedLU['nC3crpsI']= theGenomesFPU.isin(['4']).sum(axis=1)
assignedLU['nC4']= theGenomesFPU.isin(['5']).sum(axis=1)
assignedLU['nC4I']= theGenomesFPU.isin(['6']).sum(axis=1)
assignedLU['nRice']= theGenomesFPU.isin(['7']).sum(axis=1)
assignedLU['nRiceI']= theGenomesFPU.isin(['8']).sum(axis=1)
assignedLU['nPasture']= theGenomesFPU.isin(['9']).sum(axis=1)

assignedLU.describe()

assignedLU.to_csv(path_or_buf= outPath + "fpuFrontAssignedLandUse.csv", sep=",", index=False)

#in total 1242 solution

sns.set_palette(sns.cubehelix_palette(5, start=.5, rot=-.75))

# percent pnv assigned across solutions
assignedLU['Labels_nPnv'] = pd.cut(assignedLU['nPnv'], bins=[-1,0, 360, 720, 1080, 1242], labels= [0,1,2,3,4])
sns.pairplot(x_vars=['Lon'], y_vars=['Lat'], data= assignedLU, hue="Labels_nPnv", plot_kws= {"s": 8})

# percent c3 cereals assigned across solutions
assignedLU['Labels_nC3cer'] = pd.cut(assignedLU['nC3cer'], bins=[-1,0, 360, 720, 1080, 1242], labels= [0,1,2,3,4])
sns.pairplot(x_vars=['Lon'], y_vars=['Lat'], data= assignedLU, hue="Labels_nC3cer", plot_kws= {"s": 8})

# percent c3 non cereals assigned across solutions
assignedLU['Labels_nC3crps'] = pd.cut(assignedLU['nC3crps'], bins=[-1,0, 360, 720, 1080, 1242], labels= [0,1,2,3,4])
sns.pairplot(x_vars=['Lon'], y_vars=['Lat'], data= assignedLU, hue="Labels_nC3crps", plot_kws= {"s": 8})

# percent c4 crops assigned across solutions
assignedLU['Labels_nC4'] = pd.cut(assignedLU['nC4'], bins=[-1,0, 360, 720, 1080, 1242], labels= [0,1,2,3,4])
sns.pairplot(x_vars=['Lon'], y_vars=['Lat'], data= assignedLU, hue="Labels_nC4", plot_kws= {"s": 8})

# percentpasture assigned across solutions
assignedLU['Labels_nPasture'] = pd.cut(assignedLU['nPasture'], bins=[-1,0, 360, 720, 1080, 1242], labels= [0,1,2,3,4])
sns.pairplot(x_vars=['Lon'], y_vars=['Lat'], data= assignedLU, hue="Labels_nPasture", plot_kws= {"s": 8})