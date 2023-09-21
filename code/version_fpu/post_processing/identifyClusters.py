#################
## calculate the distance between solutions
## and create k-means clustzer based on the distance
################

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

basePath = "c:/sven/Nutzer/projekte_neu/ess_tradeoffs/LPJ_pareto/"
path2GenomeTable= basePath + 'results/version_july_2019/fpu/genome4cellLevel/'

runs= [5]
rcps= [26, 60]
timePeriods= [2042, 2090]

for run in runs:
    for rcp in rcps:
        for timePeriod in timePeriods:
            runRcpTime = "run_" + str(run) + "_" + str(rcp) + "_" + str(timePeriod)
            filesPrefix = runRcpTime + "_"
            fn = path2GenomeTable + runRcpTime + "/" + runRcpTime + "fpuFrontGenomeTable.csv"
            fpuTableGridBiome = pd.read_csv(fn, sep=",")
            # drop unnecessary columns
            theGenomesFPU = fpuTableGridBiome.drop(
                labels=['Lon', 'Lat', 'ATbiome_x', 'FPU_code', 'newFPU', 'ID', 'LONLAT', 'FPU', 'ATbiome_y'], axis=1)
            n = theGenomesFPU.shape[1]

            # create distance matrix
            distM = np.zeros((n, n))

            for i in range(0, n):
                print str(i)
                for j in range(i + 1, n):
                    # diagonal has already distance 0 and the matrix is symetric
                    # i.e. start j at i+1
                    cellsDifferent = theGenomesFPU.apply(lambda x: 0 if x.iloc[i] == x.iloc[j] else 1, axis=1)
                    distM[i, j] = len(cellsDifferent[cellsDifferent == 1].index)
                    distM[j, i] = distM[i, j]  # symetric matrix

            fn = path2GenomeTable + runRcpTime + "/" + runRcpTime + "distanceMatrix.csv"
            np.savetxt(fn, distM)
            fn = path2GenomeTable + runRcpTime + "/" + runRcpTime + "distanceMatrix.npz"
            np.savez_compressed(fn, distM=distM)

plt.matshow(distM)
plt.colorbar()

theBins= [0,1,10,100,1000, 10**4]
np.histogram(distM, bins='auto')
plt.hist(distM[np.triu_indices(n)], bins='auto')