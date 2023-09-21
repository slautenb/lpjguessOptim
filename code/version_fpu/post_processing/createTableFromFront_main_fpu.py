####################
# Create one big genome table from front archive
# here for the fpu level archive
# this is the main function calling the respective function
# version for updated data december 2019
###################

import pandas as pd
import sys
# import h5py

fpuCodePostPath = sys.path[0]
fpuCodePath = fpuCodePostPath + "/../"
basePath = fpuCodePath + "/../../"
gaBasePath = fpuCodePath + "/../ecspy-0.4/"

sys.path.append(gaBasePath)
sys.path.append(fpuCodePath)
sys.path.append(fpuCodePath)
sys.path.append(fpuCodePostPath)

import postProcessingFunctions_fpu as postprocessing

outPath = basePath + 'results/fpu/genomeCellLevel/combinedRuns/genomeCellLevel/'
inPathData = basePath + "data/"
path2ArchiveBase = basePath + 'results/fpu/'

# constraints
constraints = pd.read_csv(inPathData + "constraints.csv", delim_whitespace=True)
# constraints
# adjust units!
# and apply correction factor
corFac = 1
constraints['carbon'] = constraints['carbon'] * 10 ** 3 * corFac
constraints['food'] = constraints['food'] * 10 ** 9 * corFac
constraints['water'] = constraints['water'] * 10 ** 3 * corFac
constraints['forage'] = constraints['forage'] * 10 ** 3 * corFac

foodProv = pd.read_csv(inPathData + 'ESquant/ESquant_rcp26_2033-2042_INT_9412_v8b_4gcms_CStor_cor4ProtArea.txt',
                       delim_whitespace=True)
# foodProv
landuseDf = pd.DataFrame({'Lon': foodProv['Lon'], 'Lat': foodProv['Lat']})

fpuGridFn = basePath + "data/fpu2biomeGrid.csv"
fpuGenomeFn = basePath + "results/biomes/fpuGenomes/fpuGenomes_26_2033-2042/fpuGenome_1.csv"
# fpuGrid is used to bring the fpu level genome to the cell level
fpuGrid = pd.read_csv(fpuGridFn, sep=";")
aFpuGenome = pd.read_csv(fpuGenomeFn, sep=";")

runList = [3]
rcpList = [26]  # [26,60]
timePeriodList = [2042]  # [2042, 2090]

resList = []
for aRcp in rcpList:
    for aTimePeriod in timePeriodList:
        title = "Rcp: " + str(aRcp) + ", " + str(aTimePeriod)
        print("Processing now " + title)
        # call function
        res = postprocessing.createTablefromFronts(rcp=aRcp, timePeriodShort=aTimePeriod, runList=runList,
                                                   constraints=constraints,
                                                   landuseDf=landuseDf, basePath=path2ArchiveBase, outPath=outPath,
                                                   fpuGrid=fpuGrid, refFPU=aFpuGenome, title=title)
        resList.append(res)
