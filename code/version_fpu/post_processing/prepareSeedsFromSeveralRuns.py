import sys
import pandas as pd

fpuCodePostPath = sys.path[0]
fpuCodePath = fpuCodePostPath + "/../"
basePath = fpuCodePath + "/../../"
gaBasePath = fpuCodePath + "/../ecspy-0.4/"

sys.path.append(gaBasePath)
sys.path.append(fpuCodePath)
sys.path.append(fpuCodePostPath)
sys.path.append(fpuCodePath)

outPath= basePath + 'results/fpu/genome4cellLevel/'
inPathData = basePath + "data/"
path2Archive = basePath + 'results/fpu/'


import postProcessingFunctions_fpu as postprocessing

fpuGridFn= basePath + "data/fpu2biomeGrid.csv"
biomeFn= inPathData + "biomes_lai_tslice2008-2017_rcp26_INT_4gcms_0p4_FINAL_ATbiomes_9426.txt"
fpuGenomeFn=   basePath + "results/biomes/fpuGenomes/fpuGenomes_26_2033-2042/fpuGenome_1.csv"

# constraints
constraints = pd.read_csv(inPathData + "constraints.csv", delim_whitespace=True )
#constraints
# adjust units!
# add adjustment factor
adjFact= 0.95
constraints['carbon']= constraints['carbon'] * 10**3
constraints['food']= constraints['food'] * 10**9
constraints['water']= constraints['water'] * 10**3
constraints['forage']= constraints['forage'] * 10**3

constraintsAdj= constraints
constraintsAdj['carbon']= constraintsAdj['carbon']*adjFact
constraintsAdj['food']= constraintsAdj['food']*adjFact
constraintsAdj['water']= constraintsAdj['water']*adjFact
constraintsAdj['forage']= constraintsAdj['forage']*adjFact

foodProv = pd.read_csv(inPathData + 'ESquant_final_v2/ESquant_rcp26_2033-2042_INT_9426_v8b_4gcms_CStor_cor4ProtArea.txt', delim_whitespace=True)
#foodProv
landuseDf = pd.DataFrame({'Lon': foodProv['Lon'], 'Lat': foodProv['Lat']})

fpuGridFn= basePath + "data/fpu2biomeGrid.csv"
biomeFn= inPathData + "biomes_lai_tslice2008-2017_rcp26_INT_4gcms_0p4_FINAL_ATbiomes_9426.txt"
fpuGenomeFn= basePath + "results/biomes/fpuGenomes/fpuGenomes_26_2033-2042/fpuGenome_1.csv"

fpuGrid = pd.read_csv(fpuGridFn, sep=";")
refFPUGenome= pd.read_csv(fpuGenomeFn, sep=";")
refFPU= refFPUGenome.newFPU

runList= [5,7,8]
rcps= [26, 60]
timePeriods= [2042, 2090]


for rcp in rcps:
    for timePeriod in timePeriods:
        print(str(rcp) + ", " + str(timePeriod) + "," + str(run))
        theConst = constraints[(constraints['rcp'] == int(rcp)) & (constraints['timePeriod_short'] == timePeriod)]
        theConstAdj=  constraintsAdj[(constraintsAdj['rcp'] == int(rcp)) & (constraintsAdj['timePeriod_short'] == timePeriod)]
        combinedPop = postprocessing.combinePops(rcp=rcp, timePeriod=timePeriod, runList=runList, path=path2Archive)
        csFitBetter = postprocessing.subsetPop(pop=combinedPop, constraints=theConstAdj)
        csFitPareto = postprocessing.getParetoSolutions(csFitBetter)
        # create cell level representation from fpu level representation
        cellLevelRepresentation= postprocessing.createFpuGenomeTable4Pop(fpuLevelPop=csFitPareto['pop'], fpuGrid=fpuGrid, refFPU=refFPU,
                                                                         outPath=outPath, rcp=rcp, timePeriod=timePeriod)

        postprocessing.plotFront(csFitPareto, path=outPath, rcp=rcp, timePeriod=timePeriod, constraints=theConst,
                                 title="FPU")
        # create summary table
        summaryTable = postprocessing.createSummaryTable(landuseDf=cellLevelRepresentation, outPath=outPath, rcp=rcp,
                                          timePeriodShort=timePeriod)
        # plot maps
        postprocessing.plotMaps(landuseDf=cellLevelRepresentation, summaryTable=summaryTable, outPath=outPath, rcp=rcp,
                 timePeriodShort=timePeriod, title="FPU")
