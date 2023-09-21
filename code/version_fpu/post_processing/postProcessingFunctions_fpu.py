###################################################
## function to combine different runs
## for the same rcp, time-horizon combination
## export them to pickle and csv
## calculate the percentage of solution for which
## a cell become assigned a specific land use
##
## author: sven.lautenbach@uni-heidelberg.de
## date: 22.11.2019
###################################################
import os
import pickle
import matplotlib

matplotlib.use('PDF')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.backends.backend_pdf import PdfPages
import seaborn as sns
import fnmatch
import re
import itertools
import paretoSort
import pandas as pd


def assure_path_exists(path):
    dir = os.path.dirname(path)
    if not os.path.exists(dir):
        os.makedirs(dir)


# get pop file
# if optimization finished this will be nsGA_1.pkl
# if not present, take the latest intermediate
# archive file

def getPopFile(folderName):
    listOfFiles = os.listdir(folderName)
    finalFile = fnmatch.filter(listOfFiles, "nsGA_1.pkl")
    if (len(finalFile) == 0):
        finalFile = fnmatch.filter(listOfFiles, "interm_pop_*.pkl")
        if (len(finalFile) == 0):
            print("No pop archive file present in " + folderName)
            return (-1)
        if (len(finalFile) > 1):
            # get largest number
            theIDs = []
            for aFn in finalFile:
                theIDs.append(int(re.findall(r'\d+', aFn)[0]))
            theID = max(theIDs)
            finalFile = fnmatch.filter(listOfFiles, "interm_pop_*" + str(theID) + ".pkl")

    # print("Files used in " + folderName + ": " + str(finalFile) )
    return (finalFile[0])


def getPop(rcp, timePeriod, run, path):
    folderName = path + "run_" + str(run) + "_" + str(rcp) + "_" + str(timePeriod) + "/"
    fn = folderName + getPopFile(folderName)
    picklefile = open(fn, 'r')
    # final_pop = pd.read_pickle(picklefile)
    final_pop = pickle.load(picklefile)
    picklefile.close()
    return (final_pop)


def combinePops(rcp, timePeriod, runList, path):
    combinedPop = []
    for aRun in runList:
        thePop = getPop(rcp, timePeriod, aRun, path=path)
        combinedPop.extend(thePop[:])
    return combinedPop


def subsetPop(pop, constraints):
    res_fitnessC = []
    res_fitnessF = []
    res_fitnessW = []

    for ind in pop:
        # counter += 1
        cs = ind.candidate
        # save genome to file
        # numpy.savetxt(outPath  + "finalPop_" + str(counter) + ".txt", cs, fmt='%i',  delimiter=';')
        res_fitnessC.append(ind.fitness[0])
        res_fitnessF.append(ind.fitness[1])
        res_fitnessW.append(ind.fitness[2])

    # create Dataframe
    fitnessDF = pd.DataFrame(
        {'carbonStorage': res_fitnessC, 'foodProvisioning': res_fitnessF, 'waterSupply': res_fitnessW})
    # fitnessDF.shape
    # identify duplicates
    whoIsDuplicated = fitnessDF.duplicated(subset=('carbonStorage', 'foodProvisioning', 'waterSupply'), keep='first')
    fitnessDFNoDup = fitnessDF[~whoIsDuplicated]

    selectedSubsetBool = (fitnessDFNoDup['carbonStorage'] >= float(constraints['carbon'])) & (
            fitnessDFNoDup['foodProvisioning'] >= float(constraints['food'])) & (
                                 fitnessDFNoDup['waterSupply'] >= float(constraints['water']))
    selectedSubsetBool.sum()

    selectedSubset = fitnessDFNoDup[selectedSubsetBool].copy()
    selectedSubset.describe()

    # get candidates based on the index of the selected fitness rows
    idx = selectedSubset.index.get_values().tolist()
    # type(idx)
    # idx[0]
    selected_solutions = [pop[i] for i in idx]

    return {'pop': selected_solutions, 'fitness': selectedSubset}


def getParetoSolutions(popFitness):
    # get the pareto optimal solutions from the set of solutions which are better than cuirrent land use/reference situation
    pop = popFitness['pop']
    fitness = popFitness['fitness']
    # combine into one ndarray
    fitness_np = fitness.values

    # get index of pareto optimal solutions
    isParetoArray = paretoSort.is_pareto(fitness_np, maximise=True)
    fitnessPareto = fitness[isParetoArray].copy()
    # subset selected population to those pareto optimal
    popPareto = list(itertools.compress(pop, isParetoArray))

    return {'pop': popPareto, 'fitness': fitnessPareto}


def createFpuGenomeTable4Pop(fpuLevelPop, fpuGrid, refFPU, outPath, rcp, timePeriod):
    i = 0
    for ind in fpuLevelPop:
        colName = 'sol_' + str(i)
        cs = ind.candidate
        fpuGenome = {'newFPU': refFPU['newFPU'],
                     'PFT': cs}
        fpuGenome = pd.DataFrame(fpuGenome)
        fpuGridWithGenome = pd.merge(fpuGrid, fpuGenome, on=('newFPU'))
        # need to order accordingly to the objectives
        fpuGridWithGenome.sort_values(by='ID', inplace=True)
        fpuGridWithGenome.reset_index(drop=True, inplace=True)

        luList = fpuGridWithGenome.PFT
        if i == 0:
            cellLevelRepresentation = pd.DataFrame({'Lon': fpuGridWithGenome['Lon'], 'Lat': fpuGridWithGenome['Lat'],
                                                    'newFPU': fpuGridWithGenome['newFPU'], colName: luList})
        cellLevelRepresentation = cellLevelRepresentation.assign(**{colName: luList})
        i = i + 1
    # serialize and export as csv
    outFn = outPath + "fpuGenomes4CellLevel_" + str(rcp) + "_" + str(timePeriod)
    cellLevelRepresentation.drop(labels=['Lon', 'Lat', 'newFPU'], axis=1).to_pickle(path=outFn + ".pkl")
    cellLevelRepresentation.drop(labels=['Lon', 'Lat', 'newFPU'], axis=1).to_csv(path_or_buf=outFn + ".csv", sep=";",
                                                                                 index=False)
    return (cellLevelRepresentation)


def plotFront(popFitness, path, rcp, timePeriod, constraints, title):
    assure_path_exists(path)
    fitness = popFitness['fitness']

    fig = plt.figure()
    fig.set_size_inches(18.5, 10.5)
    ax = fig.add_subplot(111, projection='3d')
    ax.scatter(fitness['carbonStorage'], fitness['foodProvisioning'],
               fitness['waterSupply'], marker='o', c='b')
    ax.scatter(constraints['carbon'], constraints['food'], constraints['water'], marker='o', c='r')

    ax.set_xlabel('carbon storage')
    ax.set_ylabel('food provisioning')
    ax.set_zlabel('water supply')
    ax.set_title(title + ", " + str(rcp) + ", " + str(timePeriod))
    fig.savefig(path + 'frontierPlusCurrentSelection_fpu_' + str(rcp) + '_' + str(timePeriod) + '.pdf', dpi=100)
    # plt.show()


def createSummaryTable(landuseDf, outPath, rcp, timePeriodShort):
    summaryTable = pd.DataFrame({'Lat': landuseDf['Lat'], 'Lon': landuseDf['Lon']})

    summaryTable['nPnv'] = landuseDf.isin(['0']).sum(axis=1)
    summaryTable['nC3cer'] = landuseDf.isin(['1']).sum(axis=1)
    summaryTable['nC3cerI'] = landuseDf.isin(['2']).sum(axis=1)
    summaryTable['nC3crps'] = landuseDf.isin(['3']).sum(axis=1)
    summaryTable['nC3crpsI'] = landuseDf.isin(['4']).sum(axis=1)
    summaryTable['nC4'] = landuseDf.isin(['5']).sum(axis=1)
    summaryTable['nC4I'] = landuseDf.isin(['6']).sum(axis=1)
    summaryTable['nRice'] = landuseDf.isin(['7']).sum(axis=1)
    summaryTable['nRiceI'] = landuseDf.isin(['8']).sum(axis=1)
    summaryTable['nPasture'] = landuseDf.isin(['9']).sum(axis=1)
    summaryTable['nCropIrrigated'] = summaryTable['nC3cerI'] + summaryTable['nC3crpsI'] + summaryTable['nC4I'] + + \
    summaryTable['nRiceI']
    summaryTable['nCropRainfeed'] = summaryTable['nC3cer'] + summaryTable['nC3crps'] + summaryTable['nC4'] + + \
    summaryTable['nRice']
    summaryTable['nCrop'] = summaryTable['nCropIrrigated'] + summaryTable['nCropRainfeed']
    # calculate relative values
    n = float(landuseDf.shape[1])
    summaryTable['percPnv'] = summaryTable['nPnv'] / n
    summaryTable['percC3cer'] = summaryTable['nC3cer'] / n
    summaryTable['percC3cerI'] = summaryTable['nC3cerI'] / n
    summaryTable['perC3crps'] = summaryTable['nC3crps'] / n
    summaryTable['perC3crpsI'] = summaryTable['nC3crpsI'] / n
    summaryTable['percC4'] = summaryTable['nC4'] / n
    summaryTable['percC4I'] = summaryTable['nC4I'] / n
    summaryTable['percRice'] = summaryTable['nRice'] / n
    summaryTable['percRiceI'] = summaryTable['nRiceI'] / n
    summaryTable['percPasture'] = summaryTable['nPasture'] / n
    summaryTable['percCrop'] = summaryTable['nCrop'] / n
    summaryTable['percCropIrrigated'] = summaryTable['nCropIrrigated'] / n
    summaryTable['percCropRainfeed'] = summaryTable['nCropRainfeed'] / n

    fn = outPath + "summaryTableAssignedLU_fpu_" + str(rcp) + "_" + str(timePeriodShort) + ".csv"
    summaryTable.to_csv(path_or_buf=fn, sep=",", index=False)
    return (summaryTable)


def plotMaps(landuseDf, summaryTable, outPath, rcp, timePeriodShort, title=""):
    # title= title + ", " + str(rcp) + ", " + str(timePeriodShort)
    # in total n solution
    n = landuseDf.shape[1]
    oneFifth = round(n / 5)
    # sns.set_palette(sns.cubehelix_palette(5, start=.5, rot=-.75))

    # percent pnv assigned across solutions
    theBins = [-1, 0, oneFifth, 2 * oneFifth, 3 * oneFifth, 4 * oneFifth, n]
    theLabels = ["0%", "0-20%", "20-40%", "40-60%", "60-80%", "> 80%"]
    thePalette = sns.cubehelix_palette(6)
    theKwsDict = {"s": 10, "edgecolor": "none"}
    theHeight = 10
    theAspect = 1.85
    # Create the PdfPages object to which we will save the pages:
    # The with statement makes sure that the PdfPages object is closed properly at
    # the end of the block, even if an Exception occurs.
    fn = outPath + "summaryMaps_" + str(rcp) + "_" + str(timePeriodShort) + ".pdf"
    with PdfPages(fn) as pdf:
        summaryTable['Labels_nPnv'] = pd.cut(summaryTable['nPnv'], bins=theBins, labels=theLabels)
        thePlot = sns.pairplot(x_vars=['Lon'], y_vars=['Lat'], data=summaryTable, hue="Labels_nPnv",
                               plot_kws=theKwsDict,
                               palette=thePalette, height=theHeight, aspect=theAspect)
        thePlot.fig.suptitle(title + "PNV")
        pdf.savefig(thePlot.fig)  # saves the current figure into a pdf page

        # percent c3 cereals assigned across solutions
        summaryTable['Labels_nC3cer'] = pd.cut(summaryTable['nC3cer'], bins=theBins, labels=theLabels)
        thePlot = sns.pairplot(x_vars=['Lon'], y_vars=['Lat'], data=summaryTable, hue="Labels_nC3cer",
                               plot_kws=theKwsDict,
                               palette=thePalette, height=theHeight, aspect=theAspect)
        thePlot.fig.suptitle(title + " C3 non cereals, rainfeed")
        pdf.savefig(thePlot.fig)

        # percent c3 non cereals assigned across solutions
        summaryTable['Labels_nC3crps'] = pd.cut(summaryTable['nC3crps'], bins=theBins, labels=theLabels)
        thePlot = sns.pairplot(x_vars=['Lon'], y_vars=['Lat'], data=summaryTable, hue="Labels_nC3crps",
                               plot_kws=theKwsDict,
                               palette=thePalette, height=theHeight, aspect=theAspect)
        thePlot.fig.suptitle(title + " C3 cereals, rainfeed")
        pdf.savefig(thePlot.fig)

        # percent c4 crops assigned across solutions
        summaryTable['Labels_nC4'] = pd.cut(summaryTable['nC4'], bins=theBins, labels=theLabels)
        thePlot = sns.pairplot(x_vars=['Lon'], y_vars=['Lat'], data=summaryTable, hue="Labels_nC4", plot_kws=theKwsDict,
                               palette=thePalette, height=theHeight, aspect=theAspect)
        thePlot.fig.suptitle(title + " C4, rainfeed")
        pdf.savefig(thePlot.fig)

        # percentpasture assigned across solutions
        summaryTable['Labels_nPasture'] = pd.cut(summaryTable['nPasture'], bins=theBins, labels=theLabels)
        thePlot = sns.pairplot(x_vars=['Lon'], y_vars=['Lat'], data=summaryTable, hue="Labels_nPasture",
                               plot_kws=theKwsDict,
                               palette=thePalette, height=theHeight, aspect=theAspect)
        thePlot.fig.suptitle(title + " Pasture")
        pdf.savefig(thePlot.fig)

        # percent rice assigned across solutions
        summaryTable['Labels_nRice'] = pd.cut(summaryTable['nRice'], bins=theBins, labels=theLabels)
        thePlot = sns.pairplot(x_vars=['Lon'], y_vars=['Lat'], data=summaryTable, hue="Labels_nRice",
                               plot_kws=theKwsDict,
                               palette=thePalette, height=theHeight, aspect=theAspect)
        thePlot.fig.suptitle(title + " Rice, rainfeed")
        pdf.savefig(thePlot.fig)

        # percent rice assigned across solutions
        summaryTable['Labels_nRiceI'] = pd.cut(summaryTable['nRiceI'], bins=theBins, labels=theLabels)
        thePlot = sns.pairplot(x_vars=['Lon'], y_vars=['Lat'], data=summaryTable, hue="Labels_nRiceI",
                               plot_kws=theKwsDict,
                               palette=thePalette, height=theHeight, aspect=theAspect)
        thePlot.fig.suptitle(title + " Rice, irrigated")
        pdf.savefig(thePlot.fig)

        # percent c4 crops assigned across solutions
        summaryTable['Labels_nC4I'] = pd.cut(summaryTable['nC4I'], bins=theBins, labels=theLabels)
        thePlot = sns.pairplot(x_vars=['Lon'], y_vars=['Lat'], data=summaryTable, hue="Labels_nC4I",
                               plot_kws={"s": 8, "edgecolor": "none"}, palette=thePalette, height=theHeight,
                               aspect=theAspect)
        thePlot.fig.suptitle(title + " C4, irrigated")
        pdf.savefig(thePlot.fig)

        # percent c3 non cereals assigned across solutions
        summaryTable['Labels_nC3crpsI'] = pd.cut(summaryTable['nC3crpsI'], bins=theBins, labels=theLabels)
        thePlot = sns.pairplot(x_vars=['Lon'], y_vars=['Lat'], data=summaryTable, hue="Labels_nC3crpsI",
                               plot_kws={"s": 8, "edgecolor": "none"}, palette=thePalette, height=theHeight,
                               aspect=theAspect)
        thePlot.fig.suptitle(title + " C3 non cereals, irrigated")
        pdf.savefig(thePlot.fig)

        # percent c3  cereals assigned across solutions
        summaryTable['Labels_nC3cerI'] = pd.cut(summaryTable['nC3cerI'], bins=theBins, labels=theLabels)
        thePlot = sns.pairplot(x_vars=['Lon'], y_vars=['Lat'], data=summaryTable, hue="Labels_nC3cerI",
                               plot_kws={"s": 8, "edgecolor": "none"}, palette=thePalette, height=theHeight,
                               aspect=theAspect)
        thePlot.fig.suptitle(title + " C3 cereals, irrigated")
        pdf.savefig(thePlot.fig)


def seralizeOutput(cellLevelRepresentation, fitness, path, rcp, timePeriod, landuseDf,
                   basefnFitness="betterThanCurrent_fpu_"):
    ####################################
    # serialize fitness
    assure_path_exists(path)
    fn = path + basefnFitness + str(rcp) + "_" + str(timePeriod)
    picklefile = open(fn + ".pkl", 'wb')
    pickle.dump(fitness, picklefile)
    picklefile.close()

    # %% store in hd5 format
    # myStore = pd.HDFStore( fn + '.h5')
    # myStore['fitness'] = fitness
    # myStore
    # myStore.close()

    # write to csv for use in R
    fitness.to_csv(path_or_buf=fn + '.csv')

    ####################################
    # serialize the candidate solutions!
    picklefile = open(fn + "candidates.pkl", 'wb')
    pickle.dump(cellLevelRepresentation, picklefile)
    picklefile.close()

    # write to csv for use in R
    for i, ind in zip(xrange(0, len(cellLevelRepresentation)), cellLevelRepresentation):
        colName = 'sol_' + str(i)
        # print(colName)
        landuseDf = landuseDf.assign(**{colName: ind})

    # landuseDf.describe()

    landuseDf.to_csv(path_or_buf=fn + "FrontGenomeTable.csv", sep=",", index=False)

    return (landuseDf)


def createTablefromFronts(rcp, timePeriodShort, runList, constraints, landuseDf, basePath, outPath, fpuGrid, refFPU,
                          title=""):
    theConst = constraints[(constraints['rcp'] == int(rcp)) & (constraints['timePeriod_short'] == timePeriodShort)]
    combinedPop = combinePops(rcp=rcp, timePeriod=timePeriodShort, runList=runList, path=basePath)
    csFitBetter = subsetPop(pop=combinedPop, constraints=theConst)
    csFitPareto = getParetoSolutions(csFitBetter)

    csFitPareto['fitness'].to_csv(
        path_or_buf='{0}better_than_current_fpu_{1}_{2}.csv'.format(outPath, str(rcp), str(timePeriodShort)))

    # get the cell level representation for the fpu genomes
    cellLevelRepresentation = createFpuGenomeTable4Pop(fpuLevelPop=csFitPareto['pop'], fpuGrid=fpuGrid, refFPU=refFPU,
                                                       outPath=outPath, rcp=rcp, timePeriod=timePeriodShort)
    # landuseDf= seralizeOutput(cellLevelRepresentation, fitness=csFitBetter['fitness'], path=outPath, rcp=rcp, timePeriod=timePeriodShort, landuseDf= landuseDf)

    plotFront(csFitPareto, path=outPath, rcp=rcp, timePeriod=timePeriodShort, constraints=theConst, title=title)
    # create summary table
    summaryTable = createSummaryTable(landuseDf=cellLevelRepresentation, outPath=outPath, rcp=rcp,
                                      timePeriodShort=timePeriodShort)
    # plot maps
    plotMaps(landuseDf=cellLevelRepresentation, summaryTable=summaryTable, outPath=outPath, rcp=rcp,
             timePeriodShort=timePeriodShort, title=title)
    # return(csFitPareto)
