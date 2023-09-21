# -*- coding: utf-8 -*-
"""
Created on Fri May 15 10:48:11 2015

@author: Lautenbach
"""

from random import Random
from time import time
import numpy, numpy.random
from itertools import izip # for crossover operator
import pandas as pd 
#from ecspy import emo
#import emo_epsilon
import emo_constraint
import pandas
import six
import matplotlib.pyplot as plt

def generate_selected_candidates(args):
    inStartGenomesFolder = args['inStartGenomesFolder']
    nStartGenomesCreated = args['nStartGenomesCreated']    
    
    fn = inStartGenomesFolder + "fpuGenome_" + str(nStartGenomesCreated) + ".csv"
    try:
        fpuGenome = pandas.read_csv(fn, sep=";" )
        genom_array = fpuGenome.PFT    # if the field is named differently this needs to be changed...needs a better solution
    except (OSError, IOError) as e:
        print "generator_predef: \n" + e.message
        raise e
    
    return genom_array

def generate_initial_candidate(random, args):
    # create random list
    random.seed()
    nCells =  args['nCells']
    rangeValues = args['rangeValues']
    # for consistency I am returning a panda.Series
    # that shoule behave like a numpy.ndarray
    genom_array= pd.Series(numpy.random.choice(rangeValues, size=nCells) )

    return genom_array

#evaluator for fpu level optimization
def evaluate_lpj_fpuLevel(candidates, args):
    listOfObjectiveLookupTables = args['listOfObjectiveLookupTables']
    # listOfEpsValues= args['listOfEpsValues']
    fpuOrder = args['fpuOrder']
    fpuGrid = args['fpuGrid']
    foodMinThreshold = args['foodMinThreshold']
    fodPos = args['fodPos']
    tupleOfConstraints = args['tupleOfConstraints']
    maximumObjectives = args['maximumObjectives']
    minimumObjectives = args['minimumObjectives']

    fitness = []

    for cs in candidates:
        # print(cs.dtype)
        # combine with fpuOrder
        # all files should be in the same order!
        fpuGenome = {'newFPU': fpuOrder,
                     'PFT': cs}

        fpuGenome = pd.DataFrame(fpuGenome)
        fpuGridWithGenome = pd.merge(fpuGrid, fpuGenome, on=('newFPU'))
        # need to order accordingly to the objectives
        # need to order accordingly to the objectives
        fpuGridWithGenome.sort_values(by='ID', inplace=True)
        fpuGridWithGenome.reset_index(drop=True, inplace=True)

        luList = fpuGridWithGenome.PFT

        # call evaluation function
        indicators, constraintViolations = returnValuesSumAndConstraints(
            luList=luList,
            listOfObjectiveLookupTables=listOfObjectiveLookupTables,
            foodMinThreshold=foodMinThreshold,
            foodPos=fodPos,
            tupleOfConstraints=tupleOfConstraints,
            maximumObjectives=maximumObjectives,
            minimumObjectives=minimumObjectives)
        #print(str(indicators))
        #print(str(constraintViolations))
        #print(str(tupleOfConstraints))
        # currently all objectives need to be maximized
        # if not the ones to be minimized would have to be multiplied with -1
        #since forage is not considered as a objective, only add first 3 objectives
        fitness.append(emo_constraint.constraint_array_Pareto(values=[indicators[0], indicators[1], indicators[2]],
                                                              constraintViolations=constraintViolations))

        # new needed?
    return (fitness)

# evaluator for cell level optimization
def evaluate_lpj_cellLevel(candidates, args):
    listOfObjectiveLookupTables = args['listOfObjectiveLookupTables']
    #listOfEpsValues= args['listOfEpsValues']
    #fpuOrder = args['fpuOrder']
    #fpuGrid = args['fpuGrid']
    foodMinThreshold = args['foodMinThreshold']
    fodPos = args['fodPos']
    tupleOfConstraints = args['tupleOfConstraints']
    maximumObjectives = args['maximumObjectives']
    minimumObjectives = args['minimumObjectives']
    
    fitness = []
    
    for cs in candidates:
        indicators, constraintViolations = returnValuesSumAndConstraints( 
                luList=cs, 
                listOfObjectiveLookupTables= listOfObjectiveLookupTables,
                foodMinThreshold=foodMinThreshold,
                foodPos=fodPos,
                tupleOfConstraints=tupleOfConstraints,
                maximumObjectives=maximumObjectives,
                minimumObjectives= minimumObjectives)
        #print(indicators)
        # currently all objectives need to be maximized
        # if not the ones to be minimized would have to be multiplied with -1
        fitness.append(emo_constraint.constraint_array_Pareto(values= [indicators[0], indicators[1], indicators[2]], constraintViolations= constraintViolations) )
    
    # new needed?
    return(fitness)
# evaluator
def evaluate_lpj(candidates, args):
    
    listOfObjectiveLookupTables = args['listOfObjectiveLookupTables']
    listOfEpsValues= args['listOfEpsValues']
    fpuOrder = args['fpuOrder']
    fpuGrid = args['fpuGrid']
    foodMinThreshold = args['foodMinThreshold']
    fodPos = args['fodPos']
    
    fitness = []

    for cs in candidates:
        #print(cs.dtype)
        # combine with fpuOrder
        # all files should be in the same order!
        fpuGenome = {'newFPU' : fpuOrder,
             'PFT' : cs} 
    
        fpuGenome= pd.DataFrame(fpuGenome)
        fpuGridWithGenome = pd.merge(fpuGrid, fpuGenome, on=('newFPU'))
        # need to order accordingly to the objectives
        fpuGridWithGenome.sort_values(by='ID', inplace=True)      
        fpuGridWithGenome.reset_index(drop=True, inplace=True)
        
        luList = fpuGridWithGenome.PFT
        
               
        indicators = returnValuesSum( luList=luList, 
                                     listOfObjectiveLookupTables= listOfObjectiveLookupTables,
                                     foodMinThreshold=foodMinThreshold,
                                     foodPos=fodPos)
        #print(indicators)
        # are all values to be maximized or to be minimized?
        # now all maximized
        #fitness.append(emo.array_Pareto([-1*indicators[0], indicators[1], indicators[2]]))
        fitness.append(emo_epsilon.epsilon_array_Pareto(values= [indicators[0], indicators[1], indicators[2]], epsilons  = listOfEpsValues) )
        
        return(fitness)  

def returnValuesSum( luList, listOfObjectiveLookupTables, foodMinThreshold, offset=2, foodPos=1):
    # correct unfeasible solutions
    FoodProv = listOfObjectiveLookupTables[foodPos]
    luList= setUnfeasible2PNV(FoodProv, luList, foodMinThreshold)
    # add offset to list
    # this is needed since then we can use the PFT value directly to extract values from the numpy array
    # e.g. if PFT = PNV this is coded as 0 but since the first to columns are LAT and LON we need to get
    # the value from the position 2 (0+2)
    luList = numpy.add(luList,offset)
    results = []
    for lookupTable in listOfObjectiveLookupTables:
        #print((range(0,lookupTable.shape[0] )))
        #print(luList.shape)
        #valList = lookupTable.iloc[range(0, lookupTable.shape[0]), luList]
        #lookupTable = df_to_sarray(lookupTable.reset_index())
        valList = lookupTable[(range(0,lookupTable.shape[0] )), luList]
        #results.append(numpy.sum(valList)) 
        results.append(valList.sum())
    return(results) 

def returnValuesSumAndConstraints( luList, listOfObjectiveLookupTables,  tupleOfConstraints, maximumObjectives, minimumObjectives, foodMinThreshold, offset=2, foodPos=1):
    # to do constraint handling!!!
    
    # correct unfeasible solutions
    FoodProv = listOfObjectiveLookupTables[foodPos]
    # corre3cted 09.10.2018: offset was missing!!!!
    # this lead to lat and lon beeing included....
    # add offset to list
    # this is needed since then we can use the PFT value directly to extract values from the numpy array
    # e.g. if PFT = PNV this is coded as 0 but since the first to columns are LAT and LON we need to get
    # the value from the position 2 (0+2)
    luList = numpy.add(luList, offset) # first add the offset
    # setUnfeasibel2PNV puts 2 in if food below threshold
    #luList = setUnfeasible2PNV(FoodProv, luList, foodMinThreshold)
    # forrage has food provisioning of zero... therefore all pasture had been replaced by PNV which was not intended
    # I think we do not need this anayway since Ani has taken this already in account


   
    results = []
    constraintViolations = []
    for lookupTable, constraint, minObj, maxObj in zip(listOfObjectiveLookupTables, 
                                                       tupleOfConstraints,
                                                       maximumObjectives,
                                                       minimumObjectives):
        valList = lookupTable[(range(0, lookupTable.shape[0])), luList]
        sumValues = valList.sum()
        results.append(sumValues)
        theViolation = 0
        if(sumValues <= constraint):
            #print("Constraint violated")
            # get amount of constraint violation
            theViolation = constraint - sumValues
            # scale it to [0,1], changed 28.09.2018
            theViolation = theViolation / constraint
            
        constraintViolations.append(theViolation)
    return(results, constraintViolations) 

def getMaximumObjectiveValue(objectiveDf, startPos, endPos):
    #Get the maximum objective value
    # the columns which contain the crop types are hard coded right now
    theNames = objectiveDf.columns.values.tolist()
    theMax =objectiveDf[theNames[startPos:endPos]].max(axis=1).sum() 
    return(theMax)

def getMinimumObjectiveValue(objectiveDf, startPos, endPos):
    #Get the maximum objective value
    # the columns which contain the crop types are hard coded right now
    theNames = objectiveDf.columns.values.tolist()
    theMin = objectiveDf[theNames[startPos:endPos]].min(axis=1).sum() 
    return(theMin)

def getMaximumObjectiveDataFrame(objectiveDf, plot=False, startPos=2, endPos=12):
    #Get the single objective maximal land use dataframe
    # the columns which contain the crop types are hard coded right now
    theNames = objectiveDf.columns.values.tolist()
    maxPftPerCellStr = objectiveDf[theNames[startPos:endPos]].idxmax(axis=1)
    maxPftPerCell = maxPftPerCellStr
    maxPftPerCell[maxPftPerCellStr =="PNV"] = 0
    maxPftPerCell[maxPftPerCellStr =="C3cer_TeWW"] = 1
    maxPftPerCell[maxPftPerCellStr =="C3cer_TeWWi"] = 2
    maxPftPerCell[maxPftPerCellStr =="C3crps_TeSW"] = 3
    maxPftPerCell[maxPftPerCellStr =="C3crps_TeSWi"] = 4
    maxPftPerCell[maxPftPerCellStr =="C4crps_TeCo"] = 5
    maxPftPerCell[maxPftPerCellStr =="C4crps_TeCoi"] = 6
    maxPftPerCell[maxPftPerCellStr =="Rice_TrRi"] = 7
    maxPftPerCell[maxPftPerCellStr =="Rice_TrRii"] = 8
    maxPftPerCell[maxPftPerCellStr == "Pas"] = 9 # update 04.10.2019
    maxPftPerCell= maxPftPerCell.astype('int32')
    #maxPftPerCell.dtype
    
    
    if plot:
        fig, ax = plt.subplots()
        scatMap = plt.scatter(objectiveDf['Lon'], objectiveDf['Lat'], marker='o', c=maxPftPerCell)
        fig.colorbar(scatMap)
    return(maxPftPerCell)
            
def sigle_point_Array_crossover(random, candidates, args):
    #
    #print('crossover')
    try:
        crossover_rate = args['crossover_rate']
    except KeyError:
        crossover_rate = 1.0
        args['crossover_rate'] = crossover_rate
    cand = list(candidates)
    
    
    if len(cand) % 2 == 1:
        cand = cand[:-1]
    random.shuffle(cand)
    moms = cand[::2]
    dads = cand[1::2]
    children = []
    for mom, dad in izip(moms, dads):
        if random.random() < crossover_rate:
            # get the ndarrays            
            mom=mom.values
            dad = dad.values
            shape = mom.shape # get array size
            
                        
            mom_1d = mom.flatten()
            dad_1d = dad.flatten()
            bro = mom_1d.copy()
            sis = dad_1d.copy()
            cut_point = random.randint(0, dad_1d.shape[0])
            
            bro[cut_point:] = dad_1d[cut_point:]
            sis[:cut_point] = mom_1d[:cut_point]
            bro = bro.reshape(shape)
            sis = sis.reshape(shape)
            # backtransform to panda.Series
            bro = pandas.Series(bro)
            sis = pandas.Series(sis)
            children.append(bro)
            children.append(sis)
            
        else:
            children.append(mom)
            children.append(dad)     
    return children
    try:
        crossover_rate = args['crossover_rate']
    except KeyError:
        crossover_rate = 1.0
    LatLon = args['LatLon']
    
    # get the LAT and LON values
    LON = LatLon['Lon'].values
    LAT = LatLon['Lat'].values
    LatLon= LatLon.values
    

    #get range
    xmin = LON.min()
    xmax = LON.max()
    ymin = LAT.min()
    ymax = LAT.max()
    
    
    # list of candidate solutions
    cand = list(candidates)
    
    # select mom and dad for mating
    if len(cand) % 2 == 1:
        cand = cand[:-1]
    random.shuffle(cand)
    moms = cand[::2]
    dads = cand[1::2]
    children = []
    
    
    
    for mom, dad in izip(moms, dads):
        
        if random.random() < crossover_rate:
            # select two random points inside the range
            x1 = random.randint(xmin, xmax)
            x2 = random.randint(xmin, xmax)
            y1 = random.randint(ymin, ymax)
            y2 = random.randint(ymin, ymax)
            
            # which points are above and below the line
            # beware of extreme cases vertical and horizontal lines
            # use this line to split the genome
            isHorizontal = y1 == y2
            isVertical = x1 == x2
            while(isHorizontal & isVertical):
                # get new first point
                x1 = random.randint(xmin, xmax)
                y1 = random.randint(ymin, ymax)
                isHorizontal = x1 == x2
                isVertical = y1 == y2
            if(isHorizontal):
                # check only x
                pAbove_rows = numpy.where(LAT > y1)
                pBelow_rows = numpy.where(LAT <= y1)
            elif(isVertical):
                # check only y
                pAbove_rows = numpy.where(LON > x1)
                pBelow_rows = numpy.where(LON <= x1)
            else:
                # y = b(x-x1) + y1
                b = (y2-y1)/(x2-x1)
                yLine = y1 + (LON-x1) * (y2-y1)/(x2-x1) 
                pAbove_rows = numpy.where(LAT > yLine)
                pBelow_rows = numpy.where(LAT <= yLine)
            
            #plt.plot(LON[pAbove_rows], LAT[pAbove_rows], marker='o', color='r', ls='')
            #plt.plot(LON[pBelow_rows], LAT[pBelow_rows], marker='o', color='g', ls='')
            #plt.plot([x1,x2], [y1,y2], marker='x', color='b', ls='')
            #plt.plot([x1,x2], [y1,y2],  color='b', ls='-')
            #plt.show()
            
            # operation uses numpy functions but we have pandas series objects
            # so we need to convert them
            
            mom= mom.values
            dad= dad.values
            #shape = mom.shape # get array size
            # now combine mom and dad
            bro = mom.copy()
            sis = dad.copy()
            bro[pAbove_rows] = dad[pAbove_rows]
            sis[pAbove_rows] = mom[pAbove_rows]
            
            #plt.scatter(LON,LAT, c=mom)
            #plt.scatter(LON,LAT, c=dad)
            #plt.scatter(LON, LAT, c=bro)
            #plt.scatter(LON, LAT, c=sis)
            
            # and convert the kids to pandas series again...
            bro = pandas.Series(bro)
            sis = pandas.Series(sis)
            
            children.append(bro)
            children.append(sis)
        else:
            children.append(mom)
            children.append(dad)     
    return children


# mutation operator
def integer_mutation(random, candidates, args):
    #
    #print('mutation')
    #numpy.random.RandomState(int(time()))
    nCells =  args['nCells']
    rangeValues = args['rangeValues']
    mut_genome_array= numpy.random.choice(rangeValues, size=nCells)
    #print(mut_genome_array.dtype)
    try:
        mut_rate = args['mutation_rate']
    except KeyError:
        print("No mutation rate provided!")
        mut_rate = 0.1
        args['mutation_rate'] = mut_rate
    
    cs_copy = list(candidates)
    for i,cand in enumerate(cs_copy):
        #print "mutator"
        #print cand
        shape = cand.shape 
        test_idx = numpy.where(numpy.random.uniform(low=0,high=1,size=shape) < mut_rate)
        
        for anIdx in test_idx:
            cand[anIdx] = mut_genome_array[anIdx]
            #print(mut_genome_array[anIdx])
        cs_copy[i] = cand
    return cs_copy

def setUnfeasible2PNV(FoodProv, luList, threshold):
    # at the FPU or at the biome level it is possible that the
    # PFT is valid overall but not for single cells
    # therefore it is necessary to only rest land use at the grid cell lecel
    # for cell level optimization this could be used to replace the genome
    # as a repair strategy (or we might try a different crop for those cells)
    # get food provisioning values depending on the land use
    valList = FoodProv[(range(0,FoodProv.shape[0] )), luList]
    updatedLuList = numpy.where(valList<=threshold, 2, luList)
    # 09.10.2018
    # now we need a two since I added the offset already above...
    return(updatedLuList)


    
def df_to_sarray(df):
    """
    Convert a pandas DataFrame object to a numpy structured array.
    This is functionally equivalent to but more efficient than
    np.array(df.to_array())

    :param df: the data frame to convert
    :return: a numpy structured array representation of df
    """

    v = df.values
    cols = df.columns

    if six.PY2:  # python 2 needs .encode() but 3 does not
        types = [(cols[i].encode(), df[k].dtype.type) for (i, k) in enumerate(cols)]
    else:
        types = [(cols[i], df[k].dtype.type) for (i, k) in enumerate(cols)]
    dtype = numpy.dtype(types)
    z = numpy.zeros(v.shape[0], dtype)
    for (i, k) in enumerate(z.dtype.names):
        z[k] = v[:, i]
    return z
