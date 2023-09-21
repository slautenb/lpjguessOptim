# -*- coding: utf-8 -*-
"""
Created on Mon Feb 05 10:49:47 2018
read the selected predefined solutions and preprocess them so
that they can be used as seeds for the cell level optimization
@author: Lautenbach
"""

import pandas as pd
import fnmatch
import os
import pickle
import ec_arrayGenome_4ConstraintsAndPredefinedGenomes  as ec_arrayGenome

# Get Seeds for FPU level optimization based on the results
# of the biome level optimization. The expansion from biome to FPU has already been done.
# So I can simply read in the solutions and use them as the genom.
# takes the genomes of the selected candidates and creates a list
# that can be used as seeds
def prepareSeedsForFPU(biomesSolutionsPath, basefn,  maxNseeds=-1):
    # biomesSolutionsPath: path to the folder with the genomes stored in csv files
    # basefn: filename without the running numbeer and .csv file extension
    # maxNseeds: in case less than the total available number of seeds should be used
    #            e.g. for testing purposes
    #            if set to -1 use all

    # file numbering assumed to start with 1
    selectedGenomesList = []

    nSeeds = len(fnmatch.filter(os.listdir(biomesSolutionsPath), basefn + '*.csv'))

    for i in range(0, nSeeds):
        if (maxNseeds == -1) | (maxNseeds >= i):
            fn = biomesSolutionsPath + basefn + str(i + 1) + ".csv"
            try:
                aFpuGenome = pd.read_csv(fn, sep=";")
            except (OSError, IOError) as e:
                print "generator_predef: \n" + e.message
                raise e

            selectedGenomesList.append(aFpuGenome['PFT'])
        else:
            break
    return (selectedGenomesList)


def prepareSeedsForFPUIncorrect(biomesSolutionsPath, basefn, fpuGrid,  maxNseeds=-1):
    # biomesSolutionsPath: path to the folder with the genomes stored in csv files
    # basefn: filename without the running numbeer and .csv file extension
    # fpuGrid: the pandas data frame to align the genome to
    # maxNseeds: in case less than the total available number of seeds should be used
    #            e.g. for testing purposes
    #            if set to -1 use all

    # file numbering assumed to start with 1
    selectedGenomesList = []

    nSeeds = len(fnmatch.filter(os.listdir(biomesSolutionsPath), basefn + '*.csv'))

    for i in range(0, nSeeds):
        if (maxNseeds == -1) | (maxNseeds >= i):
            fn = biomesSolutionsPath + basefn + str(i + 1) + ".csv"
            try:
                aFpuGenome = pd.read_csv(fn, sep=";")
            except (OSError, IOError) as e:
                print "generator_predef: \n" + e.message
                raise e

            aFpuGridWithGenome = pd.merge(fpuGrid, aFpuGenome, on=('newFPU'))
            # need to order accordingly to the objectives
            aFpuGridWithGenome.sort_values(by='ID', inplace=True)
            aFpuGridWithGenome.reset_index(drop=True, inplace=True)

            selectedGenomesList.append(aFpuGridWithGenome['PFT'])
        else:
            break
    return (selectedGenomesList)


#candidates_combined10runs_betterThanCurrent.pkl
def prepareSeeds(fpuSolutionsPath, fn, fpuGrid, fpuOrder, maxNseeds=-1):
    # fpuSolutionsPath: path to the pickle file where the already selected solutions are stored
    # maxNseeds: in case less than the total available number of seeds should be used
    #            e.g. for testing purposes
    #            if set to -1 use all
    picklefile = open(fpuSolutionsPath  + fn, 'r')
    selectedFpuSolutions = pd.read_pickle(picklefile) #pickle.load(picklefile)
    picklefile.close()
    #selectedFpuSolutions = pd.read_hdf(fpuSolutionsPath  + fn, 'selectedSubset')
    selectedGenomesList =[]
    # Now create a cell level genome from the fpu level genome
    i = 0
    for aFpuCs in selectedFpuSolutions:
        i += 1
        if(maxNseeds == -1) | (maxNseeds >= i):
            #type(aFpuCs) # ec_arrayGenome_4predefinedGenomes.Individual_array
            cs = aFpuCs.candidate
            
            aFpuGenome = {'newFPU' : fpuOrder, 'PFT' : cs} 
            
            aFpuGenome= pd.DataFrame(aFpuGenome)
            
            aFpuGridWithGenome = pd.merge(fpuGrid, aFpuGenome, on=('newFPU'))
            # need to order accordingly to the objectives
            aFpuGridWithGenome.sort_values(by='ID', inplace=True)      
            aFpuGridWithGenome.reset_index(drop=True, inplace=True)
            
            selectedGenomesList.append(aFpuGridWithGenome['PFT'])
        else:
            break
    return(selectedGenomesList)

def prepareSeeds_h5(fpuSolutionsPath, fn, fpuGrid, fpuOrder, maxNseeds=-1):
    # path: path to the pickle file where the already selected solutions are stored
    # maxNseeds: in case less than the total available number of seeds should be used
    #            e.g. for testing purposes
    #            if set to -1 use all
    myStore = pd.HDFStore(fpuSolutionsPath + fn)
    selectedGenomesList =[]

    theKeys = myStore.keys()



    # Now create a cell level genome from the fpu level genome
    i = 0
    for aKey in theKeys:
        i += 1
        if(maxNseeds == -1) | (maxNseeds >= i):
            #type(aFpuCs) # ec_arrayGenome_4predefinedGenomes.Individual_array
            cs = myStore[aKey]
            
            aFpuGenome = {'newFPU' : fpuOrder, 'PFT' : cs} 
            
            aFpuGenome= pd.DataFrame(aFpuGenome)
            
            aFpuGridWithGenome = pd.merge(fpuGrid, aFpuGenome, on=('newFPU'))
            # need to order accordingly to the objectives
            aFpuGridWithGenome.sort_values(by='ID', inplace=True)      
            aFpuGridWithGenome.reset_index(drop=True, inplace=True)
            
            selectedGenomesList.append(aFpuGridWithGenome['PFT'])
        else:
            break
    myStore.close()
    return(selectedGenomesList)


def preProcessInput4ProtAreas(objDf, protAreas, PNV_pos, crop_pos):
    # objDf objective as a pandas DataFrame
    # protAreas protected areas as a pandas DataFrame
    # for each cell that is (partially) proteced
    modObjDf = objDf.copy()
    theNames = objDf.columns.values.tolist()
    PNV_name = theNames[PNV_pos]
    
    for i in crop_pos:
        aCrop = theNames[i]
        modObjDf[aCrop] = objDf[aCrop]*(1-protAreas['ProtPerc']) + objDf[PNV_name]*protAreas['ProtPerc']
        
    return(modObjDf)

    

        
    
