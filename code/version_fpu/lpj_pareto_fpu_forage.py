# -*- coding: utf-8 -*-
"""
Created on Friday 4th of October 2019
LPJ-GUESS optimization at the fpu level
now runs under the additional constraint of current forage level
@author: Sven Lautenbach
"""

# set up
import numpy, numpy.random
import pandas as pd
import sys, os, getopt, fnmatch

fpuCodePath = sys.path[0]
basePath = fpuCodePath + "/../../"
gaBasePath = fpuCodePath + "/../ecspy-0.4/"

sys.path.append(gaBasePath)
sys.path.append(fpuCodePath)

inPath = basePath + "data/" #
# path to where the output from the biomes Pareto frontier step are stored
# this should be the fpu genomes derived from the biome level solutions
inStartGenomesFolderBase = basePath + "results/biomes/fpuGenomes/"

from ecspy import terminators

import my_observers
# added for the import of predefined genomes
import ec_arrayGenome_4ConstraintsAndPredefinedGenomes  as ec_arrayGenome
import my_replacers
import my_archivers

import pickle

import my_selectors as selectors
import emo_constraint

from lpj_pareto_functions import *


import preprocessingFunctions as preprocess


def assure_path_exists(path):
        dir = os.path.dirname(path)
        if not os.path.exists(dir):
                os.makedirs(dir)

def main(argv):
    print ("started")
    run_ID = None
    mutation_rate = None
    corFac = None
    tournament_size = None
    max_generations = None
    pop_extra = None
    try:
        opts, args = getopt.getopt(argv, "r:",
                                   ['mutation_rate=', 'corFac=', 'tournament_size=', 'max_generations=', 'time_period=', 'rcp=',
                                    'pop_extra='])

    except getopt.GetoptError:
        print 'lpj_pareto_fpu_forage.py -r <runid> --mutation_rate= --corFac= --tournament_size= --max_generations= --time_period=<2042|2099> --rcp=<26|60> --pop_extra='
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-r':
            run_ID = arg
            print 'setting run_ID to ' + str(run_ID) + ' and using run_' + str(run_ID) + ' as the output folder'
        if opt == '--mutation_rate':
            mutation_rate = float(arg)
        if opt == '--corFac':
            corFac = float(arg)
        if opt == '--tournament_size':
            tournament_size = int(arg)
        if opt == '--max_generations':
            max_generations = int(arg)
        if opt == '--pop_extra':
            if int(arg) < 0:
                print("pop_size must be a positive integer")
                sys.exit(2)
            pop_extra = int(arg)
        if opt == '--time_period':
            if int(arg) == 2042:
                time_period = 2042
                time_period_long = "2033-2042"
            elif int(arg) == 2090:
                time_period = 2090
                time_period_long = "2090-2099"
            else:
                print 'time_period required. Valid values are 2042 or 2090'
                sys.exit(2)
        if opt == '--rcp':
            if int(arg) == 26:
                rcp = 26
            elif int(arg) == 60:
                rcp = 60
            else:
                print 'rcp required. Valid values are 26 or 60'
                sys.exit(2)
    if run_ID is None:
        print(
            'You need to specify a run_ID using lpj_pareto_main.py -r <runid> --mutation_rate= --corFac= --tournament_size=')
        sys.exit(2)
    if mutation_rate is None:
        mutation_rate = 0.0125
        print 'setting mutation rate to standard: 0.0125'
    if corFac is None:
        corFac = .98
        print ' Setting correction factor for constraints to default of 0.98'
    if tournament_size is None:
        tournament_size = 2
        print ' Setting tournament_size to default of 2'
    if max_generations is None:
        max_generations = 500
        print ' Setting max_generations to default of 500'
    if pop_extra is None:
        pop_extra = 180
        print ' Setting pop_extra to default of 160'
    # create output path
    outPath = basePath + "results/fpu/run_" + str(run_ID) + "_" + str(rcp) + "_" + str(time_period) + "/"
    assure_path_exists(outPath)
    ###########################################################
    # read input data
    ###########################################################
    #ESquant_rcp26_2033-2042_INT_9426_v8b_4gcms_CStor_cor4ProtArea
    # load pre-calculated LPJ-GUESS output
    fnamePart1= inPath + "/ESquant/ESquant_rcp" + str(rcp) + "_" + str(time_period_long) + "_INT_9412_v8b_4gcms_"
    fnamePart3 = "_cor4ProtArea.txt"

    Cstorage = pd.read_csv( fnamePart1 + "CStor" + fnamePart3, delim_whitespace=True)
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
    fodPos = 1  # position of food provisioning in the lookup tables list

    # get the FPU order (all genomes should be sorted following to this order!)
    inStartGenomesFolder = inStartGenomesFolderBase + "fpuGenomes_" + str(rcp) + "_" + str(time_period_long) + "/"
    try:
        aFpu = pd.read_csv(inStartGenomesFolder + "fpuGenome_1.csv", sep=";")
    except Exception as e:
        print str(e)
        print("Problem opening fupGenome file: " + inStartGenomesFolder + "fpuGenome_1.csv")
        sys.exc_clear(3)
    fpuOrder = aFpu.newFPU

    fpuGrid = pd.read_csv(basePath + "data/fpu2biomeGrid.csv", sep=";")

    # get number of fpu genomes present for the rcp/time_period combination
    #nSeeds = len(fnmatch.filter(os.listdir(inStartGenomesFolder), 'fpuGenome_*.csv'))

    # get list of solutions from the results of the biome level optimization
    theSeeds= preprocess.prepareSeedsForFPU( biomesSolutionsPath= inStartGenomesFolder, basefn= "fpuGenome_",  maxNseeds=-1)

    # todo: if required create additional seeds
    nAdditional = pop_extra / 2  # divide by two since we always get two candidates back
    # todo: update the number of seeds if additional solutions are added
    nSeeds = len(theSeeds)
    pop_size = nSeeds


    # the number of FPUs defines the size of the genome
    nCells = aFpu.shape[0]  # Cstorage.shape[0]
    rangeValues = range(0, 10)  # the number of land use options, updated to 10 (including pasture)

    # minmum calories produced used as a threshold to set PFT to PNV if below
    foodMinThreshold = 0  # has already been done in the corrected input files, 13.08.2018

    # maximum values for each objective
    # corrected for additional column (pasture), endPos=12 includes this column
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

    # constraints
    constraints = pd.read_csv(inPath + "constraints.csv", delim_whitespace=True)
    # constraints
    # adjust units!
    # and apply correction factor
    constraints['carbon'] = constraints['carbon'] * 10 ** 3 * corFac
    constraints['food'] = constraints['food'] * 10 ** 9 * corFac
    constraints['water'] = constraints['water'] * 10 ** 3 * corFac
    constraints['forage'] = constraints['forage'] * 10 ** 3 * corFac

    theConst = constraints[(constraints['rcp'] == int(rcp)) & (constraints['timePeriod_short'] == time_period)]
    tupleOfConstraints = (float(theConst['carbon']), float(theConst['food']), float(theConst['water']), float(theConst['forage']))

    # write config to file
    configreportfile = open(outPath + "config4run" + str(run_ID) + ".txt", 'wb')
    configreportfile.write("run ID: " + str(run_ID) + "\n")
    configreportfile.write("time period for simulation: " + str(time_period) + "\n")
    configreportfile.write("population size: " + str(pop_size) + "\n")
    configreportfile.write("max_generations: " + str(max_generations) + "\n")
    configreportfile.write("variator: [sigle_point_Array_crossover,  integer_mutation]\n")
    configreportfile.write("landuse range: " + str(rangeValues) + "\n")
    configreportfile.write("foodMinThreshold: " + str(foodMinThreshold) + "\n")
    configreportfile.write("tupleOfConstraints: " + str(tupleOfConstraints) + "\n")
    configreportfile.write("constraintFactor: " + str(corFac) + "\n")
    configreportfile.write("mutation_rate: " + str(mutation_rate) + "\n")
    configreportfile.write("tournament_size: " + str(tournament_size) + "\n")

    configreportfile.close()

    rand = Random()
    rand.seed(int(time()))
    my_ec = ec_arrayGenome.EvolutionaryComputation_array(rand)
    my_ec.observer = my_observers.performance2file_nsga_observer  # performance2file_observer  # observers.plot_observer  # observers.file_observer  #observers.screen_observer #
    # was ec.selectors.tournament_selection
    my_ec.selector = selectors.tournament_selection_constraints  # selectors.fitness_proportionate_selection#tournament_selection
    my_ec.variator = [sigle_point_Array_crossover, integer_mutation]
    my_ec.replacer = my_replacers.nsga_constraint_replacement  # my_replacers.nsga_replacement_array#ec.replacers.nsga_replacement #ec.replacers.truncation_replacement
    my_ec.archiver = my_archivers.best_archiver_array  # ec.archivers.best_archiver_array

    reportFile = open(outPath + 'reportNSGAII.txt', 'w')
    archiveFitness_file = open(outPath + 'reportArchiveNSGAII.txt', 'w')

    final_pop = my_ec.evolve_nsgaConstraint(evaluator=evaluate_lpj_fpuLevel,
                             generator=generate_initial_candidate,
                             generator_predef=generate_selected_candidates,
                             terminator=[terminators.generation_termination],
                             # terminators.avg_fitness_termination], #terminators.evaluation_termination,#
                             max_generations=max_generations,  # max_evaluations=1000,
                             pop_size=nSeeds,
                             seeds=theSeeds,
                             num_selected=nSeeds,
                             mutation_rate=mutation_rate,
                             crossover_rate=1.0,
                             observer_file=reportFile,
                             archiveFitness_file=archiveFitness_file,
                             tourn_size=2,
                             min_fitness_diff=0.001,
                             outFolder=outPath,
                             nCells=nCells,
                             rangeValues=rangeValues,
                             listOfObjectiveLookupTables=listOfObjectiveLookupTables,
                             tupleOfConstraints=tupleOfConstraints,
                             maximumObjectives=maximumObjectives,
                             minimumObjectives=minimumObjectives,
                             freqPopReport=2,
                             freqArchiveReport=1,
                             freqArchiveDump=10,
                             nStartGenomes=nSeeds,
                             inStartGenomesFolder=inStartGenomesFolder,
                             fpuOrder=fpuOrder,
                             fpuGrid=fpuGrid,
                             foodMinThreshold=foodMinThreshold,
                             fodPos=fodPos
                             ) # maybe I don't need inStartGenomesFolder anymore - has been/could be done before hand

    reportFile.close()
    # save results with pickle
    picklefile = open(outPath + "nsGA_1.pkl", 'wb')
    pickle.dump(final_pop, picklefile)
    picklefile.close()

    # in case we need to read the pickle file
    #picklefile = open(outPath + "nsGA_1.pkl", 'r')
    #final_pop = pickle.load(picklefile)
    #picklefile.close()

if __name__ == "__main__":
   main(sys.argv[1:])