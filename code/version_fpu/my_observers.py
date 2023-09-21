import os
import sys
import time
import numpy
import pickle
import matplotlib
matplotlib.use('PDF')
import matplotlib.pyplot as plt


def performance2file_observer(population, num_generations, num_evaluations, args):   
    
    
    try:
        observer_file = args['observer_file']
    except KeyError:
        filename = 'ecspy-observer-file-' + str(time.time()) + ".out"
        observer_file = open(filename, 'w')
        args['observer_file'] = observer_file

   
    best_fitness = population[0].fitness
    average_fitness = sum([x.fitness for x in population]) / float(len(population))
    median_fitness = population[len(population)/2].fitness
    
    data = []
    if num_generations == 0:
        observer_file.write('num_generations \t num_evaluations \t best_fitness \t average_fitness \t median_fitness \n')
        
    observer_file.write('%d \t %d \t %d \t %0.5f \t %0.5f \n'% (num_generations, num_evaluations, best_fitness,  average_fitness, median_fitness) )
        
    
    #observer_file.write('Generation: %d \n' % num_generations)
    #observer_file.write('Evaluations: %d \n' % num_evaluations)
    #avg_fit = sum([x.fitness for x in population]) / float(len(population))
    #observer_file.write('Average Fitness: %0.5f     Maximum Fitness: %0.5f \n' % (avg_fit, population[0].fitness))
    
def performance2file_nsga_observer(population, num_generations, num_evaluations, args):     
    try:
        observer_file = args['observer_file']
    except KeyError:
        filename = 'ecspy_observer_file-' + str(time.time()) + ".out"
        observer_file = open(filename, 'w')
        args['observer_file'] = observer_file
        # will not be closed....
    try:
        archiveFitness_file = args['archiveFitness_file']
    except KeyError:
        filename = 'ecspy_archiveFitness_file-' + str(time.time()) + ".out"
        archiveFitness_file = open(filename, 'w')
        args['archiveFitness_file'] = archiveFitness_file
        # will not be closed....
    try:
        outFolder = args['outFolder']
    except KeyError:
        outFolder = ""
    try:
        freqPopReport= args['freqPopReport']
    except KeyError:
        freqPopReport= freqPopReport
    try:
        freqArchiveReport= args['freqArchiveReport ']
    except KeyError:
        freqArchiveReport= freqPopReport
    try:
        freqArchiveDump= args['freqArchiveDump']
    except KeyError:
        freqArchiveDump= freqPopReport
    try:
        out_id= args['out_id']
    except KeyError:
        out_id= ""


    #print("freqPopReport:" + str(freqPopReport))
    #print("freqArchiveReport:" + str(freqArchiveReport))
    #print("freqArchiveDump:" + str(freqArchiveDump))
    
    # number of objectives
    nobj= len(population[0].fitness)
    nviol= len(population[0].constraintViolations)
   
    if freqPopReport > 0:
        # report of population fitness 
        if num_generations == 0:
            # write header
            observer_file.write('generation \t individual \t')
            for i in range(0,nobj):
                observer_file.write('obj_' + str(i) + '\t' )
            for i in range(0,nviol):
                observer_file.write('violConstr_' + str(i) + '\t' )
            observer_file.write('\n')
            # print initial population fitness and constraints
            for i, p in enumerate(population):
                observer_file.write('%d \t %d \t' % (num_generations, i) )
                for j in range(0,nobj):
                    observer_file.write('%0.5f \t'% (p.fitness[j]) )
                for j in range(0,nviol):
                    observer_file.write('%0.5f \t'% (p.constraintViolations[j]) )
                observer_file.write('\n')
            observer_file.flush()
            os.fsync(observer_file.fileno())
        # for following generations
        elif num_generations % freqPopReport == 0:
            for i, p in enumerate(population):
                observer_file.write('%d \t %d \t' % (num_generations, i) )
                for j in range(0,nobj):
                    observer_file.write('%0.5f \t'% (p.fitness[j]) )
                for j in range(0,nviol):
                    observer_file.write('%0.5f \t'% (p.constraintViolations[j]) )
                observer_file.write('\n')
            observer_file.flush()
            os.fsync(observer_file.fileno())
        
        
    
    if freqArchiveReport > 0:    
        # archive report
        if num_generations == 0:
            # write header
            archiveFitness_file.write('generation \t individual \t')
            for i in range(0,nobj):
                archiveFitness_file.write('obj_' + str(i) + '\t' )
            for i in range(0,nviol):
                archiveFitness_file.write('violConstr_' + str(i) + '\t' )
            archiveFitness_file.write('\n')
            # initial population
            for i, p in enumerate(args['_archive']):
                archiveFitness_file.write('%d \t %d \t' % (num_generations, i) )
                for j in range(0,nobj):
                    archiveFitness_file.write('%0.5f \t'% (p.fitness[j]) )
                for j in range(0,nviol):
                    archiveFitness_file.write('%0.5f \t'% (p.constraintViolations[j]) )
                archiveFitness_file.write('\n')
            archiveFitness_file.flush()
            os.fsync(archiveFitness_file.fileno())
        elif num_generations % freqArchiveReport == 0:
            # write info on archive population as well
            for i, p in enumerate(args['_archive']):
                archiveFitness_file.write('%d \t %d \t' % (num_generations, i) )
                for j in range(0,nobj):
                    archiveFitness_file.write('%0.5f \t'% (p.fitness[j]) )
                for j in range(0,nviol):
                    archiveFitness_file.write('%0.5f \t'% (p.constraintViolations[j]) )
                archiveFitness_file.write('\n')
            archiveFitness_file.flush()
            os.fsync(archiveFitness_file.fileno())
            
            
    if freqArchiveDump > 0:
        if num_generations % freqArchiveDump == 0:
            fn = outFolder + 'interm_pop_' + str(out_id) + '_' + str(num_generations) + '.pkl'
            picklefile = open(fn, 'wb')
            # new 07.04.2016: write archive not population
            pickle.dump(args['_archive'], picklefile)
            picklefile.close()
