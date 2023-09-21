# extending ecspy classes for arrays

import svens_ec_constraints as ec  # was from ecspy import ec, changed 04.01.2018

import time
import copy
import numpy
from ecspy import selectors
from ecspy import variators
from ecspy import replacers
from ecspy import migrators
from ecspy import archivers
from ecspy import terminators
from ecspy import observers

class Individual_array(object):
    """Represents an individual in an evolutionary computation.
    
    An individual is defined by its candidate solution and the
    fitness (or value) of that candidate solution.
    
    Public Attributes:
    
    - *candidate* -- the candidate solution
    - *fitness* -- the value of the candidate solution
    - *birthdate* -- the system time at which the individual was created
    - *constraintViolations* - the constraint violations of the candidate solution
    
    """
    def __init__(self, candidate = None):
        self.candidate = candidate
        self.fitness = None
        self.constraintViolations = None
        self.birthdate = time.time()
    
    def __setattr__(self, name, val):
        if name == 'candidate':
            self.__dict__[name] = val
            self.fitness = None
            self.constraintViolations = None
        else:
            self.__dict__[name] = val
    
    def __str__(self):
        return "%s : %s (constraint violations: %s)" % (str(self.candidate), str(self.fitness), str(self.constraintViolations))
        
    def __repr__(self):
        return "<Individual: candidate = %s, fitness = %s, constraint violations: %s, birthdate = %d>" % (str(self.candidate), str(self.fitness), str(self.constraintViolations), self.birthdate)
        
    def __lt__(self, other):
        if (len(self.fitness) != len(other.fitness)) or (len(self.constraintViolations) != len(other.constraintViolations)) :
            raise Exception("fitness or constraintViolations is not defined")
            
        else:
            not_worse = True
            strictly_better = False
            
            if any(self.constraintViolations) < 0:
                raise Exception("Negative constraint violations for candidate:\n" + str(self) )
            if any(other.constraintViolations) < 0:
                raise Exception("Negative constraint violations for candidate:\n" + str(other) )
            # first check case one and two    
            # number of constraints
            # first set constraints > 0 to 1 
            selfCV01 = [1 if aVal >0 else 0 for aVal in self.constraintViolations]
            otherCV01 = [1 if aVal >0 else 0 for aVal in other.constraintViolations]
            # when count them
            ncvSelf = sum(selfCV01)
            ncvOther = sum(otherCV01)
            # total size of violation: it is important that the constraints have been
            # coded relative to abs value to make them comparable
            sumcvSelf = sum(self.constraintViolations)
            sumcvOther = sum(other.constraintViolations)
            # 
            # Case 1: One solution is feasible and the other solution is infeasible.
            if (ncvSelf ==0) and (ncvOther >0):
                # self has no constraints but other has constraints
                # -> self is not worse (<) than other
                return False
            elif (ncvSelf >0) and (ncvOther == 0):
                # self has constraints but other has no constraints
                # -> self is worse (<) than other
                return True
            
            #Case 2: Both solutions are both infeasible; however,
            # one solution  has a smaller constraint violation than the other.
            # I use the sum of all constraints - however one might also first count the
            # number of objectives there violations occured...
            elif (ncvSelf > 0) and (ncvOther > 0):
                # are both equal?
                if sumcvSelf == sumcvOther:
                    # if number of constrains for self is larger then for
                    # other self is worse (<) then other
                    if ncvSelf == ncvOther:
                        # same number of dimensions violated, same sum
                        # of normalized violations
                        # return the dominant one
                        for x, y in zip(self.fitness, other.fitness):
                            # this is the same code from below
                            # not so nice to have duplicated code
                            # maybe put it in a function
                            if x > y:
                                not_worse = False
                            elif y > x:
                                strictly_better = True
                        return not_worse and strictly_better
                    else:
                        # different number of dimensions with violation
                        # if self has more dimensions violated it is worse (<)
                        # than other
                        ######
                        # check here just changed this
                        return ncvSelf > ncvOther
                # if the sum of (normalized) constaints in all dimensions
                # is higher in self it is worse (<) then other
                #if both are 
                else:
                    # if both candidates have at least on contraint violation
                    # and sum of violations is not equal
                    # return false if self has a larger sum of violations
                    return sumcvSelf > sumcvOther
                    
            #Case 3: Both solutions are both feasible
            # self = x
            # other = y
            for x, y in zip(self.fitness, other.fitness):
               
                if x > y:
                    not_worse = False
                elif y > x:
                    strictly_better = True
            return not_worse and strictly_better
        if self.fitness is not None and other.fitness is not None:
            return self.fitness < other.fitness
            

    def __le__(self, other):
        return self < other or not other < self
            
    def __gt__(self, other):
        return other < self

    def __ge__(self, other):
        return other < self or not self < other
        
    def __eq__(self, other):
        #return self.candidate == other.candidate
        return  numpy.array_equiv(self.candidate, other.candidate)
        
    def __ne__(self, other):
        return self.candidate != other.candidate



class EvolutionaryComputation_array(ec.EvolutionaryComputation):
    
    # not sure if we need the constructure explicitl here, but to be on the save side
    def __init__(self, random):
        self._random = random
        self.selector = selectors.default_selection
        self.variator = variators.default_variation
        self.replacer = replacers.default_replacement
        self.migrator = migrators.default_migration
        self.observer = observers.default_observer
        self.archiver = archivers.default_archiver
        self._kwargs = dict()
     
    def _should_terminate(self, terminator, pop, ng, ne):
        terminate = False
        try:
            for clause in terminator:
                terminate = terminate or clause(population=pop, num_generations=ng, num_evaluations=ne, args=self._kwargs)
        except TypeError:
            terminate = terminator(population=pop, num_generations=ng, num_evaluations=ne, args=self._kwargs)
        return terminate
        
    
    def evolve(self, generator, evaluator, pop_size=100, seeds=[], terminator=terminators.default_termination, **args):
        """Perform the evolution.
        
        This function creates a population and then runs it through a series
        of evolutionary epochs until the terminator is satisfied. The general
        outline of an epoch is selection, variation, evaluation, replacement,
        migration, archival, and observation. 
        
        Arguments:
        
        - *generator* -- the function to be used to generate candidate solutions 
        - *evaluator* -- the function to be used to evaluate candidate solutions
        - *pop_size* -- the number of Individuals in the population (default 100)
        - *seeds* -- an iterable collection of candidate solutions to include
          in the initial population (default [])
        - *terminator* -- the terminator (or iterable collection of terminators)
          to be used to determine whether the evolutionary process
          has finished (default terminators.default_termination)
        - *args* -- a dictionary of keyword arguments

        Note that this method will fail if the generator and evaluator 
        parameters are left with their default values. Note also that the
        *_kwargs* class variable will be initialized to the args parameter here.
        It will also be modified to include the following 'built-in' keyword
        arguments (each preceded by an underscore) unless these arguments are
        already present (which shouldn't be the case):
        
        - *_generator* -- the generator used for creating candidate solutions
        - *_evaluator* -- the evaluator used for evaluating solutions
        - *_terminator* -- the particular terminator(s) used
        - *_population_size* -- the size of the population
        - *_num_generations* -- the number of generations that have elapsed
        - *_num_evaluations* -- the number of evaluations that have been completed
        - *_population* -- the current population
        - *_archive* -- the current archive
        
        """
        self._kwargs = args
        # Add entries to the keyword arguments dictionary
        # if they're not already present.
        try:
            self._kwargs['_generator']
        except KeyError:
            self._kwargs['_generator'] = generator
        try:
            self._kwargs['_evaluator']
        except KeyError:
            self._kwargs['_evaluator'] = evaluator
        try:
            self._kwargs['_terminator']
        except KeyError:
            self._kwargs['_terminator'] = terminator
        try:
            self._kwargs['_population_size']
        except KeyError:
            self._kwargs['_population_size'] = pop_size
        # new: in case the number of starting genomes to be importet has not been set we assume it to be zero
        try:
            self._kwargs['nStartGenomes']
        except KeyError:
            self._kwargs['nStartGenomes'] = 0
        # likewise for the number of genomes read already
        try:
            self._kwargs['nStartGenomesCreated']
        except KeyError:
            self._kwargs['nStartGenomesCreated'] = 0
        
        # Create the initial population.
        try:
            iter(seeds)
        except TypeError:
            seeds = [seeds]
        initial_cs = list(seeds)
        num_generated = max(pop_size - len(seeds), 0)
        i = 0
        while i < num_generated:
            # 
            if(self._kwargs['nStartGenomesCreated'] < self._kwargs['nStartGenomes']):
                self._kwargs['nStartGenomesCreated'] += 1 
                try:
                    cs = self._kwargs['generator_predef'](args=self._kwargs)
                except (OSError, IOError) as e:
                    print "Fallback to random initialization"
                    cs = generator(random=self._random, args=self._kwargs)
                
            else:
                # if not predfined create random list
                cs = generator(random=self._random, args=self._kwargs)
            # the following 6 lines have been changed /added compared to the original ec implementation
            # before it was just #if cs not in initial_cs:
            present = False
            for ini_cs in initial_cs:
                if numpy.array_equiv(cs, ini_cs):
                    present = True
                    
            if not present: 
                initial_cs.append(cs)
                i += 1
        initial_fit = evaluator(candidates=initial_cs, args=self._kwargs)
        population = []
        archive = []
        for cs, fit in zip(initial_cs, initial_fit):
            ind = Individual_array(cs)
            ind.fitness = fit
            population.append(ind)
            
        num_evaluations = len(initial_fit)
        num_generations = 0
        self._kwargs['_num_generations'] = num_generations
        self._kwargs['_num_evaluations'] = num_evaluations
        
        population.sort(key=lambda x: x.fitness, reverse=True)
        self._kwargs['_population'] = population
        
        pop_copy = list(population)
        arc_copy = list(archive)
        archive = self.archiver(random=self._random, population=pop_copy, archive=arc_copy, args=self._kwargs)
        self._kwargs['_archive'] = archive
        
        try:
            for obs in self.observer:
                obs(population=population, num_generations=num_generations, num_evaluations=num_evaluations, args=self._kwargs)
        except TypeError:
            self.observer(population=population, num_generations=num_generations, num_evaluations=num_evaluations, args=self._kwargs)
        while not self._should_terminate(terminator, population, num_generations, num_evaluations):
            # Select individuals.
            pop_copy = list(population)
            parents = self.selector(random=self._random, population=pop_copy, args=self._kwargs)
            
            # Sort the parents just before taking out the candidate so that relative fitness
            # can be determined in the variators (e.g., differential crossover).
            parents.sort(key=lambda x: x.fitness, reverse=True)
            parent_cs = [copy.deepcopy(i.candidate) for i in parents]
            offspring_cs = parent_cs
            try:
                for op in self.variator:
                    offspring_cs = op(random=self._random, candidates=offspring_cs, args=self._kwargs)
            except TypeError:
                offspring_cs = self.variator(random=self._random, candidates=offspring_cs, args=self._kwargs)
            
             
                
            # Evaluate offspring.
            offspring_fit = evaluator(candidates=offspring_cs, args=self._kwargs)
            offspring = []
            for cs, fit in zip(offspring_cs, offspring_fit):
                off = Individual_array(cs)
                off.fitness = fit
                offspring.append(off)

            
            num_evaluations += len(offspring_fit)        
            self._kwargs['_num_evaluations'] = num_evaluations

            # Replace individuals.
            pop_copy = self.replacer(random=self._random, population=pop_copy, parents=parents, offspring=offspring, args=self._kwargs)
            
            # Migrate individuals.
            population = self.migrator(random=self._random, population=pop_copy, args=self._kwargs)
            population.sort(key=lambda x: x.fitness, reverse=True)
            self._kwargs['_population'] = population
            
            # Archive individuals.
            pop_copy = list(population)
            arc_copy = list(archive)
            archive = self.archiver(random=self._random, archive=arc_copy, population=pop_copy, args=self._kwargs)
            self._kwargs['_archive'] = archive
            
            num_generations += 1
            self._kwargs['_num_generations'] = num_generations
            try:
                for obs in self.observer:
                    obs(population=population, num_generations=num_generations, num_evaluations=num_evaluations, args=self._kwargs)
            except TypeError:
                self.observer(population=population, num_generations=num_generations, num_evaluations=num_evaluations, args=self._kwargs)        
        return archive
    
    def evolve_nsgaConstraint(self, generator, evaluator, pop_size=100, seeds=[], terminator=terminators.default_termination, **args):
        """Perform the evolution.
        
        This function creates a population and then runs it through a series
        of evolutionary epochs until the terminator is satisfied. The general
        outline of an epoch is selection, variation, evaluation, replacement,
        migration, archival, and observation. 
        
        Arguments:
        
        - *generator* -- the function to be used to generate candidate solutions 
        - *evaluator* -- the function to be used to evaluate candidate solutions
        - *pop_size* -- the number of Individuals in the population (default 100)
        - *seeds* -- an iterable collection of candidate solutions to include
          in the initial population (default [])
        - *terminator* -- the terminator (or iterable collection of terminators)
          to be used to determine whether the evolutionary process
          has finished (default terminators.default_termination)
        - *args* -- a dictionary of keyword arguments

        Note that this method will fail if the generator and evaluator 
        parameters are left with their default values. Note also that the
        *_kwargs* class variable will be initialized to the args parameter here.
        It will also be modified to include the following 'built-in' keyword
        arguments (each preceded by an underscore) unless these arguments are
        already present (which shouldn't be the case):
        
        - *_generator* -- the generator used for creating candidate solutions
        - *_evaluator* -- the evaluator used for evaluating solutions
        - *_terminator* -- the particular terminator(s) used
        - *_population_size* -- the size of the population
        - *_num_generations* -- the number of generations that have elapsed
        - *_num_evaluations* -- the number of evaluations that have been completed
        - *_population* -- the current population
        - *_archive* -- the current archive
        
        """
        self._kwargs = args
        # Add entries to the keyword arguments dictionary
        # if they're not already present.
        try:
            self._kwargs['_generator']
        except KeyError:
            self._kwargs['_generator'] = generator
        try:
            self._kwargs['_evaluator']
        except KeyError:
            self._kwargs['_evaluator'] = evaluator
        try:
            self._kwargs['_terminator']
        except KeyError:
            self._kwargs['_terminator'] = terminator
        try:
            self._kwargs['_population_size']
        except KeyError:
            self._kwargs['_population_size'] = pop_size
        # new: in case the number of starting genomes to be importet has not been set we assume it to be zero
        try:
            self._kwargs['nStartGenomes']
        except KeyError:
            self._kwargs['nStartGenomes'] = 0
        # likewise for the number of genomes read already
        try:
            self._kwargs['nStartGenomesCreated']
        except KeyError:
            self._kwargs['nStartGenomesCreated'] = 0
        
        # Create the initial population.
        try:
            iter(seeds)
        except TypeError:
            seeds = [seeds]
        initial_cs = list(seeds)
        num_generated = max(pop_size - len(seeds), 0)
        i = 0
        while i < num_generated:
            cs = generator(random=self._random, args=self._kwargs)
            # the following 6 lines have been changed /added compared to the original ec implementation
            # before it was just #if cs not in initial_cs:
            present = False
            for ini_cs in initial_cs:
                if numpy.array_equiv(cs, ini_cs):
                    present = True
                    
            if not present: 
                initial_cs.append(cs)
                i += 1
        initial_fitnesWithConstraints = evaluator(candidates=initial_cs, args=self._kwargs)
                
        # evaluator returns here object of class
        # emo_constraint.constraint_array_Pareto
        population = []
        archive = []
        for cs, fitConstr in zip(initial_cs, initial_fitnesWithConstraints):
            ind = Individual_array(cs)
            ind.fitness = fitConstr.values
            ind.constraintViolations = fitConstr.constraintViolations
            population.append(ind)
            #print(str(ind))
            
        num_evaluations = len(initial_fitnesWithConstraints)
        num_generations = 0
        self._kwargs['_num_generations'] = num_generations
        self._kwargs['_num_evaluations'] = num_evaluations
        
        #population.sort(key=lambda x: x.fitness, reverse=True)
        population.sort(reverse=True) # uses the __cmp__ function
        self._kwargs['_population'] = population
        
        pop_copy = list(population)
        arc_copy = list(archive)
        archive = self.archiver(random=self._random, population=pop_copy, archive=arc_copy, args=self._kwargs)
        self._kwargs['_archive'] = archive
        
        try:
            for obs in self.observer:
                obs(population=population, num_generations=num_generations, num_evaluations=num_evaluations, args=self._kwargs)
        except TypeError:
            self.observer(population=population, num_generations=num_generations, num_evaluations=num_evaluations, args=self._kwargs)
        while not self._should_terminate(terminator, population, num_generations, num_evaluations):
            # Select individuals.
            pop_copy = list(population)
            parents = self.selector(random=self._random, population=pop_copy, args=self._kwargs)
            
            # Sort the parents just before taking out the candidate so that relative fitness
            # can be determined in the variators (e.g., differential crossover).
            #parents.sort(key=lambda x: x.fitness, reverse=True)
            parents.sort(reverse=True)
            parent_cs = [copy.deepcopy(i.candidate) for i in parents]
            offspring_cs = parent_cs
            try:
                for op in self.variator:
                    offspring_cs = op(random=self._random, candidates=offspring_cs, args=self._kwargs)
            except TypeError:
                offspring_cs = self.variator(random=self._random, candidates=offspring_cs, args=self._kwargs)
            
             
                
            # Evaluate offspring.
            
            offspring_fitnessWithConstraints= evaluator(candidates=offspring_cs, args=self._kwargs)
            
            offspring = []
            for cs, fitConstr in zip(offspring_cs, offspring_fitnessWithConstraints):
                off = Individual_array(cs)
                off.fitness = fitConstr.values
                off.constraintViolations = fitConstr.constraintViolations
                offspring.append(off)

            
            num_evaluations += len(offspring_fitnessWithConstraints)        
            self._kwargs['_num_evaluations'] = num_evaluations

            # Replace individuals.
            pop_copy = self.replacer(random=self._random, population=pop_copy, parents=parents, offspring=offspring, args=self._kwargs)
            
            # Migrate individuals.
            population = self.migrator(random=self._random, population=pop_copy, args=self._kwargs)
            #population.sort(key=lambda x: x.fitness, reverse=True)
            population.sort(reverse=True)
            self._kwargs['_population'] = population
            
            # Archive individuals.
            pop_copy = list(population)
            arc_copy = list(archive)
            archive = self.archiver(random=self._random, archive=arc_copy, population=pop_copy, args=self._kwargs)
            self._kwargs['_archive'] = archive
            
            num_generations += 1
            print("------------------------")
            print("num_generations: " + str(num_generations))
            print("------------------------")
            self._kwargs['_num_generations'] = num_generations
            try:
                for obs in self.observer:
                    obs(population=population, num_generations=num_generations, num_evaluations=num_evaluations, args=self._kwargs)
            except TypeError:
                self.observer(population=population, num_generations=num_generations, num_evaluations=num_evaluations, args=self._kwargs)        
        return archive


class GA_array(EvolutionaryComputation_array):
    """Evolutionary computation representing a canonical genetic algorithm.
    
    This class represents a genetic algorithm which uses, by 
    default, fitness proportionate selection, n-point crossover,
    bit-flip mutation, and generational replacement. In the case
    of bit-flip mutation, it is expected that the candidate 
    solution is an iterable object of binary values. 
    
    """
    def __init__(self, random):
        EvolutionaryComputation_array.__init__(self, random)
        self.selector = selectors.fitness_proportionate_selection
        self.variator = [variators.n_point_crossover, variators.bit_flip_mutation]
        self.replacer = replacers.generational_replacement
        
    def evolve(self, generator, evaluator, pop_size=100, seeds=[], terminator=terminators.default_termination, **args):
        try:
            args['num_selected']
        except KeyError:
            args['num_selected'] = pop_size
        return EvolutionaryComputation_array.evolve(self, generator, evaluator, pop_size, seeds, terminator, **args)



