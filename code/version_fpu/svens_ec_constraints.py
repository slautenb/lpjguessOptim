"""
    This module provides the framework for making evolutionary computations.
    
    .. Copyright (C) 2009  Inspired Intelligence Initiative

    .. This program is free software: you can redistribute it and/or modify
       it under the terms of the GNU General Public License as published by
       the Free Software Foundation, either version 3 of the License, or
       (at your option) any later version.

    .. This program is distributed in the hope that it will be useful,
       but WITHOUT ANY WARRANTY; without even the implied warranty of
       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
       GNU General Public License for more details.

    .. You should have received a copy of the GNU General Public License
       along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import time
import copy
import sys
codePath = "c:/sven/Nutzer/development/python/GA/ecspy-0.4/"
sys.path.append(codePath)
from ecspy import selectors
from ecspy import variators
from ecspy import replacers
from ecspy import migrators
from ecspy import archivers
from ecspy import terminators
from ecspy import observers


class Individual(object):
    """Represents an individual in an evolutionary computation.
    
    An individual is defined by its candidate solution and the
    fitness (or value) of that candidate solution.
    
    Public Attributes:
    
    - *candidate* -- the candidate solution
    - *fitness* -- the value of the candidate solution
    - *birthdate* -- the system time at which the individual was created
    
    """
    def __init__(self, candidate = None):
        self.candidate = candidate
        self.fitness = None
        self.birthdate = time.time()
    
    def __setattr__(self, name, val):
        if name == 'candidate':
            self.__dict__[name] = val
            self.fitness = None
        else:
            self.__dict__[name] = val
    
    def __str__(self):
        return "%s : %s" % (str(self.candidate), str(self.fitness))
        
    def __repr__(self):
        return "<Individual: candidate = %s, fitness = %s, birthdate = %d>" % (str(self.candidate), str(self.fitness), self.birthdate)
        
    def __lt__(self, other):
        if self.fitness is not None and other.fitness is not None:
            return self.fitness < other.fitness
        else:
            raise Exception("fitness is not defined")

    def __le__(self, other):
        return self < other or not other < self
            
    def __gt__(self, other):
        return other < self

    def __ge__(self, other):
        return other < self or not self < other
        
    def __eq__(self, other):
        return self.candidate == other.candidate
        
    def __ne__(self, other):
        return self.candidate != other.candidate


class EvolutionaryComputation(object):
    """Represents a basic evolutionary computation.
    
    This class encapsulates the components of a generic evolutionary
    computation. These components are the selection mechanism, the
    variation operators, the replacement mechanism, the migration
    scheme, and the observers.
    
    Public Attributes:
    
    - *selector* -- the selection operator
    - *variator* -- the (possibly list of) variation operator(s)
    - *replacer* -- the replacement operator
    - *migrator* -- the migration operator
    - *archiver* -- the archival operator
    - *observer* -- the (possibly list of) observer(s)
    
    Protected Attributes:
    
    - *_random* -- the random number generator object
    - *_kwargs* -- the dictionary of keyword arguments initialized
      from the *args* parameter in the *evolve* method
    
    Public Methods:
    
    - ``evolve`` -- performs the evolution and returns the final
      population of individuals
    
    """
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
        
    
    def evolve(self, generator, evaluator,  generatorSpecial, pop_size=100, seeds=[], terminator=terminators.default_termination, **args):
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
        
        # new 07.10.2015         
        try:
            self._kwargs['_generatorSpecial']
        except KeyError:
            self._kwargs['_generatorSpecial'] = generator          
        try:
            nSpecialStartIndividuals = args['nSpecialStartIndividuals']
        except KeyError:
            nSpecialStartIndividuals = 0
            args['nSpecialStartIndividuals'] = 0
        
        # Create the initial population.
        try:
            iter(seeds)
        except TypeError:
            seeds = [seeds]
        initial_cs = list(seeds)
        num_generated = max(pop_size - len(seeds), 0)
        i = 0
        while i < num_generated:
            if(i <nSpecialStartIndividuals):
               cs = generateSpecialIndividuals(i, args=self._kwargs)
               next
            cs = generator(random=self._random, args=self._kwargs)
            if cs not in initial_cs:
                initial_cs.append(cs)
                i += 1
        initial_fit = evaluator(candidates=initial_cs, args=self._kwargs)
        
        population = []
        archive = []
        for cs, fit in zip(initial_cs, initial_fit):
            ind = Individual(cs)
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
                off = Individual(cs)
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
        

class GA(EvolutionaryComputation):
    """Evolutionary computation representing a canonical genetic algorithm.
    
    This class represents a genetic algorithm which uses, by 
    default, fitness proportionate selection, n-point crossover,
    bit-flip mutation, and generational replacement. In the case
    of bit-flip mutation, it is expected that the candidate 
    solution is an iterable object of binary values. 
    
    """
    def __init__(self, random):
        EvolutionaryComputation.__init__(self, random)
        self.selector = selectors.fitness_proportionate_selection
        self.variator = [variators.n_point_crossover, variators.bit_flip_mutation]
        self.replacer = replacers.generational_replacement
        
    def evolve(self, generator, evaluator, pop_size=100, seeds=[], terminator=terminators.default_termination, **args):
        try:
            args['num_selected']
        except KeyError:
            args['num_selected'] = pop_size
        return EvolutionaryComputation.evolve(self, generator, evaluator, pop_size, seeds, terminator, **args)


class ES(EvolutionaryComputation):
    """Evolutionary computation representing a canonical evolution strategy.
    
    This class represents an evolution strategy which uses, by 
    default, the default selection (i.e., all individuals are selected), 
    Gaussian mutation, and 'plus' replacement. It is assumed that the
    candidate solution is an iterable object of real values. 
    
    """
    def __init__(self, random):
        EvolutionaryComputation.__init__(self, random)
        self.selector = selectors.default_selection
        self.variator = variators.gaussian_mutation
        self.replacer = replacers.plus_replacement
        

class EDA(EvolutionaryComputation):
    """Evolutionary computation representing a canonical estimation of distribution algorithm.
    
    This class represents an estimation of distribution algorithm which
    uses, by default, truncation selection, estimation of distribution 
    variation, and generational replacement. It is assumed that the
    candidate solution is an iterable object of real values. 
    
    """
    def __init__(self, random):
        EvolutionaryComputation.__init__(self, random)
        self.selector = selectors.truncation_selection
        self.variator = variators.estimation_of_distribution_variation
        self.replacer = replacers.generational_replacement
        
    def evolve(self, generator, evaluator, pop_size=100, seeds=[], terminator=terminators.default_termination, **args):
        try:
            args['num_selected']
        except KeyError:
            args['num_selected'] = pop_size / 2
        try:
            args['num_offspring']
        except KeyError:
            args['num_offspring'] = pop_size
        return EvolutionaryComputation.evolve(self, generator, evaluator, pop_size, seeds, terminator, **args)


class DEA(EvolutionaryComputation):
    """Evolutionary computation representing a differential evolutionary algorithm.
    
    This class represents a differential evolutionary algorithm which uses, by 
    default, tournament selection, differential crossover, Gaussian mutation,
    and steady-state replacement. It is expected that the candidate solution 
    is an iterable object of real values. 
    
    """
    def __init__(self, random):
        EvolutionaryComputation.__init__(self, random)
        self.selector = selectors.tournament_selection
        self.variator = [variators.differential_crossover, variators.gaussian_mutation]
        self.replacer = replacers.steady_state_replacement
        
    def evolve(self, generator, evaluator, pop_size=100, seeds=[], terminator=terminators.default_termination, **args):
        try:
            args['num_selected']
        except KeyError:
            args['num_selected'] = 2
        return EvolutionaryComputation.evolve(self, generator, evaluator, pop_size, seeds, terminator, **args)


class SA(EvolutionaryComputation):
    """Evolutionary computation representing simulated annealing."""
    def __init__(self, random):
        EvolutionaryComputation.__init__(self, random)
        self.selector = selectors.default_selection
        self.variator = variators.gaussian_mutation
        self.replacer = replacers.simulated_annealing_replacement
    
    def evolve(self, generator, evaluator, pop_size=1, seeds=[], terminator=terminators.default_termination, **args):
        pop_size=1
        return EvolutionaryComputation.evolve(self, generator, evaluator, pop_size, seeds, terminator, **args)
