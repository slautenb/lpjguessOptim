"""
    This module provides the framework for making multiobjective evolutionary computations.
    
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

import numpy 

import svens_ec_constraints as ec# was from ecspy import ec,  changed 04.01.2018
from ecspy import archivers
from ecspy import selectors
from ecspy import replacers
from ecspy import terminators
import my_replacers

class constraintPareto(object):
    """Represents a Pareto multiobjective solution.
    
    With constraint handling
    
    A Pareto solution is a multiobjective value that can be compared
    to other Pareto values using Pareto preference. This means that
    a solution dominates, or is better than, another solution if it is
    better than or equal to the other solution in all objectives and
    strictly better in at least one objective.
    
    """
    def __init__(self, values=[]):
        self.values = values
        self.constraintViolations = None
    def __len__(self):
        return len(self.values)
    
    def __getitem__(self, key):
        return self.values[key]
        
    def __iter__(self):
        return iter(self.values)
    
    def __lt__(self, other):
        # Pareto comparision with constraint handling
        if (len(self.values) != len(other.values)) or (len(self.constraintViolations) != len(other.constraintViolations)) :
            raise Exception("values or constraintViolations not of equal length")
            
        else:
            not_worse = True
            strictly_better = False
            
            if any(self.constraintViolations) < 0:
                raise Exception("Negative contraint violations for candidate:\n" + str(self) ) 
            if any(other.constraintViolations) < 0:
                raise Exception("Negative contraint violations for candidate:\n" + str(other) ) 
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
                        for x, y in zip(self.values, other.values):
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
                    # return fals if self has a larger sum of violations
                    return sumcvSelf > sumcvOther
                
            #Case 3: Both solutions are both feasible
            # self = x
            # other = y
            for x, y in zip(self.values, other.values):
               
                if x > y:
                    not_worse = False
                elif y > x:
                    strictly_better = True
            return not_worse and strictly_better
        if self.values is not None and other.values is not None:
            return self.values < other.values
            



    def __le__(self, other):
        return self < other or not other < self
        
    def __gt__(self, other):
        return other < self
        
    def __ge__(self, other):
        return other < self or not self < other
    
    def __str__(self):
        return "values: " + str(self.values) + ", constraint violations " + str(self.constraintViolations)
    
    def __repr__(self):
        return "<values:  = %s, constraintViolations = %s>" % (str(self.values), str(self.constraintViolations))
    
    def __cmp__(self, other):
        # added 06.02.2018 to allow sorting of list
        if self < other:
            return(-1)
        if self > other:
            return(1)
        else:
            return(0)
 
class constraint_array_Pareto(constraintPareto):
    """Represents a Pareto multiobjective solution.
    
    A Pareto solution is a multiobjective value that can be compared
    to other Pareto values using Pareto preference. This means that
    a solution dominates, or is better than, another solution if it is
    better than or equal to the other solution in all objectives and
    strictly better in at least one objective.
    Method overwrites in addition to Pareto class the operator == 
    
    """
    def __init__(self, values=[], constraintViolations = []):
        self.values = values
        #self.epsilons = epsilons # the epsilon values in the scaled version
        self.constraintViolations= constraintViolations
    def __eq__(self, other):
        return numpy.array_equiv(self.candidate, other.candidate)
    
   
        
        
class epsilonNSGA2(ec.EvolutionaryComputation):
    """Evolutionary computation representing the nondominated sorting genetic algorithm.
    
    This class represents the nondominated sorting genetic algorithm (NSGA-II)
    of Kalyanmoy Deb et al. It uses nondominated sorting with crowding for 
    replacement, binary tournament selection to produce *population size*
    children, and a Pareto archival strategy. The remaining operators take 
    on the typical default values but they may be specified by the designer.
    
    """
    def __init__(self, random):
        ec.EvolutionaryComputation.__init__(self, random)
        self.archiver = archivers.best_archiver
        self.replacer = replacers.nsga_replacement
        self.selector = selectors.tournament_selection
    
    def evolve(self, generator, evaluator, pop_size=100, seeds=[], terminator=terminators.default_termination, **args):
        try:
            args['num_selected']
        except KeyError:
            args['num_selected'] = pop_size
        try:
            args['tourn_size']
        except KeyError:
            args['tourn_size'] = 2
        
        return ec.EvolutionaryComputation.evolve(self, generator, evaluator, pop_size, seeds, terminator, **args)

    
    

class constraintNSGA2(ec.EvolutionaryComputation):
    """Evolutionary computation representing the nondominated sorting genetic algorithm.
    
    With constraint handling
    
    This class represents the nondominated sorting genetic algorithm (NSGA-II)
    of Kalyanmoy Deb et al. It uses nondominated sorting with crowding for 
    replacement, binary tournament selection to produce *population size*
    children, and a Pareto archival strategy. The remaining operators take 
    on the typical default values but they may be specified by the designer.
    
    """
    def __init__(self, random):
        ec.EvolutionaryComputation.__init__(self, random)
        self.archiver = archivers.best_archiver
        #self.replacer = replacers.nsga_replacement
        self.replacer = myreplacers.nsga_constraint_replacement
        self.selector = selectors.tournament_selection
    
    def evolve(self, generator, evaluator, pop_size=100, seeds=[], terminator=terminators.default_termination, **args):
        try:
            args['num_selected']
        except KeyError:
            args['num_selected'] = pop_size
        try:
            args['tourn_size']
        except KeyError:
            args['tourn_size'] = 2
        
        return ec.EvolutionaryComputation.evolve(self, generator, evaluator, pop_size, seeds, terminator, **args)