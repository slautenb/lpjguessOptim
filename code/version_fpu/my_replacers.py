"""
    This module provides pre-defined replacers for evolutionary computations.
    
    All replacer functions have the following arguments:
    
    - *random* -- the random number generator object
    - *population* -- the population of Individuals
    - *parents* -- the list of parent individuals
    - *offspring* -- the list of offspring individuals
    - *args* -- a dictionary of keyword arguments

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
import math
import numpy

def default_replacement(random, population, parents, offspring, args):
    """Replaces entire population with offspring.
    
    .. Arguments:
       random -- the random number generator object
       population -- the population of Individuals
       parents -- the list of parent individuals
       offspring -- the list of offspring individuals
       args -- a dictionary of keyword arguments
    
    """
    return offspring

    
def truncation_replacement(random, population, parents, offspring, args):
    """Replaces population with the best of the population and offspring.
    
    This function performs truncation replacement, which means that
    the entire existing population is replaced by the best from among
    the current population and offspring, keeping the existing population
    size fixed. This is similar to so called 'plus' replacement in the 
    evolution strategies literature, except that 'plus' replacement 
    considers only parents and offspring for survival. However, if the
    entire population are parents (which is often the case in ES), then
    truncation replacement and plus-replacement are equivalent approaches.
    
    .. Arguments:
       random -- the random number generator object
       population -- the population of Individuals
       parents -- the list of parent individuals
       offspring -- the list of offspring individuals
       args -- a dictionary of keyword arguments
    
    """
    pool = list(population)
    pool.extend(list(offspring))
    pool.sort(key=lambda x: x.fitness, reverse=True)
    return pool[:len(population)]

    
def steady_state_replacement(random, population, parents, offspring, args):
    """Performs steady-state replacement for the offspring.
    
    This function performs steady-state replacement, which means that
    the offspring replace the least fit individuals in the existing
    population, even if those offspring are less fit than the individuals
    that they replace.
    
    .. Arguments:
       random -- the random number generator object
       population -- the population of Individuals
       parents -- the list of parent individuals
       offspring -- the list of offspring individuals
       args -- a dictionary of keyword arguments
    
    """
    off = list(offspring)
    pop = list(population)
    pop.sort(key=lambda x: x.fitness)
    num_to_replace = min(len(off), len(pop))
    pop[:num_to_replace] = off[:num_to_replace]
    return pop


def generational_replacement(random, population, parents, offspring, args):
    """Performs generational replacement with optional weak elitism.
    
    This function performs generational replacement, which means that
    the entire existing population is replaced by the offspring,
    truncating to the population size if the number of offspring is 
    larger. Weak elitism may also be specified through the num_elites
    keyword argument in args. If this is used, the best num_elites
    individuals in the current population are allowed to survive if
    they are better than the worst num_elites offspring.
    
    .. Arguments:
       random -- the random number generator object
       population -- the population of Individuals
       parents -- the list of parent individuals
       offspring -- the list of offspring individuals
       args -- a dictionary of keyword arguments

    Optional keyword arguments in args:
    
    *num_elites* -- number of elites to consider (default 0)
    
    """
    try:
        num_elites = args['num_elites']
    except KeyError:
        num_elites = 0
        args['num_elites'] = num_elites
    off = list(offspring)
    pop = list(population)
    pop.sort(key=lambda x: x.fitness, reverse=True)
    off.extend(pop[:num_elites])
    off.sort(key=lambda x: x.fitness, reverse=True)
    survivors = off[:len(population)]
    return survivors


def random_replacement(random, population, parents, offspring, args):
    """Performs random replacement with optional weak elitism.
    
    This function performs random replacement, which means that
    the offspring replace random members of the population, keeping
    the population size constant. Weak elitism may also be specified 
    through the num_elites keyword argument in args. If this is used, 
    the best num_elites individuals in the current population are 
    allowed to survive if they are better than the worst num_elites 
    offspring.
    
    .. Arguments:
       random -- the random number generator object
       population -- the population of Individuals
       parents -- the list of parent individuals
       offspring -- the list of offspring individuals
       args -- a dictionary of keyword arguments

    Optional keyword arguments in args:
    
    *num_elites* -- number of elites to consider (default 0)
    
    """
    try:
        num_elites = args['num_elites']
    except KeyError:
        num_elites = 0
        args['num_elites'] = num_elites
    off = list(offspring)
    pop = list(population)
    pop.sort(key=lambda x: x.fitness, reverse=True)
    num_to_replace = min(len(off), len(pop) - num_elites) 
    valid_indices = range(num_elites, len(pop))
    rep_index = random.sample(valid_indices, num_to_replace)
    for i in xrange(len(off)):
        pop[rep_index[i]] = off[i]
    return pop


def plus_replacement(random, population, parents, offspring, args):
    """Performs 'plus' replacement with optional adaptive mutation.
    
    This function performs 'plus' replacement, which means that
    the entire existing population is replaced by the best
    population-many elements from the combined set of parents and 
    offspring. Adaptive mutation can also be specified here through 
    the use_one_fifth_rule keyword argument in args. This adaptive
    mutation attempts to keep the rate of viable offspring near 20%.
    If the number of successful offspring is below 20%, the mutation
    rate is reduced by 20% (to allow more exploitation). If the 
    number of successful offspring is above 20%, the mutation rate
    is increased by 20% (to allow more exploration).
    
    .. Arguments:
       random -- the random number generator object
       population -- the population of Individuals
       parents -- the list of parent individuals
       offspring -- the list of offspring individuals
       args -- a dictionary of keyword arguments

    Optional keyword arguments in args:
    
    *use_one_fifth_rule* -- whether the 1/5 rule should be used (default False)
    
    """
    try:
        use_one_fifth_rule = args['use_one_fifth_rule']
    except KeyError:
        use_one_fifth_rule = False
        args['use_one_fifth_rule'] = use_one_fifth_rule
    pool = list(parents)
    pool.extend(list(offspring))
    pool.sort(key=lambda x: x.fitness, reverse=True)
    survivors = pool[:len(population)]
    if use_one_fifth_rule:
        count = len([x for x in offspring if x in survivors])
        rate = count / float(len(offspring))
        if rate < 0.2:
            try:
                args['mutation_rate'] = args['mutation_rate'] * 0.8
            except KeyError:
                args['use_one_fifth_rule'] = False
        elif rate > 0.2:
            try:
                args['mutation_rate'] = args['mutation_rate'] * 1.2
            except KeyError:
                args['use_one_fifth_rule'] = False            
    return survivors


def comma_replacement(random, population, parents, offspring, args):
    """Performs 'comma' replacement.
    
    This function performs 'comma' replacement, which means that
    the entire existing population is replaced by the best
    population-many elements from the offspring. This function
    makes the assumption that the size of the offspring is at 
    least as large as the original population. Otherwise, the
    population size will not be constant.
    
    .. Arguments:
       random -- the random number generator object
       population -- the population of Individuals
       parents -- the list of parent individuals
       offspring -- the list of offspring individuals
       args -- a dictionary of keyword arguments
    
    """
    pool = list(offspring)
    pool.sort(key=lambda x: x.fitness, reverse=True)
    survivors = pool[:len(population)]
    return survivors

    
def simulated_annealing_replacement(random, population, parents, offspring, args):
    """Replaces population using the simulated annealing schedule.
    
    This function performs simulated annealing replacement based
    on a temperature and a cooling rate. These can be specified
    by the keyword arguments 'temperature', which should be the
    initial temperature, and 'cooling_rate', which should be the
    coefficient by which the temperature is reduced. If these
    keyword arguments are not present, then the function will
    attempt to base the cooling schedule either on the ratio of 
    evaluations to the maximum allowed evaluations or on the 
    ratio of generations to the maximum allowed generations. 
    Each of these rations is of the form (max - current)/max
    so that the cooling schedule moves smoothly from 1 to 0.
    
    .. Arguments:
       random -- the random number generator object
       population -- the population of Individuals
       parents -- the list of parent individuals
       offspring -- the list of offspring individuals
       args -- a dictionary of keyword arguments
    
    Optional keyword arguments in args:    
    
    - *temperature* -- the initial temperature
    - *cooling_rate* -- a real-valued coefficient in the range (0, 1) 
                      by which the temperature should be reduced 
    
    """
    try:
        temp = args['temperature']
        cooling_rate = args['cooling_rate']
        temp = temp * cooling_rate
        args['temperature'] = temp
    except KeyError:
        try:
            num_evals = args['_num_evaluations']
            max_evals = args['max_evaluations']
            temp = float(max_evals - num_evals) / float(max_evals)
        except KeyError:
            num_gens = args['_num_generations']
            max_gens = args['max_generations']
            temp = 1 - float(max_gens - num_gens) / float(max_gens)
        
    new_pop = []
    for p, o in zip(parents, offspring):
        if o.fitness >= p.fitness:
            new_pop.append(o)
        elif random.random() < math.exp(-(p.fitness - o.fitness) / temp):
            new_pop.append(o)
        else:
            new_pop.append(p)
            
    return new_pop

    
def nsga_replacement(random, population, parents, offspring, args):
    """Replaces population using the non-dominated sorting technique from NSGA-II.
    
    .. Arguments:
       random -- the random number generator object
       population -- the population of Individuals
       parents -- the list of parent individuals
       offspring -- the list of offspring individuals
       args -- a dictionary of keyword arguments
    
    """
    survivors = []
    combined = population[:]
    combined.extend(offspring[:])
    
    # Perform the non-dominated sorting to determine the fronts.
    fronts = []
    pop = set(range(len(combined)))
    while len(pop) > 0:
        front = []
        for p in pop:
            dominated = False
            for q in pop:
                if combined[p].fitness < combined[q].fitness:
                    dominated = True
                    break
            if not dominated:
                front.append(p)
        fronts.append([dict(individual=combined[f], index=f) for f in front])
        pop = pop - set(front)
    
    # Go through each front and add all the elements until doing so
    # would put you above the population limit. At that point, fall
    # back to the crowding distance to determine who to put into the
    # next population. Individuals with higher crowding distances
    # (i.e., more distance between neighbors) are preferred.
    for i, front in enumerate(fronts):
        if len(survivors) + len(front) > len(population):
            # Determine the crowding distance.
            distance = [0 for _ in xrange(len(combined))]
            individuals = front[:]
            num_individuals = len(individuals)
            num_objectives = len(individuals[0]['individual'].fitness)
            for obj in xrange(num_objectives):
                individuals.sort(key=lambda x: x['individual'].fitness[obj])
                distance[individuals[0]['index']] = float('inf')
                distance[individuals[-1]['index']] = float('inf')
                for i in xrange(1, num_individuals-1):
                    distance[individuals[i]['index']] = distance[individuals[i]['index']] + (individuals[i+1]['individual'].fitness[obj] - individuals[i-1]['individual'].fitness[obj])
                
            crowd = [dict(dist=distance[f['index']], index=f['index']) for f in front]
            crowd.sort(key=lambda x: x['dist'], reverse=True)
            last_rank = [combined[c['index']] for c in crowd]
            r = 0
            num_added = 0
            num_left_to_add = len(population) - len(survivors)
            while r < len(last_rank) and num_added < num_left_to_add:
                if last_rank[r] not in survivors:
                    survivors.append(last_rank[r])
                    num_added += 1
                r += 1
            # If we've filled out our survivor list, then stop.
            # Otherwise, process the next front in the list.
            if len(survivors) == len(population):
                break
        else:
            for f in front:
                if f['individual'] not in survivors:
                    survivors.append(f['individual'])
    return survivors
    

##########
# to do
# don't compare the fitness but use the _lt_ method of the Individual
##########
def nsga_constraint_replacement(random, population, parents, offspring, args):
    """Replaces population using the non-dominated sorting technique from NSGA-II. Including constraint handling
    
    .. Arguments:
       random -- the random number generator object
       population -- the population of Individuals
       parents -- the list of parent individuals
       offspring -- the list of offspring individuals
       args -- a dictionary of keyword arguments
    
    """
    survivors = []
    combined = population[:]
    combined.extend(offspring[:])
    
    # Perform the non-dominated sorting to determine the fronts.
    fronts = []
    pop = set(range(len(combined)))
    while len(pop) > 0:
        front = []
        for p in pop:
            dominated = False
            for q in pop:
                if combined[p] < combined[q]:
                    # compared based on the overladed < operator not on fitness
                    # was before: 
                    #combined[p].fitness < combined[q].fitness:
                    # changed 28.09.18
                    dominated = True
                    break
            if not dominated:
                #print("individual added to front")
                front.append(p)
        fronts.append([dict(individual=combined[f], index=f) for f in front])
        pop = pop - set(front)
    
    # Go through each front and add all the elements until doing so
    # would put you above the population limit. At that point, fall
    # back to the crowding distance to determine who to put into the
    # next population. Individuals with higher crowding distances
    # (i.e., more distance between neighbors) are preferred.
    for i, front in enumerate(fronts):
        if len(survivors) + len(front) > len(population):
            # Determine the crowding distance.
            distance = [0 for _ in xrange(len(combined))]
            individuals = front[:]
            num_individuals = len(individuals)
            num_objectives = len(individuals[0]['individual'].fitness)
            for obj in xrange(num_objectives):
                # changed 28.09.18 I thought to use __cmp__ function
                # however, since the crowding distance comes only into
                # play after individuals have been assigned to fronts
                # i.e. after consideration of constraints
                # it should be fine to just use the fitness here instead of
                 #individuals.sort()
                individuals.sort(key=lambda x: x['individual'].fitness[obj])
               
                distance[individuals[0]['index']] = float('inf')
                distance[individuals[-1]['index']] = float('inf')
                for i in xrange(1, num_individuals-1):
                    distance[individuals[i]['index']] = distance[individuals[i]['index']] + (individuals[i+1]['individual'].fitness[obj] - individuals[i-1]['individual'].fitness[obj])
                
            crowd = [dict(dist=distance[f['index']], index=f['index']) for f in front]
            crowd.sort(key=lambda x: x['dist'], reverse=True)
            last_rank = [combined[c['index']] for c in crowd]
            r = 0
            num_added = 0
            num_left_to_add = len(population) - len(survivors)
            while r < len(last_rank) and num_added < num_left_to_add:
                if last_rank[r] not in survivors:
                    #print("individual added to survivors based on crowding distance")
                    survivors.append(last_rank[r])
                    num_added += 1
                r += 1
            # If we've filled out our survivor list, then stop.
            # Otherwise, process the next front in the list.
            if len(survivors) == len(population):
                break
        else:
            for f in front:
                if f['individual'] not in survivors:
                    #print("individual added to survivors based on front, no crowding distance considered")
                    survivors.append(f['individual'])
    return survivors
    
