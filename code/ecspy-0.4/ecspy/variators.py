"""
    This module provides pre-defined variators for evolutionary computations.
    
    All variator functions have the following arguments:
    
    - *random* -- the random number generator object
    - *candidates* -- the candidate solutions
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


from itertools import izip


def default_variation(random, candidates, args):
    """Return the set of candidates without variation.

    .. Arguments:
       random -- the random number generator object
       candidates -- the candidate solutions
       args -- a dictionary of keyword arguments
    
    """
    return candidates


def n_point_crossover(random, candidates, args):
    """Return the offspring of n-point crossover on the candidates.

    This function assumes that the candidate solutions are sliceable.
    It selects n random points without replacement at which to 'cut'
    the candidate solutions and recombine them.

    .. Arguments:
       random -- the random number generator object
       candidates -- the candidate solutions
       args -- a dictionary of keyword arguments

    Optional keyword arguments in args:
    
    - *crossover_rate* -- the rate at which crossover is performed 
      (default 1.0)
    - *num_crossover_points* -- the n crossover points used (default 1)
    
    """
    try:
        crossover_rate = args['crossover_rate']
    except KeyError:
        crossover_rate = 1.0
        args['crossover_rate'] = crossover_rate
    try:
        num_crossover_points = args['num_crossover_points']
    except KeyError:
        num_crossover_points = 1
        args['num_crossover_points'] = num_crossover_points
    cand = list(candidates)
    if len(cand) % 2 == 1:
        cand = cand[:-1]
    random.shuffle(cand)
    moms = cand[::2]
    dads = cand[1::2]
    children = []
    for mom, dad in izip(moms, dads):
        if random.random() < crossover_rate:
            bro = list(mom)
            sis = list(dad)
            num_cuts = min(len(mom)-1, num_crossover_points)
            cut_points = random.sample(range(1, len(mom)), num_cuts)
            cut_points.sort()
            for p in cut_points:
                bro[p:] = dad[p:]
                sis[:p] = mom[:p]
            children.append(bro)
            children.append(sis)
        else:
            children.append(mom)
            children.append(dad)            
    return children


def uniform_crossover(random, candidates, args):
    """Return the offspring of uniform crossover on the candidates.

    This function assumes that the candidate solutions are iterable.
    It chooses every odd candidate as a 'mom' and every even as a 'dad'.
    For each mom-dad pair, two offspring are produced. For each element
    of the parents, a biased coin is flipped to determine whether the 
    first offspring gets the 'mom' or the 'dad' element. An optional
    keyword argument in args, pux_bias, determines the bias.

    .. Arguments:
       random -- the random number generator object
       candidates -- the candidate solutions
       args -- a dictionary of keyword arguments

    Optional keyword arguments in args:
    
    - *crossover_rate* -- the rate at which crossover is performed 
      (default 1.0)
    - *pux_bias* -- the bias toward the first candidate in the crossover 
      (default 0.5)
    
    """
    try:
        pux_bias = args['pux_bias']
    except KeyError:
        pux_bias = 0.5
        args['pux_bias'] = pux_bias
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
            bro = []
            sis = []
            for m, d in izip(mom, dad):
                if random.random() < pux_bias:
                    bro.append(m)
                    sis.append(d)
                else:
                    bro.append(d)
                    sis.append(m)
            children.append(bro)
            children.append(sis)
        else:
            children.append(mom)
            children.append(dad)
    return children

    
def blend_crossover(random, candidates, args):
    """Return the offspring of blend crossover on the candidates.

    This function assumes that the candidate solutions are iterable
    and composed of values on which arithmetic operations are defined.
    It performs blend crossover, which is similar to a generalized 
    averaging of the candidate elements.

    .. Arguments:
       random -- the random number generator object
       candidates -- the candidate solutions
       args -- a dictionary of keyword arguments

    Optional keyword arguments in args:
    
    - *crossover_rate* -- the rate at which crossover is performed 
      (default 1.0)
    - *blx_alpha* -- the blending rate (default 0.1)
    - *lower_bound* -- the lower bounds of the chromosome elements (default 0)
    - *upper_bound* -- the upper bounds of the chromosome elements (default 1)
    
    The lower and upper bounds can either be single values, which will
    be applied to all elements of a chromosome, or lists of values of 
    the same length as the chromosome.
    
    """
    try:
        blx_alpha = args['blx_alpha']
    except KeyError:
        blx_alpha = 0.1
    try:
        crossover_rate = args['crossover_rate']
    except KeyError:
        crossover_rate = 1.0
    try:
        lower_bound = args['lower_bound']
    except KeyError:
        lower_bound = 0
    try:
        upper_bound = args['upper_bound']
    except KeyError:
        upper_bound = 1
    
    try:
        iter(lower_bound)
    except TypeError:
        clen = max([len(x) for x in candidates])
        lower_bound = [lower_bound] * clen
        
    try:
        iter(upper_bound)
    except TypeError:
        clen = max([len(x) for x in candidates])
        upper_bound = [upper_bound] * clen
        
    cand = list(candidates)
    if len(cand) % 2 == 1:
        cand = cand[:-1]
    random.shuffle(cand)
    moms = cand[::2]
    dads = cand[1::2]
    children = []
    for mom, dad in izip(moms, dads):
        if random.random() < crossover_rate:
            bro = []
            sis = []
            for index, (m, d) in enumerate(izip(mom, dad)):
                smallest = min(m, d)
                largest = max(m, d)
                delta = blx_alpha * (largest - smallest)
                bro_val = smallest - delta + random.random() * (largest - smallest + 2 * delta)
                sis_val = smallest - delta + random.random() * (largest - smallest + 2 * delta)
                bro_val = max(min(bro_val, upper_bound[index]), lower_bound[index])
                sis_val = max(min(sis_val, upper_bound[index]), lower_bound[index])
                bro.append(bro_val)
                sis.append(sis_val)
            children.append(bro)
            children.append(sis)
        else:
            children.append(mom)
            children.append(dad)
    return children

    
def differential_crossover(random, candidates, args):
    """Return the offspring of differential crossover on the candidates.

    This function assumes that the candidate solutions are iterable
    and composed of values on which arithmetic operations are defined.
    It performs differential crossover, which is similar to the update
    rule used in particle swarm optimization.

    .. Arguments:
       random -- the random number generator object
       candidates -- the candidate solutions
       args -- a dictionary of keyword arguments

    Optional keyword arguments in args:
    
    - *crossover_rate* -- the rate at which crossover is performed 
      (default 1.0)
    - *differential_phi* -- the amount of random change in the crossover 
      (default 0.1)
    - *lower_bound* -- the lower bounds of the chromosome elements (default 0)
    - *upper_bound* -- the upper bounds of the chromosome elements (default 1)
    
    The lower and upper bounds can either be single values, which will
    be applied to all elements of a chromosome, or lists of values of 
    the same length as the chromosome.
    
    """
    try:
        differential_phi = args['differential_phi']
    except KeyError:
        differential_phi = 0.1
    try:
        crossover_rate = args['crossover_rate']
    except KeyError:
        crossover_rate = 1.0
    try:
        lower_bound = args['lower_bound']
    except KeyError:
        lower_bound = 0
    try:
        upper_bound = args['upper_bound']
    except KeyError:
        upper_bound = 1
    
    try:
        iter(lower_bound)
    except TypeError:
        clen = max([len(x) for x in candidates])
        lower_bound = [lower_bound] * clen
        
    try:
        iter(upper_bound)
    except TypeError:
        clen = max([len(x) for x in candidates])
        upper_bound = [upper_bound] * clen
        
    # Be careful shuffling the candidates so that 
    # we will know which is better. 
    cand = list(candidates)
    if len(cand) % 2 == 1:
        cand = cand[:-1]
    cand_pair = [(pos, ind) for pos, ind in enumerate(cand)]
    random.shuffle(cand_pair)
    moms = cand_pair[::2]
    dads = cand_pair[1::2]
    children = []
    for mom, dad in izip(moms, dads):
        if random.random() < crossover_rate:
            bro = []
            sis = []
            for index, (m, d) in enumerate(izip(mom[1], dad[1])):
                if mom[0] > dad[0]:
                    negpos = 1
                    val = d
                else:
                    negpos = -1
                    val = m
                bro_val = val + differential_phi * random.random() * negpos * (m - d)
                sis_val = val + differential_phi * random.random() * negpos * (m - d)
                bro_val = max(min(bro_val, upper_bound[index]), lower_bound[index])
                sis_val = max(min(sis_val, upper_bound[index]), lower_bound[index])
                bro.append(bro_val)
                sis.append(sis_val)
            children.append(bro)
            children.append(sis)
        else:
            children.append(mom[1])
            children.append(dad[1])
    return children
    
    
def gaussian_mutation(random, candidates, args):
    """Return the mutants created by Gaussian mutation on the candidates.

    This function assumes that the candidate solutions are indexable
    and numeric. It performs Gaussian mutation.

    .. Arguments:
       random -- the random number generator object
       candidates -- the candidate solutions
       args -- a dictionary of keyword arguments

    Optional keyword arguments in args:
    
    - *mutation_rate* -- the rate at which mutation is performed (default 0.1)
    - *mutation_range* -- the variance used in the Gaussian function 
      (default 1.0)
    - *lower_bound* -- the lower bounds of the chromosome elements (default 0)
    - *upper_bound* -- the upper bounds of the chromosome elements (default 1)
    
    The lower and upper bounds can either be single values, which will
    be applied to all elements of a chromosome, or lists of values of 
    the same length as the chromosome.
    
    """
    try:
        mut_rate = args['mutation_rate']
    except KeyError:
        mut_rate = 0.1
        args['mutation_rate'] = mut_rate
    try:
        mut_range = args['mutation_range']
    except KeyError:
        mut_range = 1.0
        args['mutation_range'] = mut_range
    try:
        lower_bound = args['lower_bound']
    except KeyError:
        lower_bound = 0
        args['lower_bound'] = lower_bound
    try:
        upper_bound = args['upper_bound']
    except KeyError:
        upper_bound = 1
        args['upper_bound'] = upper_bound
        
    try:
        iter(lower_bound)
    except TypeError:
        clen = max([len(x) for x in candidates])
        lower_bound = [lower_bound] * clen
        
    try:
        iter(upper_bound)
    except TypeError:
        clen = max([len(x) for x in candidates])
        upper_bound = [upper_bound] * clen
        
    cs_copy = list(candidates)
    for i, cs in enumerate(cs_copy):
        for j, c in enumerate(cs):
            if random.random() < mut_rate:
                c += random.gauss(0, mut_range) * (upper_bound[j] - lower_bound[j])
                c = max(min(c, upper_bound[j]), lower_bound[j])
                cs_copy[i][j] = c
    return cs_copy


def bit_flip_mutation(random, candidates, args):
    """Return the mutants produced by bit-flip mutation on the candidates.

    This function assumes that the candidate solutions are binary values.
    It performs bit-flip mutation. 

    .. Arguments:
       random -- the random number generator object
       candidates -- the candidate solutions
       args -- a dictionary of keyword arguments

    Optional keyword arguments in args:
    
    *mutation_rate* -- the rate at which mutation is performed (default 0.1)
    
    The mutation rate is applied on a bit by bit basis.
    
    """
    try:
        rate = args['mutation_rate']
    except KeyError:
        rate = 0.1
        args['mutation_rate'] = rate
    cs_copy = list(candidates)
    for i, cs in enumerate(cs_copy):
        if len(cs) == len([x for x in cs if x in [0, 1]]):
            for j, c in enumerate(cs):
                if random.random() < rate:
                    cs_copy[i][j] = (c + 1) % 2
    return cs_copy

    
def estimation_of_distribution_variation(random, candidates, args):
    """Return the offspring produced using estimation of distribution.

    This function assumes that the candidate solutions are iterable 
    objects containing real values. It creates a statistical model 
    based on the set of candidates. The offspring are then generated 
    from this model.

    .. Arguments:
       random -- the random number generator object
       candidates -- the candidate solutions
       args -- a dictionary of keyword arguments
    
    Optional keyword arguments in args:
    
    - *num_offspring* -- the number of offspring to create (default 1)
    - *lower_bound* -- the lower bounds of the chromosome elements (default 0)
    - *upper_bound* -- the upper bounds of the chromosome elements (default 1)
    
    The lower and upper bounds can either be single values, which will
    be applied to all elements of a chromosome, or lists of values of 
    the same length as the chromosome.
    
    """
    try:
        num_offspring = args['num_offspring']
    except KeyError:
        num_offspring = 1
        args['num_offspring'] = num_offspring
    try:
        lower_bound = args['lower_bound']
    except KeyError:
        lower_bound = 0
        args['lower_bound'] = lower_bound
    try:
        upper_bound = args['upper_bound']
    except KeyError:
        upper_bound = 1
        args['upper_bound'] = upper_bound
    
    try:
        iter(lower_bound)
    except TypeError:
        clen = max([len(x) for x in candidates])
        lower_bound = [lower_bound] * clen
        
    try:
        iter(upper_bound)
    except TypeError:
        clen = max([len(x) for x in candidates])
        upper_bound = [upper_bound] * clen
        
    cs_copy = list(candidates)
    num_genes = max([len(x) for x in cs_copy])
    genes = [[x[i] for x in cs_copy] for i in xrange(num_genes)] 
    mean = [float(sum(x)) / float(len(x)) for x in genes]
    stdev = [max(sum([(x - m)**2 for x in g]) / float(len(g) - 1), 0.001) for g, m in zip(genes, mean)]
    offspring = []
    for _ in xrange(num_offspring):
        child = []
        for m, s, hi, lo in zip(mean, stdev, upper_bound, lower_bound):
            value = m + random.gauss(0, s);
            value = max(min(value, hi), lo)
            child.append(value)
        offspring.append(child)
        
    return offspring
