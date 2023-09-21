"""
    This module provides pre-defined observers for evolutionary computations.
    
    All terminator functions have the following arguments:
    
    - *population* -- the population of Individuals
    - *num_generations* -- the number of elapsed generations
    - *num_evaluations* -- the number of candidate solution evaluations
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

import time


def default_observer(population, num_generations, num_evaluations, args):
    """Do nothing."""    
    pass
    
    
def screen_observer(population, num_generations, num_evaluations, args):
    """Print the output of the EC to the screen.
    
    This function displays the results of the evolutionary computation
    to the screen. The output includes the generation number, the current
    number of evaluations, the average fitness, the maximum fitness, and 
    the full population.
    
    .. Arguments:
       population -- the population of Individuals
       num_generations -- the number of elapsed generations
       num_evaluations -- the number of candidate solution evaluations
       args -- a dictionary of keyword arguments
    
    """
    print('Generation: %d' % num_generations)
    print('Evaluations: %d' % num_evaluations)
    avg_fit = sum([x.fitness for x in population]) / float(len(population))
    print('Average Fitness: %0.5f     Maximum Fitness: %0.5f' % (avg_fit, population[0].fitness))
    for ind in population:
        print(str(ind))
    print('')

    
def file_observer(population, num_generations, num_evaluations, args):
    """Print the output of the EC to a file.
    
    This function saves the results of the evolutionary computation
    to a file. The output includes the generation number, the current
    number of evaluations, the average fitness, the maximum fitness, 
    and the full population. The default action for the file is to
    create a new file called 'ecspy-observer-file-<timestamp>' in
    which to write the information.
    
    .. Arguments:
       population -- the population of Individuals
       num_generations -- the number of elapsed generations
       num_evaluations -- the number of candidate solution evaluations
       args -- a dictionary of keyword arguments

    Optional keyword arguments in args:
    *observer_file* -- a file object (default; see text)
    
    """
    try:
        observer_file = args['observer_file']
    except KeyError:
        filename = 'ecspy-observer-file-' + str(time.time())
        observer_file = open(filename, 'w')
        args['observer_file'] = observer_file

    observer_file.write('Generation: %d \n' % num_generations)
    observer_file.write('Evaluations: %d \n' % num_evaluations)
    avg_fit = sum([x.fitness for x in population]) / float(len(population))
    observer_file.write('Average Fitness: %0.5f     Maximum Fitness: %0.5f \n' % (avg_fit, population[0].fitness))
    for ind in population:
        observer_file.write(str(ind) + '\n')
    observer_file.write('\n')
    

def archive_observer(population, num_generations, num_evaluations, args):
    """Print the current archive to the screen."""
    try:
        archive = args['_archive']
    except KeyError:
        archive = []
    for a in archive:
        print(a.candidate)

        
def plot_observer(population, num_generations, num_evaluations, args):    
    """Plot the output of the EC as a graph.
    
    This function plots the performance of the EC as a line graph 
    using the pylab library (matplotlib). The graph consists of a 
    blue line representing the best fitness, a green line representing
    the average fitness, and a red line representing the median fitness.
    It modifies the keyword arguments variable 'args' by including an
    entry called 'plot_data'.
    
    If this observer is used, the calling script should also import
    the pylab library and should end the script with 
    
    pylab.show()
    
    Otherwise, the program may generate a runtime error.
    
    .. Arguments:
       population -- the population of Individuals
       num_generations -- the number of elapsed generations
       num_evaluations -- the number of candidate solution evaluations
       args -- a dictionary of keyword arguments
    
    """
    # Import the necessary libraries here. Otherwise, they would have to be
    # installed even if this function is not called.
    import pylab
    import numpy
    
    best_fitness = population[0].fitness
    average_fitness = sum([x.fitness for x in population]) / float(len(population))
    median_fitness = population[len(population)/2].fitness
    colors = ['blue', 'green', 'red']
    labels = ['best', 'average', 'median']
    data = []
    if num_generations == 0:
        pylab.ion()
        data = [[num_evaluations], [best_fitness], [average_fitness], [median_fitness]]
        lines = []
        for i in xrange(3):
            line, = pylab.plot(data[0], data[i+1], color=colors[i], label=labels[i])
            lines.append(line)
        # Add the legend when the first data is added.
        pylab.legend(loc='lower right')
        args['plot_data'] = data
        args['plot_lines'] = lines
    else:
        data = args['plot_data']
        data[0].append(num_evaluations)
        data[1].append(best_fitness)
        data[2].append(average_fitness)
        data[3].append(median_fitness)
        lines = args['plot_lines']
        for i, line in enumerate(lines):
            line.set_xdata(numpy.array(data[0]))
            line.set_ydata(numpy.array(data[i+1]))
        args['plot_data'] = data
        args['plot_lines'] = lines
    ymin = min(min(data[1]), min(data[2]), min(data[3]))
    ymax = max(max(data[1]), max(data[2]), max(data[3]))
    yrange = ymax - ymin
    pylab.xlim((0, num_evaluations))
    pylab.ylim((ymin - 0.1*yrange, ymax + 0.1*yrange))