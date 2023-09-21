"""
    This module provides pre-defined archivers for evolutionary computations.
    
    All archiver functions have the following arguments:
    
    - *random* -- the random number generator object
    - *population* -- the population of Individuals
    - *archive* -- the current archive of individuals
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


def default_archiver(random, population, archive, args):
    """Archive the current population."""
    new_archive = []
    for ind in population:
        new_archive.append(ind)
    return new_archive
    

def best_archiver(random, population, archive, args):
    """Archive only the best individual(s).
    
    This function archives the best solutions and removes inferior ones.
    If the comparison operators have been overloaded to define Pareto
    preference (as in the ``Pareto`` class), then this archiver will form 
    a Pareto archive.
    
    """
    new_archive = archive
    for ind in population:
        if len(new_archive) == 0:
            new_archive.append(ind)
        else:
            should_remove = []
            should_add = True
            for a in new_archive:
                if ind == a:
                    should_add = False
                    break
                elif ind < a:
                    should_add = False
                elif ind > a:
                    should_remove.append(a)
            for r in should_remove:
                new_archive.remove(r)
            if should_add:
                new_archive.append(ind)
    return new_archive

    
