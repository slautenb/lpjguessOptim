"""
    This module provides pre-defined evaluators for evolutionary computations.

    All evaluator functions have the following arguments:
    
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




def parallel_evaluation(candidates, args):
    """Evaluate the candidates in parallel.

    This function allows parallel evaluation of candidate solutions.
    It uses the Parallel Python (pp) library to accomplish the 
    parallelization. This library must already be installed in order
    to use this function. The function assigns the evaluation of each
    candidate to its own job, all of which are then distributed to the
    available processing units.
    
    .. Arguments:
       candidates -- the candidate solutions
       args -- a dictionary of keyword arguments

    Required keyword arguments in args:
    
    *serial_evaluator* -- the actual evaluation function, which should take a 
    single argument representing a candidate solution (required)
    
    Optional keyword arguments in args:
    
    - *serial_dependencies* -- tuple of functional dependencies of the serial 
      evaluator (default ())
    - *serial_modules* -- tuple of modules that must be imported for the 
      functional dependencies (default ())
    - *parallel_servers* -- tuple of servers (on a cluster) that will be used 
      for parallel processing (default ("*",))
    
    """
    # Import the necessary library here. Otherwise, it would have to be
    # installed even if this function is not called.
    import pp
    
    try:
        serial_eval = args['serial_evaluator']
    except KeyError:
        return [-float('inf') for c in candidates]
    try:
        serial_depend = args['serial_dependencies']
    except KeyError:
        serial_depend = ()
    try:
        serial_mod = args['serial_modules']
    except KeyError:
        serial_mod = ()
    try:
        job_server = args['_job_server']
    except KeyError:
        try:
            parallel_servers = args['parallel_servers']
        except KeyError:
            parallel_servers = ("*",)
        job_server = pp.Server(ppservers=parallel_servers)
        args['_job_server'] = job_server
        
    func_template = pp.Template(job_server, serial_eval, serial_depend, serial_mod)
    jobs = [func_template.submit(cand) for cand in candidates]
    
    results = []
    for job in jobs:
        results.append(job())
    
    fitness = []
    for result in results:
        fitness.append(result)
    return fitness
    