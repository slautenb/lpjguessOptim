``ecspy`` -- A framework for creating evolutionary computations in Python.
--------------------------------------------------------------------------

ECsPy (pronounced "easy as pie") is a free, open source framework for 
creating evolutionary computations in Python. Additionally, ECsPy 
provides an easy-to-use canonical genetic algorithm (GA), evolution 
strategy (ES), estimation of distribution algorithm (EDA), differential 
evolution algorithm (DEA), and particle swarm optimizer (PSO) for users 
who don't need much customization.

  
Requirements
============

  * Requires at least Python 2.6.
  * Matplotlib is required if the line plot observer is used.
  * Parallel Python (pp) is required if parallel evaluation is used.



License
=======

This package is distributed under the GNU General Public License 
version 3.0 (GPLv3). This license can be found online at
http://www.opensource.org/licenses/gpl-3.0.html.

  
Package Structure
=================
  
ECsPy consists of the following modules:

  * ec.py -- provides the basic framework for the EvolutionEngine and specific ECs
  
  * evaluators.py -- defines useful evaluation schemes, such as parallel evaluation
  
  * migrators.py -- defines a basic default migration which does nothing
             
  * observers.py -- defines a few built-in observers, including screen, file, and plotting observers
  
  * replacers.py -- defines standard replacement schemes such as generational and steady-state replacement
                    
  * selectors.py -- defines standard selectors (e.g., tournament)
  
  * swarm.py -- provides a basic particle swarm optimizer
  
  * terminators.py -- defines standard terminators (e.g., exceeding a maximum number of generations)
                      
  * variators.py -- defines standard variators (crossover and mutation schemes such as n-point crossover)


Resources
=========

  * Homepage: http://ecspy.googlecode.com
  
  * Email: aaron.lee.garrett@gmail.com
  
