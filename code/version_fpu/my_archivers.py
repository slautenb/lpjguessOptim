import numpy

def best_archiver_array(random, population, archive, args):
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
                #print type(a)
                
                if ind == a:
                    #if numpy.array_equiv(ind.candidate,a.candidate):
                    #print("equal")
                    #print("a.fitness" + str(a.fitness))
                    #print("ind.fitness" + str(ind.fitness))
                    #print("a.constraintViolations" + str(a.constraintViolations))
                    #print("ind.constraintViolations" + str(ind.constraintViolations))
                    should_add = False
                    break
                elif ind < a:
                    #print("<")
                    should_add = False
                elif ind > a:
                    #print(">")
                    #print("a.fitness" + str(a.fitness))
                    #print("ind.fitness" + str(ind.fitness))
                    #print("a.constraintViolations" + str(a.constraintViolations))
                    #print("ind.constraintViolations" + str(ind.constraintViolations))
                    should_remove.append(a)
            for r in should_remove:
                #print("individual removed from archive")
                new_archive.remove(r)
            if should_add:
                #print("individual added to archive")
                new_archive.append(ind)
    return new_archive
