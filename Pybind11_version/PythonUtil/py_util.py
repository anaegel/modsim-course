# Getting UG4 Functionality 
# These libraries were written in C++ and are exported with the Pybind extension
import os
env_ug4_root = os.environ['UG4_ROOT']

import sys
sys.path.append(env_ug4_root+'/lib/')
sys.path.append(env_ug4_root+'/bin/plugins/ug4py')

import ug4py as ug4
import pylimex as limex
import pyconvectiondiffusion as cd


def ug_assert(condition, message):
    """
    Assert: Aborts program if condition is false
    ends with costum error message
    """
    if condition:
        pass
    else:
        print(message)
        sys.exit(1) # abort


def check_subsets(dom, neededSubsets):
    """
    checks if all required subsets are contained in the SubsetHandler 
    """
    # FEHLT NOCH
    return True

class UG4Util:

    def __init__(self, dim, algebra_type="CPU", order=1):
        ug4.InitUG(dim, ug4.AlgebraType(algebra_type, order))

    # ToDo: Domain Funktionalitäten in eigene Klasse auslagern => Erben
    def create_domain(self, gridName, numRefs=None, neededSubsets=None, noIntegrityCheck=False):
        
        print(f"Loading Domain {gridName} ...")
        # Create Domain Instance
        dom = ug4.Domain2d()
        # Load Domain
        ug4.LoadDomain(dom, gridName)
        print("Domain loaded.")

        if not(noIntegrityCheck):
            print("Performing integrity check on domain ...")
            
            if ug4.CheckForUnconnectedSides(dom.grid()) == True:
                print("WARNING: unconnected sides found (see above).")
                print("NOTE: You may disable this check by passing 'True' to 'noIntegrityCheck' in 'Ug4Util.create_domain()")

            print("IntegrityCheck done.")

        # Refining the Grid
        if numRefs is None:
            numRefs = 0
        
        if numRefs > 0:
            print("Refining ...")
            refiner = ug4.GlobalDomainRefiner(dom)
            for i in range(numRefs):
                ug4.TerminateAbortedRun()
                refiner.refine()
                print(f"Refining step {i} ...")

            print("Refining done")
        
        # Check whether required subsets are present
        if neededSubsets is None:
            cond = check_subsets(dom, neededSubsets)
            msg = "Something wrong with required subsets. Aborting."
            ug_assert(cond, msg)
            
        # return created domain
        return dom

