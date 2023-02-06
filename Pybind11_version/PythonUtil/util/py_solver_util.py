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


class PyUtilSolver:

    def SolveLinearTimeProblem(domainDisc,
                               solverDesc,
                               outFilePrefix,
                               timeDesc,
                               outDesc,
                               startValueCB):
        print("")
        print("util.solver: setting up linear system ...")