import os
env_ug4_root = os.environ['UG4_ROOT']

import sys
sys.path.append(env_ug4_root+'/lib/')
sys.path.append(env_ug4_root+'/bin/plugins/ug4py')
#print(sys.path)

import ug4py as ug4
import pylimex as limex
import pyconvectiondiffusion as cd

# Setup:
# defining needed subsets, gird and number of refinements
requiredSubsets = ["Inner", "Boundary"] # defining subsets
gridName= "laplace_sample_grid_2d.ugx" 
numRefs	= 4

# Choosing a domain object
# (either 1d, 2d or 3d suffix)
dom = ug4.Domain2d()

# Loading the given grid into the domain
print(f"Loading Domain {gridName} ...")
ug4.LoadDomain(dom, gridName)
print("Domain loaded.")

# Optional: Refining the grid
if numRefs > 0:
    print("Refining ...")
    refiner = ug4.GlobalDomainRefiner(dom)
    for i in range(numRefs):
        ug4.TerminateAbortedRun()
        refiner.refine()
        print(f"Refining step {i} ...")

    print("Refining done")

# Create approximation space which describes the unknowns in the equation
approxSpace = ug4.ApproximationSpace2d(dom)
approxSpace.add_fct("u", "Lagrange", 1)
approxSpace.add_fct("c", "Lagrange", 1)
approxSpace.init_levels()
#approxSpace.init_surfaces()  ## TESTEN!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
approxSpace.init_top_surface()
approxSpace.print_statistic()

# Elementdiskretisierung (FV)
M=1
Lambda=0.01
charLength = 1.0
charTime = charLength*charLength/M
A = 100.0

# Funktionen definieren
HelmholtzEnergyF0 = f"def HelmholtzEnergyF0(c):\n   return {A}*c*c*(1.0-c)*(1.0-c)"
HelmholtzEnergyF1 = f"def HelmholtzEnergyF0(c):\n   return 2.0*{A}*c*(2*c*c - 3*c + 1)"
HelmholtzEnergyF2 = f"def HelmholtzEnergyF0(c):\n   return 2.0*{A}*(6.0*(c-1)*c + 1)"

