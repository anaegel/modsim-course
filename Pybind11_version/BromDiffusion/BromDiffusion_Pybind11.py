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
requiredSubsets = ["INNER", "WALL", "IN"] # defining subsets
gridName = "brom-dose.ugx"  # Grid created with ProMesh
numRefs = 2  # Number of Refinement steps on given grid

# Choosing a domain object
# (either 1d, 2d or 3d suffix)
dom = ug4.Domain3d()

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

# checking if geometry has the needed subsets of the probelm
sh = dom.subset_handler()
for e in requiredSubsets:
    if sh.get_subset_index(e) == -1:
        print(f"Domain does not contain subset {e}.")
        sys.exit(1)

# Create approximation space which describes the unknowns in the equation
fct = "u"  # name of the function
type = "Lagrange"
order = 1  # polynom order for lagrange

approxSpace = ug4.ApproximationSpace3d(dom)
approxSpace.add_fct(fct, type, order)
approxSpace.init_levels()
approxSpace.init_surfaces()  
approxSpace.init_top_surface()
approxSpace.print_statistic()

# Create discretization for the Convection-Diffusion Equation
# perform the discretization on the actual elements
# using first order finite volumes (FV1) on the 3d grid
dif = 10

# creating instance of a convection diffusion equation
elemDisc = cd.ConvectionDiffusionFV13d("u", "INNER")
elemDisc.set_diffusion(dif)

# ug4 separates the boundary value and the discretization
# boundary conditions can be enforced through a post-process (dirichlet).
# To init at boundary, the value, function name from the Approximationspace
# and the subset name are needed
dirichletBND = ug4.DirichletBoundary3dCPU1()
dirichletBND.add(1.0, "u", "TOP")  # Zufluss
dirichletBND.add(0.0, "u", "IN")  # "offene Dose"


# create the discretization object which combines all the
# separate discretizations into one domain discretization.
domainDisc = ug4.DomainDiscretization3dCPU1(approxSpace)
domainDisc.add(elemDisc)
domainDisc.add(dirichletBND)


# Create Solver
# In this case we use an LU (Lower Upper) Solver for an exact solution
ilu=ug4.ILUCPU1()
lsolver=ug4.LinearSolverCPU1()
lsolver.set_preconditioner(ilu)

lsolver=ug4.LUCPU1()

# Solve the transient problem
# Use the Approximationspace to
# create a vector of unknowns and a vector
# which contains the right hand side
usol = ug4.GridFunction3dCPU1(approxSpace)

# Init the vector representing the unknowns with function
def InitialValue(x,y,z,t,si):
    if (z<1.75):
        return 0.0
    else:
        return 1.0

InitialValue = "def InitialValue(x, y, z, t, si):\n    return 0.0 if (z<1.75) else 1.0\n"
name = "InitialValue"
ug4.Interpolate4py(InitialValue, name, usol, "u")

# Define start time, end time and step size
startTime = 0.0
endTime = 1.0
dt = 0.0125

# create a time discretization with the theta-scheme
# takes in domain and a theta
# theta = 1.0 -> Implicit Euler, 
# 0.0 -> Explicit Euler 
timeDisc=ug4.ThetaTimeStepCPU1(domainDisc, 1.0) 

# creating Time Iterator, setting the solver and step size
timeInt = limex.LinearTimeIntegrator3dCPU1(timeDisc)
#timeInt = limex.ConstStepLinearTimeIntegrator3dCPU1(timeDisc)
timeInt.set_linear_solver(lsolver)
timeInt.set_time_step(dt)

# Solving the transient problem
try:
    timeInt.apply(usol, endTime, usol, startTime)
except Exception as inst:
    print(inst)


# Exporting the result to a vtu-file
# can be visualized in paraview or with a python extension
ug4.WriteGridFunctionToVTK(usol, "Solution_BromDif_Pybind")

# Plotting the result using pyvista
import pyvista

result = pyvista.read('Solution_BromDif_Pybind.vtu')
print()
print("Pyvista input: ")
print(result)
result.plot(scalars="u", show_edges=True, cmap='coolwarm')
