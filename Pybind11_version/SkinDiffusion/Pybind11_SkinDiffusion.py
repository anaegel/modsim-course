import os
env_ug4_root = os.environ['UG4_ROOT']

import sys
sys.path.append(env_ug4_root+'/lib/')
sys.path.append(env_ug4_root+'/bin/plugins/ug4py')
print(sys.path)

import ug4py as ug4
import pylimex as limex
import pyconvectiondiffusion as cd


# Initialize UG4
ug4.InitUG(2, ug4.AlgebraType("CPU", 1))

# Load domain.
dom = ug4.Domain2d()
ug4.LoadDomain(dom, "skin2d-aniso.ugx")

# Create approximation space.
approxSpace = ug4.ApproximationSpace2d(dom)
approxSpace.add_fct("u","Lagrange", 1)
approxSpace.init_levels()
approxSpace.init_surfaces()
approxSpace.init_top_surface()
approxSpace.print_statistic()

# Create discretization.
KDcor=1.0
KCor=1.0

elemDisc={}
elemDisc["COR"] = cd.ConvectionDiffusionFV12d("u", "COR")
elemDisc["COR"].set_diffusion(KDcor)
elemDisc["COR"].set_mass_scale(KCor)

elemDisc["LIP"] = cd.ConvectionDiffusionFV12d("u", "LIP")
elemDisc["LIP"].set_diffusion(KDcor)
elemDisc["LIP"].set_mass_scale(KCor)

dirichletBND = ug4.DirichletBoundary2dCPU1()
dirichletBND.add(1.0, "u", "TOP_SC")
dirichletBND.add(0.0, "u", "BOTTOM_SC")

domainDisc = ug4.DomainDiscretization2dCPU1(approxSpace)
domainDisc.add(elemDisc["COR"])
domainDisc.add(elemDisc["LIP"])
domainDisc.add(dirichletBND)



ilu=ug4.ILUCPU1()
lsolver=ug4.LinearSolverCPU1()
lsolver.set_preconditioner(ilu)

lsolver=ug4.LUCPU1()


usol = ug4.GridFunction2dCPU1(approxSpace)
ug4.Interpolate(0.0, usol, "u")

startTime = 0.0
endTime = 1000.0
dt = 25.0

timeDisc=ug4.ThetaTimeStepCPU1(domainDisc, 1.0)

timeInt = limex.LinearTimeIntegrator2dCPU1(timeDisc)
#timeInt = limex.ConstStepLinearTimeIntegrator2dCPU1(timeDisc)
timeInt.set_linear_solver(lsolver)
timeInt.set_time_step(dt)

try:
    timeInt.apply(usol, endTime, usol, startTime)
except Exception as inst:
    print(inst)
    
# Exporting the result to a vtk-file
ug4.vtk_export_ho(usol, domainDisc, 1, ug4.VTKOutput2d(), "skin_result")
# TODO: Demonstrator for LIMEX.
#nstages = 2    
#limex = limex.LimexTimeIntegrator2dCPU1(nstages)
#print(limex)

