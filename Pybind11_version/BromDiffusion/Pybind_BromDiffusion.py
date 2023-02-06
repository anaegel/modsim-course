import os
env_ug4_root = os.environ['UG4_ROOT']

import sys
sys.path.append(env_ug4_root+'/lib/')
sys.path.append(env_ug4_root+'/bin/plugins/ug4py')
print(sys.path)

import ug4py as ug4
import pylimex as limex
import pyconvectiondiffusion as cd

# Function Definitions
def InitialValue(x,y,z,t,si):
    if z<1.75:
        return 0.0
    else:
        return 1.0

# Initialize UG4
ug4.InitUG(3, ug4.AlgebraType("CPU", 1))

# Create the Domain
dom = ug4.Domain2d()
ug4.LoadDomain(dom, "brom-dose2.ugx")
