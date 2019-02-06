# 
#Sire script to test Sire and SOMD single point energies differences
#

import os
import re
import sys

from Sire.IO import *
from Sire.Mol import *
from Sire.CAS import *
from Sire.System import *
from Sire.Move import *
from Sire.MM import *
from Sire.FF import *
from Sire.Units import *
from Sire.Vol import *
from Sire.Maths import *
from Sire.Base import *
from Sire.Qt import *
from Sire.ID import *
from Sire.Config import *
from Sire.Analysis import *
from Sire.Tools.DCDFile import *
from Sire.Tools import Parameter, resolveParameters
import Sire.Stream
import time
import numpy as np

top = ['../io/CR_Toluene_methane_SYSTEM.top']
crd = ['../io/CR_Toluene_methane_SYSTEM.crd']
pert = ['../io/CR_Toluene_methane_MORPH.pert']


def _combR_test(top, crd, pert, combining_rules, verbose=False):
    
    if verbose:
        print("Testing whether geometric and arithmetic combining rules provide the same output")


    amber = Amber()
    (molecules, space) = amber.readCrdTop(crd, top)
    system = createSystemFreeEnergy(molecules, pert)
     
    system = setupForceFieldsFreeEnergy(system, space, combining_rules)
    if random_seed.val:
        ranseed = random_seed.val
    else:
        ranseed = RanGenerator().randInt(100000, 1000000)

    print("Setting up the simulation with random seed %s" % ranseed)

    moves = setupMovesFreeEnergy(system, ranseed, gpu.val, lambda_val.val)

    # Get energy from Sire
    nrg_sire = system.energy()

    print ("nrg_sire is %s" % nrg_sire)
    
    # Get energy from SOMD
    mdmoves = moves.moves()[0]
    integrator = mdmoves.integrator()
    nrg_somd = integrator.getPotentialEnergy(system)

    print ("nrg_somd is %s" % nrg_somd)

    diff = nrg_sire - nrg_somd

    print ("@@@@ The single point energy difference between sire and somd at lambda %s is %s " % (lambda_val.val,diff) )


def test_comb_rules(verbose=False):
    if verbose:
        print("\nTesting Combining Rules")

    _combR_test(top, crd, pert, verbose)

if __name__ == "__main__":
    test_comb_rules(True)
