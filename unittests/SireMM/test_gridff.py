
try:
    import sire as sr
    sr.use_old_api()
except ImportError:
    pass

from Sire.MM import *
from Sire.FF import *
from Sire.Vol import *
from Sire.IO import *
from Sire.System import *
from Sire.Move import *
from Sire.Maths import *
from Sire.Units import *
from Sire.Mol import *
from Sire.Qt import *

import Sire.Stream

from nose.tools import assert_almost_equal

import math

coul_cutoff = 15 * angstrom
lj_cutoff = 5 * angstrom

grid_spacing = 0.25 * angstrom

nmoves = 5000

(mols, space) = Amber().readCrdTop("../io/waterbox.crd", "../io/waterbox.top")
water_space = space

clusterbox = AABox( Vector(10), Vector(2) )
print(clusterbox)

cluster = MoleculeGroup("cluster")
waters = MoleculeGroup("waters")

for i in range(0, mols.nMolecules()):
    mol = mols[ MolIdx(i) ].molecule()

    if clusterbox.contains( mol.evaluate().center() ):
        cluster.add(mol)
    else:
        waters.add(mol)

gridff = GridFF("gridff")
gridff2 = GridFF2("gridff2")
cljff = InterGroupCLJFF("cljff")

gridff.setShiftElectrostatics(True)
gridff.setCoulombCutoff(coul_cutoff)
gridff.setLJCutoff(lj_cutoff)
gridff.setGridSpacing(grid_spacing)

gridff2.setShiftElectrostatics(True)
gridff2.setCoulombCutoff(coul_cutoff)
gridff2.setLJCutoff(lj_cutoff)
gridff2.setGridSpacing(grid_spacing)

cljff.setSwitchingFunction( HarmonicSwitchingFunction(coul_cutoff,lj_cutoff) )
cljff.setShiftElectrostatics(True)
cljff.setCombiningRules("arithmetic")

gridff.addFixedAtoms(waters)
gridff2.addFixedAtoms(waters)
cljff.add(waters, MGIdx(1))

gridff.add(cluster, MGIdx(0))
gridff2.add(cluster, MGIdx(0))
cljff.add(cluster, MGIdx(0))

system1 = System()
system2 = System()
system3 = System()

system1.add(waters)
system1.add(gridff)

system2.add(waters)
system2.add(gridff2)

system3.add(waters)
system3.add(cljff)

def test_single_point(verbose=False):

    system1.mustNowRecalculateFromScratch()
    system2.mustNowRecalculateFromScratch()
    system3.mustNowRecalculateFromScratch()

    if verbose:
        print("Comparing single point energies...")
        print("...GridFF")

    nrg1 = system1.energy( gridff.components().coulomb() ).value()

    if verbose:
        print("...GridFF2")

    nrg2 = system2.energy( gridff2.components().coulomb() ).value()

    if verbose:
        print("...InterGroupCLJFF")

    nrg3 = system3.energy( cljff.components().coulomb() ).value()

    if verbose:
        print("GridFF (%s) : GridFF2 (%s) : CLJFF(%s)" \
              % (nrg1, nrg2, nrg3))

    assert_almost_equal(nrg1, nrg2, 1)
    assert_almost_equal(nrg1, nrg3, 1)
    assert_almost_equal(nrg2, nrg3, 1)

def _pvt_test_stream(verbose, s):
    if verbose:
        print("Testing streaming %s" % s.forceFields()[FFIdx(0)])

    s = System(s)

    moves = RigidBodyMC(cluster)

    moves.move(s, 5*nmoves, False)

    oldnrg = s.energy().value()

    Sire.Stream.save(s, "tmp.s3")

    s2 = Sire.Stream.load("tmp.s3")

    newnrg = s2.energy().value()

    if verbose:
        print("%s vs %s" % (oldnrg, newnrg))

    assert_almost_equal(oldnrg, newnrg, 1)

    s2.mustNowRecalculateFromScratch()

    newnrg = s2.energy().value()

    if verbose:
        print("%s vs %s" % (oldnrg, newnrg))

    assert_almost_equal(oldnrg, newnrg, 1)

    s.mustNowRecalculateFromScratch()

    newnrg = s.energy().value()

    if verbose:
        print("%s vs %s" % (oldnrg, newnrg))

    assert_almost_equal(oldnrg, newnrg, 1)

def test_stream(verbose=False):
    _pvt_test_stream(verbose, system1)
    _pvt_test_stream(verbose, system2)
    _pvt_test_stream(verbose, system3)

def _pvt_test_simulation(verbose, s):
    if verbose:
        print("Testing simulation %s" % s.forceFields()[FFIdx(0)])

    s = System(s)

    moves = RigidBodyMC(cluster)

    moves.move(s, nmoves, False)

    oldnrg = s.energy().value()
    s.mustNowRecalculateFromScratch()
    newnrg = s.energy().value()

    if verbose:
        print("%s vs %s" % (oldnrg, newnrg))

    assert_almost_equal(oldnrg, newnrg, 1)

def test_simulation(verbose=False):
    _pvt_test_simulation(verbose, system1)
    _pvt_test_simulation(verbose, system2)
    _pvt_test_simulation(verbose, system3)

if __name__ == "__main__":
    test_single_point(True)
    test_simulation(True)
    test_stream(True)

