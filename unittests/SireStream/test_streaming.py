try:
    import sire as sr

    sr.use_old_api()
except ImportError:
    pass

from Sire.Mol import *
from Sire.IO import *
from Sire.Vol import *
from Sire.FF import *
from Sire.MM import *
from Sire.CAS import *
from Sire.Maths import *
from Sire.Qt import *
from Sire.Units import *
from Sire.System import *
from Sire.Move import *

import Sire.Stream

from nose.tools import assert_almost_equal

cljff = InterCLJFF()

mincoords = Vector(-18.3854, -18.66855, -18.4445)
maxcoords = Vector(18.3854, 18.66855, 18.4445)

vol = PeriodicBox(mincoords, maxcoords)
switchfunc = HarmonicSwitchingFunction(15 * angstrom, 14.5 * angstrom)

cljff.setSpace(vol)
cljff.setSwitchingFunction(switchfunc)

mols = PDB().read("../io/water.pdb")

i = 0

mol = mols.moleculeAt(0).molecule()

mol = (
    mol.edit()
    .atom(AtomName("O00"))
    .setProperty("LJ", LJParameter(3.15363 * angstrom, 0.1550 * kcal_per_mol))
    .molecule()
    .atom(AtomName("H01"))
    .setProperty("charge", 0.520 * mod_electron)
    .molecule()
    .atom(AtomName("H02"))
    .setProperty("charge", 0.520 * mod_electron)
    .molecule()
    .atom(AtomName("M03"))
    .setProperty("charge", -1.04 * mod_electron)
    .molecule()
    .commit()
)

charges = mol.property("charge")
ljs = mol.property("LJ")

cljff.add(mol)

for i in range(1, mols.nMolecules()):
    mol = mols.moleculeAt(i).molecule()

    mol = (
        mol.edit()
        .rename("T4P")
        .setProperty("charge", charges)
        .setProperty("LJ", ljs)
        .commit()
    )

    cljff.add(mol)

system = System()

system.add(cljff)

lam = Symbol("lambda")

system.setComponent(lam, 0.2)
system.setComponent(system.totalComponent(), lam * cljff.components().total())

mc = RigidBodyMC(cljff.group(MGIdx(0)))

moves = SameMoves(mc)


def test_stream(verbose=False):
    if verbose:
        print("saving system...")

    data = Sire.Stream.save(system)

    if verbose:
        print(("%s takes up %d bytes" % (system.what(), data.size())))

    header = Sire.Stream.getDataHeader(data)

    if verbose:
        print((header.dataType()))
        print((header.requiredLibraries()))

        print((header.createdBy()))
        print((header.createdWhere()))

        print((header.requiredMemory()))
        print((header.compressionRatio()))
        print((header.digest()))
        print((header.repository()))
        print((header.buildVersion()))
        print((header.systemInfo()))

    system2 = Sire.Stream.load(data)

    if verbose:
        print(("Re-Reading the data"))
        print(system)
        print(system2)

        print("Probing the system...")

    nrg1 = system.energy()
    nrg2 = system2.energy()

    if verbose:
        print("%s versus %s (should be equal)" % (nrg1, nrg2))

    assert_almost_equal(nrg1, nrg2, 3)


if __name__ == "__main__":
    test_stream(True)
