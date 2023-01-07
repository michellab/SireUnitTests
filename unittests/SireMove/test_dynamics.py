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
from Sire.Stream import *

import sys

from nose.tools import assert_almost_equal

cljff = InterCLJFF("cljff")

na = Molecule("Sodium")
na = na.edit().add(CGName("1")).add(AtomName("Na")).molecule().commit()

na = (
    na.edit()
    .atom(AtomName("Na"))
    .setProperty("coordinates", Vector(0, 0, 0))
    .setProperty("element", Element("Sodium"))
    .setProperty("charge", 1 * mod_electron)
    .setProperty("LJ", LJParameter(3.0522 * angstrom, 0.4598 * kcal_per_mol))
    .molecule()
    .commit()
)

cl = Molecule("Chloride")
cl = cl.edit().add(CGName("1")).add(AtomName("Cl")).molecule().commit()

cl = (
    cl.edit()
    .atom(AtomName("Cl"))
    .setProperty("coordinates", Vector(0, 0, 0))
    .setProperty("element", Element("Chlorine"))
    .setProperty("charge", -1 * mod_electron)
    .setProperty("LJ", LJParameter(4.4124 * angstrom, 0.11779 * kcal_per_mol))
    .molecule()
    .commit()
)

salt = MoleculeGroup("salt")

# create a box of salt
add_cl = False

for i in range(0, 5):
    for j in range(0, 5):
        for k in range(0, 5):
            coords = Vector(i * 3.06, j * 3.06, k * 3.06)

            if add_cl:
                cl = (
                    cl.edit()
                    .renumber()
                    .setProperty("coordinates", AtomCoords([coords]))
                    .commit()
                )

                salt.add(cl)
                add_cl = False
            else:
                na = (
                    na.edit()
                    .renumber()
                    .setProperty("coordinates", AtomCoords([coords]))
                    .commit()
                )

                salt.add(na)
                add_cl = True

cljff = InterCLJFF("salt-salt")
cljff.add(salt)

system = System()
system.add(salt)
system.add(cljff)


def test_conservation(verbose=False):

    mdmove = MolecularDynamics(salt)

    mdmove.setTimeStep(1 * femtosecond)

    old_nrg = None

    for i in range(0, 25):
        if verbose:
            print("move %d" % (i + 1))

        mdmove.move(system, 20)

        if verbose:
            print(system.energy())
            print(mdmove.kineticEnergy())
            print(system.energy() + mdmove.kineticEnergy())

        total_nrg = system.energy().value() + mdmove.kineticEnergy().value()

        if old_nrg:
            assert_almost_equal(total_nrg, old_nrg, 1)

        old_nrg = total_nrg


if __name__ == "__main__":
    test_conservation(True)
