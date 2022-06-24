
try:
    import sire as sr
    sr.use_old_api()
except ImportError:
    pass

from Sire.IO import *
from Sire.MM import *
from Sire.Mol import *
from Sire.CAS import *
from Sire.Units import *

from nose.tools import assert_equal, assert_almost_equal

OW = AtomName("OW")
HW1 = AtomName("HW1")
HW2 = AtomName("HW2")

k_oh = 0.5 * 502416.0 * (kJ_per_mol / (nanometer*nanometer))
r_oh = 0.09572 * nanometer

k_hoh = 0.5 * 628.02 * (kJ_per_mol / (radian*radian))
t_hoh = 104.52 * degrees

def test_flex(verbose=False):
    s = MoleculeParser.read("../io/waterbox.gro", "../io/waterbox.grotop",
                            {"GROMACS_PATH":"../io/gromacs"})

    water = s.search("molidx 0")[0].molecule()

    bonds = water.property("bond")
    angles = water.property("angle")

    assert_equal( bonds.nFunctions(), 2 )
    assert_equal( angles.nFunctions(), 1 )

    # Debug errors by getting the warnings from the parser...
    #g = GroTop("../io/waterbox.grotop", {"GROMACS_PATH":"../io/gromacs"})
    #
    #for warning in g.warnings():
    #    print(warning)

    boh1 = bonds.potential( OW, HW1 )
    boh2 = bonds.potential( OW, HW2 )

    if verbose:
        print("Bond OW-HW1 = %s" % boh1)
        print("Bond OW-HW2 = %s" % boh2)

    assert_equal(boh1, boh2)

    r = Symbol("r")

    amb_boh = AmberBond(boh1, r)

    if verbose:
        print("Amber Bond OW-HW1 = %s" % amb_boh)

    assert_almost_equal( amb_boh.r0(), r_oh.value() )
    assert_almost_equal( amb_boh.k(), k_oh.value() )

    hoh = angles.potential(HW1, OW, HW2)

    if verbose:
        print("Angle HW1-OW-HW2 = %s" % hoh)

    t = Symbol("theta")

    amb_hoh = AmberAngle(hoh, t)

    if verbose:
        print("Amber Angle HW1-OW-HW2 = %s" % amb_hoh)

    assert_almost_equal( amb_hoh.theta0(), t_hoh.value() )
    assert_almost_equal( amb_hoh.k(), k_hoh.value() )        

if __name__ == "__main__":
    test_flex(True)
