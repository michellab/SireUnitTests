
try:
    import sire as sr
    sr.use_old_api()
except ImportError:
    pass

from Sire.Base import *
from Sire.IO import *
from Sire.MM import *
from Sire.Mol import *

from glob import glob
from nose.tools import assert_equal, assert_almost_equal

# Check that we have support for the required parsers in this version of Sire.
has_parsers = True

try:
    p = PDB2()
    p = CharmmPSF()
    p = Supplementary()
except:
    # Parsers are missing.
    has_parsers = False

# Test multi-stage molecular parser. This should generate a fully parameterised
# CHARMM molecular system. Three files are required:
#     1) PSF        - Topology
#     2) PDB        - Coordinates
#     3) PARAMS/INP - Force-field
def test_multi_parse(verbose=False):
    if not has_parsers:
        return

    # Glob all of the NAMD example directories, removing the README.md file.
    dirs = glob('../io/namd/*')
    [dirs.remove(d) for d in dirs if d.endswith('README.md')]

    print(dirs)

    # Test in parallel and serial mode.
    for use_par in [True, False]:

        # Test each NAMD example.
        for dir in dirs:

            # Glob the input files.
            files = glob('%s/*' % dir)

            if verbose:
                print("Constructing molecular system from: %s" % files)
                print("Parallel = %s" % use_par)

            # Construct a molecular system from the three input files.
            s = MoleculeParser.read(files, {"parallel" : BooleanProperty(use_par)})

            if verbose:
                print("Passed!\n")

# Test the structure of the molecular system.
def test_molecular_structure(verbose=False):
    if not has_parsers:
        return

    # Test in parallel and serial mode.
    for use_par in [True, False]:

        # Use the tiny example, since it contains multiple molecules.
        files = ['../io/namd/tiny/tiny.psf', \
                 '../io/namd/tiny/tiny.pdb', \
                 '../io/namd/tiny/par_all22_prot.inp']

        if verbose:
            print("Constructing molecular system from: %s" % files)
            print("Parallel = %s" % use_par)

        # Construct a molecular system from the three input files.
        s = MoleculeParser.read(files, {"parallel" : BooleanProperty(use_par)})

        # Now loop over all molecules and accumulate the total number
        # of bonds, angles, dihedrals, and impropers.

        num_bonds = 0
        num_angles = 0
        num_dihedrals = 0
        num_impropers = 0

        if verbose:
            print("Checking number of bonds, angles, dihedrals, and impropers...")

        for i in range(0, s.nMolecules()):

            # Extract molecule 'i' from the system.
            m = s[MolIdx(i)]

            if m.hasProperty("bond"):
                num_bonds += m.property("bond").nFunctions()

            if m.hasProperty("angle"):
                num_angles += m.property("angle").nFunctions()

            if m.hasProperty("dihedral"):
                num_dihedrals += m.property("dihedral").nFunctions()

            if m.hasProperty("improper"):
                num_impropers += m.property("improper").nFunctions()

        # Make sure the record counts match those in the PSF file.

        assert_equal(num_bonds, 396)
        assert_equal(num_angles, 421)
        assert_equal(num_dihedrals, 354)
        assert_equal(num_impropers, 35)

        if verbose:
            print("Checking atom properties...")

        # Get the first atom from the first molecule.
        a = s[MolIdx(0)].atoms()[0]

        # Extract the atom coordinates.
        c = a.property("coordinates")

        # Validate parsed coordinates against known values.
        assert_almost_equal( c[0],   0.598 )
        assert_almost_equal( c[1], -12.496 )
        assert_almost_equal( c[2],  -1.169 )

        # Validate the atomic mass and charge.
        assert_almost_equal( a.property("mass").value(),   14.007 )
        assert_almost_equal( a.property("charge").value(), -0.300 )

        # Validate residue data.
        assert_equal( a.residue().name().value(), 'ALA' )
        assert_equal( a.residue().number().value(), 1 )

        # Get the last atom from the last molecule.
        a = s[MolIdx(s.nMolecules()-1)].atoms()[-1]

        # Extract the atom coordinates.
        c = a.property("coordinates")

        # Validate parsed coordinates against known values.
        assert_almost_equal( c[0], 1.164 )
        assert_almost_equal( c[1], 3.494 )
        assert_almost_equal( c[2], 7.532 )

        # Validate the atomic mass and charge.
        assert_almost_equal( a.property("mass").value(),   1.008 )
        assert_almost_equal( a.property("charge").value(), 0.417 )

        # Validate residue data.
        assert_equal( a.residue().name().value(), 'TIP3' )
        assert_equal( a.residue().number().value(), 9226 )

        if verbose:
            print("Passed!\n")

if __name__ == "__main__":
    test_multi_parse(True)
    test_molecular_structure(True)
