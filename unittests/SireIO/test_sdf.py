from Sire.Base import *
from Sire.IO import *
from Sire.Mol import *

from glob import glob
from nose.tools import assert_equal, assert_almost_equal

# Check that we have SDF support in this version of Sire.
has_sdf = True

try:
    sdf = SDF()
except:
    # No SDF support.
    has_sdf = False

# General test of ability to read and write Mol2 files.
# All sdf files in the "../io/" directory are parsed.
# Once the input file is parsed we then check that the parser constructs a
# Sire Molecule from the parsed data. Following this, we then check that the
# parser can convert the molecule back into the correct data format, ready to
# be written to file.
def test_read_write(verbose=False):
    if not has_sdf:
        return

    # Glob all of the SDF files in the example file directory.
    sdffiles = glob('../io/*sdf')

    # Loop over all test files.
    for file in sdffiles:

        # Test in parallel and serial mode.
        for use_par in [True, False]:

            if verbose:
                print("Reading SDF file: %s" % file)
                print("Parallel = %s" % use_par)

            # Parse the file into a Mol2 object.
            p = SDF(file, {"parallel" : wrap(use_par)})

            if verbose:
                print("Constructing molecular system...")

            # Construct a Sire molecular system.
            s = p.toSystem()

            if verbose:
                print("Reconstructing SDF data from molecular system...")

            # Now re-parse the molecular system.
            p = SDF(s, {"parallel" : wrap(use_par)})

            s2 = p.toSystem()

            assert s.nMolecules() == s2.nMolecules()

            m1 = s[0]
            m2 = s2[0]

            assert m1.nAtoms() == m2.nAtoms()

            for key in m1.properties().keys():
                v1 = m1.property(key)
                v2 = m2.property(key)

                assert type(v1) == type(v2)

            for i in range(0, m1.nAtoms()):
                a1 = m1.atoms()[i]
                a2 = m2.atoms()[i]
                assert a1.property("coordinates") == a2.property("coordinates")
                assert a1.property("radical") == a2.property("radical")
                assert a1.property("element") == a2.property("element")
                assert a1.name() == a2.name()

            c1 = m1.property("connectivity")
            c2 = m2.property("connectivity")

            for bond in c1.getBonds():
                p1 = c1.properties(bond)
                p2 = c2.properties(bond)

                for k in p1.keys():
                    assert p1[k] == p2[k]

            if m1.hasProperty("sdf_fields"):
                assert m1.property("sdf_fields") == m2.property("sdf_fields")

            assert m1.property("sdf_counts") == m2.property("sdf_counts")

            if m1.hasProperty("sdf_properties"):
                assert m1.property("sdf_properties") == m2.property("sdf_properties")

            if verbose:
                print("Passed!\n")


# Specific atom coordinate data validation test for file "../io/challenge.sdf".
def test_atom_coords(verbose=False):
    if not has_sdf:
        return

    # Test coordinates.
    coords = [[ 0.5369,    0.9749,    0.0000], 
              [ 1.4030,    0.4749,    0.0000],
              [ 2.2690,    0.9749,    0.0000],
              [ 1.8015,    0.0000,    0.0000],
              [ 1.0044,    0.0000,    0.0000],
              [ 1.9590,    1.5118,    0.0000],
              [ 2.8059,    1.2849,    0.0000],
              [ 2.5790,    0.4380,    0.0000], 
              [ 0.0000,    0.6649,    0.0000]]
 
    # Test in parallel and serial mode.
    for use_par in [True, False]:

        if verbose:
            print("Reading SDF file: ../io/challenge.sdf")
            print("Parallel = %s" % use_par)

        # Parse the SDF file.
        p = SDF('../io/challenge.sdf', {"parallel" : wrap(use_par)})

        if verbose:
            print("Constructing molecular system...")

        # Create a molecular system.
        s = p.toSystem()

        # Get the first molecule.
        m = s[MolIdx(0)]

        if verbose:
            print("Checking atomic coordinates...")

        # Loop over all of the atoms.
        for i in range(0, len(coords)):

            # Extract the atom.
            a = m.atoms()[i]

            # Extract the atom coordinates.
            c = a.property("coordinates")

            # Validate parsed coordinates against known values.
            assert_almost_equal( c[0], coords[i][0] )
            assert_almost_equal( c[1], coords[i][1] )
            assert_almost_equal( c[2], coords[i][2] )

        if verbose:
            print("Passed!\n")

if __name__ == "__main__":
    test_read_write(True)
    test_atom_coords(True)
