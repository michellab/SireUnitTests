try:
    import sire as sr

    sr.use_old_api()
except ImportError:
    pass

from Sire.Base import *
from Sire.IO import *
from Sire.Mol import *

from glob import glob
from nose.tools import assert_equal, assert_almost_equal

# Check that we have PDB2 support in this version of Sire.
has_pdb2 = True

try:
    p = PDB2()
except:
    # No PDB2 support.
    has_pdb2 = False

# General test of ability to read and write PDB files.
# All PDB files in the "../io/" directory are parsed.
# Once the input file is parsed we then check that the parser constructs a
# Sire Molecule from the parsed data. Following this, we then check that the
# parser can convert the molecule back into the correct data format, ready to
# be written to file.
def test_read_write(verbose=False):
    if not has_pdb2:
        return

    # Glob all of the PDB files in the example file directory.
    pdbfiles = glob("../io/*pdb")

    # Loop over all files.
    for file in pdbfiles:

        # Test in parallel and serial mode.
        for use_par in [True, False]:

            if verbose:
                print("Reading PDB file: %s" % file)
                print("Parallel = %s" % use_par)

            # Parse the file into a PDB2 object.
            # Errors should be thrown if the record data in a file
            # doesn't match the PDB format.
            p = PDB2(file, {"parallel": wrap(use_par)})

            if verbose:
                print("Constructing molecular system...")

            # Construct a Sire molecular system.
            s = p.toSystem()

            if verbose:
                print("Reconstructing PDB data from molecular system...")

            # Now re-parse the molecular system.
            p = PDB2(s, {"parallel": wrap(use_par)})

            if verbose:
                print("Passed!\n")


# Specific atom coordinate data validation test for file "../io/ntrc.pdb".
def test_atom_coords(verbose=False):
    if not has_pdb2:
        return

    # Test atoms.
    atoms = ["CA", "CB", "N", "O", "HB"]

    # Test coordinates.
    coords = [
        [-13.721, -3.484, 14.690],
        [-10.695, -0.294, 14.759],
        [-8.536, -2.557, 13.277],
        [-7.037, -1.615, 9.350],
        [-5.045, 2.118, 8.812],
    ]

    # Test in parallel and serial mode.
    for use_par in [True, False]:

        if verbose:
            print("Reading PDB file: ../io/ntrc.pdb")
            print("Parallel = %s" % use_par)

        # Parse the PDB file.
        p = PDB2("../io/ntrc.pdb", {"parallel": wrap(use_par)})

        if verbose:
            print("Constructing molecular system...")

        # Create a molecular system.
        s = p.toSystem()

        # Get the first molecule.
        m = s[MolIdx(0)]

        if verbose:
            print("Checking atomic coordinates...")

        # Loop over all of the atoms.
        for i in range(0, len(atoms)):

            # Extract the atom from the residue "i + 1".
            a = m.atom(AtomName(atoms[i]) + ResNum(i + 1))

            # Extract the atom coordinates.
            c = a.property("coordinates")

            # Validate parsed coordinates against known values.
            assert_almost_equal(c[0], coords[i][0])
            assert_almost_equal(c[1], coords[i][1])
            assert_almost_equal(c[2], coords[i][2])

        if verbose:
            print("Passed!\n")


# Test that elements are inferred correctly when the PDB file is a PSF companion
# used for a NAMD simulation. Here the element data is often missing, and this
# section of the record line is used to label the residue to which each atom
# belongs.
def test_psf_companion(verbose=False):

    # Load the PDB file.
    p = PDB2("../io/psf_companion.pdb")

    # Create the molecular system.
    s = p.toSystem()

    # Create the list of elements.
    elements = [
        Element("C"),
        Element("C"),
        Element("O"),
        Element("N"),
        Element("H"),
        Element("C"),
        Element("C"),
        Element("C"),
        Element("O"),
        Element("N"),
    ]

    # Now assert that the elements are correct.
    for i, atom in enumerate(s.molecule(MolIdx(0)).atoms()):
        assert atom.property("element") == elements[i]


if __name__ == "__main__":
    test_read_write(True)
    test_atom_coords(True)
    test_psf_companion(True)
