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

# Helper function for calling arbitrary parsers.
# Add in support as needed.
def parse(format, x=None, use_par=False):
    if format == "PDB":
        if x == None:
            return PDB2()
        else:
            return PDB2(x, {"parallel": wrap(use_par)})
    elif format == "Mol2":
        if x == None:
            return Mol2()
        else:
            return Mol2(x, {"parallel": wrap(use_par)})
    else:
        return False


# Test for round-trip file conversion,
# i.e. format1 --> format2 --> format1
# Eventually this should check that a minimal set of properties that
# are essential for simulation are preserved during interconversion.
def _test_round_trip(file, format1, format2, verbose=False):
    if not parse(format1):
        print("ERROR: No support for format: %s" % format1)
        return

    if not parse(format2):
        print("ERROR: No support for format: %s" % format2)
        return

    print(
        "Testing round-trip file "
        "conversion: %s --> %s --> %s\n" % (format1, format2, format1)
    )

    # Test in parallel and serial mode.
    for use_par in [True, False]:

        if verbose:
            print("Reading file: %s" % file)
            print("Parallel = %s" % use_par)

        # Parse the file into the parser object.
        p = parse(format1, file, use_par)

        if verbose:
            print("Constructing first molecular system...")

        # Construct the original Sire molecular system.
        s1 = p.toSystem()

        if verbose:
            print(
                "Converting first molecular system to %s format..." % format2
            )

        # Now construct data records in the second format using the molecular system.
        p = parse(format2, s1, use_par)

        p.writeToFile("test.pdb")

        if verbose:
            print("Constructing second molecular system...")

        # Construct second Sire molecular system.
        s2 = p.toSystem()

        if verbose:
            print(
                "Converting second molecular system back to %s format..."
                % format1
            )

        # Now construct original data records from the molecular system.
        p = parse(format1, s2, use_par)

        if verbose:
            print("Constructing final molecular system...")

        # Construct final Sire molecular system.
        s3 = p.toSystem()

        if verbose:
            print("Checking that molecular systems are the same...")

        # Check that the molecular systems contain the number of molecules.
        assert_equal(s1.nMolecules(), s2.nMolecules(), s3.nMolecules())

        # Now check that each molecule is the same.
        for i in range(0, s1.nMolecules()):

            # Extract the molecule from each system.
            m1 = s1[MolIdx(i)]
            m2 = s2[MolIdx(i)]
            m3 = s3[MolIdx(i)]

            # Number of atoms.
            assert_equal(m1.nAtoms(), m2.nAtoms(), m3.nAtoms())

            # Number of chains.
            assert_equal(m1.nChains(), m2.nChains(), m3.nChains())

            # Number of residues.
            assert_equal(m1.nResidues(), m2.nResidues(), m3.nResidues())

            # Loop over all of the atoms.
            for j in range(0, m1.nAtoms()):

                # Extract the atom from the residue "i + 1".
                a1 = m1.atom(AtomIdx(j))
                a2 = m2.atom(AtomIdx(j))
                a3 = m3.atom(AtomIdx(j))

                # Check that the atom names match.
                assert_equal(a1.name().value(), a2.name().value())
                assert_equal(a1.name().value(), a2.name().value())

                # Extract the atom coordinates.
                c1 = a1.property("coordinates")
                c2 = a2.property("coordinates")
                c3 = a3.property("coordinates")

                # Check that coordinates match.
                # We need to round these since assert_almost_equal is useless,
                # i.e. the places flag doesn't always work as expected.
                assert_almost_equal(round(c1[0], 3), round(c2[0], 3))
                assert_almost_equal(round(c1[0], 3), round(c3[0], 3))
                assert_almost_equal(round(c1[1], 3), round(c2[1], 3))
                assert_almost_equal(round(c1[1], 3), round(c3[1], 3))
                assert_almost_equal(round(c1[2], 3), round(c2[2], 3))
                assert_almost_equal(round(c1[2], 3), round(c3[2], 3))

                # Check that residue data matches.
                if a1.isWithinResidue():
                    assert_equal(
                        a1.residue().name().value()[0:2],
                        a2.residue().name().value()[0:2],
                    )
                    assert_equal(
                        a1.residue().name().value()[0:2],
                        a3.residue().name().value()[0:2],
                    )

                    assert_equal(
                        a1.residue().number().value(),
                        a2.residue().number().value(),
                    )
                    assert_equal(
                        a1.residue().number().value(),
                        a3.residue().number().value(),
                    )

                # Check that chain data matches.
                if a1.isWithinChain():
                    assert_equal(
                        a1.chain().name().value(), a2.chain().name().value()
                    )
                    assert_equal(
                        a1.chain().name().value(), a3.chain().name().value()
                    )

                # Check the element symbol matches.
                if a1.hasProperty("element"):
                    assert_equal(
                        a1.property("element").symbol(),
                        a2.property("element").symbol(),
                    )
                    assert_equal(
                        a1.property("element").symbol(),
                        a3.property("element").symbol(),
                    )

        if verbose:
            print("Passed!\n")


def test_round_trip1(verbose=False):
    _test_round_trip("../io/ntrc.pdb", "PDB", "Mol2", verbose)


def test_round_trip2(verbose=False):
    _test_round_trip("../io/complex.mol2", "Mol2", "PDB", verbose)


if __name__ == "__main__":
    test_round_trip1(True)
    test_round_trip2(True)
