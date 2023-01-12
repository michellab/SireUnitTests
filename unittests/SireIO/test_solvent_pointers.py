try:
    import sire as sr

    sr.use_old_api()
except ImportError:
    pass

from Sire.IO import *

from nose.tools import assert_equal

# testing that the "SOLVENT_POINTERS" section is written correctly


def _getSolventPointers(f):
    FILE = open(f, "r")

    line = FILE.readline()

    while line:
        if line.find("SOLVENT_POINTERS") != -1:
            line = FILE.readline()
            line = FILE.readline()
            words = line.split()
            return (int(words[0]), int(words[1]), int(words[2]))

        line = FILE.readline()


def _test(topfile, verbose):
    if verbose:
        print("\nTesting '%s'" % topfile)
        print("Reading the file...")

    s = MoleculeParser.read(topfile)

    if verbose:
        print("Writing the file...")

    a = AmberPrm(s)
    a.writeToFile("test.top")

    if verbose:
        print("Comparing SOLVENT_POINTERS...")

    orig_pointers = _getSolventPointers(topfile)
    new_pointers = _getSolventPointers("test.top")

    if verbose:
        print("Original SOLVENT_POINTERS == %s" % str(orig_pointers))
        print("     New SOLVENT_POINTERS == %s" % str(new_pointers))

    assert_equal(orig_pointers, new_pointers)


def test_solvent_pointers(verbose=False):
    _test("../io/ala.top", verbose)
    _test("../io/ethane.top", verbose)
    _test("../io/ose.top", verbose)
    _test("../io/waterbox.top", verbose)
    _test("../io/NA16.top", verbose)


if __name__ == "__main__":
    test_solvent_pointers(True)
