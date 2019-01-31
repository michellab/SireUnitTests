
from Sire.IO import *
from Sire.Mol import *
from Sire.MM import *

from nose.tools import assert_equal

def test_extract(verbose=False):
    if verbose:
        print("Loading molecules")

    s = MoleculeParser.load("../io/thrombin.rst7", "../io/thrombin.top")

    m = s[MolIdx(1)]

    ala = m.search("resname /ala/i").join()[0]

    if verbose:
        print("Extracting alanine residues...")
        print(ala)

    ala_mol = ala.extract()

    assert(ala_mol.selectedAll())

    assert_equal(ala_mol.nAtoms(), ala.nAtoms())
    assert_equal(ala_mol.nResidues(), ala.nResidues())

    if verbose:
        print("Comparing energies...")

    ff = InternalFF()
    ff2 = InternalFF()

    ff.setStrict(True)
    ff2.setStrict(True)

    ff.add(ala_mol)
    ff2.add(ala)

    nrg = ff.energy().value()
    nrg2 = ff2.energy().value()

    if verbose:
        print("Energies: %s versus %s" % (nrg, nrg2))

    assert( abs(nrg-nrg2) < 0.001 )

if __name__ == "__main__":
    test_extract(True)

