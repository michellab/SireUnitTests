
from Sire.IO import *

def test_netcdf(verbose=False):

    r = AmberRst("../io/NA16.rst")

    print(r)

    print(r.creatorApplication())

    s = MoleculeParser.read("../io/NA16.rst", "../io/NA16.top")

    PDB().write(s.molecules(), "test.pdb")

if __name__ == "__main__":
    test_netcdf(True)

