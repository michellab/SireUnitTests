
from Sire.IO import *

def test_netcdf(verbose=False):

    r = AmberRst("../io/NA16.rst")

    print(r)

    print(r.creatorApplication())
    print(r.time())

    s = MoleculeParser.read("../io/NA16.rst", "../io/NA16.top")

    r7 = AmberRst7(s)
    r7.writeToFile("test.rst7")

    print(r7)

    PDB().write(s.molecules(), "test.pdb")

if __name__ == "__main__":
    test_netcdf(True)

