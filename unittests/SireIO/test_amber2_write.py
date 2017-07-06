
from Sire.IO import *
from Sire.Mol import *
from Sire.MM import *

def test_write(verbose=False):

    try:
        s = MoleculeParser.read
    except:
        return

    if verbose:
        print("Reading...")

    s = MoleculeParser.read("../io/ose.top", "../io/ose.crd")
    #s = MoleculeParser.read("../io/proteinbox.top", "../io/proteinbox.crd")

    #m = s[MolIdx(0)].molecule()
    #print(m.property("amberparameters").bornRadii())

    if verbose:
        print("Extracting...")

    p = AmberPrm(s)

    if verbose:
        print("Writing...")

    p.write("test.top")

if __name__ == "__main__":
    test_write(True)

