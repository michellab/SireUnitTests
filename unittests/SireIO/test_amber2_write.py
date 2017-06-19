
from Sire.IO import *

def test_write(verbose=False):

    try:
        s = MoleculeParser.read
    except:
        return

    #s = MoleculeParser.read("../io/ose.top", "../io/ose.crd")

    if verbose:
        print("Reading...")

    s = MoleculeParser.read("../io/proteinbox.top", "../io/proteinbox.crd")

    if verbose:
        print("Extracting...")

    p = AmberParm(s)

    if verbose:
        print("Writing...")

    p.write("test.top")

if __name__ == "__main__":
    test_write(True)

