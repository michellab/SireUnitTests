
import pytest

def test_index_search():

    import Sire as sr

    try:
        mols = sr.load("../io/P38.pdb")
        mol = mols[0]
    except Exception:
        print("Skipping as older version of Sire")
        return

    print(mol.residues())

    res = mol["residx 1"]

    print(res)

    res2 = mol.residue(1)

    print(res2)

    assert res == res2

    print(res, res.atoms("atomname C, CA"))

    for atom in res.atoms("atomname C, CA"):
        print(atom.residue(), atom)

        resname = res.name().value()

        print(atom.residues(f"resname {resname}"))

        with pytest.raises(KeyError):
            print(atom.residues("resname ALA"))

    print(mol.residues("resname ALA"))


if __name__ == "__main__":
    test_index_search()

