
def test_names_numbers():
    try:
        import Sire as sr
        mols = sr.load("../io/P38.pdb")
        mol = mols[0]
    except Exception:
        print("Skipping as older version of Sire")

    assert mol.nResidues() == len(mol.residues().names())

    for i, name in enumerate(mol.residues().names()):
        assert name == mol.residue(i).name()

    assert mol.nResidues() == len(mol.residues().numbers())

    for i, number in enumerate(mol.residues().numbers()):
        assert number == mol.residue(i).number()


if __name__ == "__main__":
    test_names_numbers()

