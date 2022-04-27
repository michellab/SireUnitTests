
def test_regexp_search():
    try:
        import Sire as sr
        mols = sr.load("../io/p38.pdb")
        mol = mols[0]
    except Exception:
        print("Skipping as older version of Sire")
        return

    for atom in mol["atomname /N?[0-9]/"]:
        n = atom.name().value()

        assert n.startswith("N")
        assert int(n[-1]) <= 9

    for atom in mol["atomname /n*/i"]:
       assert atom.name().value().startswith("N")


if __name__ == "__main__":
    test_regexp_search()

