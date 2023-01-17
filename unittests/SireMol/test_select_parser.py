
try:
    import sire as sr
    sr.use_old_api()
except ImportError:
    pass

from Sire.Mol import *
from Sire.IO import *
from Sire.System import *

from nose.tools import assert_equal

try:
    s = Select()
    have_select = True
except:
    have_select = False

def _writeToFile(result, filename):
    s = System()

    m = MoleculeGroup("all")
    for r in result:
        m.add(r)

    s.add(m)

    p = PDB2(s)
    p.writeToFile(filename)

def test_selections(verbose=False):
    if not have_select:
        return

    if verbose:
        print("\nBase name selections\n")

    s = Select("atomnam hello")

    if verbose:
        print(s)

    s = Select("resname goodbye")

    if verbose:
        print(s)

    s = Select("chainname world")

    if verbose:
        print(s)

    s = Select("(atomnam bracket)")

    if verbose:
        print(s)

    s = Select("((atomnam doublebracket))")

    if verbose:
        print(s)

    s = Select("(atomnam one and resnam two) and (resnam three)")

    if verbose:
        print(s)

    s = Select("(atomnam one and resnam two) and (resnam three and atomname four)")

    if verbose:
        print(s)

    s = Select("atomnam hello and resname goodbye")

    if verbose:
        print(s)

    s = Select("atomnam hello and resname goodbye and chainname world")

    if verbose:
        print(s)

    s = Select("atomnam hello and resname goodbye and (chainname world and resname one)")

    if verbose:
        print(s)

    #s = Select("atomnam hello and (resname goodbye and (chainname world and (resname one and atomname two)))")

    #if verbose:
    #    print(s)

    if verbose:
        print("\nMulti name\n")

    s = Select("atomname C,CA , N,O")

    if verbose:
        print(s)

    s = Select("atomname 'C 1',O")

    if verbose:
        print(s)

    s = Select("atomname /hello/")

    if verbose:
        print(s)

    s = Select("atomname /goodbye/i")

    if verbose:
        print(s)

    s = Select("atomname C,/hello/,'C 1',/goodbye/i")

    if verbose:
        print(s)

    s = Select("(atomname /CA/i or atomname /C/i) and resname /ALA/i")

    if verbose:
        print(s)

    if verbose:
        print("\nNumbers...\n")

    s = Select("atomnum 5")

    if verbose:
        print(s)

    s = Select("resnum 1:10")

    if verbose:
        print(s)

    s = Select("cgidx 1:10:2, 5, 2:10")

    if verbose:
        print(s)

    s = Select("resnum 100 and chainname A")

    if verbose:
        print(s)

    s = Select("resnum > 5 and resnum <= 100")

    if verbose:
        print(s)

    if verbose:
        print("\nWith/In/not...\n")

    #s = Select("molecules with (resname /ala/i,/asp/i,/gly/i and atomname /ca/i)")

    #if verbose:
    #    print(s)

    s = Select("atoms in residx 1:10")

    if verbose:
        print(s)

    s = Select("not resname ALA")

    if verbose:
        print(s)

    s = Select("chainname B and not atomname CA,C,N,O")

    if verbose:
        print(s)

    if verbose:
        print(s)

    if verbose:
        print("\nComments\n")

    s = Select("resname /* This is a comment */ ALA")

    if verbose:
        print(s)

    if verbose:
        print("\nWithin distance\n")

    s = Select("molecules within 5.0 of resname ALA")

    if verbose:
        print(s)

    s = Select("atoms within 10 nm of 1,2,3")

    if verbose:
        print(s)

    s = Select("residues within 30 A of 1A,2nm,3pm")

    if verbose:
        print(s)

    if verbose:
        print("\nelement parsing\n")

    s = Select("element aluminium,iron,As,ca,iridium")

    if verbose:
        print(s)

    s = Select("element dummy")

    if verbose:
        print(s)

def test_engines(verbose=False):

    if verbose:
        print("\nReading molecules...")

    mols = MoleculeParser.read("../io/kigaki.gro", "../io/kigaki.top",
                               {"GROMACS_PATH":"../io/gromacs"})

    if verbose:
        print("Testing selections...")

    s = Select("atomname /HW\.*/")
    r = s(mols)

    if verbose:
        print(s)
        print(r)

    s = Select("atomname CA")
    r = s(mols)

    if verbose:
        print(s)
        print(r)

    s = Select("resname /ala/i")
    r = s(mols)

    if verbose:
        print(s)
        print(r)

    s = Select("resnum 1,3,5,7")
    r = s(mols)

    if verbose:
        print(s)
        print(r)

    s = Select("residx > -3")
    r = s(mols)

    if verbose:
        print(s)
        print(r)

    s = Select("atomidx 5:-1:3")
    r = s(mols)

    if verbose:
        print(s)
        print(r)

    s = Select("resname /ala/i,/lys/i,/leu/i and atomnam /ca/i")
    r = s(mols)

    if verbose:
        print(s)
        print(r)

    for atom in r:
        assert_equal( atom.name().value(), "CA" )
        resnam = atom.residue().name().value()

        assert( resnam == "ALA" or resnam == "LYS" or resnam == "LEU" )

    s = Select("atomnam /ca/i or atomnam /n\.*/i")
    r = s(mols)

    if verbose:
        print(s)
        print(r)

    s = Select("resnam /ala/i or atomnam /ca/i")
    r = s(mols)

    if verbose:
       print(s)
       print(r)

    s = Select("atomnam /ca/i or resnam /ala/i")
    r = s(mols)

    if verbose:
       print(s)
       print(r)

    s = Select("not resname /ala/i")
    r = s(mols)

    if verbose:
        print(s)
        print(r)

    assert_equal( Select("resname /ala/i")(r).isEmpty(), True )

    s = Select("not (resname /glu/i or atomname /ca/i)")
    r = s(mols)

    if verbose:
        print(s)
        print(r)

    assert_equal( Select("resname /glu/i or atomname /ca/i")(r).isEmpty(), True )

    s = Select("join (resname /glu/i or atomname /h\.*/i)")
    r = s(mols)

    if verbose:
        print(s)
        print(r)

    #s = Select("{atomname /h\.*/i}[-12:-1:2]")
    #r = s(mols)

    #if verbose:
    #    print(s)
    #    print(r)

    s = Select("resname /ala/i[0]")
    r = s(mols)

    if verbose:
        print(s)
        print(r)

    s = Select("atoms in resname /ala/i")
    r = s(mols)

    if verbose:
        print(s)
        print(r)

    #s = Select("residues in (resname /glu/i or atomname /c\.*/i)")
    #r = s(mols)

    #if verbose:
    #    print(s)
    #    print(r)

    s = Select("molecules with atomname /hw\.*/i")
    r = s(mols)

    if verbose:
        print(s)
        print(r)

    s = Select("resname /ala/i and element carbon,H,o")
    r = s(mols)

    if verbose:
        print(s)
        print(r)

    for atom in r:
        assert_equal( atom.residue().name().value(), "ALA" )

        element = atom.property("element")

        ok = (element == Element("C") or element == Element("H") or \
              element == Element("O"))

        assert_equal( ok, True )

    s = Select("atoms within 2 of resname /ala/i")
    r = s(mols)

    if verbose:
        print(s)
        print(r)

    _writeToFile(r,"test.pdb")

if __name__ == "__main__":
    test_selections(True)
    test_engines(True)
