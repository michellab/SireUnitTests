try:
    import sire as sr

    sr.use_old_api()
except ImportError:
    pass


from Sire.Tools import resolveParameters, readParams, Parameter
from Sire.Units import *

distance_restraints_dict = Parameter(
    "distance restraints dictionary",
    {},
    """Dictionary of pair of atoms whose distance is restrained, and restraint
                                     parameters. Syntax is {(atom0,atom1):(reql, kl, Dl)} where atom0, atom1 are atomic
                                     indices. reql the equilibrium distance. Kl the force constant of the restraint.
                                     D the flat bottom radius. WARNING: PBC distance checks not implemented, avoid
                                     restraining pair of atoms that may diffuse out of the box.""",
)

params = readParams("configfile")

params["test1"] = "5 angstrom"
params["test2"] = "0.5 * gram/(centimeter*centimeter*centimeter)"
params["test3"] = 500


@resolveParameters
def _test_resolveparams():

    print(params)

    # this should work
    dic_items = list(distance_restraints_dict.val.items())

    print(dic_items)


def test_resolveparams(verbose=False):
    _test_resolveparams(params)


if __name__ == "__main__":
    test_resolveparams()
