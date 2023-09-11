try:
    import sire as sr

    sr.use_old_api()
except ImportError:
    pass

from math import sqrt
from nose.tools import assert_almost_equal
from random import random, randint

from Sire.IO import *
from Sire.Maths import *
from Sire.MM import *
from Sire.Units import *
from Sire.Vol import *

# The number of coordinate groups.
num_groups = 10

# The number of coordinates within each group.
num_coords = 10

# Generate non-zero box lengths.
x = 0
while x == 0:
    x = randint(-10, 10)
y = 0
while y == 0:
    y = randint(-10, 10)
z = 0
while z == 0:
    z = randint(-10, 10)

# Create the box and store the absolute values of the lengths.
box = [x, y, z]
abs_box = [abs(x), abs(y), abs(z)]

# Create the box, sorted by absolute value.
box = [i for _, i in sorted(zip(abs_box, box))]

# Create a periodic box.
p = PeriodicBox(Vector(box))

# Create a triclinic box.
t = TriclinicBox(Vector(box[0], 0, 0), Vector(0, box[1], 0), Vector(0, 0, box[2]))

# Rotate and reduce the box. The rotation will ensure that the lattice vector
# with the shorest length is aligned with the x axis, the next shortest lies in
# the x-y plane, and the largest has a positive z component.
t.rotate()
t.reduce()

# Generate a vector of random coordinates.
cg_array = CoordGroupArray()

# Create num_groups coordinate groups.
for x in range(0, num_groups):
    coords = []

    # Randomise the center of the coordinate group.
    center = randint(-50, 50) * Vector(random(), random(), random())

    # Append num_coords random coordinates to each group.
    for y in range(0, num_coords):
        coords.append(
            Vector(1 - 2 * random(), 1 - 2 * random(), 1 - 2 * random())
            + center
        )
    cg_array.append(coords)


def test_distances(verbose=False):
    # Assert that the distance between all pairs of vectors in cg_array is the
    # same when computed in the periodic and triclinic space. For completeness,
    # we test distance, distance squared, and the distance vector.
    for a in range(0, num_groups):
        for b in range(0, num_coords):
            v0 = cg_array[a][b]
            for c in range(0, num_groups):
                for d in range(0, num_coords):
                    v1 = cg_array[c][d]

                    dp = p.calcDist(v0, v1)
                    dt = t.calcDist(v0, v1)
                    assert_almost_equal(dp, dt)

                    d2p = p.calcDist2(v0, v1)
                    d2t = t.calcDist2(v0, v1)
                    assert_almost_equal(d2p, d2t)

                    vp = p.calcDistVector(v0, v1)
                    vt = t.calcDistVector(v0, v1)
                    assert_almost_equal(vp.x(), vt.x())
                    assert_almost_equal(vp.y(), vt.y())
                    assert_almost_equal(vp.z(), vt.z())


def test_minimum_distance(verbose=False):
    # Assert that the minimum distance between coordinate groups and coordinate
    # group axis-aligned boxes are the same in both spaces.
    for a in range(0, num_groups):
        for b in range(0, num_groups):
            dp_array = []
            dt_array = []
            mp = p.minimumDistance(cg_array[a], cg_array[b])
            mt = t.minimumDistance(cg_array[a], cg_array[b])
            assert_almost_equal(mp, mt)
            mp = p.minimumDistance(cg_array[a].aaBox(), cg_array[b].aaBox())
            mt = t.minimumDistance(cg_array[a].aaBox(), cg_array[b].aaBox())
            assert_almost_equal(mp, mt)


def test_minimum_image(verbose=False):
    # Assert that the minimum image of a point with respect to center is the
    # same in both spaces.
    for a in range(0, num_groups):
        for b in range(0, num_coords):
            for c in range(0, 10):
                center = Vector(
                    50 - 100 * random(),
                    50 - 100 * random(),
                    50 - 100 * random(),
                )
                mp = p.getMinimumImage(cg_array[a][b], center)
                mt = t.getMinimumImage(cg_array[a][b], center)
                assert_almost_equal(mp.x(), mt.x())
                assert_almost_equal(mp.y(), mt.y())
                assert_almost_equal(mp.z(), mt.z())

    # Assert that the minimum image of a coordinate group with respect to center
    # is the same in both spaces.
    for a in range(0, num_groups):
        for b in range(0, 10):
            center = Vector(
                50 - 100 * random(), 50 - 100 * random(), 50 - 100 * random()
            )
            mp = p.getMinimumImage(cg_array[a], center)
            mt = t.getMinimumImage(cg_array[a], center)

            assert len(mp) == len(mt)

            for c in range(0, len(mp)):
                assert_almost_equal(mp[c].x(), mt[c].x())
                assert_almost_equal(mp[c].y(), mt[c].y())
                assert_almost_equal(mp[c].z(), mt[c].z())

    # Assert that the minimum image of a coordinate group array with respect to
    # center is the same in both spaces.
    for a in range(0, 10):
        center = Vector(
            50 - 100 * random(), 50 - 100 * random(), 50 - 100 * random()
        )

        # Whether to shift the entire array as one.
        for b in [False, True]:
            mp = p.getMinimumImage(cg_array, center, b)
            mt = t.getMinimumImage(cg_array, center, b)

            assert len(mp) == len(mt)

            for c in range(0, len(mp)):
                assert len(mp[c]) == len(mt[c])
                for d in range(0, len(mp[c])):
                    assert_almost_equal(mp[c][d].x(), mt[c][d].x())
                    assert_almost_equal(mp[c][d].y(), mt[c][d].y())
                    assert_almost_equal(mp[c][d].z(), mt[c][d].z())

    # Assert that the minimum image of a coordinate group axis-algined boxes
    # with respect to center is the same in both spaces.
    for a in range(0, num_groups):
        for b in range(0, 10):
            center = Vector(
                50 - 100 * random(), 50 - 100 * random(), 50 - 100 * random()
            )
            mp = p.getMinimumImage(cg_array[a].aaBox(), center)
            mt = t.getMinimumImage(cg_array[a].aaBox(), center)
            mp_min = mp.minCoords()
            mp_max = mp.maxCoords()
            mt_min = mt.minCoords()
            mt_max = mt.maxCoords()
            assert_almost_equal(mp_min.x(), mt_min.x())
            assert_almost_equal(mp_min.y(), mt_min.y())
            assert_almost_equal(mp_min.z(), mt_min.z())
            assert_almost_equal(mp_max.x(), mt_max.x())
            assert_almost_equal(mp_max.y(), mt_max.y())
            assert_almost_equal(mp_max.z(), mt_max.z())


def test_angles(verbose=False):
    # Assert that angles between points computed in both spaces are equivalent.
    for a in range(0, num_groups):
        for b in range(0, num_coords):
            v0 = cg_array[a][b]
            for c in range(0, num_coords):
                v1 = cg_array[a][c]
                for d in range(0, num_coords):
                    if b != c and b != d and c != d:
                        v2 = cg_array[a][d]
                        ap = p.calcAngle(v0, v1, v2)
                        at = t.calcAngle(v0, v1, v2)
                        assert_almost_equal(ap.value(), at.value())


def test_dihedrals(verbose=False):
    # Assert that dihedrals between points computed in both spaces are equivalent.
    for a in range(0, num_groups):
        for b in range(0, num_coords):
            v0 = cg_array[a][b]
            for c in range(0, num_coords):
                v1 = cg_array[a][c]
                for d in range(0, num_coords):
                    v2 = cg_array[a][d]
                    for e in range(0, num_coords):
                        if (
                            b != c
                            and b != d
                            and b != e
                            and c != d
                            and c != e
                            and d != e
                        ):
                            v3 = cg_array[a][e]
                            dp = p.calcDihedral(v0, v1, v2, v3)
                            dt = t.calcDihedral(v0, v1, v2, v3)
                            assert_almost_equal(dp.value(), dt.value())


def test_images_within(verbose=False):
    # Assert that the periodic images of point with respect to center within
    # dist are the same in both spaces.
    max_dist = max([abs(x), abs(y), abs(z)])
    for a in range(0, 10):
        point = Vector(
            x * (0.5 - random()), y * (0.5 - random()), z * (0.5 - random())
        )
        for b in range(0, 10):
            center = Vector(
                x * (0.5 - random()),
                y * (0.5 - random()),
                z * (0.5 - random()),
            )
            for c in range(0, 10):
                dist = max_dist * random()
                ip = p.getImagesWithin(point, center, dist)
                it = t.getImagesWithin(point, center, dist)
                assert len(ip) == len(it)
                for d in range(0, len(ip)):
                    assert_almost_equal(ip[d].x(), it[d].x())
                    assert_almost_equal(ip[d].y(), it[d].y())
                    assert_almost_equal(ip[d].z(), it[d].z())


def test_energy(verbose=False):
    # Assert that intermolecular energies calculated in both spaces agree.
    # Note that, other than using InterCLJFF, this comparison can only be
    # made between Cartesian triclinic spaces and a periodic box.

    (mols, periodic_space) = Amber().readCrdTop(
        "../io/waterbox.crd", "../io/waterbox.top"
    )

    dimensions = periodic_space.dimensions()

    triclinic_space = TriclinicBox(
        Vector(dimensions.x(), 0, 0),
        Vector(0, dimensions.y(), 0),
        Vector(0, 0, dimensions.z()),
    )

    long_coul_cutoff = 25 * angstrom
    coul_cutoff = 15 * angstrom
    lj_cutoff = 10 * angstrom

    switchfunc = HarmonicSwitchingFunction(
        coul_cutoff, coul_cutoff, lj_cutoff, lj_cutoff
    )

    p_oldff = InterCLJFF("p_oldff")
    p_oldff.setSwitchingFunction(switchfunc)
    p_oldff.setSpace(periodic_space)
    p_oldff.setShiftElectrostatics(True)
    p_oldff.add(mols)

    t_oldff = InterCLJFF("t_oldff")
    t_oldff.setSwitchingFunction(switchfunc)
    t_oldff.setSpace(triclinic_space)
    t_oldff.setShiftElectrostatics(True)
    t_oldff.add(mols)

    p_newff = InterFF("p_newff")
    p_newff.setProperty("switchingFunction", switchfunc)
    p_newff.setProperty("space", periodic_space)
    p_newff.add(mols)

    t_newff = InterFF("t_newff")
    t_newff.setProperty("switchingFunction", switchfunc)
    t_newff.setProperty("space", triclinic_space)
    t_newff.add(mols)

    p_cljfunc = CLJShiftFunction(coul_cutoff, lj_cutoff)
    p_cljfunc.setSpace(periodic_space)

    t_cljfunc = CLJShiftFunction(coul_cutoff, lj_cutoff)
    t_cljfunc.setSpace(triclinic_space)

    cljmols = CLJBoxes(CLJAtoms(mols.molecules()))
    cljcalc = CLJCalculator()
    (p_cnrg, p_ljnrg) = cljcalc.calculate(p_cljfunc, cljmols)
    (t_cnrg, t_ljnrg) = cljcalc.calculate(t_cljfunc, cljmols)

    p_oldnrgs = p_oldff.energies()
    p_oldcnrg = p_oldff.energy(p_oldff.components().coulomb()).value()
    p_oldljnrg = p_oldff.energy(p_oldff.components().lj()).value()

    t_oldnrgs = t_oldff.energies()
    t_oldcnrg = t_oldff.energy(t_oldff.components().coulomb()).value()
    t_oldljnrg = t_oldff.energy(t_oldff.components().lj()).value()

    p_newnrgs = p_newff.energies()
    p_newcnrg = p_newff.energy(p_newff.components().coulomb()).value()
    p_newljnrg = p_newff.energy(p_newff.components().lj()).value()

    t_newnrgs = t_newff.energies()
    t_newcnrg = t_newff.energy(t_newff.components().coulomb()).value()
    t_newljnrg = t_newff.energy(t_newff.components().lj()).value()

    assert_almost_equal(p_cnrg, t_cnrg)
    assert_almost_equal(p_ljnrg, t_ljnrg)
    assert_almost_equal(p_newcnrg, t_newcnrg)
    assert_almost_equal(p_newljnrg, t_newljnrg)
    assert_almost_equal(p_oldcnrg, t_oldcnrg)
    assert_almost_equal(p_oldljnrg, t_oldljnrg)


def test_equivalence(verbose=False):
    # Assert that TriclinicBox objects constructed using lattice vectors,
    # or lattice vector magnitudes and angles are equivalent. (Note that
    # other equivalent angles could be supplied, which would only change
    # the sign lattice vector components, not their magnitudes, e.g. an
    # equivalent left-tilting box.)

    # Cubic.
    t0 = TriclinicBox.cubic(1)
    t1 = TriclinicBox(1, 1, 1, 90 * degrees, 90 * degrees, 90 * degrees)
    assert_almost_equal(t0.vector0().x(), t1.vector0().x())
    assert_almost_equal(t0.vector0().y(), t1.vector0().y())
    assert_almost_equal(t0.vector0().z(), t1.vector0().z())
    assert_almost_equal(t0.vector1().x(), t1.vector1().x())
    assert_almost_equal(t0.vector1().y(), t1.vector1().y())
    assert_almost_equal(t0.vector1().z(), t1.vector1().z())
    assert_almost_equal(t0.vector2().x(), t1.vector2().x())
    assert_almost_equal(t0.vector2().y(), t1.vector2().y())
    assert_almost_equal(t0.vector2().z(), t1.vector2().z())

    # Rhombic-dodecahedron (square).
    t0 = TriclinicBox.rhombicDodecahedronSquare(1)
    t1 = TriclinicBox(1, 1, 1, 120 * degrees, 120 * degrees, 90 * degrees)
    assert_almost_equal(t0.vector0().x(), t1.vector0().x())
    assert_almost_equal(t0.vector0().y(), t1.vector0().y())
    assert_almost_equal(t0.vector0().z(), t1.vector0().z())
    assert_almost_equal(t0.vector1().x(), t1.vector1().x())
    assert_almost_equal(t0.vector1().y(), t1.vector1().y())
    assert_almost_equal(t0.vector1().z(), t1.vector1().z())
    assert_almost_equal(t0.vector2().x(), t1.vector2().x())
    assert_almost_equal(t0.vector2().y(), t1.vector2().y())
    assert_almost_equal(t0.vector2().z(), t1.vector2().z())

    # Rhombic-dodecahedron (hexagon).
    t0 = TriclinicBox.rhombicDodecahedronHexagon(1)
    t1 = TriclinicBox(1, 1, 1, 120 * degrees, 60 * degrees, 120 * degrees)
    assert_almost_equal(t0.vector0().x(), t1.vector0().x())
    assert_almost_equal(t0.vector0().y(), t1.vector0().y())
    assert_almost_equal(t0.vector0().z(), t1.vector0().z())
    assert_almost_equal(t0.vector1().x(), t1.vector1().x())
    assert_almost_equal(t0.vector1().y(), t1.vector1().y())
    assert_almost_equal(t0.vector1().z(), t1.vector1().z())
    assert_almost_equal(t0.vector2().x(), t1.vector2().x())
    assert_almost_equal(t0.vector2().y(), t1.vector2().y())
    assert_almost_equal(t0.vector2().z(), t1.vector2().z())


if __name__ == "__main__":
    test_distances(True)
    test_minimum_distance(True)
    test_minimum_image(True)
    test_angles(True)
    test_dihedrals(True)
    test_images_within(True)
    test_energy(True)
    test_equivalence(True)
