try:
    import sire as sr

    sr.use_old_api()
except ImportError:
    pass

from Sire.Base import *
from Sire.Maths import *
from Sire.Qt import *


def test_linearap(verbose=False):

    costs = NMatrix(10, 10)

    rand = RanGenerator()

    if verbose:
        print("Generating the costs...")

    t = QTime()
    t.start()

    for i in range(0, costs.nRows()):
        for j in range(0, costs.nColumns()):
            d = rand.rand(0, 5)
            costs.set(i, j, d * d)

    ms = t.elapsed()

    if verbose:
        print("...took %d ms" % ms)
        print("\nFinding the combination with the lowest total cost...")

    t.start()

    rows_to_columns = solve_linear_assignment(costs)

    ms = t.elapsed()

    t.start()

    solve_linear_assignment(costs, True)

    ms2 = t.elapsed()

    if verbose:
        print("\nSolution:")
        print(rows_to_columns)

        print("\nSolution took %d ms (%d ms with checking)" % (ms, ms2))

        print("Total cost = %f" % calculate_total_cost(costs, rows_to_columns))

    t.start()
    rows_to_columns2 = brute_force_linear_assignment(costs)
    ms = t.elapsed()

    if verbose:
        print("\nBrute force solution:")
        print(rows_to_columns2)

        print(
            "Total cost = %f" % calculate_total_cost(costs, rows_to_columns2)
        )

        print("\nSolution took %d ms" % ms)

    assert rows_to_columns == rows_to_columns2


if __name__ == "__main__":
    test_linearap(True)
