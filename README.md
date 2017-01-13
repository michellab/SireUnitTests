# SireUnitTests

This is the repository that contains all of the unit tests for [Sire](http://siremol.org). This is held
separately from the main [Sire GitHub repository](https://github.com/michellab/Sire) as these tests are
updated on a seperate release cycle to Sire, and are downloaded dynamically by the sire_tests executable.
This uses `git pull` to pull the latest version of all of the unit tests every time it runs, thereby
enabling us to add new tests that will be performed by all currently installed versions of the code.

## Layout

The main directory is `unittests`. Inside this are subdirectories of all of the tests to perform. The
`sire_test` program (which is actually the [sire_test.py](https://github.com/michellab/Sire/blob/devel/wrapper/python/scripts/sire_test.py)),
will loop over each of these directories, and run all of the tests that are contained within each one.

Each test in each directory should be called `test_XXX.py`, where `XXX` is the thing you are testing. For example,
[unittests/SireMM/test_boxing.py](https://github.com/michellab/SireUnitTests/blob/devel/unittests/SireMM/test_boxing.py) is
testing the boxing algorithm within the SireMM library.

The tests in each `test_XXX.py` python file should have the following format;

```python

# Import the needed Sire libraries
from Sire.MM import *

# Import any needed utilities from nose (the testing framework)
from nose.tools import assert_almost_equal

# Each test should be in its own function, called test_YYY where YYY is the
# thing to be tested

def test_something(verbose=False):
    """Test function. If 'verbose' is False, then this test must not print to the screen"""    

    # the code to run the test. This should raise an exception
    # if the test fails. For example, 
    x = 3.5
    y = 9 * 0.5

    assert_almost_equal( x, y )

    # (this obviously fails!)

def test_something_else(verbose=False):
    """Another test function. This one should pass"""

    x = 5.2
    y = 2 * 2.6

    if verbose:
        print("Does %s equal %s?" % (x, y))

    assert_almost_equal( x, y )

# Finally, allow the test to be run manually from the command line.
# This should run the tests with 'verbose' set to True
if __name__ == "__main__":
    test_something(True)
    test_something_else(True)
```

Writing the test in this way will allow both for it to be automatically run by `sire_test`, and also for
it to be run from the command line, e.g. via

```
$SIRE/bin/python test_XXX.py
````

## How does sire_test connect to this repository?

The `sire_test` script should be able to detect the branch of Sire that you are using, and will automatically
clone that branch from the SireUnitTests repository. If you create a feature branch of SireUnitTests that has
the same name as your feature branch in Sire, then the feature branch copy of `sire_test` should automatically
get the right tests.

When you submit a pull request for your feature branch to Sire, then you should also submit a pull request for
SireUnitTests so that your tests can be merged in as well.

## What about formal releases of Sire?

Formal releases of Sire will still look for the tests using a `git clone`. Where `git` is not available on the
local machine, then `sire_test` will attempt to download a packaged version of the tests from the download
part of the [Sire website](http://siremol.org).

To create such a release package, outside the SireUnitTests directory type

```
tar -jcvf unittests_YEAR_RELEASE_PATH.tar.bz2 --exclude ".git/*" --exclude ".git" --exclude ".gitignore" SireUnitTests
```

replacing `YEAR_RELEASE_PATCH` with the respective year, release and patch number of the release (i.e. 2016.3.1).


This will create a `unittests_YEAR_RELEASE_PATCH.tar.bz2` file that can be uploaded by the release managers
to the Sire website.

## What are the C++ tests?

You may see that `sire_test` also runs some C++ tests. These are tests that are compiled directly into Sire, and
that are part of the C++ libraries. These are based on [SireBase::UnitTest](https://github.com/michellab/Sire/blob/devel/corelib/src/libs/SireBase/unittest.h). This provides a mechanism to easily create C++ unit tests that can be automatically located
and run by the C++ test harness.

To create a C++ test, add the following C++ code anywhere within Sire;

```c++

// include the unittest.h header file
#include "SireBase/unittest.h"

/** Function to run a test. The function must have the signature

    void function_name(bool verbose)

    where if 'verbose' is false, the test should not print anything to the screen 
*/
void test_something(bool verbose)
{
    // write your tests here

    // you can use SireBase::assert_equal. This compare the first
    // two arguments for equality, with CODELOC passed so that the
    // location in the code of this test can be printed by the raised
    // exception if this test fails
    SireBase::assert_equal( 3.0, 2 * 1.5, CODELOC );


    // There is also SireBase::assert_not_equal, SireBase::assert_nearly_equal,
    // SireBase::assert_true and SireBase::assert_false

    // This is also SireBase::assert_throws, which can be used to check that
    // an exception of the supplied type is thrown in the contained code, i.e.
    SireBase::assert_throws( [](){ throw SireError::program_bug(); }, SireError::program_bug(), CODELOC );
}

// Once you have written your test, you need to register it with the unit testing
// system using the macro SIRE_UNITTEST, i.e.
SIRE_UNITTEST( test_something )
```

Once you have added your test, it will be automatically registered with the registry in `SireBase::UnitTest`, and can
be queried and ran using the static functions attached to this class.
