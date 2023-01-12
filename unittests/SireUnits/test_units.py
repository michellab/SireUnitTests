
try:
    import sire as sr
    sr.use_old_api()
except ImportError:
    pass

from Sire.Units import *
from Sire.Move import *

from Sire.Maths import pi

def _pvt_close( val0, val1 ):
   
   if not (abs(float(val1)-float(val0)) < 0.000001):
        print("FAIL! %f != %f" % (val0, val1))
        return False
   else:
        return True

def test_units(verbose=False):
    temp = convert(100, fahrenheit, celsius)

    if verbose:   
        print("100 F == %f C" % (100*fahrenheit).to(celsius))

    assert( _pvt_close(temp, (100-32)/1.8) )
    assert( _pvt_close(temp, (100*fahrenheit).to(celsius)) )

    mc = RigidBodyMC()

    mc.setTemperature( 100 * fahrenheit )

    if verbose:
        print(mc.temperature().to(celsius))

    assert( _pvt_close(mc.temperature().to(fahrenheit), 100) )

    k = (4 * pi * 8.854187817e-12 * farad / meter)

    if verbose:
        print(k)

    k = 1 / k

    if verbose:
        print(k)

    assert( _pvt_close(k.value(), 332.063710) )

if __name__ == "__main__":
    test_units(True)
