
import pytest

def test_gil_exception_default_args():
    """Test that exceptions raised from functions that 
       have default arguments are caught correctly,
       without crashing the interpreter
    """
    try:
        from Sire import thumbs_up
    except Exception:
        print("Skipping as not supported in this version.")
        return

    from Sire.IO import AmberPrm

    with pytest.raises(IOError):
        AmberPrm.parse("file_does_not_exist")


if __name__ == "__main__":
    test_gil_exception_default_args()
