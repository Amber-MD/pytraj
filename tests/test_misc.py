import unittest
from pytraj.base import *
from pytraj import io as mdio
from pytraj.utils.check_and_assert import assert_almost_equal
from pytraj.misc import action_dict

def main():
    import sys
    try:
        key = sys.argv[1]
        print (key)
        action_dict[key]().help()
    except:
        print(list(action_dict.keys()))

if __name__ == "__main__":
    main()
