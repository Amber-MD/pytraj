import unittest
from pytraj.actions.GridAction import GridAction

# TODO : more tests
class Test(unittest.TestCase):
    def test_1(self):
        ga = GridAction()
        print(ga)

if __name__ == "__main__":
    unittest.main()

