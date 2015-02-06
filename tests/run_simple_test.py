import unittest

class Test(unittest.TestCase):
    def test_0(self):
        from pytraj import run_tests
        run_tests()

if __name__ == "__main__":
    unittest.main()
