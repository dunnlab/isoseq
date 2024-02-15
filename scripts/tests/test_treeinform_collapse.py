import unittest
import sys
sys.path.append("../")  # Add the parent directory to the Python path
from treeinform_collapse import *  # import the functions you want to test

class TestTreeinformCollapse(unittest.TestCase):

	def setUp(self):
		# This method will be called before each test. You can add setup code here.
		pass

	def tearDown(self):
		# This method will be called after each test. You can add cleanup code here.
		pass

	def test_function1(self):
		# Replace 'function1' with the name of the function you want to test
		# Call the function with some inputs, and assert the expected output
		pass

	def test_function2(self):
		# Replace 'function2' with the name of the function you want to test
		# Call the function with some inputs, and assert the expected output
		pass

# This allows the tests to be run if the file is run as a script, e.g. with `python -m unittest treeinform_collapse.py`
if __name__ == '__main__':
	unittest.main()