from __future__ import with_statement, absolute_import
import sys
from os import name as __name
available_tests = ['unit_tests']

del __name

from os.path import abspath as __abspath
from os.path import join as __join
from os.path import split as __split
from os.path import sep as __sep

nampy_directory = __abspath(__join(__split(__abspath(__file__))[0], ".."))
nampy_location = __abspath(__join(nampy_directory, ".."))
data_directory = nampy_directory + "/data/"
nampy_directory += '/core/'

humannet_filename = __join(data_directory, "HumanNet_v1_join_networkonly.txt")
recon_1_filename = __join(data_directory, "recon_1_nampy_v01")

del __abspath, __join, __split, __sep

def create_test_model(test_network_file = humannet_filename):
    """Returns a nampy model for testing.  

    test_network_file: a two-column network file to be read
    
    """
    # from os import name as __name
    from nampy.networkio import networkio
    
    the_nampy_model = networkio.create_network_model_from_textfile('humannet', test_network_file, verbose = False)
    
    return the_nampy_model


def create_test_suite():
    """create a unittest.TestSuite with available tests"""
    from unittest import TestLoader, TestSuite
    loader = TestLoader()
    suite = TestSuite()
    for test_name in available_tests:
        exec("from . import " + test_name)
        suite.addTests(loader.loadTestsFromModule(eval(test_name)))
    return suite

suite = create_test_suite()

def test_all():
    """###running unit tests on nampy###"""
    from unittest import TextTestRunner
    TextTestRunner(verbosity=2).run(create_test_suite())
