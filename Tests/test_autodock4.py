#
# 
#
# $Id: test_autodock4.py,v 1.2 2005/08/25 15:08:05 rhuey Exp $
#
"""

"""


import os
import unittest
from AutoDockTools.Docking import Docking

computed_dlg = False
expected_value = -7.13


class Autodock4_1pgp_test(unittest.TestCase):
    """
    Test that autodock4 executes using an extremely short run
"""    
    
    def setUp(self):
        """Set up for autodock4 tests.
        Locate the autodock binary now during setUp.
        """
        global computed_dlg
        self.autodock = "../autodock4"
        dpf_filename = "1pgp.dpf"
        self.dlg_filename = "test_1pgp.dlg"
        command = "rm -f " + self.dlg_filename
        os.system(command)
        cmd_str = "%s -p %s -l %s" % \
                  (self.autodock, dpf_filename, self.dlg_filename)
        #print "\ncomputing new " + self.dlg_filename + ":\n"
        (i,o,e) = os.popen3(cmd_str) # trap all the outputs
        #print 'waiting...'
        os.wait() # for the child process to finish
        #print "after wait\n"
        computed_dlg = True


    def tearDown(self):
        #??do something here??
        computed_dlg = False


    def test_check_result_exists(self):
        """ check that run finished and a new dlg has been computed """
        self.assertEqual(computed_dlg, True)


    def test_check_result_energy(self):
        """ check the final energy is expected value """
        global expected_value
        d = Docking()
        d.readDlg(self.dlg_filename)
        c = d.ch.conformations[0]
        self.assertAlmostEqual(c.intermol_energy, expected_value)
        #self.assertAlmostEqual(c.intermol_energy, -7.13)


class Autodock4_1pgp_no_parameter_file_test(unittest.TestCase):
    """
    Test that autodock4 works using default parameter library
"""    
    
    def setUp(self):
        """Set up for autodock4 tests.
        Locate the autodock binary now during setUp.
        """
        global computed_dlg_no_parameter_library
        self.autodock = "../autodock4"
        dpf_filename = "1pgp_no_parameter_file.dpf"
        self.dlg_filename = "test_1pgp_no_parameter_file.dlg"
        command = "rm -f " + self.dlg_filename
        os.system(command)
        cmd_str = "%s -p %s -l %s" % \
                  (self.autodock, dpf_filename, self.dlg_filename)
        #print "\ncomputing new " + self.dlg_filename + ":\n"
        (i,o,e) = os.popen3(cmd_str) # trap all the outputs
        #print 'waiting...'
        os.wait() # for the child process to finish
        #print "after wait\n"
        computed_dlg_no_parameter_library = True


    def tearDown(self):
        #??do something here??
        computed_dlg_no_parameter_library = False


    def test_check_result_exists_default_parameter_file(self):
        """ using default parameter file: check that a run finished .... """
        self.assertEqual(computed_dlg_no_parameter_library, True)


    def test_check_result_energy_default_parameter_file(self):
        """ using default parameter file: check the final energy is expected value """
        global expected_value
        d = Docking()
        d.readDlg(self.dlg_filename)
        c = d.ch.conformations[0]
        self.assertAlmostEqual(c.intermol_energy, expected_value)
        #self.assertAlmostEqual(c.intermol_energy, -7.13)


if __name__ == '__main__':
    #this syntax allows you to run all tests 
    #or conveniently comment out tests you're not interested in
    test_cases = [
        'Autodock4_1pgp_test',
        'Autodock4_1pgp_no_parameter_file_test',
    ]
    unittest.main( argv=([__name__ ,] + test_cases) )
    #for verbose output use this:
    #unittest.main(argv=([__name__, '-v'] + test_cases))
