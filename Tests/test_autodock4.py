#
# $Id: test_autodock4.py,v 1.6 2006/06/01 00:06:55 garrett Exp $
#
"""
Unit Tests for AutoDock 4.
"""

import os
import unittest
from DlgParser import DlgParser

computed_dlg = False
computed_dlg_no_parameter_library = False
expected_value = -6.17
expected_internal_energy_value = -1.58


class Autodock4_1pgp_test(unittest.TestCase):
    """Test that autodock4 executes using an extremely short run."""
    def setUp(self):
        """Set up for autodock4 tests. Locate the autodock binary now during setUp."""
        global computed_dlg
        self.dlg_filename = "test_1pgp.dlg"
        if computed_dlg is False:
            self.autodock = "../autodock4"
            dpf_filename = "1pgp.dpf"
            command = "rm -f " + self.dlg_filename
            os.system(command)
            cmd_str = "%s -p %s -l %s" % \
                      (self.autodock, dpf_filename, self.dlg_filename)
            print "\ncomputing new " + self.dlg_filename + ":"
            (i,o,e) = os.popen3(cmd_str) # trap all the outputs
            os.wait() # for the child process to finish
            computed_dlg = True

    def test_check_result_exists(self):
        """Check that run finished and a new dlg has been computed."""
        self.assertEqual(computed_dlg, True)

    def test_check_result_energy(self):
        """Check the final energy is expected value."""
        global expected_value, expected_internal_energy_value
        parser = DlgParser()
        parser.parse(self.dlg_filename)
        docked = parser.clist[0]  #dictionary of results
        intermol_energy = docked['intermol_energy']  #-6.17
        #d = Docking()
        #d.readDlg(self.dlg_filename)
        #c = d.ch.conformations[0]
        #self.assertAlmostEqual(c.intermol_energy, -6.17, places=6)
        #self.assertAlmostEqual(c.intermol_energy, expected_value, places=6)
        self.assertAlmostEqual(intermol_energy, expected_value, places=6)
        internal_energy = docked['internal_energy']  # -1.58
        self.assertAlmostEqual(internal_energy, expected_internal_energy_value, places=6)


class Autodock4_1pgp_no_parameter_file_test(unittest.TestCase):
    """Test that autodock4 works using default parameter library."""
    def setUp(self):
        """Set up for autodock4 tests. Locate the autodock binary now during setUp."""
        global computed_dlg_no_parameter_library
        self.dlg_filename = "test_1pgp_no_parameter_file.dlg"
        if computed_dlg_no_parameter_library is False:
            self.autodock = "../autodock4"
            dpf_filename = "1pgp_no_parameter_file.dpf"
            command = "rm -f " + self.dlg_filename
            os.system(command)
            cmd_str = "%s -p %s -l %s" % \
                      (self.autodock, dpf_filename, self.dlg_filename)
            print "\ncomputing new " + self.dlg_filename + ":"
            (i,o,e) = os.popen3(cmd_str) # trap all the outputs
            os.wait() # for the child process to finish
            computed_dlg_no_parameter_library = True

    def test_check_result_exists_default_parameter_file(self):
        """Using default parameter file: check that a run finished... """
        self.assertEqual(computed_dlg_no_parameter_library, True)

    def test_check_result_energy_default_parameter_file(self):
        """ check the final energy is expected value """
        global expected_value, expected_internal_energy_value
        parser = DlgParser()
        parser.parse(self.dlg_filename)
        docked = parser.clist[0]  #dictionary of results
        intermol_energy = docked['intermol_energy']  #-6.17
        #d = Docking()
        #d.readDlg(self.dlg_filename)
        #c = d.ch.conformations[0]
        #self.assertAlmostEqual(c.intermol_energy, -6.17, places=6)
        #self.assertAlmostEqual(c.intermol_energy, expected_value, places=6)
        self.assertAlmostEqual(intermol_energy, expected_value, places=6)
        internal_energy = docked['internal_energy']  # -1.58
        self.assertAlmostEqual(internal_energy, expected_internal_energy_value, places=6)


if __name__ == '__main__':
    #this syntax allows you to run all tests or conveniently comment out tests you're not interested in
    test_cases = [
        'Autodock4_1pgp_test',
        'Autodock4_1pgp_no_parameter_file_test',
    ]
    unittest.main( argv=([__name__ ,] + test_cases) )
    #for verbose output use this:
    #unittest.main(argv=([__name__, '-v'] + test_cases))
