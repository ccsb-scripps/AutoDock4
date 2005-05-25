#
# 
#
# $Id: test_autodock4.py,v 1.1 2005/05/25 19:03:18 rhuey Exp $
#
"""

"""


import os
import unittest
from AutoDockTools.Docking import Docking

computed_dlg = False

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
        """ check that a new dlg has been computed """
        self.assertEqual(computed_dlg, True)


    def test_check_result_energy(self):
        """ check the final energy """
        d = Docking()
        d.readDlg(self.dlg_filename)
        c = d.ch.conformations[0]
        self.assertAlmostEqual(c.intermol_energy, -7.13)


if __name__ == '__main__':
    unittest.main()
