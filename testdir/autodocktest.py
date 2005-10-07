#
# Last modified on Wed Feb 26 17:10:16 PST 2003 by lindy
#
# $Id: autodocktest.py,v 1.4 2003/02/27 01:18:14 lindy Exp $
#
"""
We assume here that the test is being run from the autodock
root directory that contains  the autodock source files.

USAGE: python2.2 testdir/autdocktest.py [-v]

"""

__version__ = "$Revision: 1.4 $"
#__test_for__ = '../autodock3"

import os
import types
import unittest

from AutoDockTools import Docking

class AutdockTestError(Exception):
    pass



class AutodockTestCase(unittest.TestCase):
    def setUp(self):
        """Set up for autodock3 tests.
        
        The current working directory will be different for each
        individual test. Locate the autodock binary now during setUp.
        """
        self.saved_cwd = os.getcwd()
        self.autodock = self.saved_cwd + '/autodock3'


    def tearDown(self):
        """Restore the current working directory"""
        os.chdir(self.saved_cwd)


    def run_cmd(self, cmd_str):
        """Fork the autodock3 command and wait for completion"""
        (i,o,e) = os.popen3(cmd_str) # trap all the outputs
        os.wait() # for the child process to finish


    def compare_dlgs(self, dlg1_filename, dlg2_filename):
        """Compare the conformations in two dlgs attribute by attribute
        """
        # read in the old and new docking logs
        d1 = Docking.Docking()
        d1.readDlg(dlg1_filename)
        d2 = Docking.Docking()
        d2.readDlg(dlg2_filename)
        # compare dlg's conformation-by-conformation
        self.assertEqual(len(d1.ch.conformations), len(d2.ch.conformations))
        for c1, c2 in zip(d1.ch.conformations, d2.ch.conformations):
            # compare one attribute,
            ##self.assertEqual(c1.binding_energy, c2.binding_energy)
            # or introspect and compare all attributes
            for attr in c1.__dict__.keys():
                #print "comparing %s (%s) %s" % (attr, type(getattr(c1, attr)),
                #                                repr(getattr(c1,attr)))
                if type(getattr(c1, attr)) != types.InstanceType:
                    self.assertEqual(getattr(c1, attr), getattr(c2, attr))


    def remove_previous_dlg(self, dlg_filename):
        """remove the previous docking log (if any)

        If autodock fails to create a new docking log, we don't want
        the tests to find and compare a previously written log.
        """
        try:
            os.remove(dlg_filename)
        except OSError, (error_number, message):
            # don't complain if the file does not exist
            import errno
            if error_number != errno.ENOENT:
                # pass onward the error
                raise OSError, (error_number, message)

# AutodockTestCase



class BasicTestCase(AutodockTestCase):

    def test_1hvr(self):
        """1hvr test case"""
        os.chdir('testdir/1hvr')
        dpf_filename = 'xk2A.1hvr.GASW.0.dpf'
        dlg_filename = 'xk2A.1hvr.GASW.0.dlg'
        old_dlg_filename = 'xk2A.1hvr.GASW.0.dlg.saved'
        # run autodock3
        self.remove_previous_dlg(dlg_filename)
        cmd_str = "%s -p %s -l %s" % \
                  (self.autodock, dpf_filename, dlg_filename)
        self.run_cmd(cmd_str)
        # compare resulting docking log with saved dlg
        self.compare_dlgs(dlg_filename, old_dlg_filename)


#
# These tests aren't ready yet...
#
##     def test_1phd(self):
##         """1phd test case"""
##         self.assertEqual(0, 1)
##     def test_1stp(self):
##         """1stp test case"""
##         self.assertEqual(0, 1)


    def test_2cpp(self):
        """2cpp test case"""
        # locate testdir and files
        os.chdir('testdir/2cpp')
        dpf_filename = 'cam.2cpp.GASW.dpf'
        dlg_filename = 'cam.2cpp.GASW.dlg'
        old_dlg_filename = 'cam.2cpp.GASW.dlg.saved'
        # run autodock3
        self.remove_previous_dlg(dlg_filename)
        cmd_str = "%s -p %s -l %s" % \
                  (self.autodock, dpf_filename, dlg_filename)
        self.run_cmd(cmd_str)
        # compare resulting docking log with saved dlg
        self.compare_dlgs(dlg_filename, old_dlg_filename)


    def test_2mcp(self):
        """2mcp test case"""
        # locate testdir and files
        os.chdir('testdir/2mcp')
        dpf_filename = 'pc.mcp2.GASW.dpf'
        dlg_filename = 'pc.mcp2.GASW.dlg'
        old_dlg_filename = 'pc.mcp2.GASW.dlg.saved'
        # run autodock3
        self.remove_previous_dlg(dlg_filename)
        cmd_str = "%s -p %s -l %s" % \
                  (self.autodock, dpf_filename, dlg_filename)
        self.run_cmd(cmd_str)
        # compare resulting docking log with saved dlg
        self.compare_dlgs(dlg_filename, old_dlg_filename)


    def test_3ptb(self):
        """3ptb test case"""
        # locate testdir and files
        os.chdir('testdir/3ptb')
        dpf_filename = 'benA.3ptb.GASW.dpf'
        dlg_filename = 'benA.3ptb.GASW.dlg'
        old_dlg_filename = 'benA.3ptb.GASW.dlg.saved'
        # run autodock3
        self.remove_previous_dlg(dlg_filename)
        cmd_str = "%s -p %s -l %s" % \
                  (self.autodock, dpf_filename, dlg_filename)
        self.run_cmd(cmd_str)
        # compare resulting docking log with saved dlg
        self.compare_dlgs(dlg_filename, old_dlg_filename)


    def test_4dfr(self):
        """4dfr test case"""
        # locate testdir and files
        os.chdir('testdir/4dfr')
        dpf_filename = 'mtxA.4dfr.GASW.dpf'
        dlg_filename = 'mtxA.4dfr.GASW.dlg'
        old_dlg_filename = 'mtxA.4dfr.GASW.dlg.saved'
        # run autodock3
        self.remove_previous_dlg(dlg_filename)
        cmd_str = "%s -p %s -l %s" % \
                  (self.autodock, dpf_filename, dlg_filename)
        self.run_cmd(cmd_str)
        # compare resulting docking log with saved dlg
        self.compare_dlgs(dlg_filename, old_dlg_filename)


    def test_4hmg(self):
        """4hmg test case"""
        # locate testdir and files
        os.chdir('testdir/4hmg')
        dpf_filename = 'sip.hmg4.GASW.dpf'
        dlg_filename = 'sip.hmg4.GASW.dlg'
        old_dlg_filename = 'sip.hmg4.GASW.dlg.saved'
        # run autodock3
        self.remove_previous_dlg(dlg_filename)
        cmd_str = "%s -p %s -l %s" % \
                  (self.autodock, dpf_filename, dlg_filename)
        self.run_cmd(cmd_str)
        # compare resulting docking log with saved dlg
        self.compare_dlgs(dlg_filename, old_dlg_filename)

# BasicTestCase


 
if __name__ == '__main__':
    unittest.main()
    # FYI
    test_ids = ['1hvr', '1phd', '1stp', '2cpp', '2mcp', '3ptb', '4dfr', '4hmg']
    seeds = ('863546106 2863')
