#
# $Id: test_autodock4.py,v 1.9 2006/11/25 10:00:39 garrett Exp $
#

"""
Unit Tests for AutoDock 4.
"""

#______________________________________________________________________________

import os
import unittest
from DlgParser import DlgParser

#______________________________________________________________________________
#
# Global variables

autodock_to_test = "../autodock4"

computed_dlg = False
computed_dlg_no_parameter_library = False
computed_dlg_no_elecmap = False
computed_dlg_no_desolvmap = False
computed_dlg_no_elec_desolv_maps = False

expected_intermol_energy = -6.17
expected_internal_energy = -1.58

#______________________________________________________________________________

def run_AutoDock( autodock_executable, dpf_filename, dlg_filename ):
    """Launch AutoDock, using the specified AutoDock executable and DPF,
    create the specified DLG, and trap all the outputs from standard output
    and standard error."""
    command = "rm -f " + dlg_filename
    os.system( command )
    command = "%s -p %s -l %s" % ( autodock_executable, dpf_filename, dlg_filename )
    print "\nRunning " + autodock_executable + ", saving results in \"" + dlg_filename + "\":"
    try:
        ( i, o, e ) = os.popen3( command ) # trap all the outputs
        os.wait() # for the child process to finish
        return True
    except:
        print "\nUnable to run " + autodock_executable + "."
        return False

#______________________________________________________________________________

def parse_energy_from_DLG( dlg_filename ):
    """Parse the AutoDock DLG, and return the intermolecular and internal
    energies as a tuple."""
    parser = DlgParser()
    parser.parse( dlg_filename )
    docked = parser.clist[0]  #dictionary of results
    intermol_energy = docked['intermol_energy']  #-6.17
    internal_energy = docked['total_internal']  # -1.58
    return ( intermol_energy, internal_energy )

#______________________________________________________________________________

def find_success_in_DLG( dlg_filename ):
    """Open the AutoDock DLG, and look for the string "Successful Completion"
    in the last 10 lines of the file."""
    fptr = open( dlg_filename )
    lines = fptr.readlines()
    fptr.close()
    success = False
    for l in lines[-10:]:
        if l.find( "Successful Completion" ) > -1:
            success = True
    return success

#______________________________________________________________________________

class AutoDock4_1pgp_test( unittest.TestCase ):
    """Test that autodock4 executes using an extremely short run."""
    def setUp( self ):
        """Set up for autodock4 tests. Locate the autodock binary now during setUp."""
        global computed_dlg
        self.dlg_filename = "test_1pgp.dlg"
        if computed_dlg is False:
            computed_dlg = run_AutoDock( autodock_to_test, "1pgp.dpf", self.dlg_filename )

    def test_check_result_exists( self ):
        """Check that run finished and a new dlg has been computed."""
        self.assertEqual( computed_dlg, True )

    def test_check_result_energy( self ):
        """Check the final energy is expected value."""
        global expected_intermol_energy, expected_internal_energy
        (intermol_energy, internal_energy) = parse_energy_from_DLG( self.dlg_filename )
        self.assertAlmostEqual( intermol_energy, expected_intermol_energy, places=6 )
        self.assertAlmostEqual( internal_energy, expected_internal_energy, places=6 )

#______________________________________________________________________________

class AutoDock4_1pgp_no_parameter_file_test( unittest.TestCase ):
    """Test that autodock4 works using default parameter library."""
    def setUp( self ):
        """Set up for autodock4 tests. Locate the autodock binary now during setUp."""
        global computed_dlg_no_parameter_library
        self.dlg_filename = "test_1pgp_no_parameter_file.dlg"
        if computed_dlg_no_parameter_library is False:
            computed_dlg_no_parameter_library = run_AutoDock( autodock_to_test, "1pgp_no_parameter_file.dpf", self.dlg_filename )

    def test_check_result_exists_default_parameter_file( self ):
        """Using default parameter file: check that a run finished... """
        self.assertEqual( computed_dlg_no_parameter_library, True )

    def test_check_result_energy_default_parameter_file( self ):
        """ check the final energy is expected value """
        global expected_intermol_energy, expected_internal_energy
        (intermol_energy, internal_energy) = parse_energy_from_DLG( self.dlg_filename )
        self.assertAlmostEqual( intermol_energy, expected_intermol_energy, places=6 )
        self.assertAlmostEqual( internal_energy, expected_internal_energy, places=6 )

#______________________________________________________________________________

class AutoDock4_1pgp_no_elecmap_test( unittest.TestCase ):
    """Test that autodock4 stops early if no "elecmap" keyword is specified."""
    def setUp( self ):
        """Set up for autodock4 tests. Locate the autodock binary now during setUp."""
        global computed_dlg_no_elecmap
        self.dlg_filename = "test_1pgp_no_elecmap.dlg"
        if computed_dlg_no_elecmap is False:
            computed_dlg_no_elecmap = run_AutoDock( autodock_to_test, "1pgp_no_elecmap.dpf", self.dlg_filename )

    def test_check_result_not_successful( self ):
        """Using parameter file with no 'elecmap': check that run does not reach Successful Completion... """
        success = find_success_in_DLG( self.dlg_filename )
        self.assertEqual( success, False )
        
#______________________________________________________________________________

class AutoDock4_1pgp_no_desolvmap_test( unittest.TestCase ):
    """Test that autodock4 stops early if no "desolvmap" keyword is specified."""
    def setUp( self ):
        """Set up for autodock4 tests. Locate the autodock binary now during setUp."""
        global computed_dlg_no_desolvmap
        self.dlg_filename = "test_1pgp_no_desolvmap.dlg"
        if computed_dlg_no_desolvmap is False:
            computed_dlg_no_desolvmap = run_AutoDock( autodock_to_test, "1pgp_no_desolvmap.dpf", self.dlg_filename )

    def test_check_result_not_successful( self ):
        """Using parameter file with no "desolvmap": 
        check that run does not reach Successful Completion... """
        success = find_success_in_DLG( self.dlg_filename )
        self.assertEqual( success, False )

#______________________________________________________________________________

class AutoDock4_1pgp_no_elec_desolv_maps_test( unittest.TestCase ):
    """Test that autodock4 stops early if no elecmap and no desolvmap 
    keywords are specified."""
    def setUp( self ):
        """Set up for autodock4 tests. Locate the autodock binary now during setUp."""
        global computed_dlg_no_elec_desolv_maps
        self.dlg_filename = "test_1pgp_no_elec_desolv_maps.dlg"
        if computed_dlg_no_elec_desolv_maps is False:
            computed_dlg_no_elec_desolv_maps = run_AutoDock( autodock_to_test, "1pgp_no_elec_desolv_maps.dpf", self.dlg_filename )

    def test_check_result_not_successful( self ):
        """Using parameter file with no 'elecmap' and no "desolvmap": 
        check that run does not reach Successful Completion..."""
        success = find_success_in_DLG( self.dlg_filename )
        self.assertEqual( success, False )

#______________________________________________________________________________

if __name__ == '__main__':
    #  This syntax lets us run all the tests,
    #  or conveniently comment out tests we're not interested in.
    #  NOTE:  Remember to add new TestCase class names to the list "test_cases"
    test_cases = [
        'AutoDock4_1pgp_test',
        'AutoDock4_1pgp_no_parameter_file_test',
        'AutoDock4_1pgp_no_elecmap_test',
        'AutoDock4_1pgp_no_desolvmap_test',
        'AutoDock4_1pgp_no_elec_desolv_maps_test',
    ]
    unittest.main( argv=( [__name__ ,] + test_cases ) )
    #  The call "unittest.main()" automatically runs all the TestCase classes in
    #  alphabetical order; calling with argv=([]), lets us specify the order.
    #  NOTE: "unittest.main()" saves us having to remember to add new tests to the 
    #  list of test cases.
    #unittest.main()
    #  For verbose output, use this:
    #unittest.main( argv=( [__name__, '-v'] + test_cases ) )
