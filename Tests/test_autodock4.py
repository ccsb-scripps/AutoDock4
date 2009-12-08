#! /usr/bin/env python
#
# $Id: test_autodock4.py,v 1.29 2009/12/08 18:06:20 rhuey Exp $
#

"""
Unit Tests for AutoDock 4.
"""

#______________________________________________________________________________

import sys
import os
import unittest
import getopt
import subprocess
from DlgParser import DlgParser

#______________________________________________________________________________
#
# Global variables

autodock_executable = "../autodock4" # where the AutoDock executable resides
dpf_directory = '.' # where the input DPF files reside
test_output_directory = '.' # where the DLG files will be written

try:
    opts, argv = getopt.getopt(sys.argv[1:], "d:e:o:",
    ["dpf-directory=","executable=","test-output-directory="]) 
except getopt.GetoptError, v:
    usage() 
    sys.exit(2)

for o,a in opts:
    if o in ("-d", "--dpf-directory"):
        dpf_directory = a
    if o in ("-e", "--executable"):
        autodock_executable = a
    if o in ("-o","--test-output-directory"):
        test_output_directory = a


#______________________________________________________________________________

def usage():
    """Print out the usage of this command."""
    print """Usage:  python test_autodock4.py [-d <string>] [-e <string>] [-o <string>]

where:
    -d, --dpf-directory
        specifies the directory containing the DPFs to be tested;
        this flag is optional; default is '.'
    -e, --executable
        specifies the path to the AutoDock executable to be tested;
        this flag is optional; default is '../autodock4'
    -o, --test-output-directory
        specifies the directory where the output DLGs will be written;
        this flag is optional; default is '.'

NOTE:  these may be relative to the directory where this script was invoked.
"""

#______________________________________________________________________________

def run_AutoDock( dpf_filename, dlg_filename ):
    """Launch AutoDock, using the specified AutoDock executable and DPF,
    create the specified DLG, and trap all the outputs from standard output
    and standard error."""
    dpf = dpf_directory + os.sep + dpf_filename
    dlg = test_output_directory + os.sep + dlg_filename
    command = "rm -f " + dlg
    os.system( command )
    command =   [autodock_executable, '-p', dpf, '-l', dlg] 
    print '\nRunning ' + autodock_executable + ' using DPF "'+dpf+'", saving results in "'+dlg+'":'
    try:
        #( i, o, e ) = os.popen3( command ) # trap all the outputs
        subprocess.call( command )
        # TODO os.wait seems to return (pid, exit_status)
        #os.wait() # for the child process to finish
        # return True # this should really be os.wait()'s exit_status
        return find_success_in_DLG( dlg_filename )
    except:
        print "\nUnable to run " + autodock_executable + "."
        return False

#______________________________________________________________________________

def parse_energy_from_DLG( dlg_filename, energy_list):
    """Parse the AutoDock DLG, and return the intermolecular and internal
    energies as a tuple."""
    parser = DlgParser()
    dlg = test_output_directory + os.sep + dlg_filename
    parser.parse( dlg )
    docked = parser.clist[0]  #dictionary of results
    result = []
    for energy_type in energy_list:
        newVal = docked.get(energy_type, 'ERROR')
        print energy_type, ' is now ', newVal
        result.append(docked.get(energy_type, 'ERROR'))
    #intermol_energy = docked['intermol_energy']  #-6.17
    #internal_energy = docked['total_internal']  # -1.58
    #print "docked[binding_energy]=", docked['binding_energy']
    #print "docked[electrostatic_energy]=", docked['electrostatic_energy']
    #print "docked[intermol_energy]=", docked['intermol_energy']
    #print "docked[total_internal]=", docked['total_internal']
    #unbound_energy = docked['unbound_energy']
    #print "unbound_energy=", unbound_energy
    #return ( intermol_energy, internal_energy )
    return result

#______________________________________________________________________________

def find_success_in_DLG( dlg_filename ):
    """Open the AutoDock DLG, and look for the string "Successful Completion"
    in the last 10 lines of the file."""
    dlg = test_output_directory + os.sep + dlg_filename
    try:
        fptr = open( dlg )
        lines = fptr.readlines()
        fptr.close()
        success = False
        for l in lines[-10:]:
            if l.find( "Successful Completion" ) > -1:
                success = True
        return success
    except:
        return False

#______________________________________________________________________________

class AutoDock_base_test( unittest.TestCase ):
    """Base Class for AutoDock testing."""
    dpf_stem = "BaseClass"
    computed = False
    def setUp( self ):
        """Set up for autodock4 tests. Locate the autodock binary now during setUp."""
        self.dlg_filename = "test_" + self.dpf_stem + ".dlg"
        self.computed = run_AutoDock( self.dpf_stem + ".dpf", self.dlg_filename )

    #def test_dlg_exists( self ):
    #    """Check that run finished and a new DLG has been computed."""
    #    # Check that run finished and a new DLG has been computed.
    #    if (self.expected_outcome == True ):
    #        print "Testing that DLG exists and AutoDock successfully completed."
    #    else:
    #        print "Testing that DLG exists and AutoDock did not complete."
    #    self.assertEqual( self.computed, self.expected_outcome )


#______________________________________________________________________________

class AutoDock_simple_test( unittest.TestCase ):
    """Base Class for AutoDock testing."""
    dpf_stem = "BaseClass"
    computed = False
    def setUp( self ):
        """Set up for autodock4 tests. Locate the autodock binary now during setUp."""
        self.dlg_filename = "test_" + self.dpf_stem + ".dlg"
        self.computed = run_AutoDock( self.dpf_stem + ".dpf", self.dlg_filename )

    def test_dlg_exists( self ):
        """Check that run finished and a new DLG has been computed."""
        # Check that run finished and a new DLG has been computed.
        if (self.expected_outcome == True ):
            print "Testing that DLG exists and AutoDock successfully completed."
        else:
            print "Testing that DLG exists and AutoDock did not complete."
        self.assertEqual( self.computed, self.expected_outcome )
#______________________________________________________________________________

class AutoDock4_1pgp_no_extension( AutoDock_simple_test ):
    """Test that autodock4 stops early if .dpf extension is missing
    keywords are specified."""
    dpf_stem = "1pgp_no_extension"
    print "in 1pgp_no_extension"
    expected_outcome = False # True means Successful Completion!
    def setUp( self ):
        """Set up for autodock4 tests. Locate the autodock binary now during setUp."""
        self.dlg_filename = "test_" + self.dpf_stem + ".dlg"
        self.computed = run_AutoDock( self.dpf_stem , self.dlg_filename )

#______________________________________________________________________________

class AutoDock4_1pgp_wrong_extension( AutoDock_simple_test ):
    """Test that autodock4 stops early if extension is not '.dpf'
    keywords are specified."""
    dpf_stem = "1pgp.fpd"
    print "in 1pgp_wrong_extension"
    expected_outcome = False # True means Successful Completion!
    def setUp( self ):
        """Set up for autodock4 tests. Locate the autodock binary now during setUp."""
        self.dlg_filename = "test_" + self.dpf_stem + ".dlg"
        self.computed = run_AutoDock( self.dpf_stem , self.dlg_filename )

#______________________________________________________________________________


class AutoDock4_1pgp_two_extensions( AutoDock_simple_test ):
    """Test that autodock4 stops early if dpf name includes two .dpf
    keywords are specified."""
    dpf_stem = "1pgp.dpf"
    print "in 1pgp_two_extensions"
    expected_outcome = False # True means Successful Completion!
    def setUp( self ):
        """Set up for autodock4 tests. Locate the autodock binary now during setUp."""
        self.dlg_filename = "test_" + self.dpf_stem + ".dlg"
        self.computed = run_AutoDock( self.dpf_stem +'.dpf', self.dlg_filename )


#______________________________________________________________________________


class AutoDock4_1pgp_ligand_types_map_mismatch( AutoDock_simple_test ):
    """Test that autodock4 stops early if number of maps do not equal number
    of ligand types"""
    dpf_stem = "1pgp_ligand_types_map_mismatch"
    expected_outcome = False # True means Successful Completion!
#______________________________________________________________________________

class AutoDock4_1pgp_illegal_keyword_test( AutoDock_simple_test ):
    """Test that autodock4 stops early if it finds an illegal keyword 
    in dpf """
    dpf_stem = "1pgp_illegal_keyword"
    expected_outcome = False # True means Successful Completion!
#______________________________________________________________________________

class AutoDock4_1pgp_no_elecmap_test( AutoDock_simple_test ):
    """Test that autodock4 stops early if no "elecmap" keyword is specified."""
    dpf_stem = "1pgp_no_elecmap"
    expected_outcome = False # True means Successful Completion!
#______________________________________________________________________________

class AutoDock4_1pgp_no_desolvmap_test( AutoDock_simple_test ):
    """Test that autodock4 stops early if no "desolvmap" keyword is specified."""
    dpf_stem = "1pgp_no_desolvmap"
    expected_outcome = False # True means Successful Completion!
#______________________________________________________________________________

class AutoDock4_1pgp_no_elec_desolv_maps_test( AutoDock_simple_test ):
    """Test that autodock4 stops early if no elecmap and no desolvmap 
    keywords are specified."""
    dpf_stem = "1pgp_no_elec_desolv_maps"
    expected_outcome = False # True means Successful Completion!
#______________________________________________________________________________

class AutoDock4_1pgp_too_many_torsions( AutoDock_simple_test ):
    """Test that autodock4 stops early if too many torsions 
    are specified. (current limit is 32)"""
    dpf_stem = "1pgp_too_many_torsions"
    expected_outcome = False # True means Successful Completion!
#______________________________________________________________________________

class AutoDock4_1pgp_just_right_number_torsions( AutoDock_simple_test ):
    """Test that autodock4 completes with current limit of number of  torsions 
    are specified. (current limit is 32)"""
    dpf_stem = "1pgp_just_right_number_torsions"
    expected_outcome = True # True means Successful Completion!
#______________________________________________________________________________


class AutoDock4_1pgp_too_many_ligand_types_test( AutoDock_simple_test ):
    """Test that autodock4 stops early if too many ligand types 
    are specified. (current limit is 14)"""
    dpf_stem = "1pgp_too_many_ligand_types"
    expected_outcome = False # True means Successful Completion!
#______________________________________________________________________________

class AutoDock4_1pgp_two_ligands_test( AutoDock_simple_test ):
    """Test that autodock4 can run dpf specifying two ligands."""
    dpf_stem = "1pgp_two_ligands"
    expected_outcome = True # True means Successful Completion!
#______________________________________________________________________________

class AutoDock4_1pgp_two_mapsets_test( AutoDock_simple_test ):
    """Test that autodock4 can run dpf specifying two sets of maps and one ligand."""
    dpf_stem = "1pgp_two_mapsets"
    expected_outcome = True # True means Successful Completion!
#______________________________________________________________________________

class AutoDock4_1pgp_unbound_set_illegal_test( AutoDock_simple_test ):
    """Test that autodock 4.1 works when unbound is set to 'foo' in the DPF."""
    dpf_stem = "1pgp_unbound_set_illegal"
    expected_outcome = False # True means Successful Completion!
#______________________________________________________________________________

class AutoDock4_1pgp_unbound_model_illegal_test( AutoDock_simple_test ):
    """Test that autodock4 stops early if it finds an illegal unbound_model 
    in dpf """
    dpf_stem = "1pgp_unbound_model_illegal"
    expected_outcome = False # True means Successful Completion!
#______________________________________________________________________________

class AutoDock4_1pgp_ga_select_tournament_test( AutoDock_simple_test ):
    """Test that autodock 4.1 works when ga_select_tournament is set in the DPF."""
    dpf_stem = "1pgp_ga_select_tournament"
    expected_outcome = True # True means Successful Completion!

#______________________________________________________________________________

class AutoDock_test( AutoDock_base_test ):
    """Class for AutoDock testing."""

    def test_dlg_exists_and_test_energy( self ):
        """Check that run finished and a new DLG has been computed.
        Also check the final energy is the expected value."""
        # Check that run finished and a new DLG has been computed.
        if (self.expected_outcome == True ):
            print "Testing that DLG exists and AutoDock successfully completed."
        else:
            print "Testing that DLG exists and AutoDock did not complete."
        self.assertEqual( self.computed, self.expected_outcome )
        # Check the final energy is expected value.
        expected_intermol_energy = -6.17
        expected_internal_energy = -1.80
        (intermol_energy, internal_energy) = parse_energy_from_DLG( self.dlg_filename, ['intermol_energy','total_internal'] )
        print "Testing that intermolecular energy = %.2f kcal/mol." % (expected_intermol_energy,)
        self.assertEqual( round(intermol_energy,6), round(expected_intermol_energy,6))
        print "Testing that internal energy = %.2f kcal/mol." % (expected_internal_energy,)
        self.assertEqual( round(internal_energy,6), round(expected_internal_energy,6))
#______________________________________________________________________________

class AutoDock4_1pgp_test( AutoDock_test ):
    """Test that autodock4 executes using an extremely short run."""
    dpf_stem = "1pgp"
    expected_outcome = True # True means Successful Completion!
#______________________________________________________________________________

class AutoDock4_1pgp_smaller_test( AutoDock_test ):
    """Test that autodock4 executes using fewer parameters and an extremely short run."""
    dpf_stem = "1pgp_smaller"
    expected_outcome = True # True means Successful Completion!
#______________________________________________________________________________

class AutoDock4_1pgp_no_parameter_file_test( AutoDock_test ):
    """Test that autodock4 works using default parameter library."""
    dpf_stem = "1pgp_no_parameter_file"
    expected_outcome = True # True means Successful Completion!
#______________________________________________________________________________

class AutoDock4_unbound_test( AutoDock_base_test ):
    """Class for AutoDock testing unbound energy."""
    expected_unbound_energy = None

    def test_dlg_exists_and_test_energy( self):
        """Check that run finished and a new DLG has been computed.
        Also check the final energy is the expected value."""
        # Check that run finished and a new DLG has been computed.
        if (self.expected_outcome == True ):
            print "Testing that DLG exists and AutoDock successfully completed."
        else:
            print "Testing that DLG exists and AutoDock did not complete."
        self.assertEqual( self.computed, self.expected_outcome )
        # Check the final energy is expected value.
        #expected_unbound_energy = -1.80
        (unbound_energy) = parse_energy_from_DLG( self.dlg_filename, ['unbound_energy'])[0]
        print "Testing that unbound energy = %.2f kcal/mol." % (self.expected_unbound_energy,)
        print "unbound_energy=", unbound_energy
        self.assertEqual( round(unbound_energy,6), round(self.expected_unbound_energy,6))
#______________________________________________________________________________

class AutoDock4_1pgp_unbound_default_test( AutoDock4_unbound_test ):
    """Test that autodock 4.1 works when unbound is NOT set in the DPF."""
    dpf_stem = "1pgp_unbound_default"
    expected_unbound_energy = -1.80
    expected_outcome = True # True means Successful Completion!
#______________________________________________________________________________

class AutoDock4_1pgp_unbound_model_extended( AutoDock4_unbound_test ):
    """Test that autodock 4.1 works when unbound_model is set to extended."""
    dpf_stem = "1pgp_unbound_model_extended"
    expected_unbound_energy = -0.28
    #expected_unbound_energy = -0.66 #prior to 4/2009
    expected_outcome = True # True means Successful Completion!
#______________________________________________________________________________

class AutoDock4_1pgp_unbound_compute_unbound_extended( AutoDock4_unbound_test ):
    """Test that autodock 4.1 works when unbound_model is set to extended."""
    dpf_stem = "1pgp_unbound_compute_unbound_extended"
    expected_unbound_energy = -0.28
    #expected_unbound_energy = -0.66 #prior to 4/2009
    expected_outcome = True # True means Successful Completion!
#______________________________________________________________________________

class AutoDock4_1pgp_unbound_model_value( AutoDock4_unbound_test ):
    """Test that autodock 4.1 works when unbound_model is set to a value in the DPF."""
    dpf_stem = "1pgp_unbound_model_value"
    expected_unbound_energy = -3.12 
    expected_outcome = True # True means Successful Completion!
#______________________________________________________________________________


class AutoDock4_1pgp_unbound_model_compact( AutoDock4_unbound_test ):
    """Test that autodock 4.1 works when unbound_model is set to compact."""
    dpf_stem = "1pgp_unbound_model_compact"
    expected_unbound_energy =  0.00 #@FixMe 3/2009 do not know how to calc this
    expected_outcome = True # True means Successful Completion!
#______________________________________________________________________________

class AutoDock4_1pgp_unbound_model_bound( AutoDock4_unbound_test ):
    """Test that autodock 4.1 works when unbound_model is set to bound."""
    dpf_stem = "1pgp_unbound_model_bound"
    expected_unbound_energy = -1.80
    expected_outcome = True # True means Successful Completion!
#______________________________________________________________________________

class AutoDock4_1pgp_unbound_set0_test( AutoDock4_unbound_test ):
    """Test that autodock 4.1 works when unbound is set to 0 in the DPF."""
    dpf_stem = "1pgp_unbound_set0"
    expected_unbound_energy = 0.
    expected_outcome = True # True means Successful Completion!
#______________________________________________________________________________
class AutoDock4_1pgp_unbound_set10_test( AutoDock4_unbound_test ):
    """Test that autodock 4.1 works when unbound is set to 10 in the DPF."""
    dpf_stem = "1pgp_unbound_set10"
    expected_unbound_energy = 10.
    expected_outcome = True # True means Successful Completion!
#______________________________________________________________________________

if __name__ == '__main__':
    #  This syntax lets us run all the tests,
    #  or conveniently comment out tests we're not interested in.
    #  NOTE:  Remember to add new TestCase class names to the list "test_cases"
    test_cases = [
        # tests for .dpf extension:
        'AutoDock4_1pgp_no_extension',
        'AutoDock4_1pgp_wrong_extension',
        'AutoDock4_1pgp_two_extensions',
        # simple tests:
        'AutoDock4_1pgp_ligand_types_map_mismatch',
        'AutoDock4_1pgp_illegal_keyword_test',
        'AutoDock4_1pgp_no_elecmap_test',
        'AutoDock4_1pgp_no_desolvmap_test',
        'AutoDock4_1pgp_no_elec_desolv_maps_test',
        'AutoDock4_1pgp_too_many_ligand_types_test',
        'AutoDock4_1pgp_too_many_torsions',
        'AutoDock4_1pgp_just_right_number_torsions',
        'AutoDock4_1pgp_two_ligands_test',
        'AutoDock4_1pgp_two_mapsets_test',
        'AutoDock4_1pgp_unbound_set_illegal_test',
        'AutoDock4_1pgp_unbound_model_illegal_test', #1
        'AutoDock4_1pgp_ga_select_tournament_test',
        ## tests which check for specific value
        'AutoDock4_1pgp_test',
        'AutoDock4_1pgp_smaller_test',
        'AutoDock4_1pgp_no_parameter_file_test',
        ## tests for unbound values 
        'AutoDock4_1pgp_unbound_default_test',
        'AutoDock4_1pgp_unbound_set0_test',
        'AutoDock4_1pgp_unbound_set10_test',
        # tests for unbound_model choices
        'AutoDock4_1pgp_unbound_model_compact',
        'AutoDock4_1pgp_unbound_model_bound',
        'AutoDock4_1pgp_unbound_model_extended',
        'AutoDock4_1pgp_unbound_compute_unbound_extended',
        'AutoDock4_1pgp_unbound_model_value',
    ]
    unittest.main( argv=( [__name__ ,] + test_cases ) )
    #  The call "unittest.main()" automatically runs all the TestCase classes in
    #  alphabetical order; calling with argv=([]), lets us specify the order.
    #  NOTE: "unittest.main()" saves us having to remember to add new tests to the 
    #  list of test cases.
    #unittest.main()
    #  For verbose output, use this:
    #unittest.main( argv=( [__name__, '-v'] + test_cases ) )
