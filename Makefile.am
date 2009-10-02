## Process this file with automake to produce Makefile.in

# $Id: Makefile.am,v 1.29 2009/10/02 22:33:47 rhuey Exp $
# 
# AutoDock  
# 
# Copyright (C) 1989-2007,  Garrett M. Morris,   TSRI
# All Rights Reserved.
# 
# AutoDock is a Trade Mark of The Scripps Research Institute.
# 
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

#
# Makefile.am to build AutoDock 4
#
# NOTE: Must be run in the $(AUTODOCK_DEV) directory
#

# Autotools

EXTRA_DIST = AD4_PARM99.dat AD4_parameters.dat AD4.1_bound.dat AD4.1_extended.dat AD4.1_compact.dat paramdat2h.csh Tests

dist-hook:
	rm -rf `find $(distdir)/Tests -name CVS -o -name '*.dlg'` 

bin_PROGRAMS = autodock4 

# Define the AutoDock 4 source code files:

autodock4_SOURCES = \
    main.cc

BUILT_SOURCES = main.cc default_parameters.h

CLEANFILES = default_parameters.h

noinst_LIBRARIES = libad.a 

libad_a_SOURCES = \
    alea.cc \
    analysis.cc \
	atom_parameter_manager.cc \
    banner.cc \
    bestpdb.cc \
	calculateEnergies.cc \
    call_cpso.cc \
    call_glss.cc \
    call_gs.cc \
    call_ls.cc \
    changeState.cc \
    check_header_float.cc \
    check_header_int.cc \
    check_header_line.cc \
    cluster_analysis.cc \
    clmode.cc \
    cnv_state_to_coords.cc \
    com.cc \
    conformation_sampler.cc \
    dimLibrary.cc\
	distdepdiel.cc \
    stateLibrary.cc \
    readfield.cc \
    readmap.cc \
    readPDBQT.cc \
	read_parameter_library.cc \
    eval.cc \
    gencau.cc \
    getrms.cc \
    get_atom_type.cc \
    getInitialState.cc \
    getpdbcrds.cc \
    gs.cc \
    initautodock.cc \
    input_state.cc \
    investigate.cc \
    linpack.cc \
    ls.cc \
    mapping.cc \
    minmeanmax.cc \
    mkNewState.cc \
    mkTorTree.cc \
    mkRandomState.cc \
    nonbonds.cc \
    openfile.cc \
    output_state.cc \
    parse_dpf_line.cc \
    parse_param_line.cc \
    parse_PDBQT_line.cc \
    parse_trj_line.cc \
	parsetypes.cc \
    print_2x.cc \
    print_atomic_energies.cc \
    print_avsfld.cc \
    writePDBQT.cc \
    print_rem.cc \
    printdate.cc \
    printEnergies.cc \
    printhms.cc \
    prClusterHist.cc \
    prInitialState.cc \
    prTorConList.cc \
    qmultiply.cc \
    qtransform.cc \
    quicksort.cc \
    ranlib.cc \
    rep.cc \
    scauchy.cc \
    setflags.cc \
    simanneal.cc \
    sort_enrg.cc \
    stack.cc \
    stop.cc \
    strindex.cc \
    success.cc \
    support.cc \
    swap.cc \
    timesys.cc \
    timesyshms.cc \
    torNorVec.cc \
    torsion.cc \
    trilinterp.cc \
    usage.cc \
    weedbonds.cc \
    warn_bad_file.cc \
    coliny.cc \
    eintcal.cc \
    eintcalPrint.cc \
    intnbtable.cc \
    nbe.cc \
		PDBQT_tokens.h \
		analysis.h \
		atom_parameter_manager.h \
		alea.h \
		autocomm.h \
		autoglobal.h \
		banner.h \
		bestpdb.h \
		calculateEnergies.h \
		call_cpso.h \
		call_glss.h \
		call_gs.h \
		call_ls.h \
		changeState.h \
		check_header_float.h \
		check_header_int.h \
		check_header_line.h \
		clmode.h \
		cluster_analysis.h \
		cnv_state_to_coords.h \
		coliny.h \
		conformation_sampler.h \
		constants.h \
		default_parameters.h \
		dimLibrary.h \
		distdepdiel.h \
		dpftoken.h \
		eintcal.h \
		eintcalPrint.h \
		energy.h \
		eval.h \
		gencau.h \
		getInitialState.h \
		get_atom_type.h \
		getpdbcrds.h \
		getrms.h \
		globals.h \
		grid.h \
		gs.h \
		hybrids.h \
		initautodock.h \
		input_state.h \
		intnbtable.h \
		investigate.h \
		ls.h \
		main.h \
		mdist.h \
		mkNewState.h \
		mkRandomState.h \
		mkTorTree.h \
		molstruct.h \
		nbe.h \
		nonbonds.h \
		openfile.h \
		output_state.h \
		parameters.h \
		parse_PDBQT_line.h \
		parse_dpf_line.h \
		parse_param_line.h \
		parse_trj_line.h \
		parsetypes.h \
		partokens.h \
		prClusterHist.h \
		prInitialState.h \
		prTorConList.h \
		printEnergies.h \
		print_2x.h \
		print_atomic_energies.h \
		print_avsfld.h \
		print_rem.h \
		printdate.h \
		printhms.h \
		qmultiply.h \
		qtransform.h \
		quicksort.h \
		ranlib.h \
		readPDBQT.h \
		read_parameter_library.h \
		readfield.h \
		readmap.h \
		rep.h \
		rep_constants.h \
		setflags.h \
		simanneal.h \
		sort_enrg.h \
		stateLibrary.h \
		stack.h \
		stop.h \
		strindex.h \
		structs.h \
		success.h \
		support.h \
		swap.h \
		timesys.h \
		timesyshms.h \
		torNorVec.h \
		torsion.h \
		trilinterp.h \
		trjtokens.h \
		typedefs.h \
		usage.h \
		version.h \
		warn_bad_file.h \
		weedbonds.h \
		writePDBQT.h

AM_CXXFLAGS = -O3 # Agressive optimization

    # if gnu c++ compiler
AM_CXXFLAGS += -Wall # All warnings
AM_CXXFLAGS += -Wno-write-strings # increase compatibility with pre 4.x compilers for string handling

#AM_CXXFLAGS = -Wall -O3 # All warnings, Agressive optimization
#AM_CXXFLAGS = -Wall -O3 -ftree-vectorize # All warnings, Agressive optimization + GCC vectorization # MacOSX only?

AM_CPPFLAGS = -DNOSQRT \
              -DUSE_8A_NBCUTOFF \
              -DWRITEPDBQSTATE \
              -I$(srcdir)/../acro/include # Standard accuracy, but faster; no crossover in quaternion


#autodock4_CXXFLAGS = -Wall -O3 # All warnings, Agressive optimization
#autodock4_CXXFLAGS = -Wall -O3 -ftree-vectorize # All warnings, Agressive optimization + GCC vectorization # MacOSX only?

#autodock4_CPPFLAGS = -DUSE_8A_NBCUTOFF \
#                     -DWRITEPDBQSTATE # Standard accuracy, but faster; no crossover in quaternion


autodock4_LDADD = libad.a $(COLINY_LDFLAGS) -lm

if WANT_BIGSTACK
    autodock4_LDFLAGS = -Wl,--stack,8388608 # Cygwin, big stacksize
endif

#
# If you need to use debugging or profiling, these should also be 
# modified appropriately:
# 	DBUG & PROF
#
# If you want to use the Coliny solver library, uncomment the following:
# COLINY = yes
#

# Abbreviations:
#
# Alpha   = Compaq/Digital Equipment Corp., Alpha
# Convex  = Convex, c2
# Cygwin  = Cygwin running on Microsoft Windows
# Darwin  = Darwin
# HP      = Hewlett Packard Precision Architecture, hppa
# Linux   = Any platform that runs GNU/Linux, Linux
# MacOSX  = Apple Mac OS X 10.0 & higher, MacOSX
# SGI     = Silicon Graphics Inc., sgi4D
# Sun     = Sun Microsystems, sun4

# CFLAGS = $(CSTD) $(OPT) -DUSE_8A_NBCUTOFF # SGI, HP, Alpha, Sun, Convex, Cygwin, Linux, MacOSX: Standard accuracy, but faster; no crossover in quaternion
# CFLAGS = $(CSTD) $(OPT) -DUSE_8A_NBCUTOFF # SGI, HP, Alpha, Sun, Convex, Linux, MacOSX: Standard accuracy, but faster
# CFLAGS = $(CSTD) $(OPT) -DUSE_8A_NBCUTOFF -DQUATERNION_MUTATION # As above, but treat quaternion genes properly when doing mutations--note, this is slow!
# CFLAGS = $(CSTD) $(OPT) -DUSE_8A_NBCUTOFF -DUSE_DOUBLE # SGI, HP, Alpha, Sun, Convex, Linux, MacOSX: Standard accuracy, but faster; also use Double precision throughout
# CFLAGS = $(CSTD) $(OPT) # SGI, HP, Alpha, Sun, Convex, Cygwin, Linux, MacOSX

# CC = g++ # HP, Gnu, Cygwin, Linux, Darwin, MacOSX
# CC = g++-4.0 -arch i386 # MacOSX -- Cross-compiling for Intel on PowerPC
# CC = g++-3.3 -arch ppc # MacOSX -- Cross-compiling for PowerPC on Intel
# CC = CC # SGI
# CC = cxx # Alpha

# CSTD = $(DBUG) $(PROF) $(WARN) # SGI, Sun, Linux, MacOSX
# CSTD = $(DBUG) $(PROF) $(WARN) -DUSE_XCODE # Smaller memory footprint, good for Xcode
# CSTD = $(DBUG) $(PROF) $(WARN) # SGI, Sun, Linux, MacOSX
# CSTD = $(DBUG) $(PROF) $(WARN) -I/opt/sfw/include # Sun Soliaris 8
# CSTD = $(DBUG) $(PROF) $(WARN) -std # Convex
# CSTD = -std -verbose $(PROF) $(DBUG) $(WARN) # Alpha. Not sarah
# CSTD = -std arm -verbose $(PROF) $(DBUG) $(WARN) # Alpha. sarah
# CSTD = -DHPPA -D_HPUX_SOURCE -ansi $(PROF) $(DBUG) $(WARN) # HP

# OPTLEVEL = -O3 -ftree-vectorize # Agressive optimization + GCC vectorization
# OPTLEVEL = -O3 -ffast-math # Agressive optimization, for Intel Itanium, ia64Linux2
# OPTLEVEL = -O3 # Agressive optimization
# OPTLEVEL = -fast # Agressive optimization for the G5 on MacOSX
# OPTLEVEL = -O2 # High optimization
# OPTLEVEL = -O1 # Do optimizations that can be done quickly; default.  Recommended for unit testing
# OPTLEVEL = -O0 # Do not optimize

# OPT_SGI_IPNUM = # Alpha, HP, Sun, Convex, SGI, Cygwin, Linux, MacOSX
# OPT_SGI_IPNUM = -Ofast=ip19 # SGI, 'uname -a' says 'IP19'
# OPT_SGI_IPNUM = -Ofast=ip21 # SGI, 'uname -a' says 'IP21'
# OPT_SGI_IPNUM = -Ofast=ip25 # SGI, 'uname -a' says 'IP25' PowerChallenge is R10000, IP25
# OPT_SGI_IPNUM = -Ofast=ip27 # SGI, 'uname -a' says 'IP27'
# OPT_SGI_IPNUM = -Ofast=ip30 # SGI, 'uname -a' says 'IP30'
# OPT_SGI_IPNUM = `uname -m | sed 's/IP/-Ofast=ip/'` # SGI, dynamic

# OPT_SGI_R000 = # Alpha, HP, Sun, Convex, SGI, Cygwin, Linux, MacOSX
# OPT_SGI_R000 = -r4000 -mips2 # SGI, 'hinv' says MIPS Processor is R4000
# OPT_SGI_R000 = -r8000 -mips4 # SGI, 'hinv' says MIPS Processor is R8000
# OPT_SGI_R000 = -r10000 -mips4 # SGI, 'hinv' says MIPS Processor is R10000
# OPT_SGI_R000 = -r12000 -mips4 # SGI, 'hinv' says MIPS Processor is R12000
# OPT_SGI_R000 = -r14000 -mips4 # SGI, 'hinv' says MIPS Processor is R14000
# OPT_SGI_R000 = `hinv | grep '^CPU:' | awk '{print $3}' | sed 's/R/-r/'` -mips4 # SGI, dynamic, -mips4 (works with -r8000 to -r14000, not -r4000)

# OPT = $(OPTLEVEL) # Alpha, HP, Sun, Convex, Cygwin, Linux, MacOSX
# OPT = $(OPTLEVEL) -ffast-math # Gnu cc, fast-math is dangerous!
# OPT = $(OPTLEVEL) -n32 $(OPT_SGI_IPNUM) $(OPT_SGI_R000) -IPA $(LNO_OPT) # SGI
# OPT = $(OPTLEVEL) -n32 $(OPT_SGI_IPNUM) $(OPT_SGI_R000) -IPA $(LNO_OPT) -DUSE_INT_AS_LONG # SGI (long is 8bytes)
# OPT = $(OPTLEVEL) $(OPT_SGI_IPNUM) $(OPT_SGI_R000) $(LNO_OPT) # SGI, not new 32-bit

# LNO_OPT = # SGI, no special optimization at link time; Sun, Cygwin, Linux, MacOSX
# LNO_OPT = -LNO:auto_dist=ON:gather_scatter=2 # SGI

# LINKOPT = $(CSTD) $(OPT) # All platforms except Cygwin and Sun
# LINKOPT = $(CSTD) $(OPT) -fno-stack-limit # Cygwin, 32MB stacksize
# LINKOPT = $(CSTD) $(OPT) -Wl,--stack=0x4000000 # Cygwin, 64MB stacksize
# LINKOPT = $(CSTD) $(OPT) -L/opt/sfw/lib # Sun

# LINK = $(LINKOPT) # Linking flags
# LINK = $(LINKOPT) -cord # Procedure rearranger on SGI

# LIB = -lm # for all platforms
# LIB = -lSaturn # for profiling using MacOSX Saturn

# DBUG = # Use assert code
# DBUG = -DNDEBUG # No debugging and no assert code
# DBUG = -g # dbx, or Gnu gdb
# DBUG = -g -DDEBUG # dbx + DEBUG-specific code
# DBUG = -g3 # dbx + optimization
# DBUG = -g3 -DDEBUG # dbx + optimization, + DEBUG-specific code
# DBUG = -DDEBUG # Just DEBUG-specific code
# DBUG = -DDEBUG2 # Just DEBUG2-specific code for tracking prop.selection
# DBUG = -DDEBUG3 # Just DEBUG3-specific code for print age of individuals
# DBUG = -g -DDEBUG -DDEBUG2 -DDEBUG3 # Debug everything
# DBUG = -DDEBUG_MUTATION # Use assert code, & print out information anywhere involving mutation
# DBUG = -g3 -DDEBUG_QUAT # Use assert code, & assert quaternions are 4D-normalised
# DBUG = -g3 -DDEBUG_QUAT -DDEBUG_QUAT_PRINT # Use assert code, & assert quaternions are 4D-normalised, and print out
# DBUG = -g3 -DDEBUG_QUAT -DDEBUG_QUAT_PRINT -DDEBUG # Use assert code, & assert quaternions are 4D-normalised, and use DEBUG-specific code, and print out
# DBUG = -g3 -DDEBUG_MUTATION # dbx + optimization; Use assert code, & print out information anywhere involving mutation
# DBUG = -g3 -DDEBUG_MUTATION -DDEBUG # dbx + optimization; Use assert code, & print out information anywhere involving mutation + all debug statements
# DBUG = -g3 -DDEBUG_MUTATION -DDEBUG -DDEBUG3 # dbx + optimization; Use assert code, & print out information anywhere involving mutation + all debug statements incl. 3
# DBUG = -DASSERTQUATOK # Uses the macro "assertQuatOK()" to ensure the quaternion magnitude is near unity

# PROF = # No profiling
# PROF = -p # CC Profiling
# PROF = -pg # For GNU gcc & GNU g++ Profiling with GNU gprof

# WARN = # Default warning level
# WARN = -woff all # For no warnings
# WARN = -Wall # All warnings, gcc -- Recommended for developers
# WARN = -fullwarn -ansiE -ansiW # For full warnings during compilation

### ifeq ($(COLINY),yes)
###   ACRO_OS = -DDARWIN # Darwin, MacOSX
###   # ACRO_OS = -DLINUX # Linux
### 
###   # ACRO_LINK = -L../acro/lib -lcoliny -lcolin -lpico -lutilib -lappspack -l3po -lg2c # Linux
###   ACRO_LINK = -L../acro/lib -lcoliny -lcolin -lpico -lutilib -lappspack -l3po -L/sw/lib -lg2c # Darwin, MacOSX (Fink needed)
### 
###   ACRO_INCLUDES = -I../acro/include -DUSING_COLINY $(ACRO_FLAGS)
###   ACRO_FLAGS = -DDEBUGGING -DUNIX -DMULTITASK -DANSI_HDRS -DANSI_NAMESPACES $(ACRO_OS)
###   LIB = $(ACRO_LINK)
###   ACRO = acro
###   COLINYLIB = libcoliny.a # Using Coliny
### else
#COLINY_LDADD = # Not using Coliny
### endif

# Uncomment this section with "###" if cross-compiling on MacOSX (Intel and PowerPC)
# This section remains the same for either type of cross-compilation.
### HOST_CC="gcc-4.0"
### HOST_CXX="g++-4.0"
### RANLIB=ranlib
### AR=ar
### AS=$CC
### LD=ld
### STRIP="strip -x -S"
### CROSS_COMPILE=1

#
# Makefile targets
#

autodock4.html :
	@echo " "
	@echo "Creating a colour-coded HTML report of the AutoDock source code"
	@echo " "
	enscript -E --color -Whtml --toc -p$@ *.h *.cc

check : test

test : autodock4 $(srcdir)/Tests/test_autodock4.py
#if HAVE_PYTHON
	@echo " "
	@echo "Running the Python test_autodock4.py unit-tests on autodock4"
	@echo " "
	outputdir=`pwd` ;\
	umask 0 ;\
	test -d $$outputdir/Tests || mkdir $$outputdir/Tests ;\
    cp $(srcdir)/AD4.1_bound.dat  $$outputdir/Tests;\
    cp $(srcdir)/AD4_parameters.dat  $$outputdir/Tests;\
    cp $(srcdir)/Tests/*.py   $$outputdir/Tests;\
    cp $(srcdir)/Tests/*.dpf  $$outputdir/Tests;\
    cp $(srcdir)/Tests/*.pdbqt  $$outputdir/Tests;\
    cp $(srcdir)/Tests/*.map*  $$outputdir/Tests;\
	cd Tests ;  python ./test_autodock4.py -d . -e ../autodock4$(EXEEXT) -o .
#else
#	@echo " "
#	@echo "Python version 2.1 or higher needed for 'make check', but not found."
#	@echo "Skipping AutoDock4 unit tests."
#	@echo " "
#endif

#$(COLINYLIB) : $(OBJS) $(OBJNOSQRT) $(OBJNOMINPT)
#	@echo " "
#	@echo Making the Coliny library
#	@echo " "
#	$(AR) $(ARFLAGS) $(COLINYLIB) $(?:.cc=.o)
#	$(RANLIB) $(RANLIBFLAGS) $(COLINYLIB)
###
### competing alternative from another CVS user
###
# coliny: $(ACRO)

### acro:
### 	@echo " "
### 	@echo Making the Coliny library
### 	@echo " "
### 	@(cd ../acro;\
### 	if [ ! -e Makefile ]; then\
### 	   ./setup configure;\
### 	fi;\
### 	$(MAKE))


#
# Object dependencies:
#

gs.o : gs.cc gs.h ranlib.h eval.h rep.h support.h writePDBQT.h
	$(CXXCOMPILE) -DCHECK_ISNAN -c $(srcdir)/gs.cc

### main.o : main.cc hybrids.h ranlib.h gs.h ls.h rep.h support.h main.h constants.h autocomm.h dpftoken.h structs.h autoglobal.h  autocomm.h coliny.h parse_param_line.cc partokens.h eintcalPrint.cc eintcalPrint.h eintcal.cc atom_parameter_manager.cc atom_parameter_manager.h read_parameter_library.h calculateEnergies.cc calculateEnergies.h
### $(CXXCOMPILE) $(ACRO_INCLUDES) -c main.cc

# 2007-04-30 MP & GMM: this line is necessary to be able to compile read_parameter_library.o
read_parameter_library.o : default_parameters.h

default_parameters.h : AD4_parameters.dat AD4.1_bound.dat paramdat2h.csh
	rm -f $@
	csh $(srcdir)/paramdat2h.csh $(srcdir)/AD4_parameters.dat $(srcdir)/AD4.1_bound.dat > $@

#coliny.o : coliny.h
#	$(CXX) $(CPPFLAGS) $(CXXFLAGS) $(ACRO_INCLUDES) -c coliny.cc -o coliny.o

# see Autoconf manual 4.7.3 (p. 25) support for VPATH
# enabling multi-platform builds
srcdir = @srcdir@
VPATH  = @srcdir@

#
# EOF
#
