#
# Makefile to build AutoDock from Object files.
#
# NOTE: Must be run in the $(AUTODOCK_DEV) directory.
#
# Copyright (C) 1994-2001,  Garrett Matthew Morris,  TSRI.
#
#
# Edit this Makefile to reflect your machine architecture.
#
# Specifically, change these variables:
# LIB, CSTD, CFLAGS, OPT, OLIMIT, LINT, LINTFLAGS, LINK, WARN, & CC.
#
# If you need to use debugging or profiling, these should also be 
# modified appropriately:
# DBUG & PROF
#

#
# Edit the line defining "EXE" to your binaries directory for 
# your machine architecture.
#
# Define the destination directory for the executables:
#

EXE = . # Default destination

#
# Define the object files:
#

OBJS = \
    analysis.o \
    banner.o \
    bestpdb.o \
    call_glss.o \
    call_gs.o \
    call_ls.o \
    changeState.o \
    check_header_float.o \
    check_header_int.o \
    check_header_line.o \
    cluster_analysis.o \
    clmode.o \
    cmdmode.o \
    cnv_state_to_coords.o \
    com.o \
    stateLibrary.o \
    readfield.o \
    readmap.o \
    readPDBQ.o \
    dpftypes.o \
    eval.o \
    evaluate_energy.o \
    gencau.o \
    getrms.o \
    get_atom_type.o \
    getInitialState.o \
    getpdbcrds.o \
    gs.o \
    initautodock.o \
    input_state.o \
    investigate.o \
    linpack.o \
    ls.o \
    mapping.o \
    minmeanmax.o \
    mkNewState.o \
    mkTorTree.o \
    mkRandomState.o \
    nonbonds.o \
    openfile.o \
    output_state.o \
    parse_com_line.o \
    parse_dpf_line.o \
    parse_pdbq_line.o \
    parse_trj_line.o \
    print_2x.o \
    print_atomic_energies.o \
    print_avsfld.o \
    writeMolAsPDBQ.o \
    writePDBQ.o \
    writePDBQState.o \
    print_rem.o \
    printdate.o \
    printEnergies.o \
    printhms.o \
    prClusterHist.o \
    prInitialState.o \
    prTorConList.o \
    qmultiply.o \
    qtransform.o \
    quicksort.o \
    ranlib.o \
    rep.o \
    scauchy.o \
    set_cmd_io_std.o \
    setflags.o \
    simanneal.o \
    sort_enrg.o \
    stop.o \
    strindex.o \
    success.o \
    summarizegrids.o \
    support.o \
    swap.o \
    timesys.o \
    timesyshms.o \
    torNorVec.o \
    torsion.o \
    usage.o \
    weedbonds.o \
    warn_bad_file.o

OBJNOSQRT = \
    eintcal.o \
    eintcalPrint.o \
    intnbtable.o \
    nbe.o

OBJSQRT = \
    eintcal.sqrt.o \
    eintcalPrint.sqrt.o \
    intnbtable.sqrt.o \
    nbe.sqrt.o


OBJNOMINPT = \
    trilinterp.o

OBJMINPT = \
    trilinterp.MINPT.o


ADLIB = libad.a

ARFLAGS = r # SGI, Sun, Alpha, Linux, Mac OS X


# RANLIB = file # SGI.
RANLIB = ranlib # Linux, Mac OS X.

RANLIBFLAGS = # Linux, SGI
# RANLIBFLAGS = -s # MacOS X.

# Define lint files:

LNS = \
    analysis.ln \
    banner.ln \
    bestpdb.ln \
    changeState.ln \
    check_header_float.ln \
    check_header_int.ln \
    check_header_line.ln \
    cluster_analysis.ln \
    clmode.ln \
    cmdmode.ln \
    cnv_state_to_coords.ln \
    stateLibrary.ln \
    readfield.ln \
    readmap.ln \
    readPDBQ.ln \
    dpftypes.ln \
    evaluate_energy.ln \
    getrms.ln \
    get_atom_type.ln \
    getInitialState.ln \
    getpdbcrds.ln \
    initautodock.ln \
    input_state.ln \
    investigate.ln \
    main.ln \
    mkNewState.ln \
    mkTorTree.ln \
    mkRandomState.o \
    nonbonds.ln \
    openfile.ln \
    output_state.ln \
    parse_com_line.ln \
    parse_dpf_line.ln \
    parse_pdbq_line.ln \
    parse_trj_line.ln \
    print_2x.ln \
    print_atomic_energies.ln \
    print_avsfld.ln \
    writeMolAsPDBQ.ln \
    writePDBQ.ln \
    writePDBQState.ln \
    print_rem.ln \
    printdate.ln \
    printEnergies.ln \
    printhms.ln \
    prClusterHist.ln \
    prInitialState.ln \
    prTorConList.ln \
    qmultiply.ln \
    qtransform.ln \
    quicksort.ln \
    set_cmd_io_std.ln \
    setflags.ln \
    simanneal.ln \
    sort_enrg.ln \
    stop.ln \
    strindex.ln \
    success.ln \
    summarizegrids.ln \
    swap.ln \
    timesys.ln \
    timesyshms.ln \
    torNorVec.ln \
    torsion.ln \
    trilinterp.ln \
    usage.ln \
    weedbonds.ln \
    warn_bad_file.ln

LNSNOSQRT = \
    eintcal.ln \
    eintcalPrint.ln \
    intnbtable.ln \
    nbe.ln

LNSSQRT = \
    eintcal.sqrt.ln \
    eintcalPrint.sqrt.ln \
    intnbtable.sqrt.ln \
    nbe.sqrt.ln

#
# Abbreviations of machine architectures:
#
# SGI     = Silicon Graphics Inc., sgi4D.
# Alpha   = Compaq/Digital Equipment Corp., Alpha.
# Sun     = Sun Microsystems, sun4.
# HP      = Hewlett Packard Precision Architecture, hppa.
# Convex  = Convex, c2.
# Linux   = Any platform that runs Linux, Linux
# MacOS X = Apple Mac OS X 10.0 & higher, MacOS X
#

#
# C++ compiler
#

CC = cc # MacOS X
# CC = CC # SGI.
# CC = cxx # Alpha.
# CC = gcc # use this if you have the Gnu compiler, as on Linux, MkLinux, LinuxPPC systems.
# #  portland compiler setup #
# PGI = /usr/pgi #                       Portland compiler
# LM_LICENSE_FILE = $(PGI)/license.dat # Portland compiler
# LD_LIBRARY_PATH = $(PGI)/linux86/lib # Portland compiler
# CC = pgCC #                            Portland Compiler

LIB = -lm -lsupc++  # gcc 3.1 on MacOS X.
# LIB = -lm # SGI, Sun, Linux, MacOS X.
# LIB = -lm -lc # Alpha, Convex.
# LIB = -lm -lg++ # HP, Gnu.

CSTD = -DUSE_DOUBLE $(DBUG) $(PROF) $(WARN) # SGI, Sun, Linux, MacOS X.
# CSTD = $(DBUG) $(PROF) $(WARN) # SGI, Sun, Linux, MacOS X.
# CSTD = $(DBUG) $(PROF) $(WARN) -I/opt/sfw/include # Sun Soliaris 8
# CSTD = $(DBUG) $(PROF) $(WARN) -std # Convex.
# CSTD = -std -verbose $(PROF) $(DBUG) $(WARN) # Alpha. Not sarah
# CSTD = -std arm -verbose $(PROF) $(DBUG) $(WARN) # Alpha. sarah
# CSTD = -DHPPA -D_HPUX_SOURCE -ansi $(PROF) $(DBUG) $(WARN) # HP.

CFLAGS = $(CSTD) $(OPT) # SGI, HP, Alpha, Sun, Convex, Linux, MacOS X: Optimize the object files, too.

OLIMIT = $(CSTD) $(OPT) # SGI, Sun, HP, Convex, Linux, MacOS X.
# OLIMIT = $(CSTD) $(OPT) -OPT:Olimit=2500 # Alpha, Some SGIs.
# OLIMIT = $(CFLAGS) # Do not optimize.

# OPTLEVEL = -O3 # Agressive optimization.
# OPTLEVEL = -O2 # High optimization.
OPTLEVEL = -O1 # Do optimizations that can be done quickly; default.  Recommended for unit testing.
# OPTLEVEL = -O0 # Do not optimize.

OPT_SGI_IPNUM = # Alpha, HP, Sun, Convex, SGI, Linux, MacOS X.
# OPT_SGI_IPNUM = -Ofast=ip19 # SGI, 'uname -a' says 'IP19'
# OPT_SGI_IPNUM = -Ofast=ip21 # SGI, 'uname -a' says 'IP21'
# OPT_SGI_IPNUM = -Ofast=ip25 # SGI, 'uname -a' says 'IP25' PowerChallenge is R10000, IP25
# OPT_SGI_IPNUM = -Ofast=ip27 # SGI, 'uname -a' says 'IP27'
# OPT_SGI_IPNUM = -Ofast=ip30 # SGI, 'uname -a' says 'IP30'
# OPT_SGI_IPNUM = `uname -m | sed 's/IP/-Ofast=ip/'` # SGI, dynamic
# TSRI job = IP30
# TSRI ben = IP30
# TSRI atlas, thing1, thing2 = IP27

OPT_SGI_R000 = # Alpha, HP, Sun, Convex, SGI, Linux, MacOS X.
# OPT_SGI_R000 = -r4000 -mips2 # SGI, 'hinv' says MIPS Processor is R4000
# OPT_SGI_R000 = -r8000 -mips4 # SGI, 'hinv' says MIPS Processor is R8000
# OPT_SGI_R000 = -r10000 -mips4 # SGI, 'hinv' says MIPS Processor is R10000
# OPT_SGI_R000 = -r12000 -mips4 # SGI, 'hinv' says MIPS Processor is R12000
# OPT_SGI_R000 = -r14000 -mips4 # SGI, 'hinv' says MIPS Processor is R14000
# OPT_SGI_R000 = `hinv | grep '^CPU:' | awk '{print $3}' | sed 's/R/-r/'` -mips4 # SGI, dynamic, -mips4 (works with -r8000 to -r14000, not -r4000)
# TSRI job = R10000
# TSRI ben = R12000
# TSRI atlas = R12000
# TSRI thing1, thing2 = R14000

OPT = $(OPTLEVEL) # Alpha, HP, Sun, Convex, Linux, MacOS X.
# OPT = $(OPTLEVEL) -ffast-math # Gnu cc, fast-math is dangerous!
# OPT = $(OPTLEVEL) -n32 $(OPT_SGI_IPNUM) $(OPT_SGI_R000) -IPA $(LNO_OPT) # SGI
# OPT = $(OPTLEVEL) -n32 $(OPT_SGI_IPNUM) $(OPT_SGI_R000) -IPA $(LNO_OPT) -DUSE_INT_AS_LONG # SGI (long is 8bytes).
# OPT = $(OPTLEVEL) $(OPT_SGI_IPNUM) $(OPT_SGI_R000) $(LNO_OPT) # SGI, not new 32-bit

# LNO_OPT = # SGI, no special optimization at link time; Sun; Linux; MacOS X
LNO_OPT = -LNO:auto_dist=ON:gather_scatter=2 # SGI

# LINKOPT = $(CSTD) $(OPT) # 
LINKOPT = $(CSTD) $(OPT) -fno-stack-limit # Cygwin, 32MB stacksize
# LINKOPT = $(CSTD) $(OPT) -Wl,--stack=0x2000000 # Cygwin, 32MB stacksize
# LINKOPT = $(CSTD) $(OPT) -L/opt/sfw/lib # Sun

LINK = $(LINKOPT) # Linking flags.
# LINK = $(LINKOPT) -cord # Procedure rearranger on SGI.

LINT = lint # lint C code checking.

# LINTFLAGS = $(LIB) -c # SGI, Linux, MacOS X.
# LINTFLAGS = $(LIB) -MA -c # Alpha.
# LINTFLAGS = -DHPPA -D_HPUX_SOURCE $(LIB) -c # HP.
LINTFLAGS = -u -n -lm # Sun.

DBUG = -DNDEBUG # No debugging and no assert code.
# DBUG = # Use assert code.
# DBUG = -g # dbx, or Gnu gdb.
# DBUG = -g -DDEBUG # dbx + DEBUG-specific code.
# DBUG = -g3 # dbx + optimization.
# DBUG = -g3 -DDEBUG # dbx + optimization, + DEBUG-specific code.
# DBUG = -DDEBUG # Just DEBUG-specific code.
# DBUG = -DDEBUG2 # Just DEBUG2-specific code for tracking prop.selection.
# DBUG = -DDEBUG3 # Just DEBUG3-specific code for print age of individuals.
# DBUG = -g -DDEBUG -DDEBUG2 -DDEBUG3 # Debug everything

PROF = # No profiling.
# PROF = -p # Profiling.

WARN = # Default warning level.
# WARN = -woff all # For no warnings.
# WARN = -fullwarn -ansiE -ansiW # For full warnings during compilation.


autodock3 : main.o $(ADLIB)
	echo $(EXE)'  on  '`date`', by $(USER) using '`hostname` >> LATEST_MAKE
	echo 'Flags: '$(CC) $(LINK) -DNOSQRT -L. -lad $(LIB) >> LATEST_MAKE
	@echo " "
	@echo Making autodock3
	@echo " "
	$(CC) $(LINK) -DNOSQRT -o $@ main.o -L. -lad $(LIB)

autodock3sqrt : main.o $(ADLIB)
	$(CC) $(CFLAGS) -o $@ main.o -L. -lad $(LIB)

autodock3minpt : main.o $(ADLIB)
	echo $(EXE)'  on  '`date`', using '`hostname` >> LATEST_MAKE
	$(CC) $(CFLAGS) -DNOSQRT -o $@ main.o -L. -lad $(LIB)

autodock3alt : $(OBJS) $(OBJNOSQRT) $(OBJNOMINPT)
	echo $(EXE)'  on  '`date`', using '`hostname` >> LATEST_MAKE
	$(CC) $(LINK) -DNOSQRT -o $@ $(OBJS) $(OBJNOSQRT) $(OBJNOMINPT) $(LIB)

install :
	@echo " "
	@echo Moving autodock3 to $(AUTODOCK_BIN)
	@echo " "
	mv autodock3 $(AUTODOCK_BIN)

$(ADLIB) : $(OBJS) $(OBJNOSQRT) $(OBJNOMINPT)
	@echo " "
	@echo Making the AutoDock library
	@echo " "
	$(AR) $(ARFLAGS) $(ADLIB) $(?:.cc=.o)
	$(RANLIB) $(RANLIBFLAGS) $(ADLIB)


lcheck : $(LNS) $(LNSNOSQRT)
	$(LINT) $(LIB) $(LNS) $(LNSNOSQRT)

lchecksqrt : $(LNS) $(LNSSQRT)
	$(LINT) $(LIB) $(LNS) $(LNSSQRT)

dualmap : dualmap.c
	$(CC) $(CFLAGS) -lm dualmap.c -o $@

#
# Object dependencies:
#

analysis.o : analysis.cc analysis.h constants.h getpdbcrds.h stateLibrary.h cnv_state_to_coords.h sort_enrg.h cluster_analysis.h prClusterHist.h getrms.h eintcal.h trilinterp.h print_rem.h strindex.h print_avsfld.h
	$(CC) $(CFLAGS) -c analysis.cc

banner.o : banner.cc banner.h
	$(CC) $(CFLAGS) -c banner.cc

bestpdb.o : bestpdb.cc bestpdb.h constants.h print_rem.h strindex.h print_avsfld.h
	$(CC) $(CFLAGS) -c bestpdb.cc

call_glss.o : call_glss.cc support.h rep.h eval.h ranlib.h call_glss.h
	$(CC) $(CFLAGS) -c call_glss.cc

call_gs.o : call_gs.cc support.h rep.h eval.h ranlib.h call_gs.h autocomm.h timesyshms.h
	$(CC) $(CFLAGS) -c call_gs.cc

call_ls.o : call_ls.cc support.h rep.h eval.h ranlib.h call_ls.h
	$(CC) $(CFLAGS) -c call_ls.cc

changeState.o : changeState.cc changeState.h constants.h qmultiply.h
	$(CC) $(CFLAGS) -c changeState.cc

check_header_float.o : check_header_float.cc check_header_float.h
	$(CC) $(CFLAGS) -c check_header_float.cc

check_header_int.o : check_header_int.cc check_header_int.h constants.h print_2x.h stop.h
	$(CC) $(CFLAGS) -c check_header_int.cc

check_header_line.o : check_header_line.cc check_header_line.h constants.h
	$(CC) $(CFLAGS) -c check_header_line.cc

clmode.o : clmode.cc clmode.h constants.h openfile.h strindex.h readPDBQ.h get_atom_type.h getpdbcrds.h sort_enrg.h cluster_analysis.h prClusterHist.h bestpdb.h success.h
	$(CC) $(CFLAGS) -c clmode.cc

cluster_analysis.o : cluster_analysis.cc cluster_analysis.h constants.h getrms.h
	$(CC) $(CFLAGS) -c cluster_analysis.cc

cmdmode.o : cmdmode.cc cmdtokens.h trjtokens.h cmdmode.h constants.h set_cmd_io_std.h print_2x.h parse_com_line.h strindex.h print_avsfld.h success.h openfile.h readPDBQ.h get_atom_type.h timesys.h eintcal.h trilinterp.h qmultiply.h cnv_state_to_coords.h parse_trj_line.h input_state.h autocomm.h
	$(CC) $(CFLAGS) -c -DEINTCALPRINT cmdmode.cc

cnv_state_to_coords.o : cnv_state_to_coords.cc cnv_state_to_coords.h constants.h torsion.h qtransform.h stateLibrary.h
	$(CC) $(CFLAGS) -c cnv_state_to_coords.cc

com.o : com.cc ranlib.h
	$(CC) $(CFLAGS) -c com.cc

stateLibrary.o : stateLibrary.cc stateLibrary.h constants.h
	$(CC) $(CFLAGS) -c stateLibrary.cc

readPDBQ.o : readPDBQ.cc  readPDBQ.h constants.h openfile.h stop.h readPDBQ.h get_atom_type.h print_2x.h mkTorTree.h nonbonds.h weedbonds.h torNorVec.h success.h autocomm.h
	$(CC) $(OLIMIT) -c readPDBQ.cc

dpftypes.o : dpftypes.cc dpftypes.h constants.h dpftoken.h stop.h
	$(CC) $(CFLAGS) -c dpftypes.cc

eintcal.o : eintcal.cc eintcal.h constants.h
	$(CC) $(CFLAGS) -DNOSQRT -DBOUNDED -c eintcal.cc -o eintcal.o

eintcalPrint.o : eintcal.cc eintcal.h constants.h
	$(CC) $(CFLAGS) -DNOSQRT -DBOUNDED -DEINTCALPRINT -c eintcal.cc -o eintcalPrint.o

eval.o : eval.cc eval.h structs.h constants.h autocomm.h 
	$(CC) $(CFLAGS) -c eval.cc

evaluate_energy.o : evaluate_energy.cc evaluate_energy.h constants.h trilinterp.h eintcal.h
	$(CC) $(CFLAGS) -c evaluate_energy.cc

gencau.o : gencau.cc ranlib.h
	$(CC) $(CFLAGS) -c gencau.cc

getInitialState.o : getInitialState.cc getInitialState.h constants.h qmultiply.h stateLibrary.h initautodock.h trilinterp.h eintcal.h cnv_state_to_coords.h prInitialState.h timesys.h autocomm.h
	$(CC) $(CFLAGS) -c getInitialState.cc

get_atom_type.o : get_atom_type.cc get_atom_type.h constants.h
	$(CC) $(CFLAGS) -c get_atom_type.cc

getpdbcrds.o : getpdbcrds.cc getpdbcrds.h constants.h openfile.h
	$(CC) $(CFLAGS) -c getpdbcrds.cc

getrms.o : getrms.cc getrms.h constants.h autocomm.h
	$(CC) $(CFLAGS) -c getrms.cc

gs.o : gs.cc gs.h ranlib.h eval.h rep.h support.h writePDBQ.h
	$(CC) $(CFLAGS) -DCHECK_ISNAN -c gs.cc -o gs.o

initautodock.o : initautodock.cc initautodock.h constants.h qmultiply.h cnv_state_to_coords.h print_2x.h autocomm.h
	$(CC) $(CFLAGS) -c initautodock.cc

input_state.o : input_state.cc input_state.h constants.h qmultiply.h
	$(CC) $(CFLAGS) -c input_state.cc

investigate.o : investigate.cc investigate.h constants.h changeState.h mkRandomState.h cnv_state_to_coords.h getrms.h trilinterp.h eintcal.h getpdbcrds.h stateLibrary.h
	$(CC) $(CFLAGS) -c investigate.cc

intnbtable.o : intnbtable.cc intnbtable.h constants.h
	$(CC) $(OLIMIT) -DNOSQRT -c intnbtable.cc

linpack.o : linpack.cc
	$(CC) $(CFLAGS) -c linpack.cc

ls.o : ls.cc ls.h support.h ranlib.h
	$(CC) $(CFLAGS) -c ls.cc

main.o : main.cc hybrids.h ranlib.h gs.h ls.h rep.h support.h main.h constants.h autocomm.h dpftoken.h structs.h autoglobal.h  autocomm.h
	$(CC) $(OLIMIT) -c -DWRITEPDBQSTATE main.cc

mapping.o : mapping.cc support.h
	$(CC) $(CFLAGS) -c mapping.cc

minmeanmax.o : minmeanmax.cc  rep.h support.h
	$(CC) $(CFLAGS) -c minmeanmax.cc

mkNewState.o : mkNewState.cc mkNewState.h constants.h
	$(CC) $(CFLAGS) -c mkNewState.cc

mkTorTree.o : mkTorTree.cc pdbqtokens.h mkTorTree.h constants.h
	$(CC) $(CFLAGS) -c mkTorTree.cc

mkRandomState.o : mkRandomState.cc mkRandomState.h constants.h
	$(CC) $(CFLAGS) -c mkRandomState.cc

nbe.o : nbe.cc nbe.h constants.h
	$(CC) $(OLIMIT) -DNOSQRT -c nbe.cc

nonbonds.o : nonbonds.cc nonbonds.h constants.h
	$(CC) $(CFLAGS) -c nonbonds.cc

openfile.o : openfile.cc openfile.h constants.h autocomm.h
	$(CC) $(CFLAGS) -c openfile.cc

output_state.o : output_state.cc output_state.h constants.h
	$(CC) $(CFLAGS) -c output_state.cc

parse_com_line.o : parse_com_line.cc cmdtokens.h parse_com_line.h constants.h
	$(CC) $(CFLAGS) -c parse_com_line.cc

parse_dpf_line.o : parse_dpf_line.cc parse_dpf_line.h constants.h dpftoken.h
	$(CC) $(CFLAGS) -c parse_dpf_line.cc

parse_pdbq_line.o : parse_pdbq_line.cc pdbqtokens.h  parse_pdbq_line.h constants.h
	$(CC) $(CFLAGS) -c parse_pdbq_line.cc

parse_trj_line.o : parse_trj_line.cc trjtokens.h parse_trj_line.h constants.h
	$(CC) $(CFLAGS) -c parse_trj_line.cc

prClusterHist.o : prClusterHist.cc prClusterHist.h constants.h
	$(CC) $(CFLAGS) -c prClusterHist.cc

prInitialState.o : prInitialState.cc prInitialState.h constants.h
	$(CC) $(CFLAGS) -c prInitialState.cc

prTorConList.o : prTorConList.cc prTorConList.h constants.h autocomm.h
	$(CC) $(CFLAGS) -c prTorConList.cc

print_2x.o : print_2x.cc print_2x.h
	$(CC) $(CFLAGS) -c print_2x.cc

print_atomic_energies.o : print_atomic_energies.cc print_atomic_energies.h constants.h
	$(CC) $(CFLAGS) -c print_atomic_energies.cc

print_avsfld.o : print_avsfld.cc print_avsfld.h
	$(CC) $(CFLAGS) -c print_avsfld.cc 

writeMolAsPDBQ.o : writePDBQ.cc writePDBQ.h constants.h autocomm.h
	$(CC) $(CFLAGS) -c -DWRITEMOLASPDBQFUNC writePDBQ.cc -o writeMolAsPDBQ.o

writePDBQ.o : writePDBQ.cc writePDBQ.h constants.h autocomm.h
	$(CC) $(CFLAGS) -c writePDBQ.cc -o writePDBQ.o

writePDBQState.o : writePDBQ.cc writePDBQ.h constants.h autocomm.h
	$(CC) $(CFLAGS) -c -DWRITEPDBQSTATE writePDBQ.cc -o writePDBQState.o

print_rem.o : print_rem.cc print_rem.h
	$(CC) $(CFLAGS) -c print_rem.cc

printdate.o : printdate.cc printdate.h
	$(CC) $(CFLAGS) -c printdate.cc

printEnergies.o : printEnergies.cc printEnergies.h
	$(CC) $(CFLAGS) -c printEnergies.cc

printhms.o : printhms.cc printhms.h
	$(CC) $(CFLAGS) -c printhms.cc

qmultiply.o : qmultiply.cc qmultiply.h
	$(CC) $(CFLAGS) -c qmultiply.cc

qtransform.o : qtransform.cc qtransform.h constants.h
	$(CC) $(CFLAGS) -c qtransform.cc

quicksort.o : quicksort.cc quicksort.h
	$(CC) $(CFLAGS) -c quicksort.cc

ranlib.o : ranlib.cc ranlib.h
	$(CC) $(CFLAGS) -c ranlib.cc

readfield.o : readfield.cc readfield.h constants.h openfile.h stop.h
	$(CC) $(CFLAGS) -c readfield.cc

readmap.o : readmap.cc readmap.h constants.h openfile.h warn_bad_file.h strindex.h print_2x.h check_header_line.h warn_bad_file.h check_header_float.h check_header_int.h timesys.h autocomm.h
	$(CC) $(CFLAGS) -c readmap.cc

rep.o: rep.cc rep.h ranlib.h
	$(CC) $(CFLAGS) -c rep.cc

scauchy.o  : scauchy.cc ranlib.h
	$(CC) $(CFLAGS) -c scauchy.cc

set_cmd_io_std.o : set_cmd_io_std.cc set_cmd_io_std.h
	$(CC) $(CFLAGS) -c set_cmd_io_std.cc

setflags.o : setflags.cc setflags.h
	$(CC) $(CFLAGS) -c setflags.cc

simanneal.o : simanneal.cc simanneal.h constants.h autocomm.h
	$(CC) $(OLIMIT) -c simanneal.cc

sort_enrg.o : sort_enrg.cc sort_enrg.h constants.h
	$(CC) $(CFLAGS) -c sort_enrg.cc

stop.o : stop.cc stop.h constants.h
	$(CC) $(CFLAGS) -c stop.cc

strindex.o : strindex.cc strindex.h
	$(CC) $(CFLAGS) -c strindex.cc

success.o : success.cc success.h constants.h autocomm.h
	$(CC) $(CFLAGS) -c success.cc

summarizegrids.o : summarizegrids.cc summarizegrids.h constants.h autocomm.h
	$(CC) $(CFLAGS) -c summarizegrids.cc

support.o : support.cc support.h eval.h
	$(CC) $(CFLAGS) -c support.cc

swap.o : swap.cc swap.h
	$(CC) $(CFLAGS) -c swap.cc

timesys.o : timesys.cc timesys.h
	$(CC) $(CFLAGS) -c timesys.cc

timesyshms.o : timesyshms.cc timesyshms.h
	$(CC) $(CFLAGS) -c timesyshms.cc

torNorVec.o : torNorVec.cc torNorVec.h constants.h
	$(CC) $(CFLAGS) -c torNorVec.cc

torsion.o : torsion.cc torsion.h constants.h
	$(CC) $(CFLAGS) -c torsion.cc

trilinterp.o : trilinterp.cc trilinterp.h constants.h
	$(CC) $(CFLAGS) -c trilinterp.cc

trilinterp.MINPT.o : trilinterp.cc trilinterp.h constants.h
	$(CC) $(CFLAGS) -DMINPT -c trilinterp.cc -o trilinterp.MINPT.o

usage.o : usage.cc usage.h constants.h
	$(CC) $(CFLAGS) -c usage.cc

warn_bad_file.o : warn_bad_file.cc warn_bad_file.h constants.h
	$(CC) $(CFLAGS) -c warn_bad_file.cc

weedbonds.o : weedbonds.cc weedbonds.h constants.h
	$(CC) $(CFLAGS) -c weedbonds.cc

eintcal.sqrt.o : eintcal.cc eintcal.h constants.h
	$(CC) $(CFLAGS) -c -DBOUNDED eintcal.cc -o eintcal.sqrt.o

eintcalPrint.sqrt.o : eintcal.cc eintcal.h constants.h
	$(CC) $(CFLAGS) -c -DBOUNDED -DEINTCALPRINT eintcal.cc -o eintcalPrint.sqrt.o

intnbtable.sqrt.o : intnbtable.cc intnbtable.h constants.h
	$(CC) $(OLIMIT) -c intnbtable.cc -o intnbtable.sqrt.o

nbe.sqrt.o : nbe.cc nbe.h constants.h
	$(CC) $(CFLAGS) -c nbe.cc -o nbe.sqrt.o

#
# lcheck dependencies...
#

analysis.ln : analysis.cc
	$(LINT) $(LINTFLAGS) $?

banner.ln : banner.cc
	$(LINT) $(LINTFLAGS) $?

bestpdb.ln : bestpdb.cc
	$(LINT) $(LINTFLAGS) $?

changeState.ln : changeState.cc
	$(LINT) $(LINTFLAGS) $?

check_header_float.ln : check_header_float.cc
	$(LINT) $(LINTFLAGS) $?

check_header_int.ln : check_header_int.cc
	$(LINT) $(LINTFLAGS) $?

check_header_line.ln : check_header_line.cc
	$(LINT) $(LINTFLAGS) $?

cluster_analysis.ln : cluster_analysis.cc
	$(LINT) $(LINTFLAGS) $?

clmode.ln : clmode.cc
	$(LINT) $(LINTFLAGS) $?

cmdmode.ln : cmdmode.cc
	$(LINT) $(LINTFLAGS) $?

cnv_state_to_coords.ln : cnv_state_to_coords.cc
	$(LINT) $(LINTFLAGS) $?

stateLibrary.ln : stateLibrary.cc
	$(CC) $(CFLAGS) -c $?

readfield.ln : readfield.cc
	$(LINT) $(LINTFLAGS) $?

readmap.ln : readmap.cc
	$(LINT) $(LINTFLAGS) $?

readPDBQ.ln : readPDBQ.cc
	$(LINT) $(LINTFLAGS) $?

dpftypes.ln : dpftypes.cc
	$(LINT) $(LINTFLAGS) $?

evaluate_energy.ln : evaluate_energy.cc
	$(LINT) $(LINTFLAGS) $?

getrms.ln : getrms.cc
	$(LINT) $(LINTFLAGS) $?

get_atom_type.ln : get_atom_type.cc
	$(LINT) $(LINTFLAGS) $?

getInitialState.ln : getInitialState.cc
	$(LINT) $(LINTFLAGS) $?

getpdbcrds.ln : getpdbcrds.cc
	$(LINT) $(LINTFLAGS) $?

initautodock.ln : initautodock.cc
	$(LINT) $(LINTFLAGS) $?

input_state.ln : input_state.cc
	$(LINT) $(LINTFLAGS) $?

investigate.ln : investigate.cc
	$(LINT) $(LINTFLAGS) $?

main.ln : main.cc
	$(LINT) $(LINTFLAGS) $?

mkNewState.ln : mkNewState.cc
	$(LINT) $(LINTFLAGS) $?

mkTorTree.ln : mkTorTree.cc
	$(LINT) $(LINTFLAGS) $?

mkRandomState.ln : mkRandomState.cc
	$(LINT) $(LINTFLAGS) $?

nonbonds.ln : nonbonds.cc
	$(LINT) $(LINTFLAGS) $?

openfile.ln : openfile.cc
	$(LINT) $(LINTFLAGS) $?

output_state.ln : output_state.cc
	$(LINT) $(LINTFLAGS) $?

parse_com_line.ln : parse_com_line.cc
	$(LINT) $(LINTFLAGS) $?

parse_dpf_line.ln : parse_dpf_line.cc
	$(LINT) $(LINTFLAGS) $?

parse_pdbq_line.ln : parse_pdbq_line.cc
	$(LINT) $(LINTFLAGS) $?

parse_trj_line.ln : parse_trj_line.cc
	$(LINT) $(LINTFLAGS) $?

print_2x.ln : print_2x.cc
	$(LINT) $(LINTFLAGS) $?

print_atomic_energies.ln : print_atomic_energies.cc
	$(LINT) $(LINTFLAGS) $?

print_avsfld.ln : print_avsfld.cc
	$(LINT) $(LINTFLAGS) $?

writeMolAsPDBQ.ln : writePDBQ.cc
	$(LINT) $(LINTFLAGS) -DWRITEMOLASPDBQFUNC $?

writePDBQ.ln : writePDBQ.cc
	$(LINT) $(LINTFLAGS) $?

writePDBQState.ln : writePDBQ.cc
	$(LINT) $(LINTFLAGS) -DWRITEPDBQSTATE $?

print_rem.ln : print_rem.cc
	$(LINT) $(LINTFLAGS) $?

printdate.ln : printdate.cc
	$(LINT) $(LINTFLAGS) $?

printEnergies.ln : printEnergies.cc
	$(LINT) $(LINTFLAGS) $?

printhms.ln : printhms.cc
	$(LINT) $(LINTFLAGS) $?

prClusterHist.ln : prClusterHist.cc
	$(LINT) $(LINTFLAGS) $?

prInitialState.ln : prInitialState.cc
	$(LINT) $(LINTFLAGS) $?

prTorConList.ln : prTorConList.cc
	$(LINT) $(LINTFLAGS) $?

qmultiply.ln : qmultiply.cc
	$(LINT) $(LINTFLAGS) $?

qtransform.ln : qtransform.cc
	$(LINT) $(LINTFLAGS) $?

quicksort.ln : quicksort.cc
	$(LINT) $(LINTFLAGS) $?

set_cmd_io_std.ln : set_cmd_io_std.cc
	$(LINT) $(LINTFLAGS) $?

setflags.ln : setflags.cc
	$(LINT) $(LINTFLAGS) $?

simanneal.ln : simanneal.cc
	$(LINT) $(LINTFLAGS) $?

sort_enrg.ln : sort_enrg.cc
	$(LINT) $(LINTFLAGS) $?

stop.ln : stop.cc
	$(LINT) $(LINTFLAGS) $?

strindex.ln : strindex.cc
	$(LINT) $(LINTFLAGS) $?

success.ln : success.cc
	$(LINT) $(LINTFLAGS) $?

summarizegrids.ln : summarizegrids.cc
	$(LINT) $(LINTFLAGS) $?

swap.ln : swap.cc
	$(LINT) $(LINTFLAGS) $?

timesys.ln : timesys.cc
	$(LINT) $(LINTFLAGS) $?

timesyshms.ln : timesyshms.cc
	$(LINT) $(LINTFLAGS) $?

torNorVec.ln : torNorVec.cc
	$(LINT) $(LINTFLAGS) $?

torsion.ln : torsion.cc
	$(LINT) $(LINTFLAGS) $?

trilinterp.ln : trilinterp.cc
	$(LINT) $(LINTFLAGS) $?

usage.ln : usage.cc
	$(LINT) $(LINTFLAGS) $?

weedbonds.ln : weedbonds.cc
	$(LINT) $(LINTFLAGS) $?

warn_bad_file.ln : warn_bad_file.cc
	$(LINT) $(LINTFLAGS) $?

#
# NOSQRT conditionals activated:
#

eintcal.ln : eintcal.cc
	$(LINT) $(LINTFLAGS) -DNOSQRT $?

eintcalPrint.ln : eintcal.cc
	$(LINT) $(LINTFLAGS) -DNOSQRT -DEINTCALPRINT $?

intnbtable.ln : intnbtable.cc
	$(LINT) $(LINTFLAGS) -DNOSQRT $?

nbe.ln : nbe.cc
	$(LINT) $(LINTFLAGS) -DNOSQRT $?

#
# NOSQRT conditionals NOT activated:
#

eintcal.sqrt.ln : eintcal.cc
	$(LINT) $(LINTFLAGS) $?

eintcalPrint.sqrt.ln : eintcal.cc
	$(LINT) $(LINTFLAGS) -DEINTCALPRINT $?

intnbtable.sqrt.ln : intnbtable.cc
	$(LINT) $(LINTFLAGS) $?

nbe.sqrt.ln : nbe.cc
	$(LINT) $(LINTFLAGS) $?

#
# Clean up...
#

clean :
	/bin/rm -f *.o *.s *.ln a.out mon.out autodock3 autodock3sqrt autodock3minpt dualmap libad.a

cleanlcheck :
	/bin/rm -f *.ln

#
# Compaq/DEC Alpha
#
openalpha :
	/bin/cp -p obj.alpha/*.o .
openlintalpha :
	/bin/cp -p obj.alpha/*.ln .
#
closealpha :
	/bin/cp -p *.o obj.alpha
	/bin/rm -f *.o
closelintalpha :
	/bin/cp -p *.ln obj.alpha
	/bin/rm -f *.ln

#
# Hewlett-Packard Precision Architecture
#
openhppa :
	/bin/cp -p obj.hppa/*.o .
openlinthppa :
	/bin/cp -p obj.hppa/*.ln .
#
closehppa :
	/bin/cp -p *.o obj.hppa
	/bin/rm -f *.o
closelinthppa :
	/bin/cp -p *.ln obj.hppa
	/bin/rm -f *.ln

#
# SGI/sgi4D/Silicon Graphics
#
opensgi4D :
	/bin/cp -p obj.sgi4D/*.o .
openlintsgi4D :
	/bin/cp -p obj.sgi4D/*.ln .
#
closesgi4D :
	/bin/cp -p *.o obj.sgi4D
	/bin/rm -f *.o
closelintsgi4D :
	/bin/cp -p *.ln obj.sgi4D
	/bin/rm -f *.ln

#
# Convex
#
openc2 :
	/bin/cp -p obj.c2/*.o .
openlintc2 :
	/bin/cp -p obj.c2/*.ln .
#
closec2 :
	/bin/cp -p *.o obj.c2
	/bin/rm -f *.o
closelintc2 :
	/bin/cp -p *.ln obj.c2
	/bin/rm -f *.ln

#
# Sun
#
opensun4 :
	/bin/cp -p obj.sun4/*.o .
openlintsun4 :
	/bin/cp -p obj.sun4/*.ln .
#
closesun4 :
	/bin/cp -p *.o obj.sun4
	/bin/rm -f *.o
closelintsun4 :
	/bin/cp -p *.ln obj.sun4
	/bin/rm -f *.ln

#
# EOF
#
