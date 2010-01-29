#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
# $Id: configure.ac,v 1.15 2010/01/29 19:31:06 rhuey Exp $

# for release replace x with appropriate version (1 was last 9/09)
AC_INIT([autodock], [4.2.3], [autodock@scripps.edu])
AC_CONFIG_SRCDIR([analysis.cc])

AM_INIT_AUTOMAKE()
#AM_INIT_AUTOMAKE([foreign])
#AM_INIT_AUTOMAKE([-Wall -Werror foreign])


# Checks for programs.
AC_PROG_AWK

# don't let AC_PROG_CXX force "-g -O2"
saved_CXXFLAGS=${CXXFLAGS}
AC_PROG_CXX
CXXFLAGS=${saved_CXXFLAGS}

AC_PROG_INSTALL
AC_PROG_RANLIB

AC_CANONICAL_HOST

# Checks for header files.
AC_HEADER_STDC
#AC_CHECK_HEADERS([fcntl.h limits.h stdlib.h string.h sys/param.h sys/time.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T
AC_STRUCT_TM

# Checks for library functions.
AC_FUNC_MALLOC
AC_CHECK_FUNCS([gethostname pow sqrt strncasecmp])

case $host_os in
  *cygwin* ) CYGWIN=yes;;
         * ) CYGWIN=no;;
esac
AM_CONDITIONAL([WANT_BIGSTACK], [test "$CYGWIN" = yes])

# See "automake" manual, version 1.10, section 10.5 (Python)
#AM_PATH_PYTHON([2.1],, [:])
#AM_CONDITIONAL([HAVE_PYTHON], [test "$PYTHON" != :])

AC_CONFIG_FILES([Makefile])

AC_OUTPUT
