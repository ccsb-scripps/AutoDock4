AC_INIT([autodock], [4.0.0], [autodock@scripps.edu])
AM_INIT_AUTOMAKE()
#AM_INIT_AUTOMAKE([foreign])
#AM_INIT_AUTOMAKE([-Wall -Werror foreign])
AC_PROG_RANLIB
AC_PROG_CXX
AC_CANONICAL_HOST
case $host_os in
  *cygwin* ) CYGWIN=yes;;
         * ) CYGWIN=no;;
esac
AM_CONDITIONAL([WANT_BIGSTACK], [test "$CYGWIN" = yes])
AC_CONFIG_FILES([Makefile])
AC_OUTPUT
