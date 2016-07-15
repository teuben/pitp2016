dnl ###

AC_INIT(README)
AC_CONFIG_AUX_DIR(./config)

# try to turn off auto-invocation of autotools in the top-level Makefile
AM_MAINTAINER_MODE

# Default prefix
AC_PREFIX_DEFAULT(/usr/local/ca)



dnl ######################################################################
dnl
dnl Standard tests
dnl
dnl ######################################################################

AC_CANONICAL_HOST
hosttype=$host
AC_DEFINE_UNQUOTED(HOSTTYPE, "$host")


dnl # C, C++ compiler and some other useful tools

AC_PROG_CC
AC_PROG_CXX
AC_PROG_AWK
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET
AC_PROG_RANLIB
AC_ARG_PROGRAM

dnl # make
AC_CHECK_PROGS(MAKE, [gmake make])
if test x"${MAKE}" = x
then
    AC_MSG_ERROR([cannot find a make utility])
else
    AC_MSG_CHECKING([whether ${MAKE} is GNU Make])
    if ${MAKE} --version 2>/dev/null | head -n 1 |
       grep "GNU Make" >/dev/null
    then
        AC_MSG_RESULT([yes])
    else
        AC_MSG_RESULT([no])
        AC_MSG_ERROR([GNU Make is required, but ${MAKE} is
                      different])
    fi
fi


dnl # Fortran compiler (the f77.m4 replaced the need for AC_PROG_FC and AC_PROG_F77)

AC_F77_LIBRARY_LDFLAGS 
builtin(include,config/f77.m4)


dnl # size of things

AC_CHECK_SIZEOF(char)
AC_CHECK_SIZEOF(short)
AC_CHECK_SIZEOF(int)
AC_CHECK_SIZEOF(long)
AC_CHECK_SIZEOF(long long)
AC_CHECK_SIZEOF(float)
AC_CHECK_SIZEOF(double)
AC_CHECK_SIZEOF(long double)



builtin(include,config/hdf5.m4)

builtin(include,config/gsl.m4)

dnl builtin(include,config/fftw.m4)

builtin(include,config/python.m4)

builtin(include,config/mpicompilers.m4)

dnl builtin(include,config/cfitsio.m4)

dnl GSL

dnl GMP

builtin(include,config/ds9.m4)
builtin(include,config/cvs.m4)

dnl builtin(include,config/lapack.m4)

dnl builtin(include,config/libs.m4)

dnl builtin(include,config/qt.m4)

builtin(include,config/distro.m4)

AC_CONFIG_FILES([ Makefile ])
AC_CONFIG_FILES([ Hello_World/makedefs ])


AC_OUTPUT

echo "-------------------------------------------------------------------"
echo "Now follow some simple tests to confirm compiler sanity:"
echo "cd Hello_World ; make tests"
