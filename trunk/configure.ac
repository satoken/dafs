#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.61)
AC_INIT(dusaf, 0.0.1, no-define)
AM_INIT_AUTOMAKE([-Wall foreign])
AM_CONFIG_HEADER(config.h)
AC_CONFIG_SRCDIR([src/dusaf.cpp])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_LANG(C++)
AC_CXX_HAVE_STL
AC_PROG_RANLIB

# debug option
AC_ARG_ENABLE([debug], AS_HELP_STRING([--enable-debug], [turn on debugging (default: off)]))
if test "$enable_debug" == "yes"; then
  CFLAGS="-DDEBUG $CFLAGS"
  CXXFLAGS="-DDEBUG $CXXFLAGS"
else
  CFLAGS="-DNDEBUG $CFLAGS"
  CXXFLAGS="-DNDEBUG $CXXFLAGS"
fi

# boost libraries
AX_BOOST_BASE(1.40, , AC_MSG_ERROR("the boost library is unavailable"))
if test "x$BOOST_CPPFLAGS" != x; then
  CPPFLAGS="$CPPFLAGS $BOOST_CPPFLAGS"
fi
if test "x$BOOST_LDFLAGS" != x; then
  LDFLAGS="$LDFLAGS $BOOST_LDFLAGS"
fi
dnl dnl AC_CHECK_HEADER(boost/lambda/lambda.hpp, ,
dnl dnl                 AC_MSG_ERROR("boost/lambda/lambda.hpp is unavailable"))
dnl AC_CHECK_HEADER(boost/multi_array.hpp, ,
dnl                 AC_MSG_ERROR("boost/multi_array.hpp is unavailable"))
dnl AC_CHECK_HEADER(boost/random.hpp, ,
dnl                 AC_MSG_ERROR("boost/random.hpp is unavailable"))

dnl AX_BOOST_PROGRAM_OPTIONS
dnl if test "x$BOOST_PROGRAM_OPTIONS_LIB" != x; then
dnl   LIBS="$LIBS $BOOST_PROGRAM_OPTIONS_LIB"
dnl else
dnl   AC_MSG_ERROR("the Boost::Program_Options library is unavailable")  
dnl fi

# Checks for libraries.
AC_ARG_WITH([vienna-rna], AS_HELP_STRING([--with-vienna-rna], [support Vienna RNA package]),
  if test "$with_vienna_rna" != "yes" && test "$with_vienna_rna" != "no" ; then
     CPPFLAGS="-I${with_vienna_rna}/include $CPPFLAGS"
     LDFLAGS="-L${with_vienna_rna}/lib $LDFLAGS"
  fi
)
AC_CHECK_LIB(RNA, init_pf_fold)

# Checks for header files.
AC_CHECK_HEADER([ViennaRNA/ribo.h],
		[AC_DEFINE([HAVE_VIENNA18], [], [have the Vienna RNA package version 1.8 or above])])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE

# Checks for library functions.
AC_HEADER_STDC

AC_CONFIG_FILES([Makefile src/Makefile src/contrafold/Makefile src/probconsRNA/Makefile src/contralign/Makefile src/partalign/Makefile])
AC_OUTPUT
