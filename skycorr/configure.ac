#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.


AC_PREREQ([2.63])
AC_INIT([skycorr], [1.1.2], [usd-help@eso.org])
AC_CONFIG_SRCDIR([src/sc_basic.c])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile
                 src/Makefile
                 src/test/Makefile
                 sysdata/Makefile
                 gui/Makefile
                 reflex/Makefile
                 reflex/skycorr_reflex.xml
                 examples/Makefile
                 examples/config/sctest_xshoo_1_N4594.par
                 examples/config/sctest_sinfo_H.par
                 examples/config/XPL_d0h_NIR.par
                 examples/config/XPL_d0h_VIS.par
                 examples/config/XPL_d0h_UVB.par
                 examples/config/XPL_d2h_VIS.par
                 examples/config/XPL_d2h_UVB.par
                 examples/config/sctest_xshoo_1.par
                 examples/config/sctest_fors_1_N4594.par
                 examples/config/sctest_fors_1.par
                 examples/config/sctest_sinfo_J.par
                 examples/config/sctest_sinfo_K.par
                 examples/config/sctest_sinfo_J_N4594.par
                 examples/config/sctest_sinfo_H_N4594.par
                 examples/config/XPL_d2h_NIR.par
                 doc/Makefile])
AM_INIT_AUTOMAKE([-Wall foreign])

# Checks for programmes.
AC_PROG_CC
AM_PROG_CC_C_O

# Check for CPL presence and usability
CPL_CHECK_LIBS

AC_CHECK_LIB([m],[pow],[MATHLIB=-lm],[
    AC_MSG_ERROR([Math library (libm) not found!])
])
AC_SUBST([MATHLIB])

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h string.h stdio.h ctype.h math.h],[],[
    AC_MSG_ERROR(["A standard header file was not found!"])
])

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_FUNC_STRTOD
AC_CHECK_FUNCS([mkdir strchr strpbrk strrchr strspn strstr strtol])
AM_MAINTAINER_MODE([enable])

AC_CHECK_PROG(SVNVERSION, svnversion, found, not-found)
if test $SVNVERSION = found; then
  SVNREV=$(svnversion "$srcdir")
  echo $SVNREV > SVNREV
else
  SVNREV=$(cat SVNREV)
fi
AC_SUBST(SVNREV)

if test -z "$wkfextradir"; then
    wkfextradir='${prefix}/reflex'
fi

if test -z "$configdir"; then
    configdir='${prefix}/config'
fi

if test -z "$pyguidir"; then
    pyguidir='${prefix}/gui'
fi

if test -z "$examplesdir"; then
    examplesdir='${prefix}/examples'
fi

if test -z "$sysdatadir"; then
    sysdatadir='${prefix}/sysdata'
fi

docdir='${prefix}/doc'

AC_SUBST(configdir)
AC_SUBST(wkfextradir)
AC_SUBST(pyguidir)
AC_SUBST(examplesdir)
AC_SUBST(docdir)
AC_SUBST(sysdatadir)

AC_OUTPUT
