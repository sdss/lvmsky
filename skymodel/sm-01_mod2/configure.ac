# Configure script for Sky Model Module 2

AC_PREREQ([2.59])
AC_INIT([SkyModel_Mod2], [1.0], [cpl-help@eso.org])
AC_CONFIG_SRCDIR([src/sm_general.c])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile src/Makefile])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])

# Checks for programmes.
AC_PROG_CC
CFLAGS="-std=c99 -Wno-error $CFLAGS"
AM_PROG_CC_C_O
AC_PROG_RANLIB
m4_ifdef([AM_PROG_AR], [AM_PROG_AR])
AC_ARG_WITH([cpl],
            [AS_HELP_STRING([--with-cpl],
                            [path to the cpl base path (subdirs 'lib' and 'include' must exist)])
            ],
            [CPLDIR=$withval],
            [])
AS_IF([test "x$CPLDIR" != x],
      [ CPLLDFLAG="-L$CPLDIR/lib"
        CPLINC="-I$CPLDIR/include"
        LDFLAGS="$CPLLDFLAG $LDFLAGS"
        CFLAGS="$CPLINC $CFLAGS"
        CPPFLAGS="$CPLINC $CPPFLAGS"
      ],
      [])

# Check for CPL presence and usability
CPL_CHECK_LIBS

# Checks for libraries.
AC_CHECK_LIB([cplcore],[cpl_table_dump],
             [CPLLIB="-L$CPLDIR/lib -lcplui -lcplcore -lcext"],
             [AC_MSG_ERROR([CPL libraries not found!])
])
AC_SUBST([CPLLIB])
AC_CHECK_LIB([m], [pow], [MATHLIB="-lm"], [
  AC_MSG_ERROR([Math library (libm) not found!])
])
AC_SUBST([MATHLIB])

# Checks for header files.
AC_CHECK_HEADERS([stdlib.h string.h stdio.h ctype.h math.h],[],[
    AC_MSG_ERROR(["A standard header file was not found!"])
])
AC_CHECK_HEADER([cpl.h],[],[
    AC_MSG_ERROR(["cpl.h not found!"])
])

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_FUNC_STRTOD
AC_CHECK_FUNCS([mkdir strchr strpbrk strrchr strspn strstr strtol])
AM_MAINTAINER_MODE
AC_OUTPUT
