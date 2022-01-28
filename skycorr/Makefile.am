AUTOMAKE_OPTIONS = 1.8 foreign

ACLOCAL_AMFLAGS = -I m4macros

DISTCLEANFILES = *~ SVNREV

EXTRA_DIST = m4macros/eso.m4 m4macros/cpl.m4 SVNREV

dist_config_DATA = config/template_parameter_file.par

dist_noinst_SCRIPTS = scripts/install_package.sh \
		      scripts/skycorr_installer.run

SUBDIRS = src doc sysdata examples gui reflex

# Files to be removed with command: make maintainer-clean
MAINTAINERCLEANFILES = Makefile.in aclocal.m4 config.h.in config.log \
        config.status configure stamp-h.in configure.scan depcomp \
        install-sh missing 
