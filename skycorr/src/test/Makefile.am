
WERR = #-Werror

AM_CFLAGS = -O2 -std=c99 -fno-common -Wextra -Wall \
    -Wmissing-prototypes -Wshadow -Wpointer-arith -Wcast-align \
    -Wwrite-strings -Wnested-externs -Wundef -Wsign-compare \
    -Wmissing-declarations -Wmissing-format-attribute -Wold-style-definition \
    -Wdisabled-optimization -Wformat-y2k \
    -Winit-self -Wmissing-include-dirs -Wmissing-noreturn -Wredundant-decls \
    -Wvariadic-macros -Wvolatile-register-var -Wstrict-prototypes \
    $(WERR)

AM_CPPFLAGS = $(CPL_INCLUDES) -I$(top_srcdir)/src
AM_LDFLAGS = $(CPL_LDFLAGS)

LDADD = $(LIBCPLCORE) $(LIBCPLUI) $(LIBCPLDFS) $(LIBCPLDRS) $(MATHLIB)

TESTS = test_sc_basic test_sc_par test_sc_specdiss test_sc_interpollin \
    test_sc_contsub test_sc_fwhmest test_sc_plot test_sc_skycorr

check_PROGRAMS = test_sc_basic test_sc_par test_sc_specdiss \
    test_sc_interpollin test_sc_contsub test_sc_fwhmest \
    test_sc_plot test_sc_skycorr

EXTRA_DIST = config/sctest_sinfo_H.par \
             data/sky_sinfo_1.fits \
             data/sky_sinfo_2.fits


test_sc_basic_SOURCES = test_sc_basic.c \
    $(top_srcdir)/src/sc_basic.c

test_sc_par_SOURCES = test_sc_par.c \
    $(top_srcdir)/src/sc_par.c \
    $(top_srcdir)/src/sc_basic.c

test_sc_specdiss_SOURCES = test_sc_specdiss.c \
    $(top_srcdir)/src/sc_specdiss.c \
    $(top_srcdir)/src/sc_lines.c \
    $(top_srcdir)/src/sc_basic.c \
    $(top_srcdir)/src/sc_par.c \
    $(top_srcdir)/src/sc_readspec.c \
    $(top_srcdir)/src/sc_conv.c \
    $(top_srcdir)/src/sc_plot.c

test_sc_interpollin_SOURCES = test_sc_interpollin.c \
    $(top_srcdir)/src/sc_basic.c

test_sc_contsub_SOURCES = test_sc_contsub.c \
    $(top_srcdir)/src/sc_contsub.c \
    $(top_srcdir)/src/sc_specdiss.c \
    $(top_srcdir)/src/sc_lines.c \
    $(top_srcdir)/src/sc_basic.c \
    $(top_srcdir)/src/sc_par.c \
    $(top_srcdir)/src/sc_plot.c \
    $(top_srcdir)/src/sc_readspec.c \
    $(top_srcdir)/src/sc_conv.c

test_sc_fwhmest_SOURCES = test_sc_fwhmest.c \
    $(top_srcdir)/src/mpfit.c \
    $(top_srcdir)/src/sc_fwhmest.c \
    $(top_srcdir)/src/sc_lines.c \
    $(top_srcdir)/src/sc_contsub.c \
    $(top_srcdir)/src/sc_specdiss.c \
    $(top_srcdir)/src/sc_basic.c \
    $(top_srcdir)/src/sc_par.c \
    $(top_srcdir)/src/sc_plot.c \
    $(top_srcdir)/src/sc_readspec.c \
    $(top_srcdir)/src/sc_conv.c

test_sc_plot_SOURCES = test_sc_plot.c \
    $(top_srcdir)/src/sc_par.c \
    $(top_srcdir)/src/sc_readspec.c \
    $(top_srcdir)/src/sc_conv.c \
    $(top_srcdir)/src/sc_basic.c \
    $(top_srcdir)/src/sc_plot.c

test_sc_skycorr_SOURCES = test_sc_skycorr.c \
    $(top_srcdir)/src/sc_skycorr.c \
    $(top_srcdir)/src/sc_par.c \
    $(top_srcdir)/src/sc_basic.c \
    $(top_srcdir)/src/sc_readspec.c \
    $(top_srcdir)/src/sc_conv.c \
    $(top_srcdir)/src/sc_lines.c \
    $(top_srcdir)/src/sc_specdiss.c \
    $(top_srcdir)/src/sc_contsub.c \
    $(top_srcdir)/src/sc_fwhmest.c \
    $(top_srcdir)/src/sc_weights.c \
    $(top_srcdir)/src/sc_mpfit.c \
    $(top_srcdir)/src/mpfit.c \
    $(top_srcdir)/src/sc_modsky.c
