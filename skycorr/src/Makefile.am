
WERR = #-Werror

AM_CFLAGS = -O2 -std=c99 -fno-common -Wextra -Wall \
    -Wmissing-prototypes -Wshadow -Wpointer-arith -Wcast-align \
    -Wwrite-strings -Wnested-externs -Wundef -Wsign-compare \
    -Wmissing-declarations -Wmissing-format-attribute -Wold-style-definition \
    -Wdisabled-optimization -Wformat-y2k \
    -Winit-self -Wmissing-include-dirs -Wmissing-noreturn -Wredundant-decls \
    -Wvariadic-macros -Wvolatile-register-var -Wstrict-prototypes \
    $(WERR)

AM_CPPFLAGS = $(CPL_INCLUDES)

SUBDIRS = . test

EXTRA_DIST = DISCLAIMER

# binary programmes
bin_PROGRAMS = skycorr extract1d

skycorr_SOURCES = skycorr.c \
                  sc_skycorr.c sc_skycorr.h \
                  sc_plot.c sc_plot.h \
                  sc_par.c sc_par.h \
                  sc_basic.c sc_basic.h \
                  sc_readspec.c sc_readspec.h \
                  sc_conv.c sc_conv.h \
                  sc_lines.c sc_lines.h \
                  sc_specdiss.c sc_specdiss.h \
                  sc_contsub.c sc_contsub.h \
                  sc_fwhmest.c sc_fwhmest.h \
                  sc_weights.c sc_weights.h \
                  sc_mpfit.c sc_mpfit.h \
                  mpfit.c mpfit.h \
                  sc_modsky.c sc_modsky.h
skycorr_LDFLAGS = $(CPL_LDFLAGS)
skycorr_LDADD = $(LIBCPLCORE) $(LIBCPLUI) $(LIBCPLDFS) $(LIBCPLDRS) $(MATHLIB)

extract1d_SOURCES = extract1d.c \
                    sc_par.c sc_par.h \
                    sc_basic.c sc_basic.h \
                    sc_conv.c sc_conv.h
extract1d_LDFLAGS = $(CPL_LDFLAGS)
extract1d_LDADD = $(LIBCPLCORE) $(LIBCPLUI) $(LIBCPLDFS) $(LIBCPLDRS) $(MATHLIB)
