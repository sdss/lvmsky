LDADD = -lcplui -lcplcore $(MATHLIB)
AM_CFLAGS = -fPIC -O -Werror -pedantic -Wall -fno-common -Wextra \
    -Wmissing-prototypes -Wshadow -Wpointer-arith -Wcast-align \
    -Wwrite-strings -Wnested-externs -Wundef -Wsign-compare -std=c99 \
    -Wmissing-declarations -Wmissing-format-attribute -Wold-style-definition \
    -Wdeclaration-after-statement -Wdisabled-optimization -Wformat-y2k \
    -Winit-self -Wmissing-include-dirs -Wmissing-noreturn -Wredundant-decls \
    -Wvariadic-macros -Wvolatile-register-var -Wstrict-prototypes
AM_CFLAGS += $(CPL_INCLUDES)

SUBDIRS = .

include_HEADERS = sm_general.h sm_scatmoonlight.h sm_skyemcomp.h \
                  sm_linetrans.h sm_run.h

lib_LIBRARIES = libsm_general.a libsm_skyemcomp.a libsm_run.a
libsm_general_a_SOURCES = sm_general.c sm_general.h
libsm_skyemcomp_a_SOURCES = sm_skyemcomp.c sm_skyemcomp.h \
                            sm_scatmoonlight.c sm_scatmoonlight.h
libsm_run_a_SOURCES = sm_run.c sm_run.h

bin_PROGRAMS = preplinetrans calcskymodel testskymodel estmultiscat

preplinetrans_SOURCES = preplinetrans.c \
                        sm_general.c \
                        sm_linetrans.c \
                        sm_scatmoonlight.c \
                        sm_skyemcomp.c
preplinetrans_CFLAGS = $(AM_CFLAGS)
preplinetrans_LIBS = $(CPLLIB)

calcskymodel_SOURCES = calcskymodel.c \
                       sm_general.c \
                       sm_scatmoonlight.c \
                       sm_skyemcomp.c \
                       sm_run.c
calcskymodel_CFLAGS = $(AM_CFLAGS)
calcskymodel_LIBS = $(CPLLIB)

testskymodel_SOURCES = testskymodel.c \
                       sm_general.c \
                       sm_scatmoonlight.c \
                       sm_skyemcomp.c \
                       sm_run.c
testskymodel_CFLAGS = $(AM_CFLAGS)
testskymodel_LIBS = $(CPLLIB)

estmultiscat_SOURCES = estmultiscat.c \
                       sm_general.c \
                       sm_scatmoonlight.c \
                       sm_skyemcomp.c \
                       sm_run.c
estmultiscat_CFLAGS = $(AM_CFLAGS)
estmultiscat_LIBS = $(CPLLIB)
