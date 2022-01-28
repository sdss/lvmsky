
AUTOMAKE_OPTIONS = foreign

# temp, precompiled manual
doc_DATA = VLT-MAN-ESO-19550-5896_4.3_Skycorr_User_Manual_SM-02.pdf

# find figures/ -type f -print | sed 's/^/  /;$q;s/$/ \\/' >> subfiles.mk 
include $(srcdir)/subfiles.mk
dist_noinst_DATA = dmd-doc.sty \
                   eso-logo.eps \
                   Makefile.am \
                   SC_UM_acronyms.tex \
                   SC_UM_evaluation.tex \
                   SC_UM_installation.tex \
                   SC_UM_introduction.tex \
                   SC_UM_layout.tex \
                   SC_UM_method.tex \
                   SC_UM_running.tex \
                   SC_UM_shortcut.tex \
                   skycorr_user_manual.tex \
                   $(subdir_files)

manual: $(dist_noinst_DATA) $(subdir_files)
	pdflatex skycorr_user_manual.tex
	pdflatex skycorr_user_manual.tex

CLEANFILES = rm -f *.dvi *.aux *.toc *.out *.log *.idx *.ps *.bbl *.blg *.brf \
	     skycorr_user_manual.pdf \
	     figures/*converted-to*.pdf \
	     *converted-to*.pdf

EXTRA_DIST = $(doc_DATA) $(dist_noinst_DATA)

.NOTPARALLEL:
