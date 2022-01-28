
AUTOMAKE_OPTIONS = foreign
DISTCLEANFILES = *.xml

WORKFLOWS = skycorr_reflex.xml

#This installs the workflow and the OCA rules in reflexaux directory
wkfextra_DATA = $(WORKFLOWS)

wkfextra_SCRIPTS =  create_sc_parfile.sh

#Pack in the distribution the original workflow
EXTRA_DIST = $(WORKFLOWS).in $(wkfextra_SCRIPTS)
