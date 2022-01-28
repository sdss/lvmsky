
AUTOMAKE_OPTIONS = foreign
DISTCLEANFILES = *.par

nobase_examples_DATA = \
  data/sky_sinfo_7.fits \
  data/NGC4594-Xshooter_1.fits \
  data/sky_sinfo_2.fits \
  data/sky_xshoo_28.fits \
  data/SCI_SLIT_MERGE1D_1h_NIR.fits \
  data/NGC4594-Sinfoni_H.fits \
  data/NGC4594-Sinfoni_J.fits \
  data/sky_sinfo_5.fits \
  data/SKY_SLIT_MERGE1D_1h_UVB.fits \
  data/SKY_SLIT_MERGE1D_3h_UVB.fits \
  data/SCI_SLIT_MERGE1D_1h_UVB_ref.fits \
  data/SCI_SLIT_MERGE2D_3h_NIR.fits \
  data/SCI_SLIT_MERGE1D_1h_VIS_ref.fits \
  data/sky_sinfo_4.fits \
  data/SCI_SLIT_MERGE2D_3h_VIS.fits \
  data/sky_sinfo_6.fits \
  data/SKY_SLIT_MERGE1D_1h_VIS.fits \
  data/sky_fors_0919.fits \
  data/SKY_SLIT_MERGE1D_1h_NIR.fits \
  data/SKY_SLIT_MERGE1D_3h_NIR.fits \
  data/SCI_SLIT_MERGE2D_3h_UVB.fits \
  data/SCI_SLIT_MERGE1D_1h_VIS.fits \
  data/sky_sinfo_1.fits \
  data/sky_fors_0918.fits \
  data/SCI_SLIT_MERGE1D_1h_NIR_ref.fits \
  data/SCI_SLIT_MERGE2D_1h_NIR.fits \
  data/SCI_SLIT_MERGE1D_1h_UVB.fits \
  data/SCI_SLIT_MERGE2D_1h_UVB.fits \
  data/SKY_SLIT_MERGE1D_3h_VIS.fits \
  data/SCI_SLIT_MERGE2D_1h_VIS.fits \
  data/NGC4594-FORS_1.fits \
  data/sky_xshoo_29.fits \
  config/XPL_d0h_NIR.par \
  config/XPL_d0h_VIS.par \
  config/XPL_d0h_UVB.par \
  config/XPL_d2h_VIS.par \
  config/sctest_xshoo_1.par \
  config/sctest_xshoo_1_N4594.par \
  config/sctest_fors_1_N4594.par \
  config/sctest_sinfo_H.par \
  config/sctest_sinfo_J.par \
  config/sctest_sinfo_K.par \
  config/sctest_sinfo_J_N4594.par \
  config/sctest_sinfo_H_N4594.par \
  config/sctest_fors_1.par \
  config/XPL_d2h_UVB.par \
  config/XPL_d2h_NIR.par

EXTRA_DIST = $(nobase_examples_DATA)

install-data-hook:
	ln -sf ../sysdata $(examplesdir)/sysdata

