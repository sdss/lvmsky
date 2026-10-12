"""Check normalization of the nightly profile visualization."""
import sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
import numpy as np
import pandas as pd
from publish_lunation_report import profile_summary,plot_fiber_metrics,SCATTER_METRICS,ABSOLUTE_SCATTER_METRICS
from measurement_view import exposure_positions,measurement_markup

rows=pd.DataFrame({f'lsf_mass_{i:02d}':[1/11,1/11] for i in range(11)})
rows['offset_half_width_angstrom']=2.5
x,y=profile_summary(rows)
assert y.min()>=0 and np.isclose(np.trapz(y,x),1,atol=1e-6)
assert np.isclose(np.trapz(y*x,x),0,atol=1e-6)
print('Nightly profile normalization checks passed')

coverage=pd.DataFrame([dict(expnum=expnum,fiberid=fiber,spectrographid=fiber,telescope=telescope,channel=channel,
                           accepted=(fiber!=2),gaia_verified=(fiber!=2))
    for fiber,telescope in [(1,'Sci'),(2,'SkyE'),(3,'SkyW')]
    for channel in ['B','R','Z'] for expnum in range(3)])
repeat=pd.DataFrame([dict(fiberid=fiber,telescope=telescope,channel=channel,
    **{metric:(70. if fiber==1 and channel=='B' else 0.) for metric,_,_ in SCATTER_METRICS})
    for fiber,telescope in [(1,'Sci'),(3,'SkyW')] for channel in ['B','R','Z']])
fig=plot_fiber_metrics(coverage,repeat,3,'2025-01-01',60000)
assert len(fig.axes)==15
assert fig.axes[3].get_yscale()=='symlog' and fig.axes[3].get_ylim()[1]>70
assert fig.axes[3].collections[0].get_offsets()[0].tolist()==[1.,70.]
assert fig.axes[5].collections[0].get_offsets()[0].tolist()==[3.,0.]
assert 'Sky East' in fig.axes[1].get_title()
assert len(fig.axes[1].collections[1].get_offsets())==1
assert fig.axes[0].collections[0].get_offsets()[0].tolist()==[1.,3.]
widths=[fig.axes[i].get_position().width for i in range(3)]
assert np.isclose(widths[0]/widths[1],2) and np.isclose(widths[1]/widths[2],1)
assert [float(line.get_xdata()[0]) for line in fig.axes[0].lines]==[1.5,2.5]
import matplotlib.pyplot as plt
plt.close(fig)
print('Grouped fiber plots retain zero scatter, outliers and unknown-coverage markers')

for metric,_,_ in ABSOLUTE_SCATTER_METRICS:
    repeat[metric]=np.where(repeat.fiberid.eq(1)&repeat.channel.eq('B'),.125,0.)
fig=plot_fiber_metrics(coverage,repeat,3,'2025-01-01',60000,absolute=True)
assert fig.axes[3].collections[0].get_offsets()[0].tolist()==[1.,.125]
assert fig.axes[3].get_ylabel()=='FWHM scatter (Å)'
assert fig.axes[6].get_ylabel()=='W50 scatter (Å)'
assert fig.axes[9].get_ylabel()=='W90 scatter (Å)'
assert fig.axes[3].get_ylim()[1]>.125
plt.close(fig)
print('Absolute width scatter plots use stored Angstrom spreads directly')

levels=pd.DataFrame([dict(fiberid=fiber,telescope=telescope,spectrographid=sp,channel=channel,
    **{f'{metric}_{statistic}':value for metric in ['fwhm_angstrom','w50_angstrom','w90_angstrom','centroid_angstrom']
       for statistic,value in [('mean',fiber+.1),('median',float(fiber))]})
    for fiber,telescope,sp in [(1,'Sci',1),(4,'Sci',2),(7,'SkyW',3)] for channel in ['B','R','Z']])
coverage_extra=coverage[coverage.fiberid.eq(1)].copy();coverage_extra['fiberid']=4;coverage_extra['spectrographid']=2
fixture_coverage=coverage.copy();fixture_coverage.loc[fixture_coverage.fiberid.eq(3),'fiberid']=7
fixture_repeat=repeat.copy();fixture_repeat.loc[fixture_repeat.fiberid.eq(3),'fiberid']=7
fig=plot_fiber_metrics(pd.concat([fixture_coverage,coverage_extra]),fixture_repeat,3,'2025-01-01',60000,absolute=True,levels=levels)
assert len(fig.axes)==27
assert len({tuple(np.round(ax.get_position().bounds,6)) for ax in fig.axes})==27
assert fig.axes[3].get_yscale()=='linear'
assert fig.axes[6].get_yscale()=='symlog'
value_lines=[line for line in fig.axes[3].lines if line.get_linestyle() in ['-','--'] and line.get_color() in ['#0072B2','#D55E00','#009E73']]
assert len(value_lines)==12
for line in value_lines:
    ids=np.asarray(line.get_xdata());assert len(ids)==1
    assert line.get_marker()=='_'
    expected=ids.astype(float)+(0.1 if line.get_linestyle()=='--' else 0.)
    assert np.allclose(line.get_ydata(),expected)
assert fig.axes[3].get_ylabel()=='FWHM value (Å)'
assert fig.axes[6].get_ylabel()=='FWHM scatter (Å)'
assert fig.get_size_inches()[1]==21
plt.close(fig)
print('Mean/median values and scatter occupy separate panels, each with one Y-axis; lines stop at spectrograph boundaries')

positions=exposure_positions(pd.DataFrame({'expnum':[3,1,2],
    'obstime':['2025-02-01','2025-01-01','2025-01-02'],
    'full_moon_utc':['2025-02-01','2025-01-01','2025-01-01']}))
assert positions.expnum.tolist()==[1,2,3]
assert positions.lunation_number.tolist()==[1,1,2]
assert np.allclose(positions.plot_x,[.8,1.2,2.])
print('Observed-lunation blocks preserve exposure ordering and singleton periods')

markup=measurement_markup({'periods':[dict(lunation_number=1,full_moon_utc='2025-01-01',selected_exposures=3,epoch_start_mjd='60000')],
    'fibers':[dict(fiberid=1,telescope='SkyW',spectrographid=1,exposures=3,measurements=9),
              dict(fiberid=3,telescope='Sci',spectrographid=1,exposures=3,measurements=9)]})
assert markup.count('class="measurement-figure"')==2
assert 'measurement_previews/fiber-0001.png' in markup and 'measurement_previews/fiber-0003.png' in markup
assert markup.count('<select class="measurement-fiber-control">')==2
assert 'measurement-telescope' not in markup and 'measurement-spectrograph' not in markup
assert '<details' not in markup and '<summary' not in markup
assert markup.index('data-fiber="3"')<markup.index('data-fiber="1"')
print('Telescope/spectrograph groups stay visible; controls choose only a physical fiber within a group')
