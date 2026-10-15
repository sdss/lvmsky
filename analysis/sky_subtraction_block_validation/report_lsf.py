"""Full-wavelength block maps and unrecentered individual Science profile comparisons."""
from pathlib import Path
import json
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from skysub.sky_decomp.lsf_spline2d import evaluate_lsf_diagnostics,mspline_basis
from skysub.sky_decomp.lsf_surface_iterative import evaluate_bspline_basis
from skysub.sky_subtraction_block.lsf import evaluate_block

HERE=Path(__file__).resolve().parent
CORPUS=Path('/data/sas/sdsswork/users/u6058164/lvm_lsf_epoch_coverage_drp132')
METRICS={'fwhm_angstrom':'FWHM (Å)','w50_angstrom':'W50 (Å)','w90_angstrom':'W90 (Å)','centroid_angstrom':'Numerical centering error (Å)','asymmetry':'Quantile asymmetry'}


def surface_map(r,expnum):
    template=r.state('MED_',r.labels.index('Sci'));slit=r.h['SLITMAP'].data
    fig,axes=plt.subplots(5,2,figsize=(12,15),layout='constrained')
    for j,ch in enumerate('RZ'):
        knots=template.knot_vectors[ch];degree=template.degrees[ch]
        wave=np.linspace(knots[degree]+.001,knots[-degree-1]-.001,33)
        rows=[];positions=[];boundaries=[]
        for sp in [1,2,3]:
            domain=slit['fiberid'][slit['spectrographid']==sp]
            for position in np.linspace(domain.min(),domain.max(),49):
                state=evaluate_block(template,r.blocks,sp,position)
                rows.append(evaluate_lsf_diagnostics(state,wave));positions.append(position)
            boundaries.append(float(domain.max())+.5)
        for i,(key,label) in enumerate(METRICS.items()):
            values=np.stack([row[key] for row in rows]);signed=key in ('centroid_angstrom','asymmetry')
            limit=max(float(np.nanmax(abs(values))),1e-6) if signed else None
            ax=axes[i,j];im=ax.pcolormesh(wave,positions,values,shading='nearest',cmap='RdBu_r' if signed else 'viridis',vmin=-limit if signed else float(np.nanmin(values)),vmax=limit if signed else float(np.nanmax(values)))
            ax.set_ylabel('Physical fiber ID');ax.set_xlabel('Air wavelength (Å)')
            for boundary in boundaries[:-1]:ax.axhline(boundary,lw=.6,color='black')
            if i==0:ax.set_title(ch+' · full fitted domain')
            fig.colorbar(im,ax=ax,label=label,shrink=.83)
    fig.suptitle(f'Exposure {expnum}: centered block LSF surfaces, R and Z\n49 positions per spectrograph × 33 wavelengths per channel',fontsize=16)
    fig.savefig(HERE/f'global_lsf_maps_{expnum}.png',dpi=120);plt.close(fig)


def parameter_files():
    root=Path('/data/sas/sdsswork/users/u6058164/lvm_solarblue_all_fibers_drp132_moonalt60_fli095/parameters')
    return sorted(root.glob('solarblue-shard-*.fits'))+sorted((CORPUS/'additional_clean_fits/parameters').glob('solarblue-shard-*.fits'))


def profile_comparison(r,expnum,measured):
    selected=measured[(measured.expnum==expnum)&(measured.telescope=='Sci')]
    wanted=[]
    for sp in [1,2,3]:
        group=selected[selected.spectrographid==sp]
        counts=group[group.channel.isin(['R','Z'])].groupby('fiberid').channel.nunique()
        complete=counts[counts==2].index.to_numpy()
        if len(complete):wanted.append(int(complete[len(complete)//2]))
    if not wanted:return ''
    profiles={};sources=[]
    for path in parameter_files():
        with fits.open(path,memmap=True) as h:
            meta=h['INPUT_META'].data
            match=np.flatnonzero((meta['expnum']==expnum)&np.isin(meta['fiberid'],wanted))
            if not len(match):continue
            config={str(row['channel']):row for row in h['LSF_CONFIG'].data}
            for index in match:
                f=int(meta['fiberid'][index]);profiles[f]={}
                for ch in 'RZ':
                    j='BRZ'.index(ch);conf=config[ch];profiles[f][ch]=(np.asarray(h['LSF_COEF'].data[index,j],float).copy(),
                        np.asarray(h['LSF_KNOTS'].data[index,j],float).copy(),int(conf['degree']),float(conf['offset_half_width_angstrom']))
            sources.append(str(path))
        if len(profiles)==len(wanted):break
    if len(profiles)!=len(wanted):raise ValueError('Accepted Science profile has no matching parameter FITS')
    fig,axes=plt.subplots(2*len(wanted),2,figsize=(11,4.5*len(wanted)),layout='constrained')
    checks=[]
    for n,f in enumerate(wanted):
        state=r.global_state(f)
        for j,ch in enumerate('RZ'):
            point=selected[(selected.fiberid==f)&(selected.channel==ch)].iloc[0];wave=float(point.wavelength_angstrom)
            coef,knots,degree,width=profiles[f][ch]
            support=max(width,state.config['offset_half_width_angstrom'])
            delta=np.linspace(-support,support,1601)
            reference=evaluate_bspline_basis(np.array([wave]),knots,degree)@coef.T@mspline_basis(delta,n_basis=coef.shape[0],half_width=width).T
            block=evaluate_bspline_basis(np.array([wave]),state.knot_vectors[ch],state.degrees[ch])@state.coefficients[ch].T@mspline_basis(delta,n_basis=coef.shape[0],half_width=state.config['offset_half_width_angstrom']).T
            ax=axes[2*n,j];ax.plot(delta,reference[0],lw=1,label='Individual Science fit');ax.plot(delta,block[0],lw=1,label='Sky-only block')
            ax.set_title(f'{ch} · fiber {f} · {wave:g} Å');ax.set_ylabel('LSF density (Å⁻¹)');ax.set_xlabel('Offset from input wavelength (Å)')
            ax=axes[2*n+1,j];ax.plot(delta,block[0]-reference[0],lw=.8,color='#a34242');ax.axhline(0,color='#666666',lw=.6)
            ax.set_ylabel('Block − individual fit (Å⁻¹)');ax.set_xlabel('Offset from input wavelength (Å)')
            checks.append(dict(fiberid=f,channel=ch,wavelength_A=wave,integrated_abs_difference=float(np.trapz(abs(block-reference)[0],delta))))
    handles,labels=axes[0,0].get_legend_handles_labels()
    fig.legend(handles,labels,loc='lower center',bbox_to_anchor=(.5,.002),ncol=2,fontsize=10)
    fig.get_layout_engine().set(rect=(0,.035,1,.965))
    fig.suptitle(f'Exposure {expnum}: full Science LSF profiles and signed differences\nMiddle physical fiber with accepted R and Z fits in each spectrograph; profiles retain their centroids',fontsize=15)
    fig.savefig(HERE/f'lsf_profiles_{expnum}.png',dpi=125);plt.close(fig)
    (HERE/f'lsf_profile_sources_{expnum}.json').write_text(json.dumps(dict(parameter_files=sources,comparisons=checks,profiles_recentered=False),indent=2))
    return f'<img src="lsf_profiles_{expnum}.png" alt="Full individual and global LSF profiles with signed differences for exposure {expnum}"><p>Profiles evaluated from the original individual parameter FITS and stored block coefficients, at the same wavelengths and on the same offset grid; centroids are retained. Examples use the middle physical fiber with accepted R and Z fits in each spectrograph. <a href="lsf_profile_sources_{expnum}.json">Parameter sources and integrated profile differences</a>.</p>'
