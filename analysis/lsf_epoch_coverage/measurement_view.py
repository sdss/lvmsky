"""Publish individual-fiber exposure measurements with observed-lunation boundaries."""
import json
import hashlib
import inspect
import multiprocessing
from concurrent.futures import ProcessPoolExecutor,as_completed
from pathlib import Path

import numpy as np
import pandas as pd

METRICS=['fwhm_angstrom','w50_angstrom','w90_angstrom','centroid_angstrom']


def exposure_positions(exposures):
    result=exposures.sort_values('obstime').copy()
    moons=sorted(result.full_moon_utc.unique())
    result['lunation_number']=result.full_moon_utc.map({moon:i+1 for i,moon in enumerate(moons)})
    groups=result.groupby('full_moon_utc',sort=False)
    result['plot_x']=result.lunation_number+.8*((groups.cumcount()+.5)/groups.expnum.transform('size')-.5)
    assert result.expnum.is_unique
    assert ((result.plot_x-result.lunation_number).abs()<.5).all()
    return result


def publish_measurements(central,exposures,out):
    from plotly.offline import get_plotlyjs
    positions=exposure_positions(exposures)
    data=central.merge(positions[['expnum','plot_x','lunation_number','obstime']],on='expnum',validate='many_to_one')
    assert not data.duplicated(['fiberid','channel','expnum']).any()
    assert np.isfinite(data[METRICS].to_numpy()).all()
    directory=out/'measurement_fibers';directory.mkdir(exist_ok=True)
    identity=data[['fiberid','telescope','spectrographid']].drop_duplicates()
    assert identity.fiberid.is_unique
    fibers=[]
    for fiber,group in data.groupby('fiberid',sort=True):
        meta=group.iloc[0]
        channels={}
        for channel,points in group.groupby('channel',sort=True):
            points=points.sort_values('observation_mjd_utc')
            channels[channel]={column:points[column].tolist() for column in ['expnum','plot_x','lunation_number',
                'obstime','mjd','calibration_epoch',*METRICS]}
        item=dict(fiberid=int(fiber),telescope=str(meta.telescope),spectrographid=int(meta.spectrographid),channels=channels)
        (directory/f'fiber-{int(fiber):04d}.json').write_text(json.dumps(item,separators=(',',':'),allow_nan=False))
        fibers.append(dict(fiberid=int(fiber),telescope=str(meta.telescope),spectrographid=int(meta.spectrographid),
                           exposures=int(group.expnum.nunique()),measurements=len(group)))
    pd.DataFrame(fibers).to_csv(out/'measurement_fiber_inventory.csv',index=False)
    data[['fiberid','telescope','spectrographid','channel','expnum','obstime','mjd','calibration_epoch',
          'full_moon_utc','lunation_number','plot_x','wavelength_angstrom',*METRICS]].to_parquet(
              out/'individual_lsf_measurements.parquet',index=False)
    periods=positions.groupby('lunation_number').agg(full_moon_utc=('full_moon_utc','first'),
        selected_exposures=('expnum','size'),epoch_start_mjd=('nominal_epoch',lambda values:','.join(map(str,sorted(values.unique()))))).reset_index()
    configuration=dict(fibers=fibers,periods=periods.to_dict('records'))
    (out/'lsf_measurements_config.json').write_text(json.dumps(configuration,allow_nan=False))
    (out/'plotly.min.js').write_text(get_plotlyjs())
    (out/'lsf_measurements.js').write_text(VIEW_SCRIPT)
    render_previews(out,configuration)
    return measurement_markup(configuration)


def measurement_markup(configuration):
    periods=pd.DataFrame(configuration['periods'])
    sections=[];navigation=[]
    for telescope,title in [('Sci','Science'),('SkyE','Sky East'),('SkyW','Sky West')]:
        for spectrograph in [1,2,3]:
            group=sorted((fiber for fiber in configuration['fibers'] if fiber['telescope']==telescope and fiber['spectrographid']==spectrograph),key=lambda fiber:fiber['fiberid'])
            identifier=f'measurements-{telescope}-sp{spectrograph}'
            navigation.append(f'<a href="#{identifier}">{title} · SP{spectrograph} ({len(group)} fibers)</a>')
            sections.append(f'<section id="{identifier}"><h3>{title} · spectrograph {spectrograph}: individual measurements</h3>')
            if group:
                default=max(group,key=lambda fiber:fiber['exposures'])
                options=''.join(f'<option value="{fiber["fiberid"]}"'+(' selected' if fiber['fiberid']==default['fiberid'] else '')+f'>{fiber["fiberid"]} ({fiber["exposures"]} exposures)</option>' for fiber in group)
                fiber=default
                number=fiber['fiberid']
                sections.append(f'<figure class="measurement-figure" data-telescope="{telescope}" data-spectrograph="{spectrograph}" data-fiber="{number}">'
                    f'<p><label>Physical fiber <select class="measurement-fiber-control">{options}</select></label></p>'
                    f'<figcaption class="measurement-caption"><b>{title} · SP{spectrograph} · physical fiber {number}</b>: '
                    f'{fiber["exposures"]} accepted exposures; {fiber["measurements"]} channel measurements.</figcaption>'
                    f'<div class="measurement-canvas"><img loading="lazy" width="1400" height="1000" src="measurement_previews/fiber-{number:04d}.png" '
                    f'alt="Physical fiber {number}: individual FWHM, W50, W90 and centroid values across all observed lunations">'
                    '<div class="measurement-interactive"></div></div><p class="measurement-status" role="status"></p></figure>')
            sections.append('</section>')
    return f'''<h2 id="individual-measurements">Individual LSF measurements: all telescope and spectrograph groups</h2>
<p>All nine telescope/spectrograph groups have separate, always-visible plots. Each group's control chooses only its physical fiber and changes only that plot.
The initial fiber has the highest accepted exposure count in its group. All four metrics show <b>individual recovered values in Å</b>, with B/R/Z columns.
The X axis numbers the 20 <b>observed</b> lunations chronologically; positions within each block follow order among all Moon/Sun-selected exposures. A fiber's accepted points may leave gaps.
Blocks have equal width, so this is an ordered sampling view, not an elapsed-time axis. Gray lines separate lunations; darker lines mark changes of calibration epoch.
Static plots are always available; the same plot becomes interactive automatically as it enters the viewport. Point hover gives exposure ID, exact UTC time, epoch and the measured value. Drag to zoom; double-click restores the full range.</p>
<p>Each measurement panel uses its own Y-axis range; compare numerical axis values when comparing channels or fibers.
<a href="#lunation-01">Go to per-lunation coverage and scatter plots</a>.</p>
<p>Points include accepted, finite diagnostics; each figure contains one physical fiber. The plots retain every value and use linear Å axes.
Median reference values and uncertainties are not substituted for individual measurements.
<a href="individual_lsf_measurements.parquet">All plotted measurements</a> · <a href="measurement_fiber_inventory.csv">Fiber inventory</a></p>
<p><a href="measurement_previews/">Static images for every physical fiber</a>.</p>
<div class="scroll">{periods.rename(columns={'lunation_number':'Observed lunation number','full_moon_utc':'Nearest full Moon (UTC)',
    'selected_exposures':'Moon/Sun-selected exposures before Gaia','epoch_start_mjd':'Calibration epoch start MJD'}).to_html(index=False)}</div>
<p>{' · '.join(navigation)}</p>
<style>.measurement-figure{{margin:35px 0}}.measurement-canvas{{position:relative}}.measurement-canvas img{{display:block;width:100%;height:auto}}.measurement-interactive{{position:absolute;inset:0}}.measurement-status{{font-size:13px}}</style>
{''.join(sections)}
<script src="plotly.min.js"></script><script src="lsf_measurements.js"></script>'''


def render_preview(task):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    filename,periods,output=task
    data=json.loads(Path(filename).read_text())
    fig,axes=plt.subplots(4,3,figsize=(14,10),sharex=True)
    for row,(metric,label) in enumerate(zip(METRICS,['FWHM (Å)','W50 (Å)','W90 (Å)','Centroid offset (Å)'])):
        for col,(channel,color) in enumerate(zip(['B','R','Z'],['#0072B2','#D55E00','#009E73'])):
            ax=axes[row,col];points=data['channels'].get(channel)
            if points:ax.plot(points['plot_x'],points[metric],'.',markersize=3,color=color)
            else:ax.text(.5,.5,'No accepted channel measurements',ha='center',transform=ax.transAxes,fontsize=9)
            for index in range(1,len(periods)):
                changed=periods[index]['epoch_start_mjd']!=periods[index-1]['epoch_start_mjd']
                ax.axvline(index+.5,color='#444' if changed else '#ccc',linewidth=1 if changed else .4,zorder=0)
            ax.set_xlim(.5,len(periods)+.5);ax.set_xticks(range(1,len(periods)+1));ax.tick_params(labelsize=8)
            ax.set_ylabel(label);ax.grid(axis='y',alpha=.2)
            if row==0:ax.set_title(f'{channel}: {[4700,6600,8500][col]} Å')
            if row==3:ax.set_xlabel('Observed lunation number')
    fig.suptitle(f'{data["telescope"]} · spectrograph {data["spectrographid"]} · physical fiber {data["fiberid"]}: individual measurements in Å')
    fig.text(.5,.01,'Every accepted value is shown. Equal-width lunation blocks order exposures; exact timestamps are retained in the data and interactive hover.',ha='center',fontsize=9)
    fig.tight_layout(rect=[0,.03,1,.96])
    target=Path(output);temporary=target.with_suffix('.tmp')
    fig.savefig(temporary,format='png',dpi=100);plt.close(fig);temporary.replace(target)
    return data['fiberid']


def render_previews(out,configuration):
    directory=out/'measurement_previews';directory.mkdir(exist_ok=True)
    manifest=out/'measurement_preview_manifest.json'
    previous=json.loads(manifest.read_text()) if manifest.exists() else {}
    code=inspect.getsource(render_preview)+json.dumps(configuration['periods'])
    hashes={};tasks=[]
    for fiber in configuration['fibers']:
        number=fiber['fiberid'];source=out/'measurement_fibers'/f'fiber-{number:04d}.json'
        digest=hashlib.sha256(code.encode()+source.read_bytes()).hexdigest();hashes[str(number)]=digest
        target=directory/f'fiber-{number:04d}.png'
        if previous.get(str(number))!=digest or not target.exists():tasks.append((str(source),configuration['periods'],str(target)))
    if tasks:
        print(f'Rendering {len(tasks)} always-visible individual-fiber plots',flush=True)
        with ProcessPoolExecutor(max_workers=16,mp_context=multiprocessing.get_context('fork')) as pool:
            futures=[pool.submit(render_preview,task) for task in tasks]
            for index,future in enumerate(as_completed(futures),1):
                future.result()
                if index%100==0 or index==len(tasks):print(f'Individual-fiber plots: {index}/{len(tasks)}',flush=True)
    assert all((directory/f'fiber-{fiber["fiberid"]:04d}.png').exists() for fiber in configuration['fibers'])
    temporary=manifest.with_suffix('.tmp');temporary.write_text(json.dumps(hashes));temporary.replace(manifest)


VIEW_SCRIPT=r'''(() => {
  const base = new URL('.', document.currentScript.src);
  const figures = Array.from(document.querySelectorAll('.measurement-figure'));
  const channels = ['B', 'R', 'Z'];
  const colors = ['#0072B2', '#D55E00', '#009E73'];
  const metrics = ['fwhm_angstrom', 'w50_angstrom', 'w90_angstrom', 'centroid_angstrom'];
  const labels = ['FWHM (Å)', 'W50 (Å)', 'W90 (Å)', 'Centroid offset (Å)'];
  let config;
  const states = new WeakMap();
  async function read(name) {
    const response = await fetch(new URL(name, base));
    if (!response.ok) throw new Error(`HTTP ${response.status}: ${name}`);
    return response.json();
  }
  async function draw(figure) {
    const state=states.get(figure), selected=Number(figure.dataset.fiber);
    const plot=figure.querySelector('.measurement-interactive'), status=figure.querySelector('.measurement-status');
    state.loading=true;
    try {
      const data = await read(`measurement_fibers/fiber-${String(selected).padStart(4, '0')}.json`);
      if (!state.visible||Number(figure.dataset.fiber)!==selected) return;
      const traces = [], shapes = [], annotations = [];
      const layout = {height:plot.clientHeight||800, margin:{l:70,r:25,t:60,b:75}, hovermode:'closest', showlegend:false,
        title:{text:`${data.telescope} · spectrograph ${data.spectrographid} · physical fiber ${selected}`},
        shapes, annotations};
      for (let row=0; row<4; row++) for (let col=0; col<3; col++) {
        const channel=channels[col], points=data.channels[channel], n=row*3+col+1;
        const xref=n===1?'x':`x${n}`, yref=n===1?'y':`y${n}`;
        const xkey=n===1?'xaxis':`xaxis${n}`, ykey=n===1?'yaxis':`yaxis${n}`;
        const bottom=1-(row+1)/4+.035, top=1-row/4-.015;
        layout[xkey]={domain:[col/3+.035,(col+1)/3-.025],anchor:yref,range:[.5,config.periods.length+.5],
          tickmode:'array',tickvals:config.periods.map(p=>p.lunation_number),
          ticktext:config.periods.map(p=>String(p.lunation_number)),tickfont:{size:10},
          showticklabels:row===3,showspikes:false,matches:n===1?undefined:'x',
          title:row===3?{text:'Observed lunation number',font:{size:12}}:undefined};
        layout[ykey]={domain:[bottom,top],anchor:xref,type:'linear',showspikes:false,
          title:{text:labels[row],font:{size:12}},automargin:true};
        if (points) traces.push({type:'scatter',mode:'markers',x:points.plot_x,y:points[metrics[row]],xaxis:xref,yaxis:yref,
          marker:{size:4,color:colors[col]},customdata:points.expnum.map((id,i)=>[id,points.obstime[i],points.calibration_epoch[i],points.lunation_number[i]]),
          hovertemplate:'%{y:.6f} Å<br>exp %{customdata[0]} · L%{customdata[3]}<br>%{customdata[1]} UTC<br>epoch %{customdata[2]}<extra></extra>'});
        for (let i=1;i<config.periods.length;i++) {
          const changed=config.periods[i].epoch_start_mjd!==config.periods[i-1].epoch_start_mjd;
          shapes.push({type:'line',xref,yref:`${yref} domain`,x0:i+.5,x1:i+.5,y0:0,y1:1,
            line:{color:changed?'#444':'#ccc',width:changed?1.5:.5},layer:'below'});
        }
        if (row===0) annotations.push({xref:'paper',yref:'paper',x:(col+.5)/3,y:1.045,text:`${channel}: ${[4700,6600,8500][col]} Å`,showarrow:false});
      }
      await Plotly.react(plot,traces,layout,{responsive:true,displaylogo:false});
      if (!state.visible||Number(figure.dataset.fiber)!==selected) { Plotly.purge(plot);return; }
      state.rendered=true;
      const all=Object.values(data.channels), exposures=new Set(all.flatMap(p=>p.expnum)), periods=new Set(all.flatMap(p=>p.lunation_number));
      status.textContent=`Fiber ${selected}: ${exposures.size} accepted exposures in ${periods.size} observed lunations; ${all.reduce((n,p)=>n+p.expnum.length,0)} channel measurements.`;
    } catch(error) { Plotly.purge(plot);status.textContent=`Static plot available; interactive measurements could not be loaded: ${error.message}`; }
    finally { state.loading=false;if(state.visible&&Number(figure.dataset.fiber)!==selected)draw(figure); }
  }
  read('lsf_measurements_config.json').then(value=>{
    config=value;
    const observer=new IntersectionObserver(entries=>{
      for (const entry of entries) {
        const figure=entry.target,state=states.get(figure);state.visible=entry.isIntersecting;
        if (state.visible&&!state.rendered&&!state.loading) draw(figure);
        else if (!state.visible&&state.rendered) { Plotly.purge(figure.querySelector('.measurement-interactive'));state.rendered=false; }
      }
    },{rootMargin:'100px'});
    for (const figure of figures) {
      states.set(figure,{visible:false,loading:false,rendered:false});
      figure.querySelector('.measurement-fiber-control').addEventListener('change',event=>{
        const selected=Number(event.target.value),state=states.get(figure);
        const meta=config.fibers.find(item=>item.fiberid===selected);
        if(!meta||meta.telescope!==figure.dataset.telescope||meta.spectrographid!==Number(figure.dataset.spectrograph))return;
        figure.dataset.fiber=String(selected);
        figure.querySelector('img').src=new URL(`measurement_previews/fiber-${String(selected).padStart(4,'0')}.png`,base).href;
        figure.querySelector('img').alt=`Physical fiber ${selected}: individual measurements across observed lunations`;
        figure.querySelector('.measurement-caption').textContent=`${meta.telescope} · SP${meta.spectrographid} · physical fiber ${selected}: ${meta.exposures} accepted exposures; ${meta.measurements} channel measurements.`;
        Plotly.purge(figure.querySelector('.measurement-interactive'));state.rendered=false;
        if(state.visible&&!state.loading)draw(figure);
      });
      observer.observe(figure);
    }
  }).catch(error=>{for(const figure of figures) figure.querySelector('.measurement-status').textContent=`Static plot available; ${error.message}`;});
})();'''
