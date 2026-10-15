"""Survey-LSF counts and recorded primary reasons from the current FITS QA."""
from collections import Counter
from html import escape
import numpy as np
from astropy.io import fits

REASONS={'line_information':'Line coverage / contrast / S/N',
         'profile_guard':'Profile bounds / departure from reference',
         'not_converged':'Convergence', 'insufficient_lines':'Numerical line/basis support',
         'continuum_solver':'Continuum solver', 'survey_only':'Survey-only run'}


def summarize_fallbacks(h,expnum):
    slit=h['SLITMAP'].data;qa=h['LSF_QA'].data
    science=np.flatnonzero(np.char.strip(slit['telescope'].astype(str))=='Sci')
    processed=np.intersect1d(science,h['LSF_FIBERS'].data['input_row'])
    output=dict(expnum=int(expnum),science=len(science),processed=len(processed),channels={})
    survey_sets=[]
    for channel in 'RZ':
        rows=qa[(qa['channel']==channel)&np.isin(qa['input_row'],processed)]
        if len(rows)!=len(processed) or set(rows['input_row'])!=set(processed):
            raise ValueError('Expected one LSF_QA record per processed Science fiber/channel')
        invalid=rows['reason']=='unusable_survey_fallback'
        accepted=rows['accepted'].astype(bool)
        if np.any(invalid&accepted):raise ValueError('An invalid channel cannot have an adopted in-situ LSF')
        survey=~accepted&~invalid
        reasons=dict(Counter(str(x) for x in rows['reason'][survey]))
        output['channels'][channel]=dict(insitu=int(accepted.sum()),survey=int(survey.sum()),
            invalid=int(invalid.sum()),reasons=reasons)
        survey_sets.append(set(map(int,rows['input_row'][survey])))
    output['survey_either_channel']=len(survey_sets[0]|survey_sets[1])
    return output


def fallback_summaries(exposures):
    rows=[]
    for e in exposures:
        with fits.open(e['output'],memmap=True) as h:rows.append(summarize_fallbacks(h,e['expnum']))
    return rows


def fallback_table(rows):
    body='<section id="survey-lsf-summary"><h2>Science-fiber LSF choices</h2><p>Counts from current FITS LSF_QA. Percentages use processed Science fibers. B uses survey by design. Survey reasons are recorded primary outcomes; unlisted reasons have zero counts. <a href="#lsf-selection">Thresholds</a>.</p>'
    for channel in 'RZ':
        body+=f'<h3>{channel}-channel</h3><table><tr><th>Exposure</th><th>Processed / all Science</th><th>In situ · fibers (%)</th><th>Survey · fibers (%)</th><th>Reasons for survey · fibers</th></tr>'
        for row in rows:
            total=row['processed'];data=row['channels'][channel]
            counts=[f'{data[key]} ({data[key]/total:.1%})' if total else '0 (—)' for key in ('insitu','survey')]
            reasons=[]
            for code,count in sorted(data['reasons'].items()):
                label=REASONS.get(code,'LSF solver: '+code[len('lsf_solver_'):] if code.startswith('lsf_solver_') else code.replace('_',' '))
                reasons.append(f'<div>{escape(label)}: <b>{count}</b></div>')
            invalid=f'<br>Invalid channel (not survey): {data["invalid"]}' if data['invalid'] else ''
            body+=f'<tr><td><a href="#exp{row["expnum"]}">{row["expnum"]}</a></td><td>{total} / {row["science"]}{invalid}</td><td>{counts[0]}</td><td>{counts[1]}</td><td>{"".join(reasons) or "0"}</td></tr>'
        body+='</table>'
    return body+'<p><a href="lsf_fallback_summary.json">Exact counts, reason codes and unique fibers using survey in R or Z</a></p></section>'
