"""Direct SFrame flux comparisons and inspectable individual-fiber LSF recovery."""
from pathlib import Path
import sys
import json
import html
import os
os.environ.setdefault('MPLBACKEND','Agg')
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from plotly.offline import get_plotlyjs

ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT.parents[1]))
from report_overview import overview

CSS='body{font:17px/1.5 system-ui;color:#182333;margin:30px auto;max-width:1250px;padding:0 20px}h1,h2{line-height:1.2}table{border-collapse:collapse;width:100%;font-size:15px}td,th{text-align:left;padding:8px;border-bottom:1px solid #ccd5df}.card{padding:18px;border:1px solid #ccd5df;border-radius:10px;margin:22px 0}.notice{background:#edf4fa;padding:14px}img{width:100%;height:auto}a{color:#17649d}.plot{margin:20px 0}code{font-size:14px}'
UNIT='Flux (10⁻¹⁴ erg s⁻¹ cm⁻² Å⁻¹)'


def trace(fig,x,y,name,row,col=1,color=None,dash=None,showlegend=True,xlabel="λ",xunit=" Å",width=1.3):
    fig.add_trace(go.Scatter(x=np.asarray(x,dtype=float),y=np.asarray(y,dtype=float),mode='lines',name=name,showlegend=showlegend,
        line=dict(color=color,width=width,dash=dash),connectgaps=False,
        hovertemplate=xlabel+' %{x:.2f}'+xunit+'<br>%{y:.4g}<extra>'+name+'</extra>'),row=row,col=col)


def finish(fig,title,height):
    fig.update_layout(title=title,height=height,template='plotly_white',hovermode='closest',
        legend=dict(orientation='h',yref='container',y=.015,yanchor='bottom',x=0),
        margin=dict(l=90,r=35,t=110,b=210),font=dict(size=14))
    fig.update_xaxes(showspikes=False)
    fig.update_yaxes(showspikes=False)
    return fig.to_html(full_html=False,include_plotlyjs=False,config=dict(responsive=True,displaylogo=False,scrollZoom=False))


def page(title,body):
    return '<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1"><title>'+html.escape(title)+'</title><style>'+CSS+'</style><script src="plotly.min.js"></script></head><body>'+body+'</body></html>'



def main():
    plt.rcParams.update({'font.size':10})
    (ROOT/'plotly.min.js').write_text(get_plotlyjs())
    exposures=json.loads((ROOT/'run_manifest_v2.json').read_text())
    overview(exposures,[],page)



if __name__=='__main__':main()
