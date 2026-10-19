"""Reconstruct block LSFs and individual Sky-fiber spectral fits in memory."""
import numpy as np
from skysub.sky_subtraction.reconstruction import SkyMLReconstructor
from skysub.sky_subtraction.pipeline import render,FACTOR
from .lsf import evaluate_block


class SkyBlockReconstructor(SkyMLReconstructor):
    def __init__(self,path,cframe=None):
        super().__init__(path,cframe)
        self.blocks={}
        for i,row in enumerate(self.h['BLOCK_LSF_META'].data):
            self.blocks[(int(row['spectrograph']),str(row['channel']))]=dict(
                coef=np.asarray(self.h['BLOCK_LSF_COEF'].data[i],float),
                slit_knots=np.asarray(self.h['BLOCK_SLIT_KNOTS'].data[i],float),
                wave_knots=np.asarray(self.h['BLOCK_WAVE_KNOTS'].data[i],float),
                wave_degree=int(row['wave_degree']),qa=dict(status=str(row['status']),centered=bool(row['centered']) if 'centered' in row.array.names else True))

    def global_state(self,fiberid):
        """R/Z block surfaces at this physical fiber; B is the Science fiber's survey LSF or the Sky telescope median."""
        data=self.h['SLITMAP'].data
        row=np.flatnonzero(data['fiberid']==fiberid)
        if len(row)!=1:raise ValueError('Unknown physical fiber ID')
        row=int(row[0]);science=np.flatnonzero(self.h['LSF_FIBERS'].data['input_row']==row)
        if len(science):base=self.state('',int(science[0]))
        else:
            label=str(data['telescope'][row]).strip()
            if label not in self.labels or label=='Sci':raise ValueError('Fiber has no stored B LSF')
            base=self.state('MED_',self.labels.index(label))
        return evaluate_block(base,self.blocks,data['spectrographid'][row],fiberid)

    def sky_fiber(self,fiberid):
        rows=self.h['SKY_FIBER_COEF'].data
        rows=rows[rows['fiberid']==fiberid]
        if len(rows)!=1:raise ValueError('Fiber was not selected for the Sky-only fit')
        row=rows[0];i=int(row['input_row']);label=str(row['telescope'])
        total,components,_=render(self.model(label),row['coef'],self.global_state(fiberid))
        valid=(self.mask[i]==0)&(self.ivar[i]>0)&np.isfinite(self.observed[i])&~self.h['MED_FIT_MASK'].data.astype(bool)
        sigma=np.divide(1,np.sqrt(np.maximum(self.ivar[i],0)),out=np.full(len(self.wave),np.nan),where=self.ivar[i]>0)
        return dict(observed=self.observed[i],model=total/FACTOR,components={k:v/FACTOR for k,v in components.items()},
            sigma=sigma,valid=valid,input_row=i,telescope=label)
