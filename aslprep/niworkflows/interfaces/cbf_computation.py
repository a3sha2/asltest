#from ...pybids import BIDSLayout
import os
import numpy as np
import pandas as pd
import nibabel as nb


from nipype import logging
from nipype.utils.filemanip import fname_presuffix,split_filename,copyfiles
from nipype.interfaces.base import (
    traits, TraitedSpec, BaseInterfaceInputSpec, SimpleInterface,
    File, InputMultiPath, OutputMultiPath, isdefined,Undefined)
from nipype.interfaces.fsl.base import (FSLCommand, FSLCommandInputSpec, Info)
from nipype.interfaces import fsl

LOGGER = logging.getLogger('nipype.interface')


class _extractCBFInputSpec(BaseInterfaceInputSpec):
    in_file = File(exists=True, mandatory=True,
                              desc='preprocessed file')
    in_ASLcontext = File(exists=True, mandatory=True,
              desc='ASL conext text tsv file with label and control')
    out_file=File(exists=False,mandatory=True,desc='cbf timeries data')
    out_avg=File(exists=False,mandatory=True,desc='average control')

class _extractCBFOutputSpec(TraitedSpec):
    out_file = File(exists=True, desc='cbf timeries data')
    out_avg = File(exists=True, desc='average control')


class extractCBF(SimpleInterface):
    """
    extract  CBF timeseries
    by substracting label from control

    """

    input_spec = _extractCBFInputSpec
    output_spec =_extractCBFOutputSpec

    def _run_interface(self, runtime):
         
        aslcontext=pd.read_csv(self.inputs.in_ASLcontext,header=None)
        idasl=aslcontext[0].tolist()
        controllist= [ i for i in range(len(idasl)) if idasl[i] == 'Control' ]
        labellist=[ i for i in range(len(idasl)) if idasl[i] == 'Label' ]
        
        
        # read the nifti image 
        allasl=nb.load(self.inputs.in_file)
        dataasl=allasl.get_fdata()
        if len(dataasl.shape) == 5:
                raise RuntimeError('Input image (%s) is 5D.')
        control_img=dataasl[:,:,:,controllist]
        label_img=dataasl[:,:,:,labellist]
        cbf_data=np.subtract(control_img,label_img)
        avg_control=np.mean(control_img,axis=3)

        self._results['out_file'] = fname_presuffix(self.inputs.in_file,
                                                   suffix='_cbftimeseries', newpath=runtime.cwd)
        self._results['out_avg'] = fname_presuffix(self.inputs.in_file,
                                                   suffix='_avg_control', newpath=runtime.cwd)
        nb.Nifti1Image(
            cbf_data, allasl.affine, allasl.header).to_filename(
            self._results['out_file'])
        nb.Nifti1Image(
            avg_control, allasl.affine, allasl.header).to_filename(
            self._results['out_avg'])

        return runtime




class _computeCBFInputSpec(BaseInterfaceInputSpec):
    #in_file = File(exists=True, mandatory=True,
                              #desc='asl raw')
    in_cbf = File(exists=True,mandatory=True,desc= 'cbf nifti')
    in_metadata = traits.Dict(exists=True, mandatory=True,
              desc='metadata for CBF ')
    in_m0file = File(exists=True, mandatory=False,
              desc='M0 nifti file')
    out_cbf=File(exists=False,mandatory=True,desc='cbf timeries data')
    out_mean=File(exists=False,mandatory=True,desc='average control')
    out_att=File(exists=False,mandatory=False,desc='Arterial Transit Time')
    

class _computeCBFOutputSpec(TraitedSpec):
    out_cbf=File(exists=False,mandatory=True,desc='cbf timeries data')
    out_mean=File(exists=False,mandatory=True,desc='average control')
    out_att=File(exists=False,mandatory=False,desc='Arterial Transit Time')



class computeCBF(SimpleInterface):

    """
    compute cbf pASL or pCASL

    """

    input_spec = _computeCBFInputSpec
    output_spec =_computeCBFOutputSpec
    
    def _run_interface(self, runtime):

        labeltype=self.inputs.in_metadata['LabelingType']
        tau=self.inputs.in_metadata['LabelingDuration']
        plds=np.array(self.inputs.in_metadata['InitialPostLabelDelay'])
        m0scale=self.inputs.in_metadata['M0']
        magstrength=self.inputs.in_metadata['MagneticFieldStrength']
        t1blood=110*int(magstrength[:-1])+1316
        inverstiontime=np.add(tau,plds)
        if self.inputs.in_metadata['LabelingEfficiency']:
            labeleff=self.inputs.in_metadata['LabelingEfficiency']
        elif 'CASL' in labeltype: 
             labeleff=0.72
        elif 'PASL' in labeltype:
            labeleff=0.8
        else:
            print( 'no labelelling effiecieny')

        part_coeff=0.9 # brain partition coefficient
    
        if 'CASL' in  labeltype:
            pf1=(6000*part_coeff)/(2*labeleff*t1blood*(1-np.exp(-(tau/t1blood))))
            perfusion_factor=pf1*np.exp(plds/t1blood)
        elif 'PASL' in labeltype: 
            pf1=(6000*part_coeff)/(2*labeleff)
            perfusion_factor=(pf1*np.exp(inverstiontime/t1blood))/inverstiontime
        perfusion_factor=np.array([perfusion_factor])
        # get control  now 
        avg_control=[]
        
        mzero=nb.load(self.inputs.in_m0file).get_fdata()
        if len(mzero.shape) > 3:
            avg_control=np.mean(mzero,axis=3)
        else:
            avg_control=mzero
        if not m0scale:
            m0scale=1


        cbf_data=nb.load(self.inputs.in_cbf).get_fdata()
        m0_data=np.zeros(cbf_data.shape)
        for i in range(cbf_data.shape[3]):
            m0_data[:,:,:,i]=avg_control
        
        cbf1=(cbf_data/m0_data)/m0scale

        # for compute cbf for each PLD and TI 
        att=None
        if len(perfusion_factor) > 1: 
            cbf_data_ts=np.zeros(np.concatenate(cbf.shape,len(perfusion_factor)))
            dm1factor=(2*labeleff*1.5*(1-np.exp(tau/1.5)))*avg_control
            deltaM=np.zeros(np.concatenate(avg_control.shape,len(perfusion_factor)))
            for i in range(len(perfusion_factor)):
                cbf_data_ts[:,:,:,:,i]=cbf1*perfusion_factor[i]
                deltaM[:,:,:,i]=dm1factor*(np.exp(-plds[i]/t1blood))
            cbf=np.mean(cbf_data_ts,axis=4)
            # compute  arterial transisttime
            deltaM2=np.zeros(np.concatenate(avg_control.shape,len(perfusion_factor)))
            for i in range(len(perfusion_factor)):
                deltaM2[:,:,:,i]=deltaM[:,:,:,i]*plds[i]
            att=np.sum(deltaM2,axis=4)/np.sum(deltaM,axis=4)
        else:
            cbf=cbf1*perfusion_factor
        ## cbf is timeseries
        meancbf=np.mean(cbf,axis=3)
        self._results['out_cbf'] = fname_presuffix(self.inputs.in_file,
                                                   suffix='_cbf', newpath=runtime.cwd)
        self._results['out_mean'] = fname_presuffix(self.inputs.in_file,
                                                   suffix='_meancbf', newpath=runtime.cwd)
        samplecbf=nb.load(self.inputs.in_cbf)
        nb.Nifti1Image(
            cbf, samplecbf.affine, samplecbf.header).to_filename(
            self._results['out_cbf'])
        nb.Nifti1Image(
            meancbf, samplecbf.affine, samplecbf.header).to_filename(
            self._results['out_mean'])
        if att is not None:
            self._results['out_att'] = fname_presuffix(self.inputs.in_file,
                                                   suffix='_att', newpath=runtime.cwd)
            nb.Nifti1Image(
            att, samplecbf.affine, samplecbf.header).to_filename(
            self._results['out_att'])
        
        return runtime


#score and scrub 
class _scorescrubCBFInputSpec(BaseInterfaceInputSpec):
    in_file = File(exists=True, mandatory=True,
                              desc='computed CBF from computeCBF')
    in_greyM = File(exists=True, mandatory=True,desc='grey  matter')
    in_whiteM = File(exists=True, mandatory=True,desc='white  matter')
    in_mask = File(exists=True, mandatory=True,desc='mask')
    in_csf = File(exists=True, mandatory=True,desc='csf')
    in_thresh=traits.Float(default_value=0.7,exists=True,mandatory=False,desc='threshold of propbaility matter')
    in_wfun=traits.Str(exists=True,mandatory=False,default_value='huber',
              option=['bisquare','andrews','cauchy','fair','logistics','ols','talwar','welsch'],
               desc='wavelet fun ')
    out_score=File(exists=False,mandatory=True,desc='score timeseries data')
    out_avgscore=File(exists=False,mandatory=True,desc='average score')
    out_scrub=File(exists=False,mandatory=True,desc='average scrub')
    out_scoreindex=File(exists=False,mandatory=True,desc='index of volume remove or leave by score')

class _scorescrubCBFOutputSpec(TraitedSpec):
    out_score=File(exists=True,mandatory=True,desc='score timeseries data')
    out_avgscore=File(exists=True,mandatory=True,desc='average score')
    out_scrub=File(exists=True,mandatory=True,desc='average scrub')
    out_scoreindex=File(exists=True,mandatory=True,desc='index of volume remove or leave by score')

class scorescrubCBF(SimpleInterface):

    """
    compute score and scrub 

    """
    input_spec = _scorescrubCBFInputSpec
    output_spec =_scorescrubCBFOutputSpec

    def _run_interface(self, runtime):
        cbf_ts=nb.load(self.inputs.in_file).get_fdata()
        mask=nb.load(self.inputs.in_mask).get_fdata()
        greym=nb.load(self.inputs.in_greyM).get_fdata()
        whitem=nb.load(self.inputs.in_whiteM).get_fdata()
        csf=nb.load(self.inputs.in_csf).get_fdata()
        cbf_scorets,index_score=_getcbfscore(cbfts=cbf_ts,wm=whitem,gm=greym,csf=csf,
                       thresh=self.inputs.in_thresh)
        cbfscrub=_scrubcbf(cbf_ts=cbf_scorets,gm=greym,wm=whitem,csf=csf,mask=mask,
                          wfun=self.inputs.in_wfun,thresh=self.inputs.in_thresh)
        
        self._results['out_score'] = fname_presuffix(self.inputs.out_file,
                                                   suffix='_cbfscorets', newpath=runtime.cwd)
        self._results['out_avgscore'] = fname_presuffix(self.inputs.out_file,
                                                   suffix='_meancbfscore', newpath=runtime.cwd)
        self._results['out_scrub'] = fname_presuffix(self.inputs.out_file,
                                                   suffix='_cbfscrub', newpath=runtime.cwd)
        self._results['out_scoreindex'] =runtime.cwd+'/scorescrub.txt'
                                                   
        samplecbf=nb.load(self.inputs.in_mask)
        nb.Nifti1Image(
            cbf_scorets, samplecbf.affine, samplecbf.header).to_filename(
            self._results['out_score'])
        nb.Nifti1Image(
            np.mean(cbf_scorets,axis=3), samplecbf.affine, samplecbf.header).to_filename(
            self._results['out_avgscore'])
        nb.Nifti1Image(
            cbfscrub, samplecbf.affine, samplecbf.header).to_filename(
            self._results['out_scrub'])
        np.savetxt(self._results['out_scoreindex'],index_score, delimiter=',')
        return runtime

def _weightfun(x,wfun='huber'):
    """"
    get weight fun and tuner

    """
    if wfun == 'andrews':
        tuner=1.339
        weight=(np.abs(x)<np.pi)*np.sin(x)
    elif wfun== 'bisquare':
        tuner=4.685
        weight=(np.abs(x)<1)*np.power((1-np.power(x,2)),2)
    elif wfun == 'cauchy':
        tuner=2.385
        weight=1/(1+np.power(x,2))
    elif wfun == 'logistic':
        tuner=1.205
        weight == np.tanh(x)/x
    elif wfun == 'ols':
        tuner=1
        weight=np.repeat(1,len(x))
    elif wfun== 'talwar':
        tuner=2.795
        weight=1*(np.abs(x)<1)
    elif wfun == 'welsch':
        tuner=2.985
        weight=np.exp(-(np.power(x,2)))
    else:
        tuner=1.345
        weight=1/np.abs(x)
    return weight,tuner

def _tune(wfun='huber'):
    """"
    get weight fun and tuner

    """
    if wfun == 'andrews':
        tuner=1.339
    elif wfun== 'bisquare':
        tuner=4.685
    elif wfun == 'cauchy':
        tuner=2.385
    elif wfun == 'logistic':
        tuner=1.205
    elif wfun == 'ols':
        tuner=1
    elif wfun== 'talwar':
        tuner=2.795
    elif wfun == 'welsch':
        tuner=2.985
    else:
        tuner=1.345
    return tuner


def _getchisquare(n):
    a=[0.000000, 15.484663, 8.886835, 7.224733, 5.901333, 5.126189, 4.683238, 4.272937, 4.079918, 
      3.731612, 3.515615, 3.459711, 3.280471, 3.078046, 3.037280, 2.990761, 2.837119, 2.795526, 2.785189, 
      2.649955, 2.637642, 2.532700, 2.505253, 2.469810, 2.496135, 2.342210, 2.384975, 2.275019, 2.244482, 
      2.249109, 2.271968, 2.210340, 2.179537, 2.133762, 2.174928, 2.150072, 2.142526, 2.071512, 2.091061, 
      2.039329, 2.053183, 2.066396, 1.998564, 1.993568, 1.991905, 1.981837, 1.950225, 1.938580, 1.937753, 
      1.882911, 1.892665, 1.960767, 1.915530, 1.847124, 1.947374, 1.872383, 1.852023, 1.861169, 1.843109, 
      1.823870, 1.809643, 1.815038, 1.848064, 1.791687, 1.768343, 1.778231, 1.779046, 1.759597, 1.774383, 
      1.774876, 1.751232, 1.755293, 1.757028, 1.751388, 1.739384, 1.716395, 1.730631, 1.718389, 1.693839, 
      1.696862, 1.691245, 1.682541, 1.702515, 1.700991, 1.674607, 1.669986, 1.688864, 1.653713, 1.641309, 
      1.648462, 1.630380, 1.634156, 1.660821, 1.625298, 1.643779, 1.631554, 1.643987, 1.624604, 1.606314, 
      1.609462]
    b=[0, 2.177715, 1.446966, 1.272340, 1.190646, 1.151953, 1.122953, 1.103451, 1.089395, 1.079783, 
      1.071751, 1.063096, 1.058524, 1.054137, 1.049783, 1.046265, 1.043192, 1.039536, 1.038500, 1.037296, 
      1.033765, 1.032317, 1.031334, 1.029551, 1.028829, 1.027734, 1.024896, 1.024860, 1.025207, 1.024154, 
      1.022032, 1.021962, 1.021514, 1.020388, 1.019238, 1.020381, 1.019068, 1.018729, 1.018395, 1.017134, 
      1.016539, 1.015676, 1.015641, 1.015398, 1.015481, 1.015566, 1.014620, 1.014342, 1.013901, 1.013867, 
      1.013838, 1.013602, 1.013322, 1.012083, 1.013168, 1.012667, 1.011087, 1.011959, 1.011670, 1.011494, 
      1.010463, 1.010269, 1.010393, 1.010004, 1.010775, 1.009399, 1.011000, 1.010364, 1.009831, 1.009563, 
      1.010085, 1.009149, 1.008444, 1.009455, 1.009705, 1.008597, 1.008644, 1.008051, 1.008085, 1.008550, 
      1.008265, 1.009141, 1.008235, 1.008002, 1.008007, 1.007660, 1.007993, 1.007184, 1.008093, 1.007816, 
      1.007770, 1.007932, 1.007819, 1.007063, 1.006712, 1.006752, 1.006703, 1.006650, 1.006743, 1.007087]
    return a[n],b[n]

def _getcbfscore(cbfts,wm,gm,csf,thresh=0.7):
    gm[gm<thresh]=0; gm[gm>0]=1
    wm[wm<thresh]=0; wm[wm>0]=1
    csf[csf<thresh]=0;csf[csf>0]=1

    # get the total number of voxle within csf,gm and wm 
    nogm =np.sum(gm==1)-1;  nowm=np.sum(wm==1)-1;  nocf=np.sum(csf==0)-1  
    mask=gm+wm+csf  
    #msk=sum(mask>0)

    # mean  of times series cbf within greymatter
    mgmts=np.squeeze(np.mean(cbfts[gm==1,:],axis=0))
    # robiust mean and meadian
    medmngm = np.median(mgmts); sdmngm=np.mean(np.abs(mgmts - np.mean(mgmts)))/0.675
    indx=1*(np.abs(mgmts-medmngm)>(2.5*sdmngm))
    R=np.mean(cbfts[:,:,:,indx==0],axis=3)
    V=nogm*np.var(R[gm==1]) + nowm*np.var(R[wm==1]) + nocf*np.var(R[csf==1])
    V1=V+1
    while V < V1:
        V1=V;CC =(-2*np.repeat(1,cbfts.shape[3]))*1
        for s in range(cbfts.shape[3]):
              if indx[s] != 0 :
                  break 
              else:
                  tmp1 = cbfts[:,:,:,s]
                  CC[s]=np.corrcoef(R[mask>0],tmp1[mask>0])[0][1]
        inx=np.argmax(CC); indx[inx]=2
        R=np.mean(cbfts[:,:,:,indx==0],axis=3)
        V=nogm*np.var(R[gm==1]) + nowm*np.var(R[wm==1]) + nocf*np.var(R[csf==1])
    cbfts_recon=cbfts[:,:,:,indx==0]

    return cbfts_recon,indx


def _roubustfit(Y,mu,Globalprior,modrobprior,lmd=0,localprior=0,wfun='huber',tune=1.345,flagstd=1,flagmodrobust=1,flagprior=1,thresh=0.7):
    dimcbf=Y.shape
    priow=np.ones([dimcbf[0],dimcbf[1]]);sw=1
    X=priow
    b=(np.sum(X*Y,axis=0)+mu*Globalprior+lmd*localprior)/(np.sum(X*X,axis=0)+mu+lmd)
    b0=np.repeat(0,len(b))
    h1=X/np.power(np.tile(np.sqrt(np.sum(X*X,axis=0)),(dimcbf[0],1)),2)
    h0=0.9999*np.ones([dimcbf[0],dimcbf[1]])
    h=np.minimum(h0,h1)
    adjfactor=1/(np.sqrt(1-h/priow))
    tiny_s=(1e-6)*(np.std(h,axis=0));tiny_s[tiny_s==0]=1
    D=np.sqrt(np.finfo(float).eps)
    iter =0; interlim=100
  
    while iter<interlim:
        print('iteration  ', iter,"\n")
        iter=iter + 1
        if np.sum(np.abs(b-b0)) < np.sum(D*np.minimum(np.abs(b),np.abs(b0))):
            print(' \n converged after ', iter,"iterations\n")
            break
        r = Y - X*(np.tile(b,(dimcbf[0],1))) 
        radj = r * adjfactor/sw
        if flagstd == 1 :
            s=np.sqrt(np.mean(np.power(radj,2),axis=0)) 
        else:
            rs=np.sort(np.abs(radj),axis=0); s=np.median(rs,axis=0)/0.6745 
        r1=radj*(1-flagmodrobust*np.exp(-np.tile(modrobprior,(dimcbf[0],1))))/np.tile(np.maximum(s,tiny_s)*tune,(dimcbf[0],1))
        w,_=_weightfun(r1,wfun)
        b0=b; z=np.sqrt(w); x = X*z; yz = Y*z
        b=(np.sum(x*yz,axis=0)+mu*Globalprior+lmd*localprior)/(np.sum(x*x,axis=0)+mu+lmd)

    return b
        

def _scrubcbf(cbf_ts,gm,wm,csf,mask,wfun='huber',thresh=0.7):
    gm=mask*gm;wm=mask*wm; csf=csf*mask
    gmidx=gm[mask==1]; gmidx[gmidx<thresh]=0; gmidx[gmidx>0] = 1
    wmidx=wm[mask==1]; wmidx[wmidx<thresh]=0; wmidx[wmidx>0] = 1
    csfidx = csf[mask==1]; csfidx[csfidx<thresh] = 0; csfidx[csfidx>0] =1
    #midx = mask[mask==1]
    meancbf=np.mean(cbf_ts,axis=3)
    y=np.transpose(cbf_ts[mask==1,:,])
    VV=np.var(y,axis=0)
    thresh1,thresh3=_getchisquare(y.shape[0])
    mu1=VV/(np.median(VV[gmidx==1])*thresh3)
    mu =((mu1>thresh1)&(mu1<10*thresh1))*(mu1-thresh1) +(mu1 >=10*thresh1)*(1/(2*thresh1*10)*np.power(mu1,2))+(thresh1*10/2 - thresh1)
    M=meancbf*mask; M[mask==1]=mu; modrobprior = mu/10
    gmidx2 = 1*([gm.flatten()>thresh] and [M.flatten()==0] and [wm.flatten() > csf.flatten()])[0]
    wmidx2 =1*([wm.flatten()>thresh] and [M.flatten()==0] and [gm.flatten() > csf.flatten()])[0]
    if np.sum(gmidx2)==0 or np.sum(wmidx2)==0:
        gmidx2 =1*(gm.flatten()>thresh); wmidx2 = 1*(wm.flatten()>thresh)
    idxx =gmidx2 + wmidx2; idxx[idxx>0]=1
    X = np.zeros([len(idxx),2])
    X[:,0] = gm.flatten()[gm.flatten()>=(0)]*idxx
    X[:,1] = wm.flatten()[wm.flatten()>=(0)]*idxx
    A=(meancbf.flatten()[idxx >= 0])*idxx 
    c=np.linalg.lstsq(X,A)[0]
    Globalpriorfull=c[0]*gm.flatten() +c[1]*wm.flatten()
    Globalprior =Globalpriorfull[mask.flatten()==1]
    localprior=0;lmd=0
    tune=_tune(wfun=wfun)
    
    bb=_roubustfit(Y=y,mu=mu,Globalprior=Globalprior,modrobprior=modrobprior,lmd=lmd,
           localprior=localprior,wfun=wfun,tune=tune,flagstd=1,flagmodrobust=1,flagprior=1,thresh=0.7)
    newcbf=meancbf*mask
    newcbf[mask==1]=bb
    
    return newcbf



# basil and pvcorr
class _BASILCBFInputSpec(FSLCommandInputSpec):
    # We use position args here as list indices - so a negative number
    # will put something on the end
    in_file = File(
        exists=True,
        desc="input file cbf after substracting tag-control or control-tag",
        argstr=" -i %s",
        position=0,
        mandatory=True,
    )
    mask= File(exists=True,argstr=" -m %s ",desc="mask in the same space as in_infile",mandatory=True,)
    mzero=File(exists=True,argstr=" -c %s ",desc='m0 scan',mandatory=True,)
    m0scale=traits.Float(desc='calibration of asl',argstr=" --cgain %.2f ",mandatory=True,)
    m0tr=traits.Float(desc='Mzero TR',argstr=" --tr %.2f ",mandatory=True,)
    tis=traits.Float(desc='invertion recovery time =plds+bolus',argstr=" --tis %.2f ",mandatory=True,)
    pcasl=traits.Bool(desc='label type:defualt is PASL',argstr=" --casl ",mandatory=False,default_value=False)
    bolus=traits.Float(desc='bolus or tau: label duraytion',argstr=" --bolus %.2f ",mandatory=True,)
    pvc=traits.Bool(desc='calibration of asl',mandatory=False,argstr=" --pvcorr ",default_value=True)
    pvgm=File(exists=True,mandatory=False,desc='grey matter probablity matter ',argstr=" --pvgm %s ",)
    pvwm=File(exists=True,mandatory=False,desc='white matter probablity matter ',argstr=" --pvwm %s ",)
    out_basename=File(desc="base name of output files", argstr=" -o %s ",mandatory=True,) 
    #environ=traits.Str('FSLOUTPUTTYPE': 'NIFTI_GZ'}

class _BASILCBFOutputSpec(TraitedSpec):
    out_cbfb=File(exists=False,mandatory=True,desc='cbf with spatial correction')
    out_cbfpv=File(exists=False,mandatory=False,desc='cbf with spatial correction')
    out_attb=File(exists=False,mandatory=False,desc='aretrial transist time')
   
class BASILCBF(FSLCommand):
    _cmd = " oxford_asl "
    input_spec = _BASILCBFInputSpec
    output_spec = _BASILCBFOutputSpec

    def _run_interface(self, runtime):
        import shutil
        if os.path.isdir(self.inputs.out_basename+'/native_space'):
            shutil.rmtree(self.inputs.out_basename+'/native_space')
            shutil.rmtree(self.inputs.out_basename+'/calib')    
        runtime = super(BASILCBF, self)._run_interface(runtime)
        outputs = self.input_spec().get()
        outputs["out_cbfb"]=self.inputs.out_basename+'/basilcbf.nii.gz'
        from shutil import copyfile
        copyfile(self.inputs.out_basename+'/native_space/perfusion_calib.nii.gz',outputs["out_cbf"])
        if len(np.array([self.inputs.tis])) > 1:
            outputs["out_att"]=self.inputs.out_basename+'/arrivaltime.nii.gz'
            copyfile(self.inputs.out_basename+'/native_space/arrival.nii.gz',outputs["out_att"])
        else:
            outputs["out_attb"]=Undefined
        if self.inputs.pvc:
            outputs["out_cbfpv"]=self.inputs.out_basename+'/basilcbfpv.nii.gz'
            copyfile(self.inputs.out_basename+'/native_space/pvcorr/perfusion_calib.nii.gz',outputs["out_cbfpv"])
        else:
            outputs["out_cbfpv"]=Undefined  
        if runtime.stderr:
             self.raise_exception(runtime)
        shutil.rmtree(self.inputs.out_basename+'/native_space')
        shutil.rmtree(self.inputs.out_basename+'/calib')
        return runtime

    
