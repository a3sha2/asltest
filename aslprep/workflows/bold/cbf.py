from nipype.pipeline import engine as pe
from nipype.interfaces import utility as niu, fsl
from ...niworkflows.engine.workflows import LiterateWorkflow as Workflow
from ...niworkflows.interfaces import NormalizeMotionParams
from ...niworkflows.interfaces.fixes import FixHeaderApplyTransforms as ApplyTransforms
from ...niworkflows.interfaces.itk import MCFLIRT2ITK
from ...niworkflows.interfaces.cbf_computation import (extractCBF,computeCBF
       ,scorescrubCBF,BASILCBF,refinemask)
import nibabel as nb 
import numpy as np
import os,sys
from ...config import DEFAULT_MEMORY_MIN_GB


def init_cbf_compt_wf(mem_gb,metadata,aslcontext,pcasl,omp_nthreads, name='cbf_compt_wf'):
    workflow = Workflow(name=name)
    workflow.__desc__ = """\
        it is coming 
        """
    


    inputnode = pe.Node(niu.IdentityInterface(
        fields=['bold', 'bold_mask','t1w_tpms','t1w_mask','t1_bold_xform']),
        name='inputnode')
    outputnode = pe.Node(niu.IdentityInterface(
        fields=['out_cbf', 'out_mean','out_score','out_avgscore','out_scrub',
             'out_scoreindex','out_cbfb','out_cbfpv']),
        name='outputnode')


    
    # convert tmps to bold_space
    csf_tfm = pe.Node(ApplyTransforms(interpolation='NearestNeighbor', float=True),
                      name='csf_tfm', mem_gb=0.1)
    wm_tfm = pe.Node(ApplyTransforms(interpolation='NearestNeighbor', float=True),
                     name='wm_tfm', mem_gb=0.1)
    gm_tfm = pe.Node(ApplyTransforms(interpolation='NearestNeighbor', float=True),
                     name='gm_tfm', mem_gb=0.1)
    
     
    labeltype=metadata['LabelingType']
    if 'CASL' in labeltype: 
        pcasl=True
    elif 'PASL' in labeltype:
        pcasl=False
    else:
        print('unknown label type')

        
    extractcbf = pe.Node(extractCBF(in_ASLcontext=aslcontext),mem_gb=0.2,run_without_submitting=True,name="extractcbf") 
    computecbf = pe.Node(computeCBF(in_metadata=metadata),mem_gb=0.2,
              run_without_submitting=True,name="computecbf")
    scorescrub= pe.Node(scorescrubCBF(in_thresh=0.7,in_wfun='huber'),
              name='scorescrub',run_without_submitting=True,mem_gb=0.2)
    basilcbf= pe.Node(BASILCBF(m0scale=metadata["M0"],
               bolus=metadata["InitialPostLabelDelay"],m0tr=metadata['RepetitionTime'],pvc=True,
               tis=np.add(metadata["InitialPostLabelDelay"],metadata["LabelingDuration"]),
               pcasl=pcasl,out_basename=os.getcwd()),
              name='basilcbf',run_without_submitting=True,mem_gb=0.2) 

    refinemaskj=pe.Node(refinemask(),mem_gb=0.2,run_without_submitting=True,name="refinemask")

    
    #def _getTR(file):
        #import nibabel as nb
        #motr=nb.load(file).header.get_zooms()[3]
        #return motr
    
    def _pick_csf(files):
        return files[0]
    
    def _pick_gm(files):
        return files[1]

    def _pick_wm(files):
        return files[-1]

    
    workflow.connect([
        # extract CBF data and compute cbf
        (inputnode,  extractcbf, [('bold','in_file')]),
        (extractcbf, computecbf, [('out_file','in_cbf'),('out_avg','in_m0file')]),
        #(inputnode,computecbf,[('bold_mask','in_mask')]),
        (inputnode,refinemaskj,[('t1w_mask','in_t1mask'),('bold_mask','in_boldmask'),
                                ('t1_bold_xform','transforms')]),
        (inputnode,computecbf,[('bold_mask','in_mask')]),
        (inputnode,scorescrub,[('bold_mask','in_mask')]),
        (inputnode,basilcbf,[('bold_mask','mask')]),

        # extract probability maps
        (inputnode, csf_tfm, [('bold_mask', 'reference_image'),
                              ('t1_bold_xform', 'transforms')]),
        (inputnode, csf_tfm, [(('t1w_tpms', _pick_csf), 'input_image')]),
        (inputnode, wm_tfm, [('bold_mask', 'reference_image'),
                              ('t1_bold_xform', 'transforms')]),
        (inputnode, wm_tfm, [(('t1w_tpms', _pick_wm), 'input_image')]),
        (inputnode, gm_tfm, [('bold_mask', 'reference_image'),
                              ('t1_bold_xform', 'transforms')]),
        (inputnode, gm_tfm, [(('t1w_tpms', _pick_gm), 'input_image')]),
        (computecbf,scorescrub,[('out_cbf','in_file')]),
        (gm_tfm,scorescrub,[('output_image','in_greyM')]),
        (wm_tfm,scorescrub,[('output_image','in_whiteM')]),
        (csf_tfm,scorescrub,[('output_image','in_csf')]),
        #(inputnode,scorescrub,[('bold_mask','in_mask')]),
        (extractcbf,basilcbf,[('out_file','in_file')]),
        (gm_tfm,basilcbf,[('output_image','pvgm')]),
        (wm_tfm,basilcbf,[('output_image','pvwm')]),
        #(inputnode,basilcbf,[('bold_mask','mask')]),
        (extractcbf,basilcbf,[('out_avg','mzero')]),
        (basilcbf,outputnode,[('out_cbfb','out_cbfb'),
                ('out_cbfpv','out_cbfpv')]),
        (computecbf,outputnode,[('out_cbf','out_cbf'),
                     ('out_mean','out_mean')]),
        (scorescrub,outputnode,[('out_score','out_score'),('out_scoreindex','out_scoreindex'),
                    ('out_avgscore','out_avgscore'),('out_scrub','out_scrub')]),
         ])

    return workflow


        
                  
        


    


    


